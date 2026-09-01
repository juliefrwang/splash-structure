"""Core simulation primitives: plant a stem-loop, generate SVE/SVP/null targets,
build abundance vectors. All randomness flows through a numpy ``Generator``.
"""
from __future__ import annotations

import numpy as np
from typing import Iterable

from struct_rna.src.non_wcf import V_EXT, BASES
from struct_rna.src.process_targets import find_stem_ind_wobble

TARGET_LEN = 27
STEM_L_RANGE = (5, 8)  # uniform inclusive
LOOP_MIN = 2

# Pair p (0-indexed from outermost) sits at:
#   left position  = stem_start + p
#   right position = rc_end - p     (matches find_mutation_ext)
# So pair 0 is the outermost (closest to the flanks) and pair L-1 is innermost
# (closest to the loop). This convention is fixed by find_mutation_ext.


# ---------------------------------------------------------------------
# Geometry + base-target construction
# ---------------------------------------------------------------------

def random_geometry(rng: np.random.Generator) -> tuple[int, int, int, int, int, int]:
    """Draw (L, loop_len, stem_start, stem_end, rc_start, rc_end) per the design."""
    L = int(rng.integers(STEM_L_RANGE[0], STEM_L_RANGE[1] + 1))
    max_loop = TARGET_LEN - 2 * L
    if max_loop < LOOP_MIN:
        raise ValueError(f"L={L} leaves no room for a loop ≥ {LOOP_MIN}")
    loop_len = int(rng.integers(LOOP_MIN, max_loop + 1))
    max_start = TARGET_LEN - 2 * L - loop_len
    stem_start = int(rng.integers(0, max_start + 1))
    stem_end = stem_start + L - 1
    rc_start = stem_end + 1 + loop_len
    rc_end = rc_start + L - 1
    return L, loop_len, stem_start, stem_end, rc_start, rc_end


def _plant_once(rng: np.random.Generator, valid=V_EXT) -> dict:
    """One attempt at planting: draw geometry, fill stem from ``valid``
    (the data-side V(N); default V_EXT = V_WCF ∪ G·U), randomize
    flank/loop. Returns the dict; caller checks recovery."""
    L, loop_len, stem_start, stem_end, rc_start, rc_end = random_geometry(rng)
    seq = list(rng.choice(BASES, size=TARGET_LEN))
    # sorted(), not list(): `valid` is a frozenset, and set iteration order
    # for string tuples varies with PYTHONHASHSEED across processes, which
    # made the simulator irreproducible at fixed seed.
    pair_list = sorted(valid)
    pair_idx = rng.integers(0, len(pair_list), size=L)
    n_wobble = 0
    for p in range(L):
        bL, bR = pair_list[int(pair_idx[p])]
        seq[stem_start + p] = bL
        seq[rc_end - p] = bR
        if (bL, bR) not in {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}:
            n_wobble += 1
    return {
        "sequence": "".join(seq),
        "L": L,
        "loop_len": loop_len,
        "stem_start": stem_start,
        "stem_end": stem_end,
        "rc_start": rc_start,
        "rc_end": rc_end,
        "n_wobble": n_wobble,
    }


def make_base_target(rng: np.random.Generator, max_tries: int = 200,
                     valid=V_EXT) -> dict:
    """Generate a 27-nt base target with a planted stem-loop, verifying that
    ``find_stem_ind_wobble`` recovers the planted indices.

    ``valid`` is the data-side V(N): stem pairs are drawn from it and the
    recovery check uses it, so a planted N-pair stem is detectable.
    Default ``V_EXT`` ⇒ byte-identical to the prior G·U behaviour.

    The wobble stem-finder returns the LONGEST stem with the FEWEST wobble
    pairs. Random flanks/loop can occasionally form an alternative stem that
    outranks ours; in that case we regenerate the flanks/loop until recovery
    succeeds. n_wobble emerges naturally from the pair-uniform-V_EXT draw.
    """
    last_planted = None
    last_found = None
    for _ in range(max_tries):
        geom = _plant_once(rng, valid)
        seq = geom["sequence"]
        found = find_stem_ind_wobble(seq, stem_L=5, valid=valid)
        # find_stem_ind_wobble returns (s_start, s_end, rc_start, rc_end, L)
        if (
            found[0] == geom["stem_start"]
            and found[1] == geom["stem_end"]
            and found[2] == geom["rc_start"]
            and found[3] == geom["rc_end"]
            and found[4] == geom["L"]
        ):
            return geom
        last_planted = (geom["stem_start"], geom["stem_end"], geom["rc_start"], geom["rc_end"], geom["L"])
        last_found = tuple(found)
    raise RuntimeError(
        f"failed to plant a recoverable stem after {max_tries} tries "
        f"(last planted={last_planted}, last found={last_found})"
    )


# ---------------------------------------------------------------------
# Mutation primitives
# ---------------------------------------------------------------------

def _stem_positions(geom: dict) -> set[int]:
    return set(range(geom["stem_start"], geom["stem_end"] + 1)) | set(
        range(geom["rc_start"], geom["rc_end"] + 1)
    )

# ---------------------------------------------------------------------
# Ti/Tv-aware mutation identity sampling
# ---------------------------------------------------------------------
#
# Ti/Tv ratio R = #Ti events / #Tv events under the per-mutation null.
# Per-mutation identity weights (each base has 1 Ti and 2 Tv alternatives):
#   P(Ti)       = R / (R + 1)
#   P(each Tv)  = 1 / (2 * (R + 1))
# At R = 0.5: P(Ti) = 1/3, P(each Tv) = 1/3 → uniform.
# At R = 2  : P(Ti) = 2/3, P(each Tv) = 1/6 → realistic biology.

_TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def _mut_weight(base: str, new_base: str, ti_tv_ratio: float) -> float:
    """Weight (unnormalised) of mutation base → new_base under Ti/Tv ratio."""
    if new_base == base:
        return 0.0
    if new_base == _TRANSITION[base]:
        return ti_tv_ratio
    return 1.0 / 2.0  # each of the two Tv alternatives shares mass 1, so 1/2 each


def _pick_alt_titv(
    rng: np.random.Generator, base: str, ti_tv_ratio: float
) -> str:
    """Sample a mutation identity ≠ base, weighted by Ti/Tv ratio."""
    alts = [b for b in BASES if b != base]
    weights = np.array([_mut_weight(base, a, ti_tv_ratio) for a in alts])
    weights = weights / weights.sum()
    return str(rng.choice(alts, p=weights))


def _pick_among_titv(
    rng: np.random.Generator,
    base: str,
    candidates: list[str],
    ti_tv_ratio: float,
) -> str:
    """Sample from a constrained list of mutation identities (each ≠ base),
    weighted by Ti/Tv ratio."""
    weights = np.array([_mut_weight(base, c, ti_tv_ratio) for c in candidates])
    weights = weights / weights.sum()
    return str(rng.choice(candidates, p=weights))


def _pick_pair_titv(
    rng: np.random.Generator,
    base_L: str,
    base_R: str,
    candidates: list[tuple[str, str]],
    ti_tv_ratio: float,
) -> tuple[str, str]:
    """Sample a joint (new_L, new_R) from V_EXT-preserving candidates, weighted
    by joint Ti/Tv probability (independent per side)."""
    weights = np.array(
        [
            _mut_weight(base_L, c[0], ti_tv_ratio)
            * _mut_weight(base_R, c[1], ti_tv_ratio)
            for c in candidates
        ]
    )
    weights = weights / weights.sum()
    idx = int(rng.choice(len(candidates), p=weights))
    return candidates[idx]


def make_sve_target(
    rng: np.random.Generator,
    base: str,
    geom: dict,
    v: int,
    ti_tv_ratio: float = 0.5,
) -> str:
    """v mutations placed entirely outside the stem; identities Ti/Tv-weighted."""
    non_stem = [i for i in range(TARGET_LEN) if i not in _stem_positions(geom)]
    if v > len(non_stem):
        raise ValueError(f"v={v} exceeds available non-stem positions {len(non_stem)}")
    chosen = rng.choice(non_stem, size=v, replace=False)
    seq = list(base)
    for i in chosen:
        seq[i] = _pick_alt_titv(rng, base[i], ti_tv_ratio)
    return "".join(seq)


def make_svp_target(
    rng: np.random.Generator,
    base: str,
    geom: dict,
    v: int,
    ti_tv_ratio: float = 0.5,
    max_tries: int = 200,
    valid=V_EXT,
) -> str:
    """v mutations with at least one in the stem; all stem mutations preserve
    the data-side valid set ``valid`` = V(N) (default V_EXT).

    Identity choices are Ti/Tv-weighted: per-side SPC alternatives use
    `_pick_among_titv`; joint BPC alternatives use `_pick_pair_titv`. The free
    (non-stem) mutations also use Ti/Tv weighting.
    """
    stem_set = _stem_positions(geom)
    L = geom["L"]
    stem_start = geom["stem_start"]
    rc_end = geom["rc_end"]

    for _ in range(max_tries):
        positions = rng.choice(TARGET_LEN, size=v, replace=False)
        pos_set = set(int(x) for x in positions)
        stem_hit = pos_set & stem_set
        if not stem_hit:
            continue

        # Classify hits per pair: side(s) hit and base composition.
        pair_sides: dict[int, str] = {}
        for pos in stem_hit:
            if stem_start <= pos <= geom["stem_end"]:
                p = pos - stem_start
                side = "L"
            else:
                p = rc_end - pos
                side = "R"
            pair_sides[p] = pair_sides.get(p, "") + side

        seq = list(base)
        feasible = True
        for p, sides in pair_sides.items():
            bL = base[stem_start + p]
            bR = base[rc_end - p]
            if len(sides) == 1:
                if sides == "L":
                    cands = [b for b in BASES if b != bL and (b, bR) in valid]
                    if not cands:
                        feasible = False
                        break
                    seq[stem_start + p] = _pick_among_titv(rng, bL, cands, ti_tv_ratio)
                else:  # "R"
                    cands = [b for b in BASES if b != bR and (bL, b) in valid]
                    if not cands:
                        feasible = False
                        break
                    seq[rc_end - p] = _pick_among_titv(rng, bR, cands, ti_tv_ratio)
            else:  # n_p == 2
                cands = [
                    (x, y)
                    for x in BASES
                    if x != bL
                    for y in BASES
                    if y != bR and (x, y) in valid
                ]
                if not cands:
                    feasible = False
                    break
                nbL, nbR = _pick_pair_titv(rng, bL, bR, cands, ti_tv_ratio)
                seq[stem_start + p] = nbL
                seq[rc_end - p] = nbR
        if not feasible:
            continue

        for pos in sorted(pos_set - stem_set):
            seq[pos] = _pick_alt_titv(rng, base[pos], ti_tv_ratio)

        return "".join(seq)

    raise RuntimeError(f"SVP target infeasible after {max_tries} tries (v={v}, L={L})")


def make_null_target(
    rng: np.random.Generator, base: str, v: int, ti_tv_ratio: float = 0.5
) -> str:
    """v mutations placed uniformly across all 27 positions; identities
    Ti/Tv-weighted. No preservation constraint — used for null anchors."""
    chosen = rng.choice(TARGET_LEN, size=v, replace=False)
    seq = list(base)
    for i in chosen:
        seq[i] = _pick_alt_titv(rng, base[i], ti_tv_ratio)
    return "".join(seq)


# ---------------------------------------------------------------------
# Abundance
# ---------------------------------------------------------------------

def abundance_counts(K: int, kind: str, total: int = 1000) -> np.ndarray:
    """Integer counts of length K summing to ~total (most-abundant first).

    'uniform' → all equal; 'skewed' → geometric drop-off ratio 0.5 (matches
    the abundance test in the design doc).
    """
    if kind == "uniform":
        base = total // K
        out = np.full(K, base, dtype=int)
        out[0] += total - base * K  # absorb remainder into rank-1
        return out
    elif kind == "skewed":
        weights = np.array([0.5 ** i for i in range(K)])
        weights /= weights.sum()
        out = (weights * total).astype(int)
        out[0] += total - int(out.sum())
        out = np.clip(out, 1, None)  # avoid zero-count slots
        return out
    else:
        raise ValueError(f"unknown abundance kind: {kind!r}")


# ---------------------------------------------------------------------
# Anchor assembly
# ---------------------------------------------------------------------

MAX_TARGET_SLOTS = 10  # SPLASH output convention


def build_anchor_row(
    anchor: str,
    base_target: str,
    less_abundant_targets: list[str],
    counts: np.ndarray,
) -> dict:
    """One row in the SPLASH-significant-anchors TSV format.

    Slots above K (= 1 + len(less_abundant_targets)) are padded with '-' / 0.
    """
    K = 1 + len(less_abundant_targets)
    if K > MAX_TARGET_SLOTS:
        raise ValueError(f"K={K} exceeds MAX_TARGET_SLOTS={MAX_TARGET_SLOTS}")
    if len(counts) != K:
        raise ValueError(f"counts length {len(counts)} != K {K}")
    row = {"anchor": anchor, "M": int(counts.sum())}
    row["most_freq_target_1"] = base_target
    row["cnt_most_freq_target_1"] = int(counts[0])
    for i, t in enumerate(less_abundant_targets, start=2):
        row[f"most_freq_target_{i}"] = t
        row[f"cnt_most_freq_target_{i}"] = int(counts[i - 1])
    for i in range(K + 1, MAX_TARGET_SLOTS + 1):
        row[f"most_freq_target_{i}"] = "-"
        row[f"cnt_most_freq_target_{i}"] = 0
    return row


def random_anchor_seq(rng: np.random.Generator) -> str:
    return "".join(rng.choice(BASES, size=TARGET_LEN))
