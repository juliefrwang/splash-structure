import numpy as np
import pandas as pd
from math import comb
import itertools
import sys
from pandarallel import pandarallel

### 1. Target p computation ###
def target_p1_closed_form(k, v, L, c):
    """
    k: k-mer length
    v: total mutations
    L: stem length
    c: pair of compensatory mutations
    """
    p_c = 0
    for h in range(2*c, min(v, 2*L)+1):
        # print(f'h: {h}')
        p_h  = comb(2*L, h) * comb(k-2*L, v-h) / comb(k, v)
        p_c_h = 0
        for g in range(c, h//2 +1):
            p_g = 0
            for l in range(g, min(h+1, L+1)):
                l = max(l, h-l)
                if l > L:
                    continue
                sum_m = 0
                for m in range(0, min(l-g+1, h-l-g+1)):
                    sum_m += comb(l-g, m) * 2**m * comb(L-l, h-l-g-m) * 3**(h-l-g-m) 
                p_g += comb(L, l) * 3**l * comb(l, g)* sum_m / comb(2*L, h) / 3**h 
            p_c_h += p_g
        p_c += p_h * p_c_h
    return p_c

def target_p(k, stemL, totaMut, stemMut, compMut):
    """
    Return taregt p-value:
    p_1: exact p-val found using lookup table `dt` or approximate p for longer stem
    p_2: no stem mutations 
    combine multiple p: (stemMut > 0) * p_1 + (stemMut == 0) * p_2
    """

    p_1 = target_p1_closed_form(k, totaMut, stemL, compMut)
    p_2 = comb(k - 2 * stemL, totaMut) / comb(k, totaMut)
    p = (stemMut > 0) * p_1 + (stemMut == 0) * p_2
    return p

### 2. Anchor p computation ###
def target_p_outcome(k, stemL, totaMut):
    """
    Step 1: calculate each outcome of target_p
    target_p is computed on condition of k, stemL, totaMut

    Output: 
    all_possible_outcome: list of all possible outcomes of target_p
    """
    all_possible_outcome = set()
    stemMut_start = 0 if totaMut - (k - 2 * stemL) < 0 else totaMut - (k - 2 * stemL)
    for stemMut in range(stemMut_start, min(totaMut, 2*stemL)+1):
        for compMut in range((stemMut+2)//2):
            all_possible_outcome.add(target_p(k, stemL, totaMut, stemMut, compMut))
    all_possible_outcome = list(all_possible_outcome)
    all_possible_outcome.sort()
    return all_possible_outcome

def prep_for_conv(num_target, wgt_all, k, stemL_list, totaMut_list):
    """
    Step 2: prep for convolution: calculate PMF of each outcome of target_p
    
    Output: 
    wgted_target_outcomes: nested lists of (weighted) target_p for all targets of an anchor
    target_pmf: 
    """
    if num_target > 4: # cap number of targets for convolution to 4
        wgt_all = [i / sum(wgt_all[0:4]) for i in wgt_all[0:4]]
        stemL_list = stemL_list[0:4]
        totaMut_list = totaMut_list[0:4]
        num_target = 4
        
    wgted_target_outcomes = []
    target_pmf = [] 
    
    for i in range(num_target):
        targetp = target_p_outcome(k, stemL_list[i], totaMut_list[i])
        pmf = [targetp[0]]+[targetp[i+1] - targetp[i] for i in range(len(targetp)-1)]

        target_pmf.append(pmf)
        wgted_target_outcomes.append([wgt_all[i] * j for j in targetp])
    return wgted_target_outcomes, target_pmf

def pmf_anchor_score(wgted_target_outcomes, target_pmf):
    """
    Step 3: calculate the probability of each outcome of 
    
    Output: 
    anchor_p = weighted_average(target_p, anchor_p)
    """
    # all outcomes
    all_anchor_outcomes = [sum(x) for x in itertools.product(*wgted_target_outcomes)] 
    # all pmf
    anchor_pmf = [np.prod(x) for x in itertools.product(*target_pmf)] 
    return all_anchor_outcomes, anchor_pmf

def anchor_p(all_anchor_outcomes, anchor_pmf, anchor_p):
    """
    Step 4: calculate the CDF of anchor_p and find p-value of anchor_p
    """
    df = pd.DataFrame({'outcome':all_anchor_outcomes, 'prob':anchor_pmf})
    df=df.sort_values(by=['outcome']).reset_index(drop=True)
    df['cdf'] = df['prob'].cumsum()
    p_val = df[df.outcome <= anchor_p + 1e-6]['cdf'].max()
    if p_val is np.nan:
        p_val = df.iloc[0]['cdf']
    return p_val
    
def anchor_p_target_subdf(sub_df):
    """
    Step 5 (1): wrap all functions for one anchor and apply to sub-dataframe for stucture-target
    """
    p_val = sub_df['anchor_score'].iloc[0]
        
    if len(sub_df) > 1:
        all_anchor_outcomes, anchor_pmf = pmf_anchor_score(*prep_for_conv(len(sub_df),\
                                                      list(sub_df['tar_wgt_filtered']), \
                                                      len(sub_df['base_target'].iloc[0]), \
                                                      list(sub_df['stemL']), \
                                                      list(sub_df['totaMut'])))
            
        p_val = anchor_p(all_anchor_outcomes, anchor_pmf, p_val)
        
    return p_val

def anchor_p_compactor_subdf(sub_df):
    """
    Step 5 (2): wrap all functions for one anchor-split and apply to each anchor for stucture-compactor
    """
    # for compactors, we compute anchor-score for each split
    p_val = sub_df['anchor_score_per_split'].iloc[0]
    
    if len(sub_df) > 1:
        # structure evaluation length for compactor is 80 (HARDCODED)
        all_anchor_outcomes, anchor_pmf = pmf_anchor_score(*prep_for_conv(len(sub_df),\
                                                      list(sub_df['target_weight']), \
                                                      80, \
                                                      list(sub_df['stemL']), \
                                                      list(sub_df['totaMut'])))
        p_val = anchor_p(all_anchor_outcomes, anchor_pmf, p_val)
        
    return p_val

def wrap_anchor_p_target(df):
    """
    Step 6 (1): wrap all functions for one anchor and apply to the whole dataframe for stucture-target
    """
    grouped = df.groupby('anchor') # create a groupby object based on 'anchor'
    p_val_results = grouped.parallel_apply(anchor_p_target_subdf) # Apply function to each group
    # The result is a Series where the index is the group keys ('anchor' values)
    # We can now assign this back to your DataFrame, but you'll need to align the indices
    df = df.merge(p_val_results.rename('anchor_p'), left_on='anchor', right_index=True)
    return df

def wrap_anchor_p_compactor(df):
    """
    Step 6 (1): wrap all functions for one anchor and apply to the whole dataframe for stucture-compactor
    """
    grouped = df.groupby('anchor_split') # create a groupby object based on 'anchor_split'
    p_val_results = grouped.parallel_apply(anchor_p_compactor_subdf) # Apply function to each group
    # The result is a Series where the index is the group keys ('anchor_split' values)
    # We can now assign this back to your DataFrame, but you'll need to align the indices
    df = df.merge(p_val_results.rename('anchor_p'), left_on='anchor_split', right_index=True)
    return df