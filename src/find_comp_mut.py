"""
This script is to find compensatory mutations in the target sequence, along with 
all parameters that are crucial for calculating structure p-value.
The script also provides three ways for strucutre notation.
"""
import numpy as np
import pandas as pd

def rc(seq):
    """
    Take in sequence and return the reverse complement of the given sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def find_mutation(base, target, stem_start_idx, stem_end_idx, rc_start_idx, rc_end_idx):
    """
    This functions takes in base target, target and hairpin structure indices,
    and return quanties that are crucial for calculating structure p-value, which include:
    1. totaMut: total number of mutations
    2. stemMut: number of mutations in the stem
    3. compMut: number of compensatory mutations
    4. struc: structure notation
    """
    struc = ['-'] * len(base)

    totaMut = 0
    stemMut = 0
    compMut = 0                      # pair of compensatory mutations
    sumInd = stem_start_idx + rc_end_idx # sum of compensatory stem indices 
    right_checked = []
    for i in range(len(target)):
        if base[i] != target[i]:
            totaMut += 1
            if i in range(stem_start_idx, stem_end_idx+1):          # i in left stem
                comInd = sumInd - i                           # calculate the corresponding index in the right stem 
                if base[comInd] != target[comInd]:            # if the corresponding index in the right stem is mutated
                    stemMut += 1
                    if rc(target[i]) == target[comInd]:       # if the right mutation is compensatory to the left mutation. 
                        compMut += 1
                        struc[i] = target[i]
                        struc[comInd] = target[comInd]
                        right_checked.append(comInd)
                    else:
                        struc[i] = target[i].lower()
                else:                                         # if the corresponding index in the right stem is not mutated
                    stemMut += 1
                    struc[i] = target[i].lower()
            elif i in range(rc_start_idx, rc_end_idx+1):        # i in right stem, already checked right-stem mutations
                stemMut += 1
                if i not in right_checked:
                    struc[i] = target[i].lower()
            else:                                             # i outside stem                   
                stemMut += 0
                struc[i] = target[i].lower()

    struc = (struc[:stem_start_idx] + 
             ['{'] + struc[stem_start_idx: stem_end_idx+1] + 
             ['('] + struc[stem_end_idx + 1: rc_start_idx] + [')'] + 
             struc[rc_start_idx: rc_end_idx+1] + ['}'] + 
             struc[rc_end_idx+1:])

    struc = "".join(struc)
    
    # add a patch for compactor segments: if stemL == 0 (stem_end_idx == 0, indeed)
    if stem_end_idx == 0:
        struc = np.nan

    return (totaMut, stemMut, compMut, struc)

def db_notation_from_idx(stem_start_idx, stem_end_idx, rc_start_idx, rc_end_idx, target_length):
    db_notation=''
    for i in range(target_length):
        if (i >= stem_start_idx) & (i <= stem_end_idx):
            db_notation += '('
        elif (i >= rc_start_idx) & (i <= rc_end_idx):
            db_notation += ')'
        else:
            db_notation += '.'
    return db_notation

def db_notation_from_old_notaion(old_notation):
    if pd.isna(old_notation):
        return np.nan
    stem_start_idx = old_notation.find('{')
    stem_end_idx = old_notation.find('(')-2
    rc_start_idx = old_notation.find(')')-2
    rc_end_idx =  old_notation.find('}')-4
    return db_notation_from_idx(stem_start_idx, stem_end_idx, rc_start_idx, rc_end_idx, len(old_notation)-4)

def symbol_notation_from_old_notaion(old_notation, db_notation=None):
    if pd.isna(old_notation):
        return np.nan
    if db_notation is None:
        db_notation = db_notation_from_old_notaion(old_notation)
    rc_start_idx = db_notation.find(')')
    symbol_nototation = list(db_notation)
    old_notation = [item for item in old_notation if item not in ['{','}','(',')']]
    for i, item in enumerate(old_notation):
        if item.islower():
            symbol_nototation[i] = '*'
        elif item.isupper():
            if i < rc_start_idx:
                symbol_nototation[i] = '<'
            else:  
                symbol_nototation[i] = '>'
    return "".join(symbol_nototation)