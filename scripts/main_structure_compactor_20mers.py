"""
The script takes in two arguments:
1. name of structure working folder
2. path to compactors.tsv 
3. optional: -a to run element annotation on compactors
"""
import sys
import os
import argparse
import pandas as pd
import subprocess
from statsmodels.stats.multitest import multipletests
from pandarallel import pandarallel

sys.path.append('/oak/stanford/groups/horence/juliew/structure/src')
import process_targets
import find_comp_mut
import get_pval
import elem_annas

# Parse command line arguments
parser = argparse.ArgumentParser(description='Run SPLASH-structure on compactor_20.')
parser.add_argument("-a", "--element_annotation", action="store_true", help="Run element annotation on compactors. Default is False", )
parser.add_argument("compactor_file", help="Path to compactors file.")
parser.add_argument("data_handle", help="Data handle for the output folder")
args = parser.parse_args()
COMPACTOR_FILE = args.compactor_file
DATA_HANDLE = args.data_handle

def main():
    """ Step 0: Preparation """
    # Initialize parallelization. Create folder to save results
    pandarallel.initialize()
    outfolder = f'{DATA_HANDLE}_results'
    os.makedirs(outfolder, exist_ok=True)
    os.makedirs(f'{outfolder}/int_file', exist_ok=True)

    """ Step 1: Read in compactors and process dataframe using Julia script """
    # Define the command to run the Julia script
    julia_command = f"julia /oak/stanford/groups/horence/juliew/structure/src/process_compactor_20mers.jl {COMPACTOR_FILE} {outfolder}/int_file"

    # Run the Julia script using subprocess and wait for it to finish
    completed_process = subprocess.run(julia_command, shell=True)

    """ Step 2: Read in processed compactors & add number of compactors """
    # Check if the Julia script ran successfully (exit code 0)
    if completed_process.returncode == 0:
        df = pd.read_csv(f"{outfolder}/int_file/processed_compactors_20mers.tsv", sep = '\t')
    else:
        print("Julia script encountered an error or did not finish successfully.")
        sys.exit(1) 

    # exit program if no structure is found in any target
    if len(df) == 0:
        print("No structure is found for any anchor. Exiting..")
        return
        
    """
    Step 3: Find parameters that are to be used in anchor-p computation, 
    along with three types of notations
    """
    # find stem loop index
    for i in [1,2,3,4]:
        df[[f"stem_start_idx_{i}", f"stem_end_idx_{i}", f"rc_start_idx_{i}", f"rc_end_idx_{i}", f"stemL_{i}"]] = pd.DataFrame(df.parallel_apply(lambda x: process_targets.find_stem_ind(x[f"base_S{i}"], 5), axis=1).tolist())

    # Add a column of number of stem-loop structure found in the compactors
    df['num_stem_loop'] = (df['stemL_1'] != 0).astype(int) + (df['stemL_2'] != 0).astype(int) + (df['stemL_3'] != 0).astype(int) + (df['stemL_4'] != 0).astype(int)

    # drop anchors without stem using condition num_stem_loop == 0 
    df = df[df.num_stem_loop != 0].reset_index(drop = True)

    # find mutations in stem & addition two columns for structure notations
    for i in [1,2,3,4]:
        df[[f"totaMut_{i}", f"stemMut_{i}", f"compMut_{i}", f"strucNotation_{i}"]] = pd.DataFrame(df.parallel_apply(lambda x: find_comp_mut.find_mutation(x[f'base_S{i}'], x[f'S{i}'], x[f"stem_start_idx_{i}"], x[f"stem_end_idx_{i}"], x[f"rc_start_idx_{i}"], x[f"rc_end_idx_{i}"]), axis=1).tolist())
        df[f"db_strucNotation_{i}"] = df.parallel_apply(lambda x: find_comp_mut.db_notation_from_old_notaion(x[f"strucNotation_{i}"]), axis=1) 
        df[f"symbol_strucNotation_{i}"]= df.parallel_apply(lambda x: find_comp_mut.symbol_notation_from_old_notaion(x[f"strucNotation_{i}"], x[f"db_strucNotation_{i}"]), axis=1)
    
    # compute target_p using summation of four segaments
    df['stemL'] = df[['stemL_1', 'stemL_2', 'stemL_3', 'stemL_4']].sum(axis=1)
    df['totaMut'] = df[['totaMut_1', 'totaMut_2', 'totaMut_3', 'totaMut_4']].sum(axis=1)
    df['stemMut'] = df[['stemMut_1', 'stemMut_2', 'stemMut_3', 'stemMut_4']].sum(axis=1)
    df['compMut'] = df[['compMut_1', 'compMut_2', 'compMut_3', 'compMut_4']].sum(axis=1)

    """ Step 4: Calculate structure target-p """
    df["target_p"] = df.parallel_apply(lambda x: get_pval.target_p(4 * len(x['base_S1']), x['stemL'], x['totaMut'], x['stemMut'], x['compMut']), axis=1) 

    """ Step 5: Calculate anchor_score_per_split """
    df["anchor_score_per_split"] = df["target_weight"] * df["target_p"]
    df["anchor_score_per_split"] = df.groupby(["anchor", "segment_index"])["anchor_score_per_split"].transform("sum")
    df['anchor_split'] = df.apply(lambda x: x['anchor'] + '_' + x['segment_index'], axis=1)
    
    """ Step 6: Calculate anchor_p """
    df = get_pval.wrap_anchor_p_compactor(df)

    """ Step 7: BH correction on anchors with number of compactor > 2 """
    # Filter the DataFrame to keep rows with 'anchor_split' counts greater than 2
    df = df[df.groupby('anchor_split')['compactor'].transform('count') > 2].reset_index(drop=True)
    
    # BH correction on anchors with number of num_compactor_split > 2 (this filter can be added before) 
    df_temp = df[['anchor_split', 'anchor_p']].drop_duplicates()
    correction = multipletests(df_temp['anchor_p'], alpha=0.05, method='fdr_bh')
    df_temp['anchor_p_BH'] = correction[1]
    
    # Merge back
    df_temp = df_temp.drop(columns=['anchor_p'])
    df = df.merge(df_temp, on='anchor_split', how = 'left')

    # Sort dataframe
    df = df.sort_values(by = ['anchor_p_BH', 'anchor_split'], ascending=True).reset_index(drop=True)

    """ Step 8: SAVE """
    df.to_csv(f'{outfolder}/structure_on_compactors_20mers.tsv', index=False, sep='\t')

    if args.element_annotation:
        """ Step 9: elememt annotations (optional, toggle on by -a)  """
        # run element annotations
        elem_ann_folder = elem_annas.run_anns(outfolder, "compactor", 20)

        """ Step 10: merge structure results with element annotations """
        anns_file = f"{outfolder}/elem_anns/{elem_ann_folder}_work/{elem_ann_folder}/element_annotations/element_annotations_anchors.tsv"
        df_anns = pd.read_csv(anns_file, sep='\t')
        df = elem_annas.merge_anns_struc(df_anns, df, "compactor")
        df.to_csv(f'{outfolder}/structure_on_compactors_20mers.tsv', index=False, sep='\t')

main()
