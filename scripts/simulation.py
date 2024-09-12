import sys
import os
from statsmodels.stats.multitest import multipletests
import pandas as pd
from multiprocessing import Pool

sys.path.append('/oak/stanford/groups/horence/juliew/structure/src')
import simulate_target
import structure_plot

DATA_HANDLE = sys.argv[1]
SIMULATIOM_ITER = sys.argv[2]
STURCTURE_FILE = f"/oak/stanford/groups/horence/juliew/structure/0204_results/{DATA_HANDLE}_results/structure_on_targets.tsv"
SIMULATION_FILE = f"/oak/stanford/groups/horence/juliew/structure/0204_results/{DATA_HANDLE}_results/structure_simulation_on_test_targets_{SIMULATIOM_ITER}.tsv"

df = pd.read_csv(STURCTURE_FILE, sep = '\t') # read in the structure_on_targets.tsv
unique_anchors = df['anchor'].unique()
n_iter = 100
arg_generator = ((anchor, df, n_iter) for anchor in unique_anchors)
# with Pool(cpu_count()-1) as p: 
with Pool(16) as p:
    results = p.starmap(simulate_target.get_simulated_p, arg_generator)
df_anchor_result = pd.DataFrame({'anchor': unique_anchors, 'anchor_p_simulated': results})
df = df.merge(df_anchor_result, on='anchor')

#df_temp = df[['anchor', 'anchor_p_simulated']].drop_duplicates()
#df = df.merge(df_temp, on='anchor', how = 'left') # Merge back
# expand the anchor_p_simulated to multiple columns
df = pd.concat([df, df['anchor_p_simulated'].apply(pd.Series)], axis=1)
df = df.drop(columns=['anchor_p_simulated'])

#for i in range(n_iter):
#    df_temp[i] = df_temp[i].astype(float)
#    correction = multipletests(df_temp[i], alpha=0.05, method='fdr_bh')
#    df_temp[f'p_simu_{i}_BH'] = correction[1]
#    df_temp = df_temp.drop(columns=[i])

#df = df.merge(df_temp, on='anchor', how = 'left') # Merge back
# Step4 save to folder
df.to_csv(SIMULATION_FILE, index=False, sep='\t')
