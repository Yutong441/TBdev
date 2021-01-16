'''
This script performs pseudotime analysis using stream.
Install stream using `conda install -c bioconda stream`.
If unsuccessful, need to downgrade to python 3.7 then conda 4.6.14.
References:
https://github.com/pinellolab/STREAM/issues/87
https://github.com/pinellolab/STREAM/blob/master/tutorial/1.1.STREAM_scRNA-seq%20(Bifurcation).ipynb
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata
import stream as st
import utils as pu

root=''
adata = pu.from_seurat_to_anndata (root+'data/', root+'result/')
adata = pu.Epilgraph_on_GPLVM (adata, root+'result/PT_no_prior.csv')

st.plot_dimension_reduction(adata,color=['revised', 'assigned_cluster'], n_components=3,
        show_graph=True,show_text=True)
plt.show ()
stem = 'S2'

# save branch point information for B-RGPLVM
adata.obs [stem+'_pseudotime'].to_csv (root+'data/STREAM_pseudotime.csv')
branch_label = list (adata.obs.branch_id_alias)

np.unique ([str(i) for i in branch_label])
branch_dict = {'branch1': ('S2', 'S1'), 'branch2': ('S0', 'S1'), 'branch3':
        ('S3', 'S1') }
for key, val in branch_dict.items():
    branch_label = [key if i == val else i for i in branch_label]

branch_label=pd.DataFrame (branch_label, index=adata.obs.index)
branch_label.to_csv (root+'data/STREAM_branch_labels.csv')

# save expression matrix
exp_mat = pd.DataFrame (adata.X.T, columns=adata.obs.index, index=adata.var_names)
exp_mat.to_csv (root+'data/STREAM_data.csv')

# detect transition markers
# this function is rather time consuming
st.select_variable_genes(adata,loess_frac=0.01,percentile=95, n_genes=1000)
st.detect_transition_markers(adata,marker_list=adata.uns['var_genes'],
        cutoff_spearman=0.4,cutoff_logfc=0.25, root=stem,n_jobs=6)
st.plot_transition_markers(adata,fig_size=(10,5))
plt.show()

# ---------------------------------------
# obtain the trajectory for plotting in R
# ---------------------------------------
traj_list = [ ['S2', 'S1'], 
              ['S3', 'S1'],
              ['S0', 'S1']]
traj_df = pu.get_traj_order (adata, traj_list)
traj_df.to_csv (root+'result/STREAM_graph.csv')

# --------------------------------------
# compute correlation with in vivo cells
# --------------------------------------
gp_root =''
pt_mean = pd.read_csv (gp_root+'result/PT_pred_mean.csv', index_col=[0])
pt_var  = pd.read_csv (gp_root+'result/PT_pred_var.csv', index_col=[0])
pt_mean = pt_mean.loc [ adata.obs.index]
pt_var  = pt_var.loc  [ adata.obs.index]

# obtain mean and var
ori_data = pd.read_csv (root+'data/merged_.csv', index_col=[0])
ori_mean = ori_data.values.mean(1, keepdims=True)
ori_std  = ori_data.values.std(1, keepdims=True)
all_X= pd.read_csv (root+'data/merged_all.csv', index_col=[0])
all_meta= pd.read_csv (root+'data/merged_meta_all.csv', index_col=[0])

scaled_X = pd.DataFrame ((all_X.values - ori_mean)/ori_std)
scaled_X.columns = all_meta.index
scaled_X.index = all_X.index
all_data = anndata.AnnData (scaled_X.T, obs=all_meta)

all_types = np.unique (all_data.obs['assigned_cluster'])
for one_type in all_types:
    adata.obs [one_type] = pu.get_likelihood (all_data, one_type, 'assigned_cluster',
            pt_mean.values, pt_var.values)

vitro_types = ['hTSC_OKAE', 'hTSC_TURCO', 'hESC_YAN', 'hESC']
in_vitro = all_data [all_data.obs ['revised'].isin (vitro_types) ]

for one_type in vitro_types:
    adata.obs [one_type] = pu.get_likelihood (in_vitro, one_type, 'revised',
            pt_mean.values, pt_var.values)

adata.obs.to_csv (root+'result/cell_likelihood.csv')
vivo_vitro = list (all_types)
vivo_vitro.extend (vitro_types)
st.plot_dimension_reduction(adata,color=list (all_types),n_components=3,
        show_graph=False,show_text=True, fig_ncol=6, plotly=False)
plt.show()
