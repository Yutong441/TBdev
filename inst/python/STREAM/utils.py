import numpy as np
import pandas as pd
import anndata
import stream as st

def from_seurat_to_anndata (data_dir, result_dir, exp_mat_file='merged_.csv',
        meta_file='merged_meta_.csv'):
    Y = pd.read_csv(data_dir+exp_mat_file, index_col=[0])
    metadata = pd.read_csv(data_dir+meta_file, index_col=[0])

    # exclude a subcluster of CTB that complicate developmental trajectory
    select_index1 = metadata['assigned_cluster'] != 'uCTB'
    select_index2 = metadata['revised'].isin (['PE', 'EPI', 'PSA-EPI'])

    select_index = (select_index1) & (~select_index2) 
    exp_mat = Y.T.values [select_index]
    metadata = metadata [select_index]
    adata = anndata.AnnData (X=exp_mat, obs=metadata)
    adata.var_names = Y.index
    #adata.X = adata.X - adata.X.mean (0, keepdims=True)

    st.set_workdir(adata, result_dir)
    st.select_variable_genes(adata,loess_frac=0.01,percentile=95)
    st.select_top_principal_components(adata,feature=None,first_pc=True,n_pc=40)
    return adata

def Epilgraph_on_GPLVM (adata, GPLVM_dir):
    # run spectral embedding
    st.dimension_reduction(adata,method='se',feature='top_pcs',n_components=3,n_neighbors=15,n_jobs=4)
    #st.plot_dimension_reduction(adata,color=['date','revised', 'paper'],
    #        n_components=3,show_graph=False,show_text=False)
    #plt.show()

    # load GPLVM results
    DR_pt = pd.read_csv (GPLVM_dir, index_col=[0])
    adata.obsm['X_dr'] = DR_pt.loc[adata.obs.index,:].values[:, :3]
    adata.obsm['X_se'] = DR_pt.loc[adata.obs.index,:].values[:, :3]
    #st.plot_dimension_reduction(adata,color=['date','revised'],
    #        n_components=3,show_graph=False,show_text=False)
    #plt.show()

    # elastic principal graph
    np.random.seed (100)
    st.seed_elastic_principal_graph(adata,n_clusters=20)
    #st.plot_dimension_reduction(adata,color=['date','revised'],n_components=3,show_graph=True,show_text=True)
    #plt.show ()

    st.elastic_principal_graph(adata,epg_alpha=0.02,epg_mu=0.05,epg_lambda=0.01)
    #st.optimize_branching(adata,incr_n_nodes=30)
    #st.plot_dimension_reduction(adata,color=['date','revised'],n_components=3,show_graph=True,show_text=True)
    #st.plot_branches(adata,show_text=False)
    #plt.show ()
    return adata

# marker gene detection
#st.detect_leaf_markers(adata,marker_list=adata.var_names,cutoff_zscore=1.0,
#        cutoff_pvalue=0.01, root='S5',n_jobs=4)
#
#adata.uns['leaf_markers_all'].head()
#adata.uns['leaf_markers'][('S1','S2')].head()

def mean_log_likelihood (mu, var, x):
    '''This function calculates the mean log likelihood based on IID normal
    distribution:
    1/k \Sum_k log P_{i,j,k} = 1/k \Sum_k {-0.5*\frac{ ( mu_{i,j} - x_{k,j} )^2
    }{ var_{i,j}} - \sqrt { var_{i, j} } - 0.5*\log {2\pi}}

    Because broadcasting during subtraction of `mu` with `x` creates a 3D
    array, which is not memory efficient for large data. This function uses:
    \Sum_k (mu_{i,j} - x_{k,j} )^2 = mu_{i,j}^2 - 2 mu_{i,j} \Sum_k x_{k,j} +
    \Sum_k x_{k,j}^2

    Args:
        `mu`: mean of the normal distribution
        `var`: variance of the normal distribution
        `x`: sample to be compared with
    All 3 arguments should be 2D array matching in axis 1
    '''
    out = mu**2 + (x**2).mean (axis=0, keepdims=True) -2*mu*x.mean (axis=0,
            keepdims=True)
    out /= -var*2 
    out = out - np.sqrt (var) - 0.5*np.log (2*np.pi)
    return out.mean (axis=1)

def get_likelihood (x, cell_type, feature, mu, var):
    return np.exp (mean_log_likelihood (mu, var, x.X [x.obs [feature] == cell_type]))

def get_trajectory (start, end, branch_point, adj_mat_key, adj_mat_val):
    index_list = [ int (branch_point [branch_point [:,1] == start, 0]) ]
    end_point = float (branch_point [branch_point [:,1] == end, 0])

    while True:
        ind = int (np.where (np.array (adj_mat_key) == index_list[-1]) [0] )
        next_point = adj_mat_val [ ind ]
        non_include = [ i in index_list for i in next_point ]
        next_point = list (np.array (next_point) [np.array (non_include) == False]) 
        if len (next_point) > 0: next_point = next_point [0]
        else: break
        index_list.append (next_point )
        #print ('The current point is {}'.format (index_list [-1] ) )
        if next_point == end_point: break
    return index_list

def get_multi_trajectory (adata, traj_list):
    branch_point = [[key, val ['label'] ]  for key, val in adata.uns[
        'flat_tree']._node.items() ]
    # 0th column = branch index, 1s column = point name
    branch_point = np.array (branch_point)

    # obtain the adjacency matrix
    adata.uns['epg']._adj
    adj_mat_val = [ list (val.keys()) for key, val in adata.uns['epg']._adj.items() ]
    adj_mat_key = list (adata.uns ['epg']._adj.keys())

    all_dict = {}
    for traj in traj_list:
        key = ','.join (traj)
        all_dict [key] = get_trajectory (traj[0], traj[1], branch_point,
                adj_mat_key, adj_mat_val)
    return all_dict

def get_traj_order (adata, traj_list):
    graph_df = np.stack ([val ['pos'] for val in adata.uns ['epg']._node.values()], axis=0)
    graph_df = pd.DataFrame (graph_df, index=adata.uns['epg']._node.keys())
    graph_df.columns = ['x', 'y', 'z']
    all_dict = get_multi_trajectory (adata, traj_list)

    df_list = []
    blank_row = pd.DataFrame ([[np.nan, np.nan, np.nan, 'blank']])
    blank_row.columns = ['x', 'y', 'z', 'branch']
    for key, val in all_dict.items ():
        one_df = graph_df.loc [val, :]
        one_df ['branch'] = key
        df_list.append (one_df)
        df_list.append (blank_row)

    all_df = pd.concat (df_list, axis=0)
    return all_df

