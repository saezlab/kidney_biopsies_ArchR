# %%
import anndata
from scipy.sparse import csr_matrix
import h5py
import pandas as pd
import numpy as np

# %% Parse snakemake object
mat_h5 = snakemake.input['mat_h5']
embeddings_path = snakemake.input['embeddings_path']
group = snakemake.params['group'] # GeneScores
mat_h5ad = snakemake.output['mat_h5ad']

# %% Function to read h5 file and create AnnData object
def my_read_hdf(filename, key):
    with h5py.File(filename, "r") as f:
        X = f[key]
        mat = csr_matrix((X['values'], X['indices'], X['indptr']), shape = (X['dims'][0], X['dims'][1]))
        adata = anndata.AnnData(mat, dtype = int)
        adata.obs_names = X['obs_names'][()].astype(str)
        adata.var_names = X['var_names'][()].astype(str)
        # Add obs data
        obs = f['obs']
        for key in obs.keys():
            # Check if key is a factor and if so convert
            if key + '_levels' in obs.keys():
                idx = obs[key][()]
                adata.obs[key] = obs[key + '_levels'][()][idx].astype(str)#.str.decode("utf-8")
            elif '_levels' not in key and 'colnames' not in key:
                # if bytes in [x for x in obs[key][()].apply(type).unique()]:]
                #     adata.obs[key] = obs[key][()].astype(str)#.str.decode("utf-8")
                # else:
                adata.obs[key] = obs[key][()]
    return adata

# %% Read h5 file and create AnnData object
adata = my_read_hdf(filename = mat_h5, key = group)

# %% Add sample information
if 'index' in adata.obs.columns:
    adata.obs['Sample'] = adata.obs['index'].str.decode("utf-8").str.split('#',1).str[0].str.strip()

# %% Add UMAP embeddings
print("Adding UMAP embeddings to AnnData object...")
# embeddings_path = "results/coordinates_" + group + ".csv"
embeddings = pd.read_csv(embeddings_path, index_col=0)
embedding_names = np.unique(embeddings['embedding'])
for i in embedding_names:
    coordinates = embeddings[embeddings['embedding'] == i].drop(columns = 'embedding').set_index('cell')
    adata.obsm['X_' + i] = coordinates.reindex(adata.obs_names)


# %%
adata.write_h5ad(mat_h5ad)
