# %%
import anndata
from scipy.sparse import csr_matrix
import h5py
import numpy as np

# %% Parse snakemake object
if 'snakemake' in locals():
    peak_mat_h5 = snakemake.input['peak_mat_h5']
    peak_mat_h5ad = snakemake.output['peak_mat_h5ad']
else:
    peak_mat_h5 = '/Volumes/Transcend/CKD/out/archr/Atlas/ATAC/peak_mat.h5'
    peak_mat_h5ad = '/Volumes/Transcend/CKD/out/archr/Atlas/ATAC/peak_mat.h5ad'

# %%
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
                # NA values are somehow converted to -2147483648 in the hdf5 object, which sucks. 
                # I need to find a workaround to convert them back to NA in the adata object.
                # So I create a n+1'th level for the NA values
                new_idx = np.max(idx) + 1
                idx[idx == -2147483648] = new_idx
                # Now create a new value in obs[key + '_levels'] for the NA value
                level_values = obs[key + '_levels'][()].astype(str)
                level_values = np.append(level_values, 'NA')
                # Now do the mapping
                adata.obs[key] = level_values[idx]
            elif '_levels' not in key and 'colnames' not in key:
                # if bytes in [x for x in obs[key][()].apply(type).unique()]:]
                #     adata.obs[key] = obs[key][()].astype(str)#.str.decode("utf-8")
                # else:
                adata.obs[key] = obs[key][()]
    return adata

# %%
adata = my_read_hdf(peak_mat_h5, "RawCounts")

# %%
adata.obs['sample'] = adata.obs['index'].str.decode("utf-8").str.split('#',1).str[0].str.strip()

# %% Store the raw counts in the `.layers` attribute so that we can use them afterwards to generate pseudo-bulk profiles.
adata.layers['counts'] = adata.X

# %%
adata.write_h5ad(peak_mat_h5ad)


