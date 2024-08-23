#### Tool Functions for snATAC-seq with snapATAC2
import pandas as pd
import scanpy as sc
import snapatac2 as snap

#####################################
## Downsampling function - 
## adapted from Chuck's function for snRNA-seq
# standardize library sizes by removing nucs with low libs and downsampling the rest
def standardize_libs(adata, tar_lib_sz=2500, seed=123):
    """
    Standardizes library sizes in the Anndata object by removing low-count cells and downsampling.

    Parameters:
    - adata: Anndata object containing the observations data.
    - tar_lib_sz: Target library size for standardizing cell counts. Default is 2500.
    - seed: Random seed for reproducibility of the downsampling. Default is 123.

    This function performs the following steps:
    1. Identifies and prints the number of cells (nuclei) with total counts less than or equal to the target library size,
       grouped by 'library_id'.
    2. Removes cells with total counts less than the target library size.
    3. Downsamples the total counts of each remaining cell to the target library size.
    4. Resets the data type of the counts matrix to float to address a specific bug in the scanpy version used.

    Modifies:
    - The function modifies the adata object in-place by potentially removing low-count cells and adjusting the counts 
      of the remaining cells to standardize their library sizes.
    """
    # make sure counts for each nuc are current
    #sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)

    # print number of nucs to be removed
    print("Number of low lib nucs removed per run:")
    low_mask = (adata.obs['n_fragment'] <= tar_lib_sz)
    print(pd.value_counts(adata.obs['library_id'].values[low_mask]).to_string())
    
    # remove nucs with low library sizes
    #### comment out line below to keep low UMI count nuclei ####
    sc.pp.filter_cells(adata, min_counts=tar_lib_sz, copy=False)
    
    # downsample nuc total counts over target library size
    sc.pp.downsample_counts(adata, counts_per_cell=tar_lib_sz, copy=False, replace=False, random_state=seed)
    
    # due to scanpy bug in this version have to reset dtype
    adata.X = adata.X.astype(float)
    return

#####################################
## Function for spectral embedding to umap
def spec_to_umap(adata, n_comps=30, distance_metric='cosine', sample_size=1.0, n_neighbors=50, seed=123, 
                 weighted_by_sd=True, features='selected'):
    """
    Performs spectral embedding followed by UMAP for dimensionality reduction on the given Anndata object.

    Parameters:
    - adata: Anndata object containing the observations data.
    - n_comps: Number of principal components to use for spectral embedding. Default is 30.
    - distance_metric: Distance metric to use for spectral embedding. Default is 'cosine'.
    - sample_size: Proportion of the data to sample if using a distance metric other than 'cosine'. 
                   A value between 0.0 and 1.0. Default is 1.0.
    - n_neighbors: Number of neighbors to use for the k-nearest neighbors graph construction in UMAP. Default is 50.
    - seed: Random seed for reproducibility of the spectral embedding and UMAP. Default is 123.
    - weighted_by_sd: If True, weights by the standard deviation when performing spectral embedding. Default is True.

    Modifies:
    - The function modifies the adata object in-place, adding the UMAP coordinates and cluster assignments.
    """
    if distance_metric=='cosine':
        snap.tl.spectral(adata, random_state=seed, n_comps=n_comps, 
                         weighted_by_sd=weighted_by_sd, features=features)
        snap.pp.knn(adata, n_neighbors=n_neighbors, random_state=seed)
        snap.tl.leiden(adata)
        snap.tl.umap(adata, random_state=seed)
    else:
        snap.tl.spectral(adata, distance_metric=distance_metric,
                         sample_size=sample_size, random_state=seed,
                         n_comps=n_comps, weighted_by_sd=weighted_by_sd,
                         features=features)
        snap.pp.knn(adata, n_neighbors=n_neighbors, random_state=seed)
        snap.tl.leiden(adata)
        snap.tl.umap(adata, random_state=seed)
    return