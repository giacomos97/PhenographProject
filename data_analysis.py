"""
@author: Giacomo Spisni
"""

# IMPORTS
import anndata
import scanpy as sc
import scanpy.external as sce
import numpy as np
from datetime import datetime

# FUNCTIONS
def PCA(filename, temp_dir, config):
    '''
    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    temp_dir : TYPE
        DESCRIPTION.
    config : TYPE
        DESCRIPTION.

    Returns
    -------
    filename : TYPE
        DESCRIPTION.

    '''
    
    print(f"[{datetime.now()}] PCA started.")
    # Retrieve dataset
    adata = anndata.read_h5ad(f'{temp_dir}/{filename}.h5ad')
    
    # Setup
    pca_n_comps = int(config.get('parameters', 'pca_n_comps'))

    # Check number of components
    n_comps = len(adata.var_names)-1

    if (pca_n_comps <= 0):
        if (pca_n_comps == 0):
            pca_n_comps = n_comps   # No components reduction
        else:
            pca_n_comps = n_comps + pca_n_comps     # Subtract the requested amount of components
    
    # Perform PCA
    sc.tl.pca(adata,
              n_comps = pca_n_comps,
              copy = False
              ) 
    
    print(f"[{datetime.now()}] PCA completed.")
    
    adata.write(f'{temp_dir}/{filename}.h5ad')
    
    return filename


def NEIGHBORS(filename, temp_dir, config):
    '''
    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    temp_dir : TYPE
        DESCRIPTION.
    config : TYPE
        DESCRIPTION.

    Returns
    -------
    filename : TYPE
        DESCRIPTION.

    '''
    
    print(f"[{datetime.now()}] NEIGHBORS started.")
    # Retrieve dataset
    adata = anndata.read_h5ad(f'{temp_dir}/{filename}.h5ad')
    
    # Setup
    umap_n_comps = int(config.get('parameters', 'umap_n_comps'))
    n_neighbors = int(config.get('parameters', 'n_neighbors'))
    init_seed = np.random.randint(0, 100)

    # Check number of components
    n_comps = len(adata.var_names)-1

    if (umap_n_comps <= 0):
        if (umap_n_comps == 0):
            umap_n_comps = n_comps   # No components reduction
        else:
            umap_n_comps = n_comps + umap_n_comps     # Subtract the requested amount of components
   
    # Neighbors graph
    sc.pp.neighbors(adata,
                    n_neighbors = n_neighbors,
                    n_pcs = umap_n_comps,
                    random_state = init_seed,
                    method = 'umap',
                    copy = False
                   )

    print(f"[{datetime.now()}] NEIGHBORS completed.")
    
    adata.write(f'{temp_dir}/{filename}.h5ad')
    
    return filename
    

def UMAP(filename, temp_dir, config):
    '''
    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    temp_dir : TYPE
        DESCRIPTION.
    config : TYPE
        DESCRIPTION.

    Returns
    -------
    filename : TYPE
        DESCRIPTION.

    '''
    
    print(f"[{datetime.now()}] UMAP started.")
    # Retrieve dataset
    adata = anndata.read_h5ad(f'{temp_dir}/{filename}.h5ad')
    
    # Setup
    umap_n_comps = int(config.get('parameters', 'umap_n_comps'))
    init_seed = np.random.randint(0, 100)

    # Check number of components
    n_comps = len(adata.var_names)-1

    if (umap_n_comps <= 0):
        if (umap_n_comps == 0):
            umap_n_comps = n_comps   # No components reduction
        else:
            umap_n_comps = n_comps + umap_n_comps     # Subtract the requested amount of components
   
    # Perform UMAP
    sc.tl.umap(adata,
               random_state = init_seed,
               n_components = umap_n_comps,
               copy = False
              )

    print(f"[{datetime.now()}] UMAP completed.")
    
    adata.write(f'{temp_dir}/{filename}.h5ad')
    
    return filename


def PHENOGRAPH(filename, temp_dir, config):
    '''
    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    temp_dir : TYPE
        DESCRIPTION.
    config : TYPE
        DESCRIPTION.

    Returns
    -------
    filename : TYPE
        DESCRIPTION.

    '''
    
    print(f"[{datetime.now()}] PHENOGRAPH started.")
    # Retrieve dataset
    adata = anndata.read_h5ad(f'{temp_dir}/{filename}.h5ad')
    
    # Setup
    clustering_algo = config.get('parameters', 'clustering_algo')
    n_neighbors = int(config.get('parameters', 'n_neighbors'))
    n_jobs = int(config.get('parameters', 'n_jobs'))
   
    # Perform PHENOGRAPH
    '''
    communities : integer array of community assignments for each row in data.
    graph       : the graph that was used for clustering.
    Q           : the modularity score for communities on graph.
    '''
    
    communities, graph, Q = sce.tl.phenograph(adata,
                                              clustering_algo = clustering_algo,
                                              k = n_neighbors,
                                              n_jobs = n_jobs,
                                              copy = True
                                             )
    # Also store results in adata
    adata.obs['PhenoGraph_clusters'] = communities
    adata.uns['PhenoGraph_Q'] = Q
    adata.uns['PhenoGraph_k'] = n_neighbors

    print(f"[{datetime.now()}] PHENOGRAPH completed.")
    
    adata.write(f'{temp_dir}/{filename}.h5ad')
    
    return filename