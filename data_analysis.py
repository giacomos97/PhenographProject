"""
@author: Giacomo Spisni
"""

# IMPORTS
import scanpy as sc
import scanpy.external as sce
from datetime import datetime
import pandas as pd

# FUNCTIONS
def PCA(adata, temp_dir, params):
    '''
    Perform Principal Component Analysis on the dataset identified by filename in temp_dir.
    Use parameters in scanpy_PCA section of the configuration file.
    Overwrite the dataset in temp_dir with the achieved results.
    
    Parameters
    ----------
    filename : str
        Unique identifier of this dataset stored in temp_dir.
    temp_dir : path object of pathlib module.
        Directory where the temporary copy of the datasets is stored.
    config : ConfigParser object of configparser module.
        Configuration file.

    Returns
    -------
    filename : str
        Unique identifier of this dataset stored in temp_dir.

    '''
    
    print(f"[{datetime.now()}] PCA started.")
    
    # Setup
    n_comps = int(params['n_comps'])

    # Check number of components
    max_n_comps = len(adata.var_names)-1

    if (n_comps <= 0):
        if (n_comps == 0):
            n_comps = max_n_comps   # No components reduction
        else:
            n_comps = max_n_comps + n_comps     # Subtract the requested amount of components
    
    # Perform PCA
    sc.tl.pca(adata,
              n_comps = n_comps,
              copy = False
              ) 
    
    print(f"[{datetime.now()}] PCA completed.")
    
    return adata

def NEIGHBORS(adata, temp_dir, params):
    '''
    Compute a neighborhood graph on the dataset identified by filename in temp_dir.
    Use parameters in scanpy_NEIGHBORS section of the configuration file.
    Overwrite the dataset in temp_dir with the achieved results.
    
    Parameters
    ----------
    filename : str
        Unique identifier of this dataset stored in temp_dir.
    temp_dir : path object of pathlib module.
        Directory where the temporary copy of the datasets is stored.
    config : ConfigParser object of configparser module.
        Configuration file.

    Returns
    -------
    filename : str
        Unique identifier of this dataset stored in temp_dir.

    '''
    
    print(f"[{datetime.now()}] NEIGHBORS started.")
    
    # Setup
    n_comps = int(params['n_comps'])
    n_neighbors = int(params['n_neighbors'])
    init_seed = int(params['init_seed'])

    # Check number of components
    max_n_comps = len(adata.var_names)-1

    if (n_comps <= 0):
        if (n_comps == 0):
            n_comps = max_n_comps   # No components reduction
        else:
            n_comps = max_n_comps + n_comps     # Subtract the requested amount of components
   
    # Neighbors graph
    sc.pp.neighbors(adata,
                    n_neighbors = n_neighbors,
                    n_pcs = n_comps,
                    random_state = init_seed,
                    method = 'umap',
                    copy = False
                   )

    print(f"[{datetime.now()}] NEIGHBORS completed.")
    
    return adata
    
def UMAP(adata, temp_dir, params):
    '''
    Perform dimensional reduction with UMAP on the dataset identified by filename in temp_dir.
    Use parameters in scanpy_UMAP section of the configuration file.
    Overwrite the dataset in temp_dir with the achieved results.
    
    Parameters
    ----------
    filename : str
        Unique identifier of this dataset stored in temp_dir.
    temp_dir : path object of pathlib module.
        Directory where the temporary copy of the datasets is stored.
    config : ConfigParser object of configparser module.
        Configuration file.

    Returns
    -------
    filename : str
        Unique identifier of this dataset stored in temp_dir.

    '''
    
    print(f"[{datetime.now()}] UMAP started.")
    
    # Setup
    n_comps = int(params['n_comps'])
    init_seed = int(params['init_seed'])

    # Check number of components
    max_n_comps = len(adata.var_names)-1

    if (n_comps <= 0):
        if (n_comps == 0):
            n_comps = max_n_comps   # No components reduction
        else:
            n_comps = max_n_comps + n_comps     # Subtract the requested amount of components
   
    # Perform UMAP
    sc.tl.umap(adata,
               random_state = init_seed,
               n_components = n_comps,
               copy = False
              )

    print(f"[{datetime.now()}] UMAP completed.")
    
    return adata

def PHENOGRAPH(adata, temp_dir, params):
    '''
    Perform clustering with PhenoGraph on the dataset identified by filename in temp_dir.
    Use parameters in scanpy_PHENOGRAPH section of the configuration file.
    Overwrite the dataset in temp_dir with the achieved results.
    
    Parameters
    ----------
    filename : str
        Unique identifier of this dataset stored in temp_dir.
    temp_dir : path object of pathlib module.
        Directory where the temporary copy of the datasets is stored.
    config : ConfigParser object of configparser module.
        Configuration file.

    Returns
    -------
    filename : str
        Unique identifier of this dataset stored in temp_dir.

    '''
    
    print(f"[{datetime.now()}] PHENOGRAPH started.")
    
    # Setup
    clustering_algo = params['clustering_algo']
    n_neighbors = int(params['n_neighbors'])
    n_jobs = int(params['n_jobs'])
   
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
    adata.obs['PhenoGraph_clusters_Categorical'] = pd.Categorical(communities)
    adata.obs['PhenoGraph_clusters_numeric'] = communities
    adata.uns['PhenoGraph_Q'] = Q
    adata.uns['PhenoGraph_k'] = n_neighbors

    print(f"[{datetime.now()}] PHENOGRAPH completed.")
    
    return adata
