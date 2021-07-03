"""
@author: Giacomo Spisni
"""

# IMPORTS
import scanpy as sc
import pandas as pd
import plotly.graph_objects as go

# FUNCTIONS
def PCA(adata, filename, figures_dir, params):
    '''
    Plot Principal Component Analysis results from the input dataset, 
    considering parameters provided.

    Parameters
    ----------
    adata : AnnData object
        Database containing cleaned copy of work data.
    filename : str
        Unique identifier of this dataset stored in temp_dir.
    figures_dir : path object of pathlib module
        Directory where to store plots and figures.
    params : dict
        Parameters to be provided to analysis function.

    Returns
    -------
    bool
        Successful function execution.

    '''
    
    # Plot PCA results
    sc.pl.pca(adata,
              **params,
              save = f'_{filename}.png'
            )
        
    return True

def UMAP(adata, filename, figures_dir, params):
    '''
    Plot UMAP representations from the input dataset, considering parameters
    provided.

    Parameters
    ----------
    adata : AnnData object
        Database containing cleaned copy of work data.
    filename : str
        Unique identifier of this dataset stored in temp_dir.
    figures_dir : path object of pathlib module
        Directory where to store plots and figures.
    params : dict
        Parameters to be provided to analysis function.

    Returns
    -------
    bool
        Successful function execution.

    '''

    # Plot UMAP results
    sc.pl.umap(adata,
               **params,
               palette = sc.pl.palettes.vega_20_scanpy,
               color = 'PhenoGraph_clusters_Categorical',
               save = f'_{filename}.png'
            )
        
    return True

def CLUSTERS(adata, coordinates, filename, figures_dir, params):
    '''
    Create an interactive plot of the spatial distribution of clustered cells
    from the input dataset, considering parameters provided.

    Parameters
    ----------
    adata : AnnData object
        Database containing cleaned copy of work data.
    coordinates : Pandas DataFrame object
        DataFrame containing cells' coordinates.
    filename : str
        Unique identifier of this dataset stored in temp_dir.
    figures_dir : path object of pathlib module
        Directory where to store plots and figures.
    params : dict
        Parameters to be provided to analysis function.

    Returns
    -------
    bool
        Successful function execution.

    '''
    
    # Retrieve information from adata
    Q = adata.uns['PhenoGraph_Q']
    n_neighbors= adata.uns['PhenoGraph_k']
    
    communities = adata.obs['PhenoGraph_clusters_numeric']
    communities = communities.to_frame()
    communities.reset_index(inplace = True, drop = True)
    
    # Associate each coordinate to a cluster
    coordinates = pd.concat([coordinates, communities], axis=1)
    
    fig = go.Figure()
    
    # Plot each cluster with different colors
    for cell_cluster,grp in coordinates.groupby('PhenoGraph_clusters_numeric'):
        fig.add_trace(go.Scattergl(x = grp['x'],
                                   y = grp['y'],
                                   mode = 'markers',
                                   name = f'{cell_cluster}',
                                   marker_size = 5
                                )
                     )
    
    fig.update_layout(title = f'Spatial distribution of clustered cells - run {filename}',
                      **params,
                      autosize = False,
                      showlegend = True,
                      width = 800,
                      height = 800
                      )
    
    fig.show()
    fig.write_html(f"{figures_dir}/{filename}_clusters.html")
    
    return True
