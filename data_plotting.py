"""
@author: Giacomo Spisni
"""

# IMPORTS
import scanpy as sc
import pandas as pd
import plotly.graph_objects as go

# FUNCTIONS
def PCA(adata, filename, temp_dir, figures_dir, config):
    '''
    Plot Principal Component Analysis results from the dataset identified by filename in temp_dir.
    Use parameters in scanpy_PCA_plot section of the configuration file.

    Parameters
    ----------
    filename : str
        Unique identifier of this dataset stored in temp_dir.
    temp_dir : path object of pathlib module.
        Directory where the temporary copy of the datasets is stored.
    figures_dir : TYPE
        Directory where to store plots and figures.
    config : ConfigParser object of configparser module.
        Configuration file.

    Returns
    -------
    filename : str
        Unique identifier of this dataset stored in temp_dir.

    '''
    
    # Plot PCA results
    sc.pl.pca(adata,
              #**dict(config.items('scanpy_PCA_plot')),
              components = '1,2',
              projection = '2d',
              title = 'PCA results example',
              save = f'_{filename}.png'
            )
        
    return True

def UMAP(adata, filename, temp_dir, figures_dir, config):
    '''
    Plot UMAP representations from the dataset identified by filename in temp_dir.
    Use parameters in scanpy_UMAP_plot section of the configuration file.

    Parameters
    ----------
    filename : str
        Unique identifier of this dataset stored in temp_dir.
    temp_dir : path object of pathlib module.
        Directory where the temporary copy of the datasets is stored.
    figures_dir : TYPE
        Directory where to store plots and figures.
    config : ConfigParser object of configparser module.
        Configuration file.

    Returns
    -------
    filename : str
        Unique identifier of this dataset stored in temp_dir.

    '''
    
    # Plot UMAP results
    sc.pl.umap(adata,
               #**dict(config.items('scanpy_UMAP_plot')),
               components = '1,2',
               projection = '2d',
               palette = sc.pl.palettes.vega_20_scanpy,
               legend_fontsize = 10,
               color = ['PhenoGraph_clusters_Categorical'],
               save = f'_{filename}.png'
            )
        
    return True

def CLUSTERS(adata, coordinates, filename, temp_dir, figures_dir, config):
    '''
    Create an interactive plot of the spatial distribution of clustered cells
    from the dataset identified by filename in temp_dir.
    Use parameters in scanpy_CLUSTERS_plot section of the configuration file.

    Parameters
    ----------
    filename : str
        Unique identifier of this dataset stored in temp_dir.
    temp_dir : path object of pathlib module.
        Directory where the temporary copy of the datasets is stored.
    figures_dir : TYPE
        Directory where to store plots and figures.
    config : ConfigParser object of configparser module.
        Configuration file.

    Returns
    -------
    filename : str
        Unique identifier of this dataset stored in temp_dir.

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
    
    #fig.update_layout(**dict(config.items('scanpy_CLUSTERS_plot')))
    
    fig.update_layout(title = f'Spatial distribution of clustered cells (k = {n_neighbors}) - Q = {Q}',
                      showlegend = True,
                      xaxis_title = 'X',
                      yaxis_title = 'Y',
                      autosize = False,
                      width = 800,
                      height = 800
                      )
    
    fig.show()
    fig.write_html(f"{figures_dir}/{filename}_clusters.html")
    
    return True
