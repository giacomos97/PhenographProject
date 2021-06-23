"""
@author: Giacomo Spisni
"""

# IMPORTS
import anndata
import scanpy as sc
import pandas as pd
import plotly.graph_objects as go

# FUNCTIONS
def PCA(filename, temp_dir):
    '''
    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    temp_dir : TYPE
        DESCRIPTION.

    Returns
    -------
    filename : TYPE
        DESCRIPTION.

    '''
    
    # Retrieve dataset
    adata = anndata.read_h5ad(f'{temp_dir}/{filename}.h5ad')
    
    # Plot PCA results
    sc.pl.pca(adata,
              #color = 'sample',
              components = ['1,2,3','44,45,46'],
              ncols = 2,
              projection = '3d',
              title = ['PCA results example','PCA results example'],
              save = f'_{filename}.png'
            )
        
    return filename

def UMAP(filename, temp_dir):
    '''
    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    temp_dir : TYPE
        DESCRIPTION.

    Returns
    -------
    filename : TYPE
        DESCRIPTION.

    '''
    
    # Retrieve dataset
    adata = anndata.read_h5ad(f'{temp_dir}/{filename}.h5ad')
    
    # Plot PCA results
    sc.pl.umap(adata,
               color = ['PhenoGraph_clusters'],
    #           components = ['1,2','44,45'],
    #           projection = '2d',
               palette = sc.pl.palettes.vega_20_scanpy,
    #           #palette = sc.pl.palettes.godsnot_102,
               legend_fontsize = 10,
    #           ncols = 2,
               save = f'_{filename}.png'
            )
        
    return filename


def CLUSTERS(filename, temp_dir):
    '''
    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    temp_dir : TYPE
        DESCRIPTION.

    Returns
    -------
    filename : TYPE
        DESCRIPTION.

    '''
    
    # Retrieve dataset
    adata = anndata.read_h5ad(f'{temp_dir}/{filename}.h5ad')
    Q = adata.uns['PhenoGraph_Q']
    n_neighbors= adata.uns['PhenoGraph_k']
    
    communities = adata.obs['PhenoGraph_clusters']
    communities = communities.to_frame()
    communities.reset_index(inplace = True, drop = True)
    
    coordinates = pd.read_csv(f'{temp_dir}/{filename}_coordinates.csv', header = "infer")
    
    # Associate each coordinate to a cluster
    coordinates = pd.concat([coordinates, communities], axis=1)
    print (coordinates)
    
    fig = go.Figure()
    
    # Plot each cluster with different colors
    for cell_cluster,grp in coordinates.groupby('PhenoGraph_clusters'):
        fig.add_trace(go.Scattergl(x = grp['x'],
                                   y = grp['y'],
                                   mode = 'markers',
                                   name = f'{cell_cluster}',
                                   marker_size = 5
                                )
                     )
    
    fig.update_layout(title = f'Spatial distribution of clustered cells (k = {n_neighbors}) - Q = {Q}',
                      showlegend = True,
                      xaxis_title = 'X',
                      yaxis_title = 'Y',
                      #plot_bgcolor='rgba(0,0,0,0)',
                      autosize = False,
                      width = 800,
                      height = 800
                     )
    
    fig.show()
    fig.write_html(f"figures/{filename}_clusters.html")
    
    return filename
