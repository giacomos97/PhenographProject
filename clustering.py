"""
@author: Giacomo Spisni

References:
https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.pca.html
https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.pca.html
https://scanpy.readthedocs.io/en/stable/external/scanpy.external.tl.phenograph.html
https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.umap.html
https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.umap.html
"""

# %% Imports and setup
from anndata import AnnData
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import pathlib
from pathlib import Path
import plotly.graph_objects as go
from datetime import datetime

# Scanpy parameters
sc.settings.verbose = 3
sc.set_figure_params(dpi = 600)

# Prepare workspace
pathlib.Path('figures').mkdir(parents=True, exist_ok=True)
pathlib.Path('backup').mkdir(parents=True, exist_ok=True)

directory = Path.cwd()

# Files to be processed
data_source = directory/"dataset/LNBmatrix.csv" 
coordinates_source = directory/"dataset/LNBcentroids.csv"

# %% Data reading

timestamp = datetime.timestamp(datetime.now())
print(f'Current run timestamp: {timestamp} .')
print(f"[{datetime.now()}] Process started.")

# Read files 
data = pd.read_csv(data_source)
coordinates = pd.read_csv(coordinates_source, header = None, names = ['cell','x','y'])

# Clean up data and copy into AnnData
data.drop(['Unnamed: 0'], axis = 1, errors = 'ignore', inplace = True)
coordinates.drop(['cell'], axis = 1, errors = 'ignore', inplace = True)

data.astype('float', copy = False)
coordinates.astype('float', copy = False)

adata = AnnData(data)

print(f"[{datetime.now()}] Data reading completed.")

# %%  Data analysis
# Setup
clustering_algo = "louvain" # 'louvain' or 'leiden'
n_neighbors = 15
n_jobs = 4
n_comps = len(data.columns)-1 # No components reduction

# 1. Principal Component Analysis
adata = sc.tl.pca(adata,
                  n_comps = n_comps,
                  copy = True
                 ) 

# Plot PCA results
#sc.pl.pca(adata,
#          #color = 'sample',
#          components = ['1,2,3','44,45,46'],
#          ncols = 2,
#          projection = '3d',
#          title = ['PCA results example','PCA results example'],
#          save = f'_{timestamp}.png'
#         )

adata.write(f'backup/{timestamp}_PCA_backup.h5ad')
print(f'[{datetime.now()}] PCA completed.')

# 2. UMAP
# Setup
init_seed = np.random.randint(0, 100)
reduced_n_comps = 2

# Neighbors graph
sc.pp.neighbors(adata,
                n_neighbors = n_neighbors,
                n_pcs = n_comps,
                random_state = init_seed,
                method = 'umap',
                copy = False
               )
adata.write(f'backup/{timestamp}_NEIGHBORS_backup.h5ad')

# UMAP reduction
sc.tl.umap(adata,
           random_state = init_seed,
           n_components = reduced_n_comps,
           copy = False
          )
adata.write(f'backup/{timestamp}_UMAP_backup.h5ad')

# Plot UMAP results (unclustered)
#sc.pl.umap(adata,
           #color = ['PhenoGraph_clusters'],
           #components = ['1,2','44,45'],
           #projection = '2d',
#           palette = sc.pl.palettes.vega_20_scanpy,
           #palette = sc.pl.palettes.godsnot_102,
#           legend_fontsize = 10,
           #ncols = 2,
           #title = ['UMAP results','UMAP results'],
           #save = f'_{timestamp}.png'
#        )

print(f'[{datetime.now()}] UMAP completed.')

# 3. PhenoGraph clustering
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
adata.obs['PhenoGraph_clusters'] = pd.Categorical(communities)
adata.uns['PhenoGraph_Q'] = Q
adata.uns['PhenoGraph_k'] = n_neighbors

adata.write(f'backup/{timestamp}_PhenoGraph_backup.h5ad')

# Re-plot UMAP results (clustered)
sc.pl.umap(adata,
           color = ['PhenoGraph_clusters'],
#           components = ['1,2','44,45'],
#           projection = '2d',
           palette = sc.pl.palettes.vega_20_scanpy,
#           #palette = sc.pl.palettes.godsnot_102,
           legend_fontsize = 10,
#           ncols = 2,
           title = f'PhenoGraph clusters with k={n_neighbors}',
           save = f'_{timestamp}.png'
        )

print(f'[{datetime.now()}] Clustering completed.')

# %% Spatial coordinates plot
# Associate each coordinate to a cluster
coordinates['cell_group'] = communities

fig = go.Figure()

for cell_group,grp in coordinates.groupby('cell_group'):
    fig.add_trace(go.Scattergl(x = grp['x'],
                               y = grp['y'],
                               mode = 'markers',
                               name = f'{cell_group}',
                               marker_size = 5
                            )
                 )

fig.update_layout(title = f'Spatial distribution of clustered cells (k = {n_neighbors}) - Q = {Q} - Seed = {init_seed}',
                  showlegend = True,
                  xaxis_title = 'X',
                  yaxis_title = 'Y',
                  #plot_bgcolor='rgba(0,0,0,0)',
                  autosize = False,
                  width = 800,
                  height = 800
                 )

fig.show()
fig.write_html(f"figures/Clustered_{timestamp}.html")

print(f'[{datetime.now()}] Process completed.')