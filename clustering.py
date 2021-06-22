"""
@author: Giacomo Spisni
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
import configparser
import uuid

# Reading configuration file
config = configparser.ConfigParser()
config.read('configuration.txt')

# Scanpy parameters
sc.settings.verbose = 3
sc.set_figure_params(dpi = 600)

# Prepare workspace
pathlib.Path('figures').mkdir(parents=True, exist_ok=True)
pathlib.Path('backup').mkdir(parents=True, exist_ok=True)

directory = Path.cwd()

# Files to be processed
data_source = pathlib.Path(config.get('paths', 'data_source'))
coordinates_source = pathlib.Path(config.get('paths', 'coordinates_source'))

# %% Data reading

# Generate an unique identifier for this run
filename = uuid.uuid4().hex

print(f'Current run ID: {filename}.')
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
clustering_algo = config.get('parameters', 'clustering_algo')
n_neighbors = config.get('parameters', 'n_neighbors')
n_jobs = config.get('parameters', 'n_jobs')
pca_n_comps = config.get('parameters', 'pca_n_comps')
umap_n_comps = config.get('parameters', 'umap_n_comps')

# Check maximum number of components
n_comps = len(data.columns)-1

if (pca_n_comps <= 0):
    if (pca_n_comps == 0):
        pca_n_comps = n_comps   # No components reduction
    else:
        pca_n_comps = n_comps + pca_n_comps     # Subtract the requested amount of components

if (umap_n_comps <= 0):
    if (umap_n_comps == 0):
        umap_n_comps = n_comps   # No components reduction
    else:
        umap_n_comps = n_comps + umap_n_comps     # Subtract the requested amount of components
    


# 1. Principal Component Analysis
adata = sc.tl.pca(adata,
                  n_comps = pca_n_comps,
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

adata.write(f'backup/{filename}_PCA_backup.h5ad')
print(f'[{datetime.now()}] PCA completed.')

# 2. UMAP
# Setup
init_seed = np.random.randint(0, 100)

# Neighbors graph
sc.pp.neighbors(adata,
                n_neighbors = n_neighbors,
                n_pcs = umap_n_comps,
                random_state = init_seed,
                method = 'umap',
                copy = False
               )
adata.write(f'backup/{filename}_NEIGHBORS_backup.h5ad')

# UMAP reduction
sc.tl.umap(adata,
           random_state = init_seed,
           n_components = umap_n_comps,
           copy = False
          )
adata.write(f'backup/{filename}_UMAP_backup.h5ad')

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

adata.write(f'backup/{filename}_PhenoGraph_backup.h5ad')

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
           save = f'_{filename}.png'
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
fig.write_html(f"figures/Clustered_{filename}.html")

print(f'[{datetime.now()}] Process completed.')