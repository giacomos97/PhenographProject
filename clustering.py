"""
Created on May 23 2021
@author: Giacomo Spisni

https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.pca.html
https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.pca.html

https://scanpy.readthedocs.io/en/stable/external/scanpy.external.tl.phenograph.html

https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.umap.html
https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.umap.html
"""

# %% Imports e setup
from anndata import AnnData
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import plotly.graph_objects as go
from datetime import datetime

sc.settings.verbose = 3
sc.set_figure_params(dpi = 600)

directory = Path.cwd()

print(f'Process started: {datetime.now()}')

# File normalizzato da processare
#data_source = directory/"SecondoDataset/LNBmatrix.csv" 
#coordinates_source = directory/"SecondoDataset/LNBcentroidsTOP.csv"
data_source = directory/"PrimoDataset/data_full_normalized_truncated.csv" 
coordinates_source = directory/"PrimoDataset/prova_centroidi_truncated.csv"

# %% Lettura dati
data = pd.read_csv(data_source)
coordinates = pd.read_csv(coordinates_source, header = None, names = ['cell','x','y'])

data = data.drop(['Unnamed: 0'], axis = 1, errors = 'ignore')
coordinates = coordinates.drop(['cell'], axis = 1, errors = 'ignore')

data.astype('float', copy = False)
coordinates.astype('float', copy = False)

adata = AnnData(data)

# %%  Analisi dati
clustering_algo = "louvain" # 'louvain' or 'leiden'
n_neighbors = 15
n_jobs = 4
n_comps = len(data.columns)-1 # No components reduction

# 1. Principal component analysis
'''

'''
adata = sc.tl.pca(adata,
                  n_comps = n_comps,
                  copy = True
                 ) 

#sc.pl.pca(adata,
#          #color = 'sample',
#          components = ['1,2,3','44,45,46'],
#          ncols = 2,
#          projection = '3d',
#          title = ['PCA results example','PCA results example']
#         )

# 2. UMAP
'''
The UMAP implementation in SCANPY requires a neighborhood graph as the distance matrix.
'''
random_state = 0
reduced_n_comps = 1 # Reduce to 1D

sc.pp.neighbors(adata,
                n_neighbors = n_neighbors,
                n_pcs = n_comps,
                random_state = random_state,
                method = 'umap',
                copy = False
               )

sc.tl.umap(adata,
           random_state = random_state,
           n_components = reduced_n_comps,
           copy = False
          )

#sc.pl.umap(adata,
           #color = ['PhenoGraph_clusters'],
           #components = ['1,2','44,45'],
           #projection = '2d',
#           palette = sc.pl.palettes.vega_20_scanpy,
           #palette = sc.pl.palettes.godsnot_102,
#           legend_fontsize = 10,
           #ncols = 2,
           #title = ['UMAP results','UMAP results']
#        )

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

#sc.pl.umap(adata,
#           color = ['PhenoGraph_clusters'],
#           components = ['1,2','44,45'],
#           projection = '2d',
#           palette = sc.pl.palettes.vega_20_scanpy,
#           #palette = sc.pl.palettes.godsnot_102,
#           legend_fontsize = 10,
#           ncols = 2,
#           title = f'PhenoGraph clusters with k={n_neighbors}'
#        )

# %% Plot coordinate spaziali
# Associa ogni coordinata a un cluster
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
fig.write_html("Clustered.html")

print(f'Process completed: {datetime.now()}')