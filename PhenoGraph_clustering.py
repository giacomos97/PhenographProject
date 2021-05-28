"""
Created on May 23 2021
@author: Giacomo Spisni

https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.pca.html
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
sc.set_figure_params(dpi=600)

directory = Path.cwd()

print(f'Process started: {datetime.now()}')

# File normalizzato da processare
data_source = directory/"SecondoDataset/LNBmatrix.csv" 
coordinates_source = directory/"SecondoDataset/LNBcentroidsTOP.csv"

# %% Lettura dati
data = pd.read_csv(data_source)
coordinates = pd.read_csv(coordinates_source, header = None, names = ['cell','x','y'])

data = data.drop(['Unnamed: 0'], axis=1, errors='ignore')
coordinates = coordinates.drop(['cell'], axis=1, errors='ignore')

adata = AnnData(data)

# %%  Esecuzione di PhenoGraph
clustering_algo = "leiden"
n_neighbors = 15
n_jobs = 4
n_comps = len(data.columns)-1 # Determina la dimensione dello spazio

# Principal component analysis
adata = sc.tl.pca(adata, n_comps = n_comps, copy = True) 

# Clustering
communities, graph, Q = sce.tl.phenograph(adata, clustering_algo = clustering_algo, k = n_neighbors, n_jobs = n_jobs, copy = True)
'''
communities : integer array of community assignments for each row in data.
graph       : the graph that was used for clustering.
Q           : the modularity score for communities on graph.
'''

adata.obs['PhenoGraph_clusters'] = pd.Categorical(communities)
adata.uns['PhenoGraph_Q'] = Q
adata.uns['PhenoGraph_k'] = n_neighbors

sc.pp.neighbors(adata, n_neighbors = n_neighbors, n_pcs = n_comps)
sc.tl.umap(adata, random_state = 0) # Fissare random_state per riproducibilità risultati

sc.pl.umap(adata,
           color = ['PhenoGraph_clusters'],
           palette = sc.pl.palettes.vega_20_scanpy,
           #palette = sc.pl.palettes.godsnot_102,
           legend_fontsize = 10,
           title = f'PhenoGraph example with k={n_neighbors}'
        )

# Associa ogni coordinata a un cluster
coordinates['cell_group'] = communities

# %% Plot coordinate spaziali
fig = go.Figure()

for cell_group,grp in coordinates.groupby('cell_group'):
    fig.add_trace(go.Scattergl(x = grp['x'],
                               y = grp['y'],
                               mode = 'markers',
                               name = f'{cell_group}',
                               marker_size = 5
                            )
                 )

fig.update_layout(title = f"Spatial distribution of clustered cells (k = {n_neighbors}) - Q = {Q}",
                  showlegend = True,
                  xaxis_title = 'X',
                  yaxis_title = 'Y')

fig.show()
fig.write_html("Clustered.html")

print(f'Process completed: {datetime.now()}')