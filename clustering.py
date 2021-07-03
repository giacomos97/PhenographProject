"""
@author: Giacomo Spisni
"""

# IMPORTS
import scanpy as sc
import pathlib
import configparser

import data_IO
import data_analysis
import data_plotting

# SETUP
# Reading configuration file
config = configparser.ConfigParser()
config.read('configuration.txt')

# Configure directories
data_source = pathlib.Path(config.get('paths', 'data_source'))
coordinates_source = pathlib.Path(config.get('paths', 'coordinates_source'))
figures_dir = pathlib.Path(config.get('paths', 'figures_dir'))
temp_dir = pathlib.Path(config.get('paths', 'temp_dir'))

# Configure parameters
PCA_params = dict(config.items('scanpy_PCA'))
NEIGHBORS_params = dict(config.items('scanpy_NEIGHBORS'))
UMAP_params = dict(config.items('scanpy_UMAP'))
PHENOGRAPH_params = dict(config.items('scanpy_PHENOGRAPH'))

PCA_plot_params = dict(config.items('scanpy_PCA_plot'))
UMAP_plot_params = dict(config.items('scanpy_UMAP_plot'))
CLUSTERS_plot_params = dict(config.items('scanpy_CLUSTERS_plot'))

# Scanpy parameters
sc.settings.verbose = 3     # Min = 0, Max = 3
sc.set_figure_params(config.items('scanpy_figures'))
sc.settings.figdir = figures_dir

# Prepare workspace
pathlib.Path(figures_dir).mkdir(parents = True, exist_ok = True)
pathlib.Path(temp_dir).mkdir(parents = True, exist_ok = True)

# DATA READING
filename = data_IO.READ_SOURCE(data_source, coordinates_source, temp_dir)

# DATA ANALYSIS
adata = data_IO.READ_ADATA(filename, temp_dir)

adata = data_analysis.PCA(adata, temp_dir, PCA_params)
adata = data_analysis.NEIGHBORS(adata, temp_dir, NEIGHBORS_params)
adata = data_analysis.UMAP(adata, temp_dir, UMAP_params)
adata = data_analysis.PHENOGRAPH(adata, temp_dir, PHENOGRAPH_params)

filename = data_IO.WRITE_ADATA(adata, filename, temp_dir)

# DATA PLOTTING 
adata = data_IO.READ_ADATA(filename, temp_dir)
coordinates = data_IO.READ_COORDINATES(filename, temp_dir)

#data_plotting.PCA(adata, filename, temp_dir, figures_dir, PCA_plot_params)
data_plotting.UMAP(adata, filename, temp_dir, figures_dir, UMAP_plot_params)
data_plotting.CLUSTERS(adata, coordinates, filename, temp_dir, figures_dir, CLUSTERS_plot_params)
