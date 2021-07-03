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

# Select directories
data_source = pathlib.Path(config.get('paths', 'data_source'))
coordinates_source = pathlib.Path(config.get('paths', 'coordinates_source'))
figures_dir = pathlib.Path(config.get('paths', 'figures_dir'))
temp_dir = pathlib.Path(config.get('paths', 'temp_dir'))

# Scanpy parameters
sc.settings.verbose = 3     # Min = 0, Max = 3
sc.set_figure_params(config.items('scanpy_figures'))
sc.settings.figdir = figures_dir

# Prepare workspace
pathlib.Path(figures_dir).mkdir(parents = True, exist_ok = True)
pathlib.Path(temp_dir).mkdir(parents = True, exist_ok = True)

# DATA READING
filename = data_IO.READ_SOURCE(data_source, coordinates_source, temp_dir, config)

# DATA ANALYSIS
adata = data_IO.READ_ADATA(filename, temp_dir)

adata = data_analysis.PCA(adata, temp_dir, config)
adata = data_analysis.NEIGHBORS(adata, temp_dir, config)
adata = data_analysis.UMAP(adata, temp_dir, config)
adata = data_analysis.PHENOGRAPH(adata, temp_dir, config)

filename = data_IO.WRITE_ADATA(adata, filename, temp_dir)

# DATA PLOTTING 
adata = data_IO.READ_ADATA(filename, temp_dir)
coordinates = data_IO.READ_COORDINATES(filename, temp_dir)

data_plotting.UMAP(adata, filename, temp_dir, figures_dir, config)
data_plotting.CLUSTERS(adata, coordinates, filename, temp_dir, figures_dir, config)
