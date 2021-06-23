"""
@author: Giacomo Spisni
"""

# IMPORTS
import scanpy as sc
import pathlib
import configparser
import data_reading
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
sc.settings.verbose = 3
sc.set_figure_params(dpi = 600)
sc.settings.figdir = figures_dir

# Prepare workspace
pathlib.Path(figures_dir).mkdir(parents = True, exist_ok = True)
pathlib.Path(temp_dir).mkdir(parents = True, exist_ok = True)


# DATA READING
filename = data_reading.READ(data_source, coordinates_source, temp_dir)

# DATA ANALYSIS
data_analysis.PCA(filename, temp_dir, config)
data_analysis.NEIGHBORS(filename, temp_dir, config)
data_analysis.UMAP(filename, temp_dir, config)
data_analysis.PHENOGRAPH(filename, temp_dir, config)

# DATA PLOTTING 
data_plotting.UMAP(filename, temp_dir, figures_dir)
data_plotting.CLUSTERS(filename, temp_dir, figures_dir)
