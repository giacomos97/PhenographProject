"""
@author: Giacomo Spisni
"""
#IMPORTS
import os
import numpy as np
import scanpy as sc
import pathlib
import configparser

import data_IO
import data_analysis

# SETUP
# Reading configuration file
config = configparser.ConfigParser()
config.read('./testing/test_configuration.txt')

# Configure directories
data_source = pathlib.Path(config.get('paths', 'data_source'))
coordinates_source = pathlib.Path(config.get('paths', 'coordinates_source'))
temp_dir = pathlib.Path(config.get('paths', 'temp_dir'))

# Configure parameters
PCA_params = dict(config.items('scanpy_PCA'))
NEIGHBORS_params = dict(config.items('scanpy_NEIGHBORS'))
UMAP_params = dict(config.items('scanpy_UMAP'))
PHENOGRAPH_params = dict(config.items('scanpy_PHENOGRAPH'))


communities_expected = [9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
                        8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
                        7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                        6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                        5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                        3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                        2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

# TESTING FUNCTIONS
def test_source_data_reading():
    '''
    Verify that shape of tables created and read by I/O functions corrispond 
    to what expected: 150x10 for measurements dataset and 150x2 for coordinates
    dataset.
    '''
    # Clean original testing files
    filename = data_IO.READ_SOURCE(data_source, coordinates_source, temp_dir)
    
    # Read work files
    adata = data_IO.READ_ADATA(filename, temp_dir)
    coordinates = data_IO.READ_COORDINATES(filename, temp_dir)
    
    # Cleanup temp datasets
    os.remove(f'{temp_dir}/{filename}.h5ad')
    os.remove(f'{temp_dir}/{filename}_coordinates.csv')
    
    assert adata.X.shape == (150,10), 'AnnData reading incomplete.'
    assert coordinates.shape == (150,2), 'Coordinates reading incomplete.'

 
def test_clustering_reproducibility():
    # Clean original testing files
    filename = data_IO.READ_SOURCE(data_source, coordinates_source, temp_dir)
    
    # Read cleaned files
    adata = data_IO.READ_ADATA(filename, temp_dir)
    coordinates = data_IO.READ_COORDINATES(filename, temp_dir)

    # Cleanup temp datasets
    os.remove(f'{temp_dir}/{filename}.h5ad')
    os.remove(f'{temp_dir}/{filename}_coordinates.csv')    

    # Run clustering 1
    adata_1 = data_analysis.PCA(adata, PCA_params)
    adata_1 = data_analysis.NEIGHBORS(adata_1, NEIGHBORS_params)
    adata_1 = data_analysis.UMAP(adata_1, UMAP_params)
    adata_1, communities_1, graph_1, Q_1 = data_analysis.PHENOGRAPH(adata_1, PHENOGRAPH_params)
    
    # Run clustering 2
    adata_2 = data_analysis.PCA(adata, PCA_params)
    adata_2 = data_analysis.NEIGHBORS(adata_2, NEIGHBORS_params)
    adata_2 = data_analysis.UMAP(adata_2, UMAP_params)
    adata_2, communities_2, graph_2, Q_2 = data_analysis.PHENOGRAPH(adata_2, PHENOGRAPH_params)
    
    assert np.array_equal(communities_1, communities_2), 'Communitiy allocation does not match between the two executions.'
    

    

def test_clustering_reliability():
    # Clean original testing files
    filename = data_IO.READ_SOURCE(data_source, coordinates_source, temp_dir)
    
    # Read cleaned files
    adata = data_IO.READ_ADATA(filename, temp_dir)
    coordinates = data_IO.READ_COORDINATES(filename, temp_dir)

    # Cleanup temp datasets
    os.remove(f'{temp_dir}/{filename}.h5ad')
    os.remove(f'{temp_dir}/{filename}_coordinates.csv')    

    # Run clustering
    adata = data_analysis.PCA(adata, PCA_params)
    adata = data_analysis.NEIGHBORS(adata, NEIGHBORS_params)
    adata = data_analysis.UMAP(adata, UMAP_params)
    adata, communities, graph, Q = data_analysis.PHENOGRAPH(adata, PHENOGRAPH_params)
    
    assert np.array_equal(communities, communities_expected), 'Communitiy allocation not as expected.'
    

    
