"""
@author: Giacomo Spisni
"""

# IMPORTS
from anndata import AnnData
import pandas as pd
from datetime import datetime
import uuid

# FUNCTIONS
def READ(data_source, coordinates_source, temp_dir, config):
    '''
    Read original dataset containing single-cell measurements and spatial coordinates,
    then store their cleaned copy with an uniquely identifying name.

    Parameters
    ----------
    data_source : path object of pathlib module.
        Directory where to read the original dataset containing single-cell measurements.
    coordinates_source : path object of pathlib module.
        Directory where to read the original dataset containing spatial coordinates.
    temp_dir : path object of pathlib module.
        Directory where to store the temporary copy of the datasets.
    config : ConfigParser object of configparser module.
        Configuration file.

    Returns
    -------
    filename : str
        Unique identifier generated for this dataset stored in temp_dir.

    '''
    
    # Generate an unique identifier for this run
    filename = uuid.uuid4().hex
    
    print(f'Current run ID: {filename}.')
    print(f"[{datetime.now()}] Data reading started.")
    
    # Read files 
    data = pd.read_csv(data_source)
    coordinates = pd.read_csv(coordinates_source, header = None, names = ['cell','x','y'])

    # Clean up data
    data.drop(['Unnamed: 0'], axis = 1, errors = 'ignore', inplace = True)

    data.astype('float', copy = False)
    coordinates.astype('float', copy = False)
    
    # Converto to AnnData and save to temp file
    adata = AnnData(data)    
    adata.write(f'{temp_dir}/{filename}.h5ad')
    
    # Save coordinates file
    coordinates.to_csv(f'{temp_dir}/{filename}_coordinates.csv', header = True, index = False)

    print(f"[{datetime.now()}] Data reading completed.")
    
    return filename