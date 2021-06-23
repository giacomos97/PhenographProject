"""
@author: Giacomo Spisni
"""

# IMPORTS
from anndata import AnnData
import pandas as pd
from datetime import datetime
import uuid

# FUNCTIONS
def READ(data_source, coordinates_source, temp_dir):
    '''
    Parameters
    ----------
    data_source : TYPE
        DESCRIPTION.
    coordinates_source : TYPE
        DESCRIPTION.
    temp_dir : TYPE
        DESCRIPTION.

    Returns
    -------
    filename : TYPE
        DESCRIPTION.

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