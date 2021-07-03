"""
@author: Giacomo Spisni
"""

# IMPORTS
import anndata
import pandas as pd
from datetime import datetime
import uuid

# FUNCTIONS
def READ_SOURCE(data_source, coordinates_source, temp_dir):
    '''
    Read original dataset containing single-cell measurements and spatial coordinates,
    then store their cleaned copy with an uniquely identifying name.

    Parameters
    ----------
    data_source : path object of pathlib module
        Directory where to read the original dataset containing single-cell measurements.
    coordinates_source : path object of pathlib module
        Directory where to read the original dataset containing spatial coordinates.
    temp_dir : path object of pathlib module
        Directory where to store the temporary copy of the datasets.

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
    data.rename(columns = {data.columns[0]: 'cell'}, inplace = True)
    
    coordinates = pd.read_csv(coordinates_source, header = None, names = ['cell','x','y'])
    
    # Sort by cell number to ensure 1-1 correspondence of rows
    data.sort_values(by = ['cell'], inplace = True)
    coordinates.sort_values(by = ['cell'], inplace = True)

    # Clean up data
    data.drop(['cell'], axis = 1, errors = 'ignore', inplace = True)
    coordinates.drop(['cell'], axis = 1, errors = 'ignore', inplace = True)

    data.astype('float', copy = False)
    coordinates.astype('float', copy = False)
    
    # Converto to AnnData and save to temp file
    adata = anndata.AnnData(data)    
    adata.write(f'{temp_dir}/{filename}.h5ad')
    
    # Save coordinates file
    coordinates.to_csv(f'{temp_dir}/{filename}_coordinates.csv', header = True, index = False)

    print(f"[{datetime.now()}] Data reading completed.")
    
    return filename


def READ_ADATA(filename, temp_dir):
    '''
    Read and return AnnData object from temporary directory.

    Parameters
    ----------
    filename : str
        Unique identifier generated for this dataset stored in temp_dir.
    temp_dir : path object of pathlib module
        Directory where to store the temporary copy of the datasets.

    Returns
    -------
    adata : AnnData object
        Database containing work data.

    '''
    
    adata = anndata.read_h5ad(f'{temp_dir}/{filename}.h5ad')
    
    return adata


def WRITE_ADATA(adata, filename, temp_dir):
    '''
    Write AnnData object to temporary directory.

    Parameters
    ----------
    adata : AnnData object
        Database containing work data.
    filename : str
        Unique identifier generated for this dataset stored in temp_dir.
    temp_dir : path object of pathlib module
        Directory where to store the temporary copy of the datasets.

    Returns
    -------
    filename : str
        Unique identifier generated for this dataset stored in temp_dir.

    '''
    
    adata.write(f'{temp_dir}/{filename}.h5ad')
    
    return filename


def READ_COORDINATES(filename, temp_dir):
    '''
    Read and return Pandas DataFrame object from temporary directory.

    Parameters
    ----------
    filename : str
        Unique identifier generated for this dataset stored in temp_dir.
    temp_dir : path object of pathlib module
        Directory where to store the temporary copy of the datasets.

    Returns
    -------
    coordinates : Pandas DataFrame object
        DataFrame containing cells' coordinates.

    '''
    
    coordinates = pd.read_csv(f'{temp_dir}/{filename}_coordinates.csv', header = "infer")
    
    return coordinates