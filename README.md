# PhenoGraph clustering project  
Repository dedicated to my project on PhenoGraph.  

## Pre-requisites  
To run this code, install in your environment the following libraries.  
`anndata scanpy numpy pandas pathlib plotly datetime configparser uuid`  

## How to run the code  
Following these instructions you may run the code without having to change it.  
Check the *Setup* paramaters inside the [clustering.py](./clustering.py) file if you need any adjustment.  

### Structure of input data  
The program expects **two** .csv files inside *dataset* directory:  
- **LNBmatrix.csv**, containing data on single-cell measurements.  
The document should contain comma separated columns and be organized as: one row for each cell, one column for each feature.  
For example:  
|  CD14C01CD16C02HL | CD14C01CD16C02HLA | CD1aC01SynCAMC02RB | ... |
|  0.14532181994838 |  -0.7427901168535 |   0.8510099217044  |     |
| 0.009652348004606 |  0.14532181994838 |  -1.0977624202569  |     |
|        ...        |                   |                    |     |

- **LNBcentroids.csv**, containing the spatial coordinates of each cell.  
The document should contain three comma-separated columns (without header): cell number, x coordinate, y coordinate.  
For example:  
|   1  |  5  | 7.5 |
|   2  |  8  |  7  |
|   3  | 5.2 |  7  |
|  ... |     |     |  

Cell number should have exact correspondence in both files, thus they should contain the same amount of cells.  

### Structure of the project  
Once verified that your input data is structured as expected, run the [clustering.py](./clustering.py) file.  
Check messages and results appearing in the console during execution.
At each step, a backup copy of the data will be saved in *backup* directory.
All plots shown will also be saved inside *figures* directory.  

### References  
Principal component analysis:  
https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.pca.html  
PhenoGraph:  
https://scanpy.readthedocs.io/en/stable/external/scanpy.external.tl.phenograph.html  
UMAP:  
https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.umap.html  
  
Plotting:  
https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.pca.html  
https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.umap.html  