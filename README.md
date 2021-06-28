# PhenoGraph clustering project  
Repository dedicated to my project on PhenoGraph.  

## Pre-requisites  
To run this code, install in your environment the following libraries.  
`anndata scanpy numpy pandas pathlib plotly datetime configparser uuid`  

## How to run the code  
Following these instructions you may run the code without having to change it. First, you will set your own preferences into the [configuration.txt](./configuration.txt) file, such as the names of the dataset directories or the clustering parameters. Then you will have two options:  
- run the [clustering.py](./clustering.py) file in your console  
- use the interactive [demo_notebook.ipynb](./demo_notebook.ipynb)  

Both documents implement the same procedure.

### Structure of input data  
The program expects **two** .csv files:  
- **One file containing data on single-cell measurements**.  
The document should contain comma separated columns and be organized as: first column as cell number, one row for each cell, one more column for each feature.  
For example:  

|     | CD14C01CD16C02HL  | CD14C01CD16C02HLA | CD1aC01SynCAMC02RB | ... |
|:---:|:-----------------:|:-----------------:|:------------------:|:---:|
|  1  | 0.145321181994838 |  -0.7427901168535 |   0.8510099217044  |     |
|  2  | 0.009652348004606 |  0.14785165558968 |  -1.0977624202569  |     |
|  3  | 1.254862485568716 |  0.25684126698538 |  -0.0489567128639  |     |
| ... |                   |                   |                    |     |  

- **Another file containing the spatial coordinates of each cell**.  
The document should contain three comma-separated columns (without header): cell number, x coordinate, y coordinate.  
For example:  

|      |     |     |
|:----:|:---:|:---:|
|   1  |  5  | 7.5 |
|   2  |  8  |  7  |
|   3  | 5.2 |  7  |
|  ... |     |     |  
  
Cell number must have exact correspondence in both files (i.e. cell number 1 corresponds to first row in measurements file). Hence both files should contain the same amount of cells.  

### Preparation of *configuration.txt*  
Details available soon.  

### Structure of the project  
Verify that your input data is structured as expected, and tune the [configuration.txt](./configuration.txt) according to your needs. Once done, run the [clustering.py](./clustering.py) file.  
Check messages and results appearing in the console during execution.
At each step, a backup copy of the data will be saved in the *temp* directory you selected.
All plots shown will also be saved inside the *figures* directory you also selected.  

### References  
Check the following references for more details on the parameters available in each function.  

  
Principal component analysis (PCA):  
https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.pca.html  
Uniform Manifold Approximation and Projection for Dimension Reduction (UMAP):  
https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.umap.html  
PhenoGraph:  
https://scanpy.readthedocs.io/en/stable/external/scanpy.external.tl.phenograph.html  

  
Plotting:  
https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.pca.html  
https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.umap.html  