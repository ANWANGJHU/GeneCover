# GeneCover

GeneCover is a computational tool designed for scRNA-seq and spatial transcriptomics data to select pre-labeling marker genes based on gene-gene correlations.

### Installation
To install GeneCover, run:
`pip install git+https://github.com/ANWANGJHU/GeneCover.git
`

**Note: We strongly recommend installing the package in a virtual environment to avoid dependency conflicts.**

### Dependency: Gurobi
GeneCover requires Gurobi as a dependency for optimization. To install Gurobi in python, please follow [Gurobi Installation Guide](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python). Gurobi requires a license to use the interface. For information on how to obtain a license, please refer to  [this Gurobi support article](https://support.gurobi.com/hc/en-us/articles/12684663118993-How-do-I-obtain-a-Gurobi-license).

### Tutorial 
This tutorial provides a simple example of how to run GeneCover on your data.

```python 
from genecover import *
import numpy as np
"""
User Input: 
data - an array of size N * d that store the log-normalized count expression data
w - an array of size d that records cost of each gene. (By default, we let w = np.ones(data.shape[1]))
"""

# Compute the gene-gene correlation matrix
corr_mat = gene_gene_correlation(data)

# Obtain 100 GeneCover markers
markers = GeneCover(num_marker = 100, corr_mat = corr_mat, w = w, m = 3, lambdaMax = .3, lambdaMin = 0.05)

# Obtain 300 GeneCover markers through iterative selection with three iterations
iterative_markers = Iterative_GeneCover(incremental_sizes = [100, 100, 100], corr_mat = corr_mat, w = w, m = 3, lambdaMax = .3, lambdaMin = 0.05)
```


