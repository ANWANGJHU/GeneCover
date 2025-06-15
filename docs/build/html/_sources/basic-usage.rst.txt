.. _usage:

Basic Usage
===========

This section shows how to use **GeneCover** to compute gene-gene correlation and select marker genes.

Required Imports
----------------

.. code-block:: python

    from genecover import *
    import numpy as np

Data Format
-----------

- **data**: A NumPy array of shape `(N, d)` storing gene expression measurements,
  where `N` is the number of cells/spots and `d` is the number of genes.
  
  Alternatively, `data` can be a list of such arrays (e.g., one per donor or batch), as long as they all have the same number of genes.

- **w**: A 1D array of length `d` representing the cost (weight) of selecting each gene.  
  A simple choice is:

  .. code-block:: python

      w = np.ones(data.shape[1])

Gene-Gene Correlation
---------------------

Compute the gene-gene correlation matrix:

.. code-block:: python

    corr_mat = gene_gene_correlation(data)

Gene Selection (Single Run)
---------------------------

Obtain 100 marker genes using the standard GeneCover optimization:

.. code-block:: python

    markers = GeneCover(
        num_marker=100,
        corr_mat=corr_mat,
        w=w)


or obtain 100 marker genes using the greedy heuristic approach, which does not require Gurobi:

.. code-block:: python

    markers = GeneCover(
        num_marker=100,
        corr_mat=corr_mat,
        w=w,
        greedy=True)

Iterative Gene Selection
------------------------

Select 300 markers in three iterative steps (e.g., 100 markers per round):

.. code-block:: python

    iterative_markers = Iterative_GeneCover(
        incremental_sizes=[100, 100, 100],
        corr_mat=corr_mat,
        w=w)

You can also use the greedy heuristic for iterative selection:

.. code-block:: python

    iterative_markers = Iterative_GeneCover(
        incremental_sizes=[100, 100, 100],
        corr_mat=corr_mat,
        w=w,
        greedy=True)