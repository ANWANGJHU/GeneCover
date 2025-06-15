.. _introduction:

Installation
============

To install GeneCover, run:

.. code-block:: bash

    pip install git+https://github.com/ANWANGJHU/GeneCover.git 


Dependencies
------------
    - Python >= 3.8
    - numpy, scipy, gurobipy

GeneCover uses **Gurobi** for solving the minimal weight set cover problem via integer programming.  
To install Gurobi for Python, follow the official instructions here:

- `Gurobi Installation Guide <https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python>`_

Gurobi requires a valid license. For instructions on obtaining one, see:

- `How to obtain a Gurobi license <https://support.gurobi.com/hc/en-us/articles/12684663118993-How-do-I-obtain-a-Gurobi-license>`_


If you are unable to obtain a Gurobi license, GeneCover provides a **greedy heuristic** fallback, which solves the set cover problem approximately.

To use this option, pass the `greedy=True` flag when calling the algorithm.  
This option **only requires NumPy** and runs without Gurobi.

Source Code: https://github.com/ANWANGJHU/GeneCover
