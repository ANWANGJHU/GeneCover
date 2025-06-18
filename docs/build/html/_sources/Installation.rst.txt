.. _introduction:

Installation
============

To install GeneCover, run:

.. code-block:: bash

    pip install git+https://github.com/ANWANGJHU/GeneCover.git 


Dependencies
------------
    - Python >= 3.8
    - numpy, scipy, gurobipy, pyscipopt

GeneCover uses **Gurobi** for solving the minimal weight set cover problem via integer programming.  
To install Gurobi for Python, follow the official instructions here:

- `Gurobi Installation Guide <https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python>`_

Gurobi requires a valid license. For instructions on obtaining one, see:

- `How to obtain a Gurobi license <https://support.gurobi.com/hc/en-us/articles/12684663118993-How-do-I-obtain-a-Gurobi-license>`_


If you are unable to obtain a Gurobi license, GeneCover provides two alternative options:

**1. PySCIPOpt:**  
GeneCover includes a full integer‐programming implementation based on PySCIPOpt, a free Python interface to the SCIP solver with no licensing requirement. To enable this solver, pass the `solver= "SCIP"` flag when calling the algorithm.

**2. Greedy heuristic:**  
When integer programming fails to converge in reasonable amount of time (a rare occurrence in most use cases), we recommend using the built-in greedy heuristic of minimal weight set covering. It returns a set cover of weight at most :math:`\sum_{i=1}^s \frac{1}{i}` times the optimal, where :math:`s` is the size of the largest set. To use this option, pass the `solver = "Greedy"` flag when calling the algorithm. This option **only requires NumPy**.

Source Code: https://github.com/ANWANGJHU/GeneCover
