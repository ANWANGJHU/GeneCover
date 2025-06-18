import gurobipy as grb
import numpy as np
from scipy.stats import spearmanr
from pyscipopt import Model, quicksum

def gene_gene_correlation(X, method = 'spearman'):
    """
    Compute the gene-gene correlation matrix from the gene expression matrix X.

    Args:
        X (np.ndarray or List[np.ndarray]): A matrix of shape (N, d), where N is the number of cells/spots and d is the number of genes. 
            Alternatively, a list of such matrices (from different batches/samples) with consistent gene dimensions.
        method (str): Method used to compute correlation. Must be either 'spearman' or 'pearson'.

    Returns:
        np.ndarray: A gene-gene correlation matrix of shape (d, d).
    """

    if type(X) ==list:
        corr_mat_list = []
        for x in X:
            if method == 'spearman':
                corr_mat, _ = spearmanr(x)
            if method == 'pearson':
                corr_mat = np.corrcoef(x.T)
            corr_mat_list.append(corr_mat)  
        corr_mat = np.vstack(corr_mat_list)
    else :
        if method == 'spearman':
            corr_mat, _ = spearmanr(X)
        if method == 'pearson':
            corr_mat = np.corrcoef(X.T)
    return corr_mat


def covering(Z, minSize=1, alpha=0.05, weights = 1., output=None, callBack = None,
             poolSolutions=None, poolSearchMode=None, poolGap = None, timeLimit=None, LogToConsole= 1,restart=None):
    """
    Solves a weighted gene selection problem using a Gurobi-based optimization model.

    Args:
        Z (np.ndarray): A binary matrix of shape (N, d), where N is the number of samples and d is the number of genes.
        minSize (int): The minimum number of genes to select.
        alpha (float): The minimum fraction of samples that must be covered.
        weights (np.ndarray): A 1D array of weights for each gene. Higher weights indicate higher cost for selection.
        output (int): Enables or disables solver output. Set to 1 to print optimization details, 0 to suppress.
        callBack (Callable): A callback function to be invoked during optimization.
        poolSolutions (int): Number of solutions to store in the solution pool. See: https://www.gurobi.com/documentation/current/refman/poolsolutions.html
        poolSearchMode (int): Mode for exploring the MIP search tree. See: https://www.gurobi.com/documentation/current/refman/poolsearchmode.html
        poolGap (float): Relative MIP optimality gap for accepting solutions into the pool. See:https://www.gurobi.com/documentation/current/refman/poolgap.html
        timeLimit (float): Time limit (in seconds) for the optimization run.
        LogToConsole (int): Whether to print the optimization log. Set to 1 to enable.
        restart (gurobipy.Model): A Gurobi model instance to restart the optimization from.

    Returns:
        List[int]: Indices of the selected genes.
    """

    if restart is not None:
        cov = restart
        if output is not None:
            cov.Params.OutputFlag = output
        if poolSolutions is not None:
            cov.Params.PoolSolutions = poolSolutions
        if poolSearchMode is not None:
            cov.Params.PoolSearchMode = poolSearchMode
        if poolGap is not None:
            cov.Params.PoolGap = poolGap
        if timeLimit is not None:
            cov.Params.TimeLimit = timeLimit
        if LogToConsole is not None:
            cov.Params.LogToConsole = LogToConsole
        if callBack is None:
            cov.optimize()
        else:
            cov.optimize(callBack)
        return cov

    if np.isscalar(minSize):
        minSize = [minSize]
    if np.isscalar(alpha):
        alpha = [alpha]*len(minSize)
    N = Z.shape[0]
    d = Z.shape[1]
    if type(weights) == str and weights=='prob':
        w = 1 - 0.01 * np.mean(Z, axis=0)
    elif np.isscalar(weights):
        w = weights * np.ones(d)
    else:
        w = weights
    cov = grb.Model()
    if output is not None:
        cov.Params.OutputFlag=output
    if poolSolutions is not None:
        cov.Params.PoolSolutions = poolSolutions
    if poolSearchMode is not None:
        cov.Params.PoolSearchMode = poolSearchMode
    if poolGap is not None:
        cov.Params.PoolGap = poolGap
    if timeLimit is not None:
        cov.Params.TimeLimit = timeLimit
    if LogToConsole is not None:
        cov.Params.LogToConsole = LogToConsole

    nlevels = len(minSize)
    x = []
    y = []
    for l in range(nlevels):
        x.append(cov.addMVar(d, vtype=grb.GRB.BINARY))
    for l in range(nlevels):
        y.append(cov.addMVar(N, vtype=grb.GRB.BINARY))

    for l in range(nlevels):
        expr = y[l].sum()
        cov.addConstr(expr >= N*(1-alpha[l]), 'Coverage_'+str(l))

    for l in range(nlevels):
        expr = Z @ x[l] - minSize[l]*y[l]
        cov.addConstr(expr >= 0, 'covered_' + str(l))

        # if B is not None:
        #     exprB = B @ x[l] - MinMarkerPerClass
        #     cov.addConstr(exprB >= 0, 'MinMarkerPerClass_' + str(l))

    for l in range(nlevels-1):
        for j in range(d):
            cov.addConstr(x[l+1].tolist()[j] - x[l].tolist()[j] >= 0, name= 'Nesting'+str(j)+'_'+str(l))

    expr = grb.LinExpr()
    for l in range(nlevels):
        expr += (w * x[l]).sum()
    cov.setObjective(expr, grb.GRB.MINIMIZE)
    if callBack is None:
        cov.optimize()
    else:
        cov.optimize(callBack)
    return cov


def covering_scip(Z, minSize=1, weights=1.0, timeLimit=None, output=1):
    """
    Simplified weighted set cover with alpha=0 (every element must be covered).
    Raises ValueError if any element cannot be covered by any set.
    """
    N, d = Z.shape

    # Prepare weight array
    if isinstance(weights, str) and weights == 'prob':
        w = 1 - 0.01 * np.mean(Z, axis=0)
    elif np.isscalar(weights):
        w = weights * np.ones(d, dtype=float)
    else:
        w = np.asarray(weights, dtype=float)

    # Precompute non‐zero columns per row and check coverage feasibility
    cover_indices = [np.flatnonzero(Z[i]) for i in range(N)]
    for i, cols in enumerate(cover_indices):
        if cols.size == 0:
            raise ValueError(f"Element {i} has no covering sets → problem infeasible")

    # Build SCIP model
    model = Model("SetCover_simple")
    if timeLimit is not None:
        model.setParam("limits/time", timeLimit)
    if output == 0:
        model.hideOutput()

    # Variables: x[j]=1 if set j is selected
    x = [model.addVar(vtype="B", name=f"x_{j}") for j in range(d)]

    # Objective: minimize total weight
    model.setObjective(quicksum(w[j] * x[j] for j in range(d)), "minimize")

    # Cover each element i with at least minSize selected sets
    for i, cols in enumerate(cover_indices):
        model.addCons(
            quicksum(x[j] for j in cols) >= minSize,
            name=f"cover_i{i}"
        )

    # Solve
    model.optimize()

    # Extract solution
    selected = [j for j in range(d) if model.getVal(x[j]) > 0.5]
    return selected



def greedy_weighted_set_cover(Z, w) :
    """
    Greedy heuristic for the weighted set cover problem.

    Args:
        Z (np.ndarray): A binary matrix of shape (n_elements, m_sets), where `Z[i, j] == 1` 
            indicates that set `j` covers element `i`.
        w (np.ndarray): A 1D array of length `m_sets` representing the weight of each set.

    Returns:
        List[int]: Indices of the selected sets (column indices of `Z`) that form a cover.
    """

    n, m = Z.shape
    # which elements are still uncovered
    uncovered = np.ones(n, dtype=bool)
    # which sets are still available
    available = np.ones(m, dtype=bool)
    selected = []

    while uncovered.any():
        # For each set j: how many of the still-uncovered elements it would cover?
        # Z[uncovered] is an array of shape (#uncovered_elements, m)
        cover_counts = Z[uncovered].sum(axis=0)  # shape (m,)
        # zero out the ones we've already taken
        cover_counts = np.where(available, cover_counts, 0)

        # fast-path: if the best we can do is cover exactly one element per set,
        # grab them all at once and be done
        if cover_counts.max() == Z.shape[0] // Z.shape[1]:
            singletons = np.where((available) & (cover_counts == 1))[0]
            selected.extend(singletons.tolist())
            break

        # otherwise pick the set with max (covered_new_elems / weight)
        nonzero = cover_counts > 0
        if not nonzero.any():
            # nothing left can cover any new element
            break

        ratios = np.zeros(m, dtype=float)
        ratios[nonzero] = cover_counts[nonzero] / w[nonzero]
        best = int(ratios.argmax())
        selected.append(best)

        # mark its covered elements as now covered
        # Z[:, best] is the column for set "best"
        uncovered &= ~Z[:, best].astype(bool)
        # and remove that set from future consideration
        available[best] = False
    return selected



def GeneCover(num_marker, corr_mat, w, m = 3,interval = 0, lambdaMax = .3, lambdaMin = 0.05, timeLimit = 600, output = 0, solver = "Gurobi") :
    """
    Selects marker genes based on gene-gene correlation using combinatorial optimization or a greedy heuristic.

    Args:
        num_marker (int): Desired number of markers to select.
        corr_mat (np.ndarray): Gene-gene correlation matrix.
        interval (int): Allowed deviation from `num_marker`. The final number of markers may vary within this range.
        w (np.ndarray): An array of weights for the genes. Higher weights indicate higher cost for selection.
        lambdaMax (float): Maximum threshold for acceptable gene-gene correlation.
        lambdaMin (float): Minimum threshold for acceptable gene-gene correlation.
        timeLimit (float): Time limit (in seconds) for the optimization.
        ouput (int): Whether to print the optimization process. Set to 1 to enable.
        solver (str): The solver to use for the optimization. Options are "Gurobi", "SCIP", and "Greedy".
        greedy (bool): Whether to use a greedy algorithm for set cover instead of the Gurobi solver. Default: False.

    Returns:
        List[int]: Indices of the selected marker genes.
    """

    epsilon = (lambdaMax + lambdaMin)/2
    best_marker_length_gap = 1e6
    selection = np.arange(corr_mat.shape[1])
    G_v3 = corr_mat > epsilon
    if solver == "Gurobi":
        cov_sol = covering(G_v3, minSize=1, alpha=0.0, weights=w, timeLimit=timeLimit, output = output)
        cov_sol = selection[np.array(cov_sol.x)[:len(selection)] > 0.5]
    elif solver == "Greedy":
        cov_sol = greedy_weighted_set_cover(G_v3, w)
    elif solver == "SCIP":
        cov_sol = covering_scip(G_v3, minSize=1, weights=w, timeLimit=timeLimit, output=output)
    else:
        raise ValueError("Invalid solver specified. Choose from 'Gurobi', 'SCIP', or 'Greedy'.")
    markers = []
    num_batches = G_v3.shape[0] // G_v3.shape[1]
    num_genes = G_v3.shape[1]
    for i in cov_sol:
        if num_batches > 1:
            if G_v3[[i + j * num_genes for j in range(num_batches)]].sum(axis = 1).min() >= m:
                markers.append(i)
        else:
            if G_v3[i].sum() >= m:
                markers.append(i)
    n_markers = len(markers)
    current_gap = abs(n_markers - num_marker)
    best_marker_length_gap = current_gap
    best_epsilon = epsilon
    while (lambdaMax - lambdaMin) > 1e-6 and (n_markers < num_marker or n_markers > num_marker + interval):
        if n_markers< num_marker:
            lambdaMax = epsilon
        else:
            lambdaMin = epsilon
        epsilon = (lambdaMin+lambdaMax)/2
        G_v3 = corr_mat > epsilon
        if solver == "Gurobi":
            cov_sol = covering(G_v3, minSize=1, alpha=0.0, weights=w, timeLimit=timeLimit,output = output)
            cov_sol = selection[np.array(cov_sol.x)[:len(selection)] > 0.5]
        elif solver == "Greedy":
            cov_sol = greedy_weighted_set_cover(G_v3, w)
        elif solver == "SCIP":
            cov_sol = covering_scip(G_v3, minSize=1,  weights=w, timeLimit=timeLimit, output=output)
        else:
            raise ValueError("Invalid solver specified. Choose from 'Gurobi', 'SCIP', or 'Greedy'.")
        markers = []
        for i in cov_sol:
            if num_batches > 1:
                if G_v3[[i + j * num_genes for j in range(num_batches)]].sum(axis = 1).min() >= m:
                    markers.append(i)
            else:
                if G_v3[i].sum() >= m:
                    markers.append(i)
            n_markers = len(markers)

        current_gap = abs(n_markers - num_marker)
        if current_gap < best_marker_length_gap:
            best_marker_length_gap = current_gap
            best_epsilon = epsilon
            best_lambdaMin = lambdaMin
            best_lambdaMax = lambdaMax
            best_direction = n_markers < num_marker
    print("Best Gap: ", best_marker_length_gap)
    print("Best Epsilon: ", best_epsilon)
    return markers

def Iterative_GeneCover(incremental_sizes,corr_mat, w,m = 3, lambdaMin = .05,lambdaMax = .3, timeLimit = 600, output = 0, solver = "Gurobi"):
    
    """
    Performs iterative marker gene selection using the GeneCover algorithm.

    Args:
        corr_mat (np.ndarray): Gene-gene correlation matrix of shape (d, d).
        incremental_sizes (List[int]): A list indicating the number of markers to select at each iteration.
        w (np.ndarray): An array of weights for each gene. Higher weights indicate higher cost for selection.
        lambdaMax (float): Maximum threshold for gene-gene correlation.
        lambdaMin (float): Minimum threshold for gene-gene correlation.
        timeLimit (float): Time limit (in seconds) for the optimization.
        output (int): Whether to print the optimization process. Set to 1 to enable.
        solver (str): The solver to use for the optimization. Options are "Gurobi", "SCIP", and "Greedy".
    Returns:
        List[List[int]]: A list where each element is a list of indices of the selected marker genes at the corresponding iteration.
    """
    num_batches = corr_mat.shape[0] // corr_mat.shape[1]
    num_genes = corr_mat.shape[1]
    MARKERS = []
    print("Iteration 1")
    markers = GeneCover(incremental_sizes[0],corr_mat,w = w,m =m,lambdaMax=lambdaMax,lambdaMin=lambdaMin,timeLimit=timeLimit,output=output,solver=solver)
    selection = np.arange(corr_mat.shape[1])  
    MARKERS.append(markers)
    remaining_genes_idx_abs = np.setdiff1d(selection, markers)
    for t, size in enumerate(incremental_sizes[1:]): 
        print("Iteration ", t+2)
        remaining_genes_idx_abs_batches = np.array([remaining_genes_idx_abs + j * num_genes for j in range(num_batches)]).flatten()
        corr_mat_remain=corr_mat[remaining_genes_idx_abs_batches][:,remaining_genes_idx_abs]
        markers=GeneCover(size,corr_mat_remain,w=w[remaining_genes_idx_abs],m=m,lambdaMin=lambdaMin,lambdaMax=lambdaMax,timeLimit=timeLimit,output=output,solver=solver)
        MARKERS.append(remaining_genes_idx_abs[markers])
        remaining_genes_idx_abs = np.setdiff1d(remaining_genes_idx_abs, [j for i in MARKERS for j in i])
    return MARKERS