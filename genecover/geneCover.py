import gurobipy as grb
import numpy as np
from scipy.stats import spearmanr

def covering(Z, minSize=1, alpha=0.05, weights = 1., output=None, callBack = None,
             poolSolutions=None, poolSearchMode=None, poolGap = None, timeLimit=None, LogToConsole= 1,restart=None):
    """
    :param Z: a binary matrix of size N x d, where N is the number of samples and d is the number of genes
    :param minSize: the minimum number of genes to select
    :param alpha: the minimum fraction of samples that should be covered
    :param weights: an array of weights for each gene (the higher the weight, the higher cost to select the gene)
    :param output: Enables or disables solver output. Set to 1 to print the optimization process and 0 to disable it.
    :param callBack: a callback function to be called during optimization
    :param poolSolutions: the number of solutions to store in the pool: https://www.gurobi.com/documentation/current/refman/poolsolutions.html#parameter:PoolSolutions 
    :param poolSearchMode: modes for exploring the MIP search tree: https://www.gurobi.com/documentation/current/refman/poolsearchmode.html#parameter:PoolSearchMode 
    :param poolGap: the relative MIP optimality gap for the pool: https://www.gurobi.com/documentation/current/refman/poolgap.html#parameter:PoolGap 
    :param timeLimit: time limit for the optimization
    :param LogToConsole: whether to print the optimization log. (Set to 1 to print)
    :param restart: a Gurobi model to restart the optimization process

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


def GeneCover(num_marker, corr_mat, w, m = 3,interval = 0, lambdaMax = .3, lambdaMin = 0.05, timeLimit = 600, output = 0):
    """
    Args:
    :param num_marker: number of markers to select
    :param corr_mat: gene-gene correlation matrix
    :param interval: allowed deviation from num_marker
    :param w: an array of weights for the genes (the higher the weight, the higher cost to select the gene)
    :param lambdaMax: maximum threshold for gene-gene correlation
    :param lambdaMin: minimum threshold for gene-gene correlation
    :param timeLimit: time limit for the optimization
    :param ouput: whether to print the optimization process. (Set to 1 to print)

    Returns:
    :return: the list of indices of the selected markers

    """
    epsilon = (lambdaMax + lambdaMin)/2
    best_marker_length_gap = 1e6
    selection = np.arange(corr_mat.shape[1])
    G_v3 = corr_mat > epsilon
    cov_sol = covering(G_v3, minSize=1, alpha=0.0, weights=w, timeLimit=timeLimit, output = output)
    cov_sol = selection[np.array(cov_sol.x)[:len(selection)] > 0.5]
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
        cov_sol = covering(G_v3, minSize=1, alpha=0.0, weights=w, timeLimit=timeLimit,output = output)
        cov_sol = selection[np.array(cov_sol.x)[:len(selection)] > 0.5]
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

def Iterative_GeneCover(incremental_sizes,corr_mat, w,m = 3, interval = 0, lambdaMin = .05,lambdaMax = .3, timeLimit = 600, output = 0):
    
    """
    Args:
    :param corr_mat: gene-gene correlation matrix
    :param incremental_sizes: a tuple of number of markers to select at each iteration
    :param w: an array of weights for each gene (the higher the weight, the higher cost to select the gene)
    :param lambdaMax: maximum threshold for gene-gene correlation
    :param lambdaMin: minimum threshold for gene-gene correlation
    :param timeLimit: time limit for the optimization
    :param output: wheth er to print the optimization process. (Set to 1 to print)

    Returns:
    :return: a list of lists of indices of the selected markers at each iteration
    """
    num_batches = corr_mat.shape[0] // corr_mat.shape[1]
    num_genes = corr_mat.shape[1]
    MARKERS = []
    print("Iteration 1")
    markers = GeneCover(incremental_sizes[0], corr_mat, w = w, m =m, interval = interval , lambdaMax = lambdaMax, lambdaMin = lambdaMin, timeLimit = timeLimit, output=output)
    selection = np.arange(corr_mat.shape[1])  
    MARKERS.append(markers)
    remaining_genes_idx_abs = np.setdiff1d(selection, markers)
    for t, size in enumerate(incremental_sizes[1:]): 
        print("Iteration ", t+2)
        remaining_genes_idx_abs_batches = np.array([remaining_genes_idx_abs + j * num_genes for j in range(num_batches)]).flatten()
        corr_mat_remain = corr_mat[remaining_genes_idx_abs_batches][:,remaining_genes_idx_abs]
        markers = GeneCover(size, corr_mat_remain, w = w[remaining_genes_idx_abs], m =m, lambdaMin= lambdaMin,lambdaMax=lambdaMax, timeLimit = timeLimit, output=output)
        MARKERS.append(remaining_genes_idx_abs[markers])
        remaining_genes_idx_abs = np.setdiff1d(remaining_genes_idx_abs, [j for i in MARKERS for j in i])  
    return MARKERS

def gene_gene_correlation(X, method = 'spearman'):
    """
    Args:
    :param X: an numpy array of size N x d, where N is the number of cells / spots and d is the number of genes; 
    X can also be a list of numpy arrays where each array in the list is the gene expression matrix from a batch / sample (the number of genes should be the same across all batches / samples)  
    

    Returns:
    :return: gene-gene correlation matrix of dimension d x d
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