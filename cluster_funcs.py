import numpy as np

def prop_within_btwn_mod(nSpec, nModules, matrix):
    matrix_copy = np.copy(matrix)
    assert matrix_copy.shape[0] == matrix_copy.shape[1]
    np.fill_diagonal(matrix_copy, 0) # make sure there are no nonzero entries on the diagonal
    edgeCountWithin = 0
    modSize = nSpec // nModules
    for mod_i in range(nModules):
        for spec_i in range(modSize):
            i = mod_i * modSize + spec_i
            for spec_j in range(modSize):
                j = mod_i * modSize + spec_j
                if matrix_copy[i,j] != 0:
                    edgeCountWithin += 1
                    matrix_copy[i,j] = 0
    edgeCountBetween = len(np.nonzero(matrix_copy)[0]) # remaining edges are edges between modules
    possibleWithin = nModules * (modSize**2 - modSize)
    possibleBetween = (nSpec**2 - nSpec) - possibleWithin
    return edgeCountWithin, possibleWithin, edgeCountBetween, possibleBetween

def count_one_hub(matrix, modSize, hubIndex):
    ''' Assumes the first node is the hub node. Counts the number of edges at the hub node within the module and in the module not from the hub node.
    '''
    subsetMx = matrix[hubIndex:(hubIndex + modSize), hubIndex:(hubIndex + modSize)]
    np.fill_diagonal(subsetMx, 0)
    subsetMx_nz = np.not_equal(subsetMx, 0).astype(int)
    hubEdges = np.sum(subsetMx_nz[0,:]) + np.sum(subsetMx_nz[:,0])
    possibleHubEdges = 2 * (subsetMx_nz.shape[0] - 1)
    nonHubEdges = np.sum(subsetMx_nz[1:, 1:])
    possibleNonHubEdges = (subsetMx_nz.shape[0] - 1) ** 2 - (subsetMx_nz.shape[0] - 1)
    return hubEdges, possibleHubEdges, nonHubEdges, possibleNonHubEdges

def count_all_hubs(nSpec, nModules, matrix):
    totalHub = 0
    totalpossHub = 0
    totalnonHub = 0
    totalpossNonHub = 0
    modSize = nSpec // nModules
    hub_nodes = [modi * modSize for modi in range(nModules)]
    for node in hub_nodes:
        hubEdges, possibleHubEdges, nonHubEdges, possibleNonHubEdges = count_one_hub(matrix, modSize, hubIndex = node)
        totalHub += hubEdges
        totalpossHub += possibleHubEdges
        totalnonHub += nonHubEdges
        totalpossNonHub += possibleNonHubEdges
    return totalHub, totalpossHub, totalnonHub, totalpossNonHub