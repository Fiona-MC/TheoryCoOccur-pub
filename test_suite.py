import numpy as np
from importlib import reload
import networkx as nx

import metaCommunityMx
import community_sim
import cluster_funcs
reload(metaCommunityMx)
reload(community_sim)
reload(cluster_funcs)

### Test functions in metaCommunityMx module
testMx_rand_1 = metaCommunityMx.getRandomMatrix(maxInteractionStrength = 0.1, pPositive = 0.5, pInteract = 1, nSpecies = 10, pBidirectional = 0)
testMx_rand_1_pos = metaCommunityMx.getRandomMatrix(maxInteractionStrength = 0.1, pPositive = 1, pInteract = 1, nSpecies = 10, pBidirectional = 0)
testMx_rand_1_neg = metaCommunityMx.getRandomMatrix(maxInteractionStrength = 0.1, pPositive = 0, pInteract = 1, nSpecies = 10, pBidirectional = 0)
testMx_rand_0 = metaCommunityMx.getRandomMatrix(maxInteractionStrength = 0.1, pPositive = 0.5, pInteract = 0, nSpecies = 10, pBidirectional = 0)
testMx_rand_0_5_1 = metaCommunityMx.getRandomMatrix(maxInteractionStrength = 0.1, pPositive = 0.5, pInteract = 0.5, nSpecies = 10, pBidirectional = 0)
testMx_rand_0_5_2 = metaCommunityMx.getRandomMatrix(maxInteractionStrength = 0.1, pPositive = 0.5, pInteract = 0.5, nSpecies = 10, pBidirectional = 0)

testMx_mod = metaCommunityMx.generateModularMatrix(maxInteractionStrength = 0.1, pPositive = 0.5, nSpecies = 10, pWithinMod = 1, pBetweenMod = 0, nModules = 2)
#testMx3_mod = metaCommunityMx.generateModularMatrix(maxInteractionStrength = 0.1, pPositive = 0.5, nSpecies = 10, pWithinMod = 1, pBetweenMod = 0, nModules = 3)
#testMxi_mod = metaCommunityMx.generateModularMatrix(maxInteractionStrength = 0.1, pPositive = 0.5, nSpecies = 10, pWithinMod = 0, pBetweenMod = 1, nModules = 2)

testMx_hub = metaCommunityMx.generateHubMatrix(maxInteractionStrength = 0.1, pPositive = 0.5, nSpecies = 10, pWithinMod = 0, pBetweenMod = 0, nModules = 2, pHub = 1)
testMx_hub_mod = metaCommunityMx.generateHubMatrix(maxInteractionStrength = 0.1, pPositive = 0.5, nSpecies = 10, pWithinMod = 1, pBetweenMod = 0, nModules = 2, pHub = 1)
#testMxi_hub = metaCommunityMx.generateHubMatrix(maxInteractionStrength = 0.1, pPositive = 0.5, nSpecies = 10, pWithinMod = 0, pBetweenMod = 1, nModules = 2, pHub = 0)


def test_getRandomMatrix():
    assert np.all(np.diagonal(testMx_rand_1) == 0)
    assert testMx_rand_1.shape == (10,10)
    assert np.all(np.absolute(testMx_rand_1) <= 0.1)
    assert np.all(np.absolute(testMx_rand_0_5_1) <= 0.1)
    assert np.all(np.absolute(testMx_rand_1) + np.identity(10) > 0)
    assert np.all(np.absolute(testMx_rand_0) == 0)
    assert not np.all(testMx_rand_0_5_1 == testMx_rand_0_5_2)
    assert not np.all(testMx_rand_1 > 0)
    assert not np.all(testMx_rand_1 < 0)
    assert np.all(testMx_rand_1_pos + np.identity(10) > 0)
    assert np.all(testMx_rand_1_neg - np.identity(10) < 0)

def test_generateModularMatrix():
    assert np.all(np.diagonal(testMx_mod) == 0)
    assert testMx_mod.shape == (10,10)
    assert np.all(np.absolute(testMx_mod[0:5, 0:5] + np.identity(5)) > 0)
    assert np.all(np.absolute(testMx_mod[5:, 5:] + np.identity(5)) > 0)
    assert np.all(np.absolute(testMx_mod[0:5, 5:]) == 0)
    assert np.all(np.absolute(testMx_mod[5:, 0:5]) == 0)

def test_generateHubMatrix():
    assert np.all(np.diagonal(testMx_hub) == 0)
    assert testMx_hub.shape == (10,10)
    assert np.all((testMx_hub_mod == 0) == (testMx_mod == 0))
    mod1 = testMx_hub[0:5, 0:5]
    mod2 = testMx_hub[5:, 5:]
    assert np.all(np.absolute(mod1[0, 1:]) > 0)
    assert np.all(np.absolute(mod1[1:, 0]) > 0)
    assert np.all(np.absolute(mod2[0, 1:]) > 0)
    assert np.all(np.absolute(mod2[1:, 0]) > 0)

    assert np.all(np.absolute(testMx_hub[1:5, 1:5]) == 0)
    assert np.all(np.absolute(testMx_hub[1:5, 6:]) == 0)

### Test functions in community_sim module
def test_run_communityDynamics():
    assert True

def test_simulate_coOccurrence():
    assert True



### Test functions in cluster_funcs module
testMx = np.array([[0,1,1,0,0,0],
                   [1,0,1,-4,0,0],
                   [1,0,0,1,1,0],
                   [1,0,0,1,.1,0],
                   [1,0,0,1,1,0],
                   [-1,0,0,1,1,0]])
testMxDiag = testMx
np.fill_diagonal(testMxDiag, 2)

testMx3 = np.array([[0,1,1,0,0],
                   [1,0,1,-4,0],
                   [1,0,0,1,1],
                   [1,0,0,1,.1],
                   [1,0,0,1,1]])

edgeCountWithin, possibleWithin, edgeCountBetween, possibleBetween = cluster_funcs.prop_within_btwn_mod(nSpec = testMx.shape[0], nModules = 2, matrix = testMx)
edgeCountWithin2, possibleWithin2, edgeCountBetween2, possibleBetween2 = cluster_funcs.prop_within_btwn_mod(nSpec = testMxDiag.shape[0], nModules = 2, matrix = testMxDiag)
edgeCountWithin3, possibleWithin3, edgeCountBetween3, possibleBetween3 = cluster_funcs.prop_within_btwn_mod(nSpec = testMx3.shape[0], nModules = 2, matrix = testMx3)

def test_prop_within_btwn_mod():
    assert edgeCountWithin == 9
    assert edgeCountWithin2 == 9
    assert possibleWithin == 12
    assert possibleWithin2 == 12
    assert edgeCountBetween == 6
    assert edgeCountBetween2 == 6
    assert possibleBetween == 18
    assert possibleBetween2 == 18
    assert edgeCountWithin3 == 3
    assert possibleWithin3 == 4
    assert edgeCountBetween3 == 9
    assert possibleBetween3 == 16

def test_count_all_hubs():
    totalHub, totalpossHub, totalnonHub, totalpossNonHub = cluster_funcs.count_all_hubs(nSpec = testMx.shape[0], nModules = 2, matrix = testMx)
    assert totalHub == 7
    assert totalpossHub == 8
    assert totalnonHub == 2
    assert totalpossNonHub == 4
    totalHub2, totalpossHub2, totalnonHub2, totalpossNonHub2 = cluster_funcs.count_all_hubs(nSpec = testMx3.shape[0], nModules = 2, matrix = testMx3)
    assert totalHub2 == 3
    assert totalpossHub2 == 4
    assert totalnonHub2 == 0
    assert totalpossNonHub2 == 0