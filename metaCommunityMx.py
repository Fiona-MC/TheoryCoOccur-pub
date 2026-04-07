import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def getRandomMatrix(maxInteractionStrength, pPositive, pInteract, nSpecies, pBidirectional):
    '''Make a matrix of interactions with a pre-defined probabiluty of being a positive interaciton and being a bidirectional interaciton.
    
    maxInteractionStrength (float) 
    pPositive (float): Probability that an interaction is positive, given that there is an interaction
    pInteract (float): probability of an interaction between two species (Note: not exact because pBidirectional is added later)
    nSpecies (int): number of species
    pBidirectional: probability an interaction is bidirectional given that it exists

    Returns a numpy array of interactions (diagonal is zero)

    Note: Probabilities are not exact because the diagonal is zeroed at the end.
    '''
    interactionDirection = np.random.choice([-1,0,1], size = (nSpecies, nSpecies), p = [(1-pPositive) * pInteract, 1-pInteract, pPositive * pInteract])
    interactionStrength = np.random.uniform(low = 0, high = maxInteractionStrength, size = (nSpecies, nSpecies))

    # add an interaction in the other direction 
    randomMatrix = interactionDirection * interactionStrength
    addBidirectional = np.transpose(randomMatrix) != 0 # get the locations to add the other direction of interaction
    addBidirectional = addBidirectional * (randomMatrix == 0)
    addBidirectional = addBidirectional * np.random.choice([0,1], size = (nSpecies, nSpecies), p = [1-pBidirectional, pBidirectional])

    # choose actual values for interactions in the other direction
    addBidirectional = addBidirectional * np.random.choice([-1,1], size = (nSpecies, nSpecies), p = [(1-pPositive), pPositive])
    addBidirectional = addBidirectional * np.random.uniform(low = 0, high = maxInteractionStrength, size = (nSpecies, nSpecies))

    randomMatrix = randomMatrix + addBidirectional
    np.fill_diagonal(randomMatrix, 0)
    return randomMatrix

# generate one module
def generateModule(modSize, modType, maxInteractionStrength, pWithinMod, pPositive, pHub = None):
    '''Helper function for generating modular or hub matrices
    
    modSize is the dimension/number of species in the module, 
    modType is either "module" or "hub". If modType = hub, the first species in the matrix will be the hub node. If modType is "hub" you must define pHub as the 
    probability of the hub node connecting to other nodes. 
    maxInteractionStrength (float): alsolute values of interaction coefficients are drawn uniformly between 0 and maxInteractionStrength
    pWithinMod (float): probability of an interaction between two species within a module (not including hub node)
    pPositive (float): Probability that an interaction is positive, given that there is an interaction
    pHub (float): probability of an interaction between two species within modules to and from the hub species

    Returns a matrix of interactions for one module.
    '''
    if (modType == "hub"):
        modSize = modSize - 1
        assert not pHub == None

    interactionDirection = np.random.choice([-1,0,1], size = (modSize, modSize), p = [(1-pPositive) * pWithinMod, 1-pWithinMod, pPositive * pWithinMod])
    interactionStrength = np.random.uniform(low = 0, high = maxInteractionStrength, size = (modSize, modSize))

    interactionMx = interactionDirection * interactionStrength

    if (modType == "hub"):
        # get the interactions to and from the hub node and append this to the matrix as the first row and column
        hubConnectionsDirection1 = np.random.choice([-1,0,1], size = (1, modSize), p = [(1-pPositive) * pHub, 1-pHub, pPositive * pHub])
        hubConnections1 = hubConnectionsDirection1 * np.random.uniform(low = 0, high = maxInteractionStrength, size = (1, modSize))
        hubConnectionsDirection2 = np.random.choice([-1,0,1], size = (modSize + 1, 1), p = [(1-pPositive) * pHub, 1-pHub, pPositive * pHub])
        hubConnections2= hubConnectionsDirection2 * np.random.uniform(low = 0, high = maxInteractionStrength, size = (modSize + 1, 1))

        interactionMx = np.append(hubConnections1, interactionMx, axis=0)
        interactionMx = np.append(hubConnections2, interactionMx, axis=1)

    np.fill_diagonal(interactionMx, 0)
    return interactionMx

def generateModularMatrix(maxInteractionStrength, pPositive, nSpecies, pWithinMod, pBetweenMod, nModules):
    '''
    Generates interaction matrix with modules that are more connected within than between

    maxInteractionStrength (float): alsolute values of interaction coefficients are drawn uniformly between 0 and maxInteractionStrength
    pPositive (float): Probability that an interaction is positive, given that there is an interaction
    nSpecies (int): number of species
    pWithinMod (float): probability of an interaction between two species within a module 
    bBetweenMod (float): probability of an interaction between two species between modules
    nModules (int): number of modules

    returns numpy array of interactions (diagonal is 0)
    
    Notes: probabilities are not precise because I assign all the interactions and then later zero out the diagonal.
    '''
    betweenModDirection = np.random.choice([-1,0,1], size = (nSpecies, nSpecies), p = [(1-pPositive) * pBetweenMod, 1-pBetweenMod, pPositive * pBetweenMod])
    modMatrix = betweenModDirection * np.random.uniform(low = 0, high = maxInteractionStrength, size = (nSpecies, nSpecies))
    
    modSize = nSpecies // nModules
    spIndex = 0
    for ii in range(nModules):
        if ii == (nModules - 1): # we are on the last module so this is all the remaining species
            modSize = nSpecies - spIndex 
        thisModule = generateModule(modSize = modSize, modType = "module", maxInteractionStrength = maxInteractionStrength, pWithinMod = pWithinMod, pPositive = pPositive)
       
        # insert module into matrix
        modMatrix[spIndex:spIndex+modSize, spIndex:spIndex+modSize] = thisModule
        spIndex = spIndex + modSize

    np.fill_diagonal(modMatrix, 0)
    return modMatrix

def generateHubMatrix(maxInteractionStrength, pPositive, nSpecies, pWithinMod, pBetweenMod, nModules, pHub):
    '''
    Generates interaction matrix with the first species being a hub node

    maxInteractionStrength (float): alsolute values of interaction coefficients are drawn uniformly between 0 and maxInteractionStrength
    pPositive (float): Probability that an interaction is positive, given that there is an interaction
    nSpecies (int): number of species
    pWithinMod (float): probability of an interaction between two species within a module (not including hub node)
    bBetweenMod (float): probability of an interaction between two species between modules
    nModules (int): number of modules
    pHub (float): probability of an interaction between two species within modules to and from the hub species

    returns numpy array of interactions (diagonal is 0)

    Notes: probabilities are not precise because I assign all the interactions and then later zero out the diagonal.
    '''
    betweenModDirection = np.random.choice([-1,0,1], size = (nSpecies, nSpecies), p = [(1-pPositive) * pBetweenMod, 1-pBetweenMod, pPositive * pBetweenMod])
    hubMatrix = betweenModDirection * np.random.uniform(low = 0, high = maxInteractionStrength, size = (nSpecies, nSpecies))
    
    modSize = nSpecies // nModules
    spIndex = 0
    for ii in range(nModules):
        if ii == (nModules - 1): # we are on the last module so this is all the remaining species
            modSize = nSpecies - spIndex 
        thisModule = generateModule(modSize = modSize, modType = "hub", maxInteractionStrength = maxInteractionStrength, pWithinMod = pWithinMod, pPositive = pPositive, pHub = pHub)
       
        # insert module into matrix
        hubMatrix[spIndex:spIndex+modSize, spIndex:spIndex+modSize] = thisModule
        spIndex = spIndex + modSize

    np.fill_diagonal(hubMatrix, 0)
    return hubMatrix

def getNestedMutual(maxInteractionStrength, nSpecies, pInterInteract = 1, pInterInteract_off = 0, pIntraInteract = 0, pPositive = 0):
    '''Make a matrix of interactions with nested inter-group interactions (mutualisms) and a pre-defined probability of intra-group interactions.
    
    maxInteractionStrength (float)
    nSpecies (int): number of species
    pInterInteract: Probability than an inter-group interaction (mutualism) exists. Default is 100% for perfectly nested.
    pInterInteract_off: Probability that an intergroup interaction (mutualism) exists outside the perfect nested structure. Default
    is 0% for perfectly nested. 
    pIntraInteract: Probability that an intra-group interaction exists.
    pPositive (float): Probability that an intra-group interaction is positive, given that there is an interaction
    

    Returns a numpy array of interactions (diagonal is zero)

    Note: Probabilities are not exact because the diagonal is zeroed at the end.
    '''
    interactionDirection = np.zeros((nSpecies,nSpecies))
    numgroup = int(0.5*len(interactionDirection))

    for i in range(len(interactionDirection)):
        for j in range(len(interactionDirection[1])):
            if i >= numgroup and j <=numgroup and i >= j + numgroup: 
                interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract,1-pInterInteract])
            elif i >= numgroup and j <=numgroup and i < j + numgroup:
                interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract_off,1-pInterInteract_off])
            elif i <=numgroup and j >=numgroup and j>= i+numgroup:
                interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract,1-pInterInteract])
            elif i <=numgroup and j >=numgroup and j< i+numgroup:
                interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract_off,1-pInterInteract_off])
            elif (i < numgroup and j < numgroup) or (i >= numgroup and j >=numgroup):
                interactionDirection[i,j] = np.random.choice([-1,0,1], p = [(1-pPositive) * pIntraInteract, 1-pIntraInteract, pPositive * pIntraInteract])

    interactionStrength = np.random.uniform(low = 0, high = maxInteractionStrength, size = (nSpecies, nSpecies))
    
    nestedmutualMx = interactionDirection * interactionStrength

    np.fill_diagonal(nestedmutualMx, 0)
    return nestedmutualMx

def getModularMutual(maxInteractionStrength, nSpecies, pInterInteract = 1, pInterInteract_off = 0, pIntraInteract = 0, pPositive = 0, modulefrac = 1):
    '''Make a matrix of interactions with nested inter-group interactions (mutualisms) and a pre-defined probability of intra-group interactions.
    
    maxInteractionStrength (float)
    nSpecies (int): number of species
    pInterInteract: Probability than an inter-group interaction (mutualism) exists. Default is 100% for perfectly modular.
    pInterInteract_off: Probability that an intergroup interaction (mutualism) exists outside the perfect nested structure. Default
    is 0% for perfectly modular. 
    pIntraInteract: Probability that an intra-group interaction exists.
    pPositive (float): Probability that an intra-group interaction is positive, given that there is an interaction.
    modulefrac (float): Percent of a group that is in one module (homogeneous module sizes). Default is 100%

    Returns a numpy array of interactions (diagonal is zero)

    Note: Probabilities are not exact because the diagonal is zeroed at the end.'''

    interactionDirection = np.zeros((nSpecies,nSpecies))
    numgroup = int(0.5*nSpecies)
    mod_num = int(modulefrac*numgroup)
    total_mod = int(np.ceil(numgroup/mod_num))

    for i in range(len(interactionDirection)):
        for j in range(len(interactionDirection[1])):
            if i >= numgroup and j <=numgroup: 
                for k in range(total_mod):
                    if i < numgroup+(k+1)*mod_num and i>= numgroup+(k)*mod_num and j < (k+1)*mod_num and j >=(k)*mod_num: #and i >= j + (k+1)*numgroup: 
                        interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract,1-pInterInteract])
                    elif i < numgroup+(k+1)*mod_num and i>= numgroup+(k)*mod_num and (j >= (k+1)*mod_num or j <(k)*mod_num): 
                        interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract_off,1-pInterInteract_off])
            elif i <=numgroup and j >=numgroup:
                for k in range(total_mod):
                    if i < (k+1)*mod_num and i>= (k)*mod_num and j < numgroup+(k+1)*mod_num and j >=numgroup+(k)*mod_num: #and i >= j + (k+1)*numgroup: 
                        interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract,1-pInterInteract])
                    elif i < (k+1)*mod_num and i>= (k)*mod_num and (j >= numgroup+(k+1)*mod_num or j < numgroup+(k)*mod_num):
                        interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract_off,1-pInterInteract_off])
            if (i < numgroup and j < numgroup) or (i >= numgroup and j >=numgroup):
                    interactionDirection[i,j] = np.random.choice([-1,0,1], p = [(1-pPositive) * pIntraInteract, 1-pIntraInteract, pPositive * pIntraInteract])

    interactionStrength = np.random.uniform(low = 0, high = maxInteractionStrength, size = (nSpecies, nSpecies))
    
    modular_mutualMx = interactionDirection * interactionStrength

    np.fill_diagonal(modular_mutualMx, 0)
    return modular_mutualMx

def getNestedParasite(maxInteractionStrength, nSpecies, pInterInteract = 1, pInterInteract_off = 0, pIntraInteract = 0, pPositive = 0):
    '''Make a matrix of interactions with nested inter-group interactions (parasitism) and a pre-defined probability of intra-group interactions.
    
    maxInteractionStrength (float)
    nSpecies (int): number of species
    pIntraInteract: Probability that an intra-group interaction exists. 
    pInterInteract: Probability than an inter-group interaction (parasitism/antagonism) exists. Default is 100% for perfectly nested.In this construction, all hosts are listed first, and then parasites. 
                    Therefore, the upper block of the interaction matrix has negative values, and the lower block has positive.
    pInterInteract_off: Probability that an intergroup interaction (parasitism/antagonism) exists outside the perfect nested structure. Default
    is 0% for perfectly nested. 
    pPositive (float): Probability that an intra-group interaction is positive, given that there is an interaction
    

    Returns a numpy array of interactions (diagonal is zero)

    Note: Probabilities are not exact because the diagonal is zeroed at the end.
    '''
    interactionDirection = np.zeros((nSpecies,nSpecies))
    numgroup = int(0.5*len(interactionDirection))

    for i in range(len(interactionDirection)):
        for j in range(len(interactionDirection[1])):
            if i >= numgroup and j <=numgroup and i >= j + numgroup: 
                interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract,1-pInterInteract])
            elif i >= numgroup and j <=numgroup and i < j + numgroup:
                interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract_off,1-pInterInteract_off])
            elif i <=numgroup and j >=numgroup and j>= i+numgroup:
                interactionDirection[i,j] = np.random.choice([-1,0], p = [pInterInteract,1-pInterInteract])
            elif i <=numgroup and j >=numgroup and j< i+numgroup:
                interactionDirection[i,j] = np.random.choice([-1,0], p = [pInterInteract_off,1-pInterInteract_off])
            elif (i < numgroup and j < numgroup) or (i >= numgroup and j >=numgroup):
                interactionDirection[i,j] = np.random.choice([-1,0,1], p = [(1-pPositive) * pIntraInteract, 1-pIntraInteract, pPositive * pIntraInteract])

    interactionStrength = np.random.uniform(low = 0, high = maxInteractionStrength, size = (nSpecies, nSpecies))
    
    nestedparasiteMx = interactionDirection * interactionStrength

    np.fill_diagonal(nestedparasiteMx, 0)
    return  nestedparasiteMx

def getModularParasite(maxInteractionStrength, nSpecies, pInterInteract = 1, pInterInteract_off = 0, pIntraInteract = 0, pPositive = 0, modulefrac = 1):
    '''Make a matrix of interactions with nested inter-group interactions (antagonisms) and a pre-defined probability of intra-group interactions.
    
    maxInteractionStrength (float)
    nSpecies (int): number of species
    pInterInteract: Probability than an inter-group interaction (antagonism) exists. Default is 100% for perfectly modular.
    pInterInteract_off: Probability that an intergroup interaction (antagonism) exists outside the perfect nested structure. Default
    is 0% for perfectly modular. 
    pIntraInteract: Probability that an intra-group interaction exists.
    pPositive (float): Probability that an intra-group interaction is positive, given that there is an interaction.
    modulefrac (float): Percent of a group that is in one module (homogeneous module sizes). Default is 100%

    Returns a numpy array of interactions (diagonal is zero)

    Note: Probabilities are not exact because the diagonal is zeroed at the end.'''

    interactionDirection = np.zeros((nSpecies,nSpecies))
    numgroup = int(0.5*nSpecies)
    mod_num = int(modulefrac*numgroup)
    total_mod = int(np.ceil(numgroup/mod_num))

    for i in range(len(interactionDirection)):
        for j in range(len(interactionDirection[1])):
            if i >= numgroup and j <=numgroup: 
                for k in range(total_mod):
                    if i < numgroup+(k+1)*mod_num and i>= numgroup+(k)*mod_num and j < (k+1)*mod_num and j >=(k)*mod_num: #and i >= j + (k+1)*numgroup: 
                        interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract,1-pInterInteract])
                    elif i < numgroup+(k+1)*mod_num and i>= numgroup+(k)*mod_num and (j >= (k+1)*mod_num or j <(k)*mod_num): 
                        interactionDirection[i,j] = np.random.choice([1,0], p = [pInterInteract_off,1-pInterInteract_off])
            elif i <=numgroup and j >=numgroup:
                for k in range(total_mod):
                    if i < (k+1)*mod_num and i>= (k)*mod_num and j < numgroup+(k+1)*mod_num and j >=numgroup+(k)*mod_num: #and i >= j + (k+1)*numgroup: 
                        interactionDirection[i,j] = np.random.choice([-1,0], p = [pInterInteract,1-pInterInteract])
                    elif i < (k+1)*mod_num and i>= (k)*mod_num and (j >= numgroup+(k+1)*mod_num or j < numgroup+(k)*mod_num):
                        interactionDirection[i,j] = np.random.choice([-1,0], p = [pInterInteract_off,1-pInterInteract_off])
            if (i < numgroup and j < numgroup) or (i >= numgroup and j >=numgroup):
                    interactionDirection[i,j] = np.random.choice([-1,0,1], p = [(1-pPositive) * pIntraInteract, 1-pIntraInteract, pPositive * pIntraInteract])

    interactionStrength = np.random.uniform(low = 0, high = maxInteractionStrength, size = (nSpecies, nSpecies))
    
    modular_parasiteMx = interactionDirection * interactionStrength

    np.fill_diagonal(modular_parasiteMx, 0)
    return modular_parasiteMx


def hubMxMatching(nSpec, nModules, pInteractRandom):
    '''
    Based on the number of species, modules, and interaction probability in a random matrix, 
    calculate the parameters for a mudular and hub matrix to have the same expected number of edges
    '''
    # calculating expected number of links in each network
    totalLinks = nSpec**2 - nSpec
    modSize = nSpec // nModules
    maxLinksWithin = nModules * (modSize**2 - modSize)

    E_edges_rand = pInteractRandom * totalLinks

    if E_edges_rand <= maxLinksWithin:
        pWithinMod = E_edges_rand / maxLinksWithin
        pBetweenMod = 0
    else:
        pWithinMod = 1
        pBetweenMod = (E_edges_rand - maxLinksWithin) / (totalLinks - maxLinksWithin)

    maxHubLinks = nModules * (modSize - 1)
    if E_edges_rand <= maxHubLinks:
        pHub = E_edges_rand / maxHubLinks
    else: 
        pHub = 1
        pWithinModHub = (E_edges_rand - nModules * (modSize - 1)) / (totalLinks - nModules * (modSize - 1))
        pBetweenModHub = (E_edges_rand - nModules * (modSize - 1)) / (totalLinks - nModules * (modSize - 1))

    return {
        'pWithinMod': pWithinMod,
        'pBetweenMod': pBetweenMod,
        'pWithinModHub': pWithinModHub,
        'pBetweenModHub': pBetweenModHub,
        'pHub': pHub
    }

def SOMmodel(nestedMatrix):

    focal = np.random.choice(len(nestedMatrix))
    binary_nested = (nestedMatrix != 0).astype(int)
    if (any(binary_nested[focal]==1)):
        old = np.random.choice(np.where(binary_nested[focal]==1)[0]) #where focal species has interactions
        newinteracts = np.where(binary_nested[focal]==0)[0] #where focal species lacks interactions
        newchoice = [num for num in newinteracts if num != focal] #removing self interaction

        if len(newchoice) > 0: #making sure there is a new interaction to form
            new = np.random.choice(newchoice)
            if (sum(binary_nested[new]) > sum(binary_nested[old])): #only connect to higher degree
                    nestedMatrix[focal,new] = nestedMatrix[focal,old]
                    nestedMatrix[focal,old] = 0

    return nestedMatrix


def generateNestedMatrix(randMx, SOMiterations):
    '''
    Generates a semi-nested interaction matrix, depending on the number of (S)elf (O)rganizing (M)odel iterations
    randMx (numpy array): Previously generated random matrix 
    SOMiterations (float): how many times to run the rewiring prescribed the SOM network model
    '''
    nestedMatrix = randMx
    for ii in range(SOMiterations):
        nestedMatrix = SOMmodel(nestedMatrix)

    return nestedMatrix


def getPowerOrExpMatrix(nSpec, interactStren, type, expected_degree, pPositive = 0.5):
    ''' type in ["random", "exponential", "power"]
    '''
    assert type in ["random", "exponential", "power"]

    if type == "exponential":
        # 1. Exponential degree distribution — build undirected then randomly orient each edge
        degrees = np.random.exponential(scale=expected_degree, size=nSpec)
        degree_sequence = np.maximum(1, np.round(degrees).astype(int))
        if sum(degree_sequence) % 2 == 1: # sum has to be even
            degree_sequence[0] += 1
        G_undirected = nx.configuration_model(degree_sequence)
        G_undirected = nx.Graph(G_undirected)
        G_undirected.remove_edges_from(nx.selfloop_edges(G_undirected))
        G = nx.DiGraph()
        G.add_nodes_from(G_undirected.nodes())
        for u, v in G_undirected.edges():
            if np.random.random() < 0.5:
                G.add_edge(u, v)
            else:
                G.add_edge(v, u)
    elif type == "power":
        # Create scale-free network using Barabási-Albert model, then randomly orient each edge
        G_undirected = nx.barabasi_albert_graph(nSpec, expected_degree // 2)
        G = nx.DiGraph()
        G.add_nodes_from(G_undirected.nodes())
        for u, v in G_undirected.edges():
            if np.random.random() < 0.5:
                G.add_edge(u, v)
            else:
                G.add_edge(v, u)
    else:
        nEdge = (nSpec * expected_degree) // 2
        # Divide this by 2 if you want to get p_ij from the paper because of it being directed
        edgeProb = nEdge / ((nSpec**2 - nSpec) / 2)
        G_undirected = nx.gnp_random_graph(nSpec, edgeProb, directed=False)
        G = nx.DiGraph()
        G.add_nodes_from(G_undirected.nodes())
        for u, v in G_undirected.edges():
            if np.random.random() < 0.5:
                G.add_edge(u, v)
            else:
                G.add_edge(v, u)

    mx = nx.adjacency_matrix(G).toarray() * np.random.uniform(low = 0, high = interactStren, size = (nSpec, nSpec)) * np.random.choice(a = [-1,1], size = (nSpec, nSpec), replace = True, p = [1-pPositive, pPositive])

    np.fill_diagonal(mx, 0)

    return G_undirected, mx

def generateInteractionTypeMx(maxInteractionStrength, nSpecies, interactionTypeFreq):
    '''
    Random matrix with controlled frequency of types
    maxInteractionStrength (float) 
    interactionTypeFreq (numpy array of shape (6)): Frequency of types of interactions (+/+, -/-, +/-, +/0, -/0, 0/0)
    '''
    assert interactionTypeFreq.shape == (6,)

    randomMatrix = np.zeros(shape = (nSpecies, nSpecies))

    # choose interaction types based on specified frequencies and then set them accordingly
    for i in range(nSpecies):
        for j in range(i + 1, nSpecies):
            interactType = np.random.choice(["+/+", "-/-", "+/-", "+/0", "-/0", "0/0"], size = 1, p = interactionTypeFreq)
            randNums = np.random.uniform(low = 0, high = maxInteractionStrength, size = 2)
            if interactType == "+/+":
                randomMatrix[i, j] = randNums[0]
                randomMatrix[j, i] = randNums[1]
            elif interactType == "-/-":
                randomMatrix[i, j] = -randNums[0]
                randomMatrix[j, i] = -randNums[1]
            elif interactType == "+/-":
                randomMatrix[i, j] = randNums[0]
                randomMatrix[j, i] = -randNums[1]
            elif interactType == "+/0":
                randomMatrix[i, j] = randNums[0]
            elif interactType == "-/0":
                randomMatrix[i, j] = -randNums[0]
    return randomMatrix

def interactionTypeFromMx(interactMx):
    '''
    interactMx (numpy array) 
    interactionTypeFreq (numpy array of shape (6)): Frequency of types of interactions (+/+, -/-, +/-, +/0, -/0, 0/0)
    '''
    interactionTypeFreq = np.zeros(shape = (6,))
    nSpecies = interactMx.shape[0]

    # count types of interactions
    for i in range(nSpecies):
        for j in range(i + 1, nSpecies):
            if interactMx[i, j] > 0 and interactMx[j, i] > 0: # interactType == "+/+"
                interactionTypeFreq[0] += 1
            elif interactMx[i, j] < 0 and interactMx[j, i] < 0: # interactType == "-/-"
                interactionTypeFreq[1] += 1
            elif (interactMx[i, j] < 0 and interactMx[j, i] > 0) or (interactMx[i, j] > 0 and interactMx[j, i] < 0): #interactType == "+/-"
                interactionTypeFreq[2] += 1
            elif (interactMx[i, j] > 0 and interactMx[j, i] == 0) or (interactMx[i, j] == 0 and interactMx[j, i] > 0): # interactType == "+/0"
                interactionTypeFreq[3] += 1
            elif (interactMx[i, j] < 0 and interactMx[j, i] == 0) or (interactMx[i, j] == 0 and interactMx[j, i] < 0): # interactType == "-/0"
                interactionTypeFreq[4] += 1
            else:
                interactionTypeFreq[5] += 1

    interactionTypeFreq = interactionTypeFreq / ((nSpecies**2 - nSpecies) / 2)
    return interactionTypeFreq