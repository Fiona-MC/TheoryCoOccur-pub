# this code makes CommunityDynamics_CoOccurrence.ipynb into some functions to run deterministic GLV sim or stochastic sim with Gaussian white noise
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import metaCommunityMx

def one_comm_dynamics(t, N, a, b, alpha):
    dNdt = N*(b*(1 - alpha*N) + (a@N))
    return dNdt
    
def one_comm_dynamics_stoch(t, N, a, b, alpha, noise):
    noise_t = noise*np.random.normal(0,1,np.size(N))
    dNdt = N*(b*(1 - alpha*N) + (a@N))+noise_t
    return dNdt
    
def run_communityDynamics(a, b, alpha, tmax = 30, y0 = None, max_step=None):
    '''
    a is the community matrix
    b is the intrinsic growth rates of the species (growth at low abundance) in the absence of other species
    alpha is the self-limitation parameter for all species
    tmax is the ending time for the differential equation
    y0 is chosen at random if not supplied
    '''
    nSpecies = a.shape[0]
    
    #random vector of birthrates with noise
    # b = np.ones(nSpecies)-np.random.uniform(low = 0, high = 1, size = (nSpecies))
    tspan = [0, tmax]

    if y0 is None:
        #random initial vector of starting densities
        y0 = np.random.uniform(low = 10, high = 20, size = (nSpecies))

    def myFun(t, N):
        return one_comm_dynamics(t, N, a, b, alpha)
    
    if not max_step is None:
        sol = solve_ivp(myFun, tspan, y0, max_step=max_step)
    else:
        sol = solve_ivp(myFun, tspan, y0)

    abds_end = sol.y[:, -1]
    abds_end_m1 = sol.y[:, -2]

    # check that the numbers are not changing too much at the end
    if not all(np.round(np.absolute(abds_end - abds_end_m1), 0) < 0.01 * np.absolute(abds_end)):
        print("System did not converge. Trying with max_step=1.0")
        sol = solve_ivp(myFun, tspan, y0, max_step = 1.0)
        abds_end = sol.y[:, -1]
        abds_end_m1 = sol.y[:, -2]
        if not all(np.round(np.absolute(abds_end - abds_end_m1), 0) < 0.01 * np.absolute(abds_end)):
            print("Did not reach equilibrium; run with higher tmax; System may not converge if interactions cause positive feedback loop.")
            print(np.max(np.round(np.absolute(abds_end - abds_end_m1), 0)))
        if not all(abds_end < 10**8):
            print("Max pop size over 10^8 -- system did not converge")
            print(np.max(abds_end))
    
    return sol
    
def run_communityDynamicsStoch(a, b, alpha, tmax = 30, y0 = None, noise = None):
    '''
    a is the community matrix
    b is the intrinsic growth rates of the species (growth at low abundance) in the absence of other species
    alpha is the self-limitation parameter for all species
    tmax is the ending time for the differential equation
    y0 is chosen at random if not supplied
    noise is the multiplier for Gaussian white noise applied to each population at each time step
    '''
    nSpecies = a.shape[0]
    
    #random vector of birthrates with noise
    # b = np.ones(nSpecies)-np.random.uniform(low = 0, high = 1, size = (nSpecies))
    tspan = [0, tmax]

    if y0 is None:
        #random initial vector of starting densities
        y0 = np.random.uniform(low = 10, high = 20, size = (nSpecies))

    def myFunStoch(t,N):
        return one_comm_dynamics_stoch(t, N, a, b, alpha, noise)
    
    sol = solve_ivp(myFunStoch, tspan, y0)
    abds_end = sol.y[:, -1]
    abds_end_m1 = sol.y[:, -2]
    print('done')
    # check that the numbers are not changing too much at the end
    if not all(np.round(np.absolute(abds_end - abds_end_m1), 0) < 0.01 * np.absolute(abds_end)):
        print("Did not reach equilibrium; run with higher tmax")
        print(np.max(np.round(np.absolute(abds_end - abds_end_m1), 0)))
    
    return sol

def simulate_coOccurrence(samplesize, numtrials, a, b, alpha, tmax = 30, y0 = None, noise = None, max_step=None):
    '''
    n is the number of samples to run
    samplesize is the number of species to sample from the community
    numtrials is the number of samples to get for co-occurrence dynamics
    a is the community matrix
    b is the birth rates (randomly chosen if not supplied)
    alpha is the is the self-limitation parameter for all species
    tmax is the ending time for the differential equation
    y0 is chosen at random if not supplied
    If noise is not provided, deterministic GLV will run; if provided, SDE will run
    '''
    nSpecies = a.shape[0]

    summarytable = np.zeros((numtrials,nSpecies))
    indexlist = np.zeros((numtrials,samplesize))

    if b is None:
        # random birthrates with noise
        b = np.ones(nSpecies)-np.random.uniform(low = 0, high = .1, size = (nSpecies))

    if noise is None:
        for i in range(numtrials):
            #get list of samplesize random integers from 0 to nSpecies (save list)
            spList = np.random.choice(range(0, nSpecies), size = (samplesize), replace = False)
            spList.sort()
            indexlist[i,:] = spList

            #take rows and columns of community interaction matrix that correspond to those integers for reduced interaction 
            #matrix
            reducedMx = a[np.ix_(spList, spList)]
            #also subsample birthrates 
            reducedb = b[spList]

            #figure out if subsampling alpha is necessary 
            if hasattr(alpha, "__len__"):
                reducedalpha = alpha[spList]
                #run ode solver for reduced interaction matrix
                sol = run_communityDynamics(a = reducedMx, b = reducedb, alpha = reducedalpha, tmax = tmax, y0 = None, max_step=max_step)
            else: 
                #run ode solver for reduced interaction matrix
                sol = run_communityDynamics(a = reducedMx, b = reducedb, alpha = alpha, tmax = tmax, y0 = None, max_step=max_step)
           
            #output: table with numtrial as rows and species populations as columns
            for x in range(samplesize):
                summarytable[i,spList[x]] = sol.y[x,-1] 

        return summarytable, indexlist
    else:
        for i in range(numtrials):
            #get list of samplesize random integers from 0 to nSpecies (save list)
            spList = np.random.choice(range(0, nSpecies), size = (samplesize), replace = False)
            spList.sort()
            indexlist[i,:] = spList

            #take rows and columns of community interaction matrix that correspond to those integers for reduced interaction 
            #matrix
            reducedMx = a[np.ix_(spList, spList)]
            #also subsample birthrates 
            reducedb = b[spList]

            if hasattr(alpha, "__len__"):
                reducedalpha = alpha[spList]
                sol = run_communityDynamicsStoch(a = reducedMx, b = reducedb, alpha = reducedalpha, tmax = tmax, y0 = None, noise = 
                                             noise)
            else:
            #run ode solver for reduced interaction matrix
                sol = run_communityDynamicsStoch(a = reducedMx, b = reducedb, alpha = alpha, tmax = tmax, y0 = None, noise = 
                                             noise)

            #output: table with numtrial as rows and species populations as columns
            for x in range(samplesize):
                summarytable[i,spList[x]] = sol.y[x,-1] 

        return summarytable, indexlist


def simulate_coOccurrence_traitmatch(samplesize, numtrials, a, b, alpha, traits, tmax = 30, y0 = None, envStrength = 0.5, environ = None):
    '''
    n is the number of samples to run
    samplesize is the number of species to sample from the community
    numtrials is the number of samples to get for co-occurrence dynamics
    a is the community matrix
    b is the initial birth rate vector (randomly chosen if not supplied)
    alpha is the is the self-limitation parameter for all species
    traits is a list of the environmental trait value optimum for all species (between 0 and 1)
    tmax is the ending time for the differential equation
    y0 is chosen at random if not supplied
    noise is the amount that the environmetn affects the species
    environ is a number or vector of numbers that are provided as the environment the species will live in
    If noise is not provided, standard GLV run; if provided, birth rates will be perturbed each trial according to trait match and noise
    '''
    nSpecies = a.shape[0]

    summarytable = np.zeros((numtrials,nSpecies))
    indexlist = np.zeros((numtrials,samplesize))

    if b is None:
        # random birthrates with noise
        b = np.ones(nSpecies)-np.random.uniform(low = 0, high = .1, size = (nSpecies))

    if envStrength is None:
        for i in range(numtrials):
            #get list of samplesize random integers from 0 to nSpecies (save list)
            spList = np.random.choice(range(0, nSpecies), size = (samplesize), replace = False)
            spList.sort()
            indexlist[i,:] = spList

            #take rows and columns of community interaction matrix that correspond to those integers for reduced interaction 
            #matrix
            reducedMx = a[np.ix_(spList, spList)]
            #also subsample birthrates 
            reducedb = b[spList]

            #run ode solver for reduced interaction matrix
            sol = run_communityDynamics(a = reducedMx, b = reducedb, alpha = alpha, tmax = tmax, y0 = None)

            #output: table with numtrial as rows and species populations as columns
            for x in range(samplesize):
                summarytable[i,spList[x]] = sol.y[x,-1] 

        return summarytable, indexlist

    else:
        for i in range(numtrials):
            #get list of samplesize random integers from 0 to nSpecies (save list)
            spList = np.random.choice(range(0, nSpecies), size = (samplesize), replace = False)
            spList.sort()
            indexlist[i,:] = spList

            #take rows and columns of community interaction matrix that correspond to those integers for reduced interaction 
            #matrix
            reducedMx = a[np.ix_(spList, spList)]
            #subsample traits
            reducedTraits = traits[spList]
           
            #pick the environment value for this sample
            if environ is None:
                print('environ was none')
                environ = np.random.uniform(low = 0, high = 1)

            base_k = 1/alpha # carrying capacity without any interactions/environment mismatch

            if type(environ) in (int, float): # this is for backwards compatibility
                environ = [environ]
            else: # environ is a vector of numbers
                assert type(environ) == list
            
            env_i = np.random.choice(environ)
            
            mismatch = abs(reducedTraits-env_i)

            #also subsample birthrates
            reducedb = b[spList]
            
            reduction_percent = envStrength * mismatch
            assert ((reduction_percent >= 0) & (reduction_percent <= 1)).all()
            k = base_k * (1 - reduction_percent) # new carrying capacities
            k = np.maximum(1, k) # if k is less than 1, round up to 1
            alpha_vect = 1/k

            #run ode solver for reduced interaction matrix
            sol = run_communityDynamics(a = reducedMx, b = reducedb, alpha = alpha_vect, tmax = tmax, y0 = None)

            #output: table with numtrial as rows and species populations as columns
            for x in range(samplesize):
                summarytable[i,spList[x]] = sol.y[x,-1] 
        
        return summarytable, indexlist