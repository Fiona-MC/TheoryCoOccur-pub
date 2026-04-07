# write degree sequence csvs for the galiana analysis in r
# locally conda activate networks
import csv
import numpy as np
from importlib import reload
from sklearn import metrics

import metaCommunityMx
import run_cooccurrence
reload(metaCommunityMx)
reload(run_cooccurrence)

import sys

# interactFactor = 1/10
# pPositive = 0.001
# nSpec = 1000

# Convert arguments to appropriate types
interactFactor = float(sys.argv[1])
pPositive = float(sys.argv[2])
nSpec = int(sys.argv[3])

save_dir = sys.argv[4] # "./degree_sequences_np/"

# parms for GLV sim
alpha = 0.001 
interactStren = alpha * interactFactor
#samplesize = int(0.5 * nSpec) # number of species to sample
samplesize = nSpec // 2 # number of species to sample
numtrials = 500
b = np.ones(nSpec) - np.random.uniform(low = 0, high = .1, size = (nSpec))
cutoffalpha = 0.01

for i in range(50):
    G_exp, mx_exp = metaCommunityMx.getPowerOrExpMatrix(nSpec, interactStren, type = "exponential", expected_degree = 10, pPositive = pPositive)
    G_power, mx_power = metaCommunityMx.getPowerOrExpMatrix(nSpec, interactStren, type = "power", expected_degree = 10, pPositive = pPositive)
    G_random, mx_random = metaCommunityMx.getPowerOrExpMatrix(nSpec, interactStren, type = "random", expected_degree = 10, pPositive = pPositive)

    ## write exponential node degrees
    degreeSeq = [d for n,d in G_exp.degree]
    spec = ["S" + str(ii) for ii in range(len(degreeSeq))]
    degreeSeqTable = [["species", "interactions"]] + list(zip(spec, degreeSeq))

    filename = save_dir + "inter/inter_nodeDegree" + str(nSpec) + "_Fmax" + str(interactFactor) + "_pPos" + str(pPositive) + "_exp_" + str(i) + ".csv"
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(degreeSeqTable)
    
    ## write power law node degrees
    degreeSeq = [d for n,d in G_power.degree]
    spec = ["S" + str(ii) for ii in range(len(degreeSeq))]
    degreeSeqTable = [["species", "interactions"]] + list(zip(spec, degreeSeq))

    filename = save_dir + "inter/inter_nodeDegree" + str(nSpec) + "_Fmax" + str(interactFactor) + "_pPos" + str(pPositive) + "_power_" + str(i) + ".csv"
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(degreeSeqTable)

    ## write random graph node degrees
    degreeSeq = [d for n,d in G_random.degree]
    spec = ["S" + str(ii) for ii in range(len(degreeSeq))]
    degreeSeqTable = [["species", "interactions"]] + list(zip(spec, degreeSeq))

    filename = save_dir + "inter/inter_nodeDegree" + str(nSpec) + "_Fmax" + str(interactFactor) + "_pPos" + str(pPositive) + "_rand_" + str(i) + ".csv"
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(degreeSeqTable)

    print("exponential max degree")
    print(max(sum(mx_exp != 0)))
    network_exp, G_coocc_exp = run_cooccurrence.get_coocc_network(interactionMx = mx_exp, cutoff = None, alpha = alpha, spsamplesize = samplesize, numtrials = numtrials, b = b, tmax=200, pval_alpha=cutoffalpha)
    print("power max degree")
    print(max(sum(mx_power != 0)))
    network_pow, G_coocc_pow = run_cooccurrence.get_coocc_network(interactionMx = mx_power, cutoff = None, alpha = alpha, spsamplesize = samplesize, numtrials = numtrials, b = b, tmax=200, pval_alpha=cutoffalpha)
    print("random max degree")
    print(max(sum(mx_random != 0)))
    network_rand, G_coocc_rand = run_cooccurrence.get_coocc_network(interactionMx = mx_random, cutoff = None, alpha = alpha, spsamplesize = samplesize, numtrials = numtrials, b = b, tmax=200, pval_alpha=cutoffalpha)
    
    ## write exponential coocurrence node degrees
    degreeSeq = [d for n,d in G_coocc_exp.degree]
    spec = ["S" + str(ii) for ii in range(len(degreeSeq))]

    degreeSeqTable = [["species", "interactions"]] + list(zip(spec, degreeSeq))

    filename = save_dir + "cooc/cooc_nodeDegree" + str(nSpec) + "_Fmax" + str(interactFactor) + "_pPos" + str(pPositive) + "_exp_cutoffalpha" + str(cutoffalpha) + "_" + str(i) + ".csv"
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(degreeSeqTable)

    ## write power law coocurrence node degrees
    degreeSeq = [d for n,d in G_coocc_pow.degree]
    spec = ["S" + str(ii) for ii in range(len(degreeSeq))]

    degreeSeqTable = [["species", "interactions"]] + list(zip(spec, degreeSeq))

    filename = save_dir + "cooc/cooc_nodeDegree" + str(nSpec) + "_Fmax" + str(interactFactor) + "_pPos" + str(pPositive) + "_power_cutoffalpha" + str(cutoffalpha) + "_" + str(i) + ".csv"
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(degreeSeqTable)
    
    ## write random graph coocurrence node degrees
    degreeSeq = [d for n,d in G_coocc_rand.degree]
    spec = ["S" + str(ii) for ii in range(len(degreeSeq))]

    degreeSeqTable = [["species", "interactions"]] + list(zip(spec, degreeSeq))

    filename = save_dir + "cooc/cooc_nodeDegree" + str(nSpec) + "_Fmax" + str(interactFactor) + "_pPos" + str(pPositive) + "_rand_cutoffalpha" + str(cutoffalpha) + "_" + str(i) + ".csv"
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(degreeSeqTable)