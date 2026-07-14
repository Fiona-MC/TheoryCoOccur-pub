import networkx as nx
import numpy as np
import community_sim
import matplotlib.pyplot as plt

def get_null_cutoff_permute(cooc_mx, pval_alpha = 0.01, shuffle_all = True, n_shuffle = 1000):
        ''' get co-occurrence pearson correlation cutoff by permuting abundances
                if shuffle_all, then flatten the matrix, permute and then reshape
                otherwise only shuffle within each species 
                then compute pairwise pearson correlations and compute the 95% confidence interval
                Then average over the n_shuffle confidence intervals
        '''
        lower_cutoffs = []
        upper_cutoffs = []

        for i in range(n_shuffle):
                if shuffle_all:
                        flat = cooc_mx.flatten()
                        np.random.shuffle(flat)
                        shuffled = flat.reshape(cooc_mx.shape)
                else:
                        shuffled = cooc_mx.copy()
                        for col in range(shuffled.shape[1]):
                                np.random.shuffle(shuffled[:, col])

                pearson_corr = np.corrcoef(shuffled, rowvar=False)
                null_flat = pearson_corr[np.triu_indices_from(pearson_corr, k=1)]

                lower_cutoffs.append(np.quantile(null_flat, pval_alpha / 2))
                upper_cutoffs.append(np.quantile(null_flat, 1 - (pval_alpha / 2)))

        lower_cutoff = np.mean(lower_cutoffs)
        upper_cutoff = np.mean(upper_cutoffs)

        return (lower_cutoff, upper_cutoff)

def get_network_pairwise_pseudoP(cooc_mx, pval_alpha = 0.01, n_shuffle = 1000):
        ''' convert co-occurrence matrix to binary network by permuting each pair
                and getting the correlation n_shuffle times and then computing the
                confidence interval (pval_alpha is the pvalue percentile cutoff) separately for each pair.
                Return three matrices, one with the upper cutoffs for each species pair,
                one with the lower cutoffs, and one with 0's 1's and -1's for the significant
                values from cooc_mx
        '''
        nSpec = cooc_mx.shape[1]
        upper_cutoff_mx = np.zeros((nSpec, nSpec))
        lower_cutoff_mx = np.zeros((nSpec, nSpec))

        actual_corr = np.corrcoef(cooc_mx, rowvar=False)

        for i in range(nSpec):
                for j in range(i + 1, nSpec):
                        null_corrs = np.empty(n_shuffle)
                        col_j = cooc_mx[:, j]
                        for k in range(n_shuffle):
                                shuffled_i = cooc_mx[:, i].copy()
                                np.random.shuffle(shuffled_i)
                                null_corrs[k] = np.corrcoef(shuffled_i, col_j)[0, 1]

                        lower = np.quantile(null_corrs, pval_alpha / 2)
                        upper = np.quantile(null_corrs, 1 - (pval_alpha / 2))

                        upper_cutoff_mx[i, j] = upper
                        upper_cutoff_mx[j, i] = upper
                        lower_cutoff_mx[i, j] = lower
                        lower_cutoff_mx[j, i] = lower

        sig_mx = (np.greater(actual_corr, upper_cutoff_mx).astype(int)
                  - np.less(actual_corr, lower_cutoff_mx).astype(int))
        np.fill_diagonal(sig_mx, 0)

        return lower_cutoff_mx, upper_cutoff_mx, sig_mx

def get_null_cutoff(nSpec, alpha, spsamplesize, numtrials, b, pval_alpha = 0.01, tmax = 100, compositional = False):
        ''' run the cooccurrence network pipeline without interactions and return the two tailed pair of cutoffs for significance
        '''
        nullInteractionMx = np.zeros(shape = (nSpec, nSpec))
        sim_res = community_sim.simulate_coOccurrence(samplesize = spsamplesize, numtrials = numtrials, a = nullInteractionMx, b = b, alpha = alpha, tmax = tmax, y0 = None)
        coAbd = np.round(sim_res[0], 0)
        if compositional:
                row_sums = coAbd.sum(axis=1, keepdims=True)
                coAbd = np.where(row_sums > 0, coAbd / row_sums, 0.0)

        pearson_corr = np.corrcoef(coAbd, rowvar = False)

        # cutoff based on the covariance of the species with no interactions
        pearson_null_flat = pearson_corr[np.triu_indices_from(pearson_corr, k = 1)]

        # compute equal tailed confidence interval
        lower_cutoff = np.quantile(pearson_null_flat, pval_alpha / 2)
        upper_cutoff = np.quantile(pearson_null_flat, 1 - (pval_alpha / 2))

        return (lower_cutoff, upper_cutoff)

def get_null_cutoff_traitmatch(nSpec, alpha, spsamplesize, numtrials, b, traits, pval_alpha = 0.01, tmax = 100, envStrength = 0.5):
        ''' run the cooccurrence network pipeline without interactions but WITH trait matching,
        and return the two tailed pair of cutoffs for significance.
        This accounts for environmental covariance so that only correlations beyond
        what is expected from shared environmental response become edges.
        '''
        nullInteractionMx = np.zeros(shape = (nSpec, nSpec))
        sim_res = community_sim.simulate_coOccurrence_traitmatch(samplesize = spsamplesize, numtrials = numtrials, a = nullInteractionMx, b = b, alpha = alpha, traits = traits, tmax = tmax, y0 = None, envStrength = envStrength, environ = None)
        coAbd = np.round(sim_res[0], 0)

        pearson_corr = np.corrcoef(coAbd, rowvar = False)

        # cutoff based on the covariance of the species with no interactions but with trait matching
        pearson_null_flat = pearson_corr[np.triu_indices_from(pearson_corr, k = 1)]

        # compute equal tailed confidence interval
        lower_cutoff = np.quantile(pearson_null_flat, pval_alpha / 2)
        upper_cutoff = np.quantile(pearson_null_flat, 1 - (pval_alpha / 2))

        return (lower_cutoff, upper_cutoff)

def get_coocc_network(interactionMx, cutoff, alpha, spsamplesize, numtrials, b, tmax = 100, plotNetworks = False, pval_alpha = 0.01, pos = None, max_step=None, permuatation_pseudoP_n = None, permutation_cutoff_n = None, shuffle_all = True, return_pearson = False, compositional = False):
        ''' run the cooccurrence network pipeline
        '''
        sim_res = community_sim.simulate_coOccurrence(samplesize = spsamplesize, numtrials = numtrials, a = interactionMx, b = b, alpha = alpha, tmax = tmax, y0 = None, max_step=max_step)
        coAbd = np.round(sim_res[0], 0)
        if compositional:
                row_sums = coAbd.sum(axis=1, keepdims=True)
                coAbd = np.where(row_sums > 0, coAbd / row_sums, 0.0)

        # Diagnostic: check for zero variance and identical vectors -- note this sometimes happens when there are too many negative interactions
        # if this happens, we will have some NaN values in the correlation matrix.
        variances = np.var(coAbd, axis=0)
        zero_var_count = np.sum(variances == 0)
        if zero_var_count > 0:
                print(f"WARNING: {zero_var_count} species have zero variance in run_cooccurrence.get_coocc_network()")
        
        # Check for identical columns (species with identical abundance patterns)
        nSpec = len(interactionMx[1,])
        identical_count = 0
        for i in range(nSpec):
                for j in range(i+1, nSpec):
                        if np.array_equal(coAbd[:, i], coAbd[:, j]):
                                identical_count += 1
        if identical_count > 0:
                print(f"WARNING: {identical_count} pairs of species have identical abundance vectors in run_cooccurrence.get_coocc_network()")
        
        pearson_corr = np.corrcoef(coAbd, rowvar = False)

        if permuatation_pseudoP_n is not None:
                _, _, network_coocc = get_network_pairwise_pseudoP(cooc_mx=pearson_corr, pval_alpha=pval_alpha, n_shuffle=permuatation_pseudoP_n)
        else:
                if cutoff == None:
                        if permutation_cutoff_n is not None:
                                lower_cutoff, upper_cutoff = get_null_cutoff_permute(cooc_mx = pearson_corr, pval_alpha = pval_alpha, shuffle_all = shuffle_all, n_shuffle = permutation_cutoff_n)
                        else:
                                lower_cutoff, upper_cutoff = get_null_cutoff(nSpec = interactionMx.shape[0], alpha = alpha, spsamplesize = spsamplesize, numtrials = numtrials, b = b, pval_alpha = pval_alpha, tmax = tmax, compositional = compositional)
                elif type(cutoff) == float:
                        lower_cutoff = -cutoff
                        upper_cutoff = cutoff
                elif len(cutoff) == 2:
                        lower_cutoff = cutoff[0]
                        upper_cutoff = cutoff[1]

                assert(lower_cutoff < upper_cutoff)
                assert(lower_cutoff < 0)
                assert(upper_cutoff > 0)

                network_coocc = np.greater(pearson_corr, upper_cutoff).astype(int) - np.less(pearson_corr, lower_cutoff).astype(int)

        np.fill_diagonal(network_coocc, 0)

        G_coocc = nx.from_numpy_array(network_coocc)
        
        if plotNetworks:
                G_inter = nx.from_numpy_array(interactionMx != 0)

                normalized_interMx = interactionMx / (interactionMx.max() - interactionMx.min())
                normalized_pearson_corr = pearson_corr / (pearson_corr.max() - pearson_corr.min())

                if pos == None:
                        pos = nx.spring_layout(G_inter, seed=42)

                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

                # Plot interaction network
                nx.draw(G_inter, pos, ax=ax1, with_labels=False, node_color='skyblue',
                        node_size=8, edge_color=['red' if normalized_interMx[u,v] < 0 else 'blue' for u,v in G_inter.edges()],
                        width=[abs(normalized_interMx[u,v]) * 3 for u,v in G_inter.edges()])
                ax1.set_title("Interaction")

                # Plot co-occurrence network
                nx.draw(G_coocc, pos, ax=ax2, with_labels=False, node_color='skyblue',
                        node_size=8, edge_color=['red' if normalized_pearson_corr[u,v] < 0 else 'blue' for u,v in G_coocc.edges()],
                        width=[abs(normalized_pearson_corr[u,v]) * 3 for u,v in G_coocc.edges()])
                ax2.set_title("Co-occurrence")

                plt.show()

        if return_pearson:
                return network_coocc, G_coocc, pearson_corr
        else:
                return network_coocc, G_coocc

# def get_coocc_network_randenv(interactionMx, cutoff, alpha, spsamplesize, numtrials, b, tmax = 150, plotNetworks = False, noise = None):
#         ''' run the cooccurrence network pipeline for random environment noise
#         '''
#         if cutoff == None:
#                 lower_cutoff, upper_cutoff = get_null_cutoff(nSpec = interactionMx.shape[0], alpha = alpha, spsamplesize = spsamplesize, numtrials = numtrials, b = b, pval_alpha = 0.01, tmax = tmax)
#         elif type(cutoff) == float:
#                 lower_cutoff = -cutoff
#                 upper_cutoff = cutoff
#         elif len(cutoff) == 2:
#                 lower_cutoff = cutoff[0]
#                 upper_cutoff = cutoff[1]

#         assert(lower_cutoff < upper_cutoff)
#         assert(lower_cutoff < 0)
#         assert(upper_cutoff > 0)

#         sim_res = community_sim.simulate_coOccurrence_randenv(samplesize = spsamplesize, numtrials = numtrials, a = interactionMx, b = b, alpha = alpha, tmax = tmax, y0 = None, noise = noise)
#         coAbd = np.round(sim_res[0], 0)
#         pearson_corr = np.corrcoef(coAbd, rowvar = False)
#         network_coocc = np.greater(pearson_corr, upper_cutoff).astype(int) - np.less(pearson_corr, lower_cutoff).astype(int)
#         np.fill_diagonal(network_coocc, 0)
#         G_coocc = nx.from_numpy_array(network_coocc)
        
#         if plotNetworks:
#                 G_inter = nx.from_numpy_array(interactionMx != 0)

#                 normalized_interMx = interactionMx / (interactionMx.max() - interactionMx.min())
#                 normalized_pearson_corr = pearson_corr / (pearson_corr.max() - pearson_corr.min())

#                 pos = nx.spring_layout(G_inter, seed=42)

#                 fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

#                 # Plot interaction network
#                 nx.draw(G_inter, pos, ax=ax1, with_labels=False, node_color='skyblue',
#                         node_size=8, edge_color=['red' if normalized_interMx[u,v] < 0 else 'blue' for u,v in G_inter.edges()],
#                         width=[abs(normalized_interMx[u,v]) * 3 for u,v in G_inter.edges()])
#                 ax1.set_title("Interaction")

#                 # Plot co-occurrence network
#                 nx.draw(G_coocc, pos, ax=ax2, with_labels=False, node_color='skyblue',
#                         node_size=8, edge_color=['red' if normalized_pearson_corr[u,v] < 0 else 'blue' for u,v in G_coocc.edges()],
#                         width=[abs(normalized_pearson_corr[u,v]) * 3 for u,v in G_coocc.edges()])
#                 ax2.set_title("Co-occurrence")

#                 plt.show()
    
#         return network_coocc, G_coocc

def get_coocc_network_traitmatch(interactionMx, cutoff, alpha, spsamplesize, numtrials, b, traits, tmax = 150, plotNetworks = False, environ = None, pval_alpha = 0.01, envStrength = 0.5):
        ''' run the cooccurrence network pipeline for trait matching environment noise
        '''
        if cutoff == None:
                lower_cutoff, upper_cutoff = get_null_cutoff_traitmatch(nSpec = interactionMx.shape[0], alpha = alpha, spsamplesize = spsamplesize, numtrials = numtrials, b = b, traits = traits, envStrength = envStrength, pval_alpha = pval_alpha, tmax = tmax)
        elif type(cutoff) == float:
                lower_cutoff = -cutoff
                upper_cutoff = cutoff
        elif len(cutoff) == 2:
                lower_cutoff = cutoff[0]
                upper_cutoff = cutoff[1]

        assert(lower_cutoff < upper_cutoff)
        assert(lower_cutoff < 0)
        assert(upper_cutoff > 0)

        sim_res = community_sim.simulate_coOccurrence_traitmatch(samplesize = spsamplesize, numtrials = numtrials, a = interactionMx, b = b, alpha = alpha, traits = traits, tmax = tmax, y0 = None, environ = environ, envStrength = envStrength)
        coAbd = np.round(sim_res[0], 0)
        pearson_corr = np.corrcoef(coAbd, rowvar = False)
        network_coocc = np.greater(pearson_corr, upper_cutoff).astype(int) - np.less(pearson_corr, lower_cutoff).astype(int)
        np.fill_diagonal(network_coocc, 0)
        G_coocc = nx.from_numpy_array(network_coocc)
        
        if plotNetworks:
                G_inter = nx.from_numpy_array(interactionMx != 0)

                normalized_interMx = interactionMx / (interactionMx.max() - interactionMx.min())
                normalized_pearson_corr = pearson_corr / (pearson_corr.max() - pearson_corr.min())

                pos = nx.spring_layout(G_inter, seed=42)

                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

                # Plot interaction network
                nx.draw(G_inter, pos, ax=ax1, with_labels=False, node_color='skyblue',
                        node_size=8, edge_color=['red' if normalized_interMx[u,v] < 0 else 'blue' for u,v in G_inter.edges()],
                        width=[abs(normalized_interMx[u,v]) * 3 for u,v in G_inter.edges()])
                ax1.set_title("Interaction")

                # Plot co-occurrence network
                nx.draw(G_coocc, pos, ax=ax2, with_labels=False, node_color='skyblue',
                        node_size=8, edge_color=['red' if normalized_pearson_corr[u,v] < 0 else 'blue' for u,v in G_coocc.edges()],
                        width=[abs(normalized_pearson_corr[u,v]) * 3 for u,v in G_coocc.edges()])
                ax2.set_title("Co-occurrence")

                plt.show()
    
        return network_coocc, G_coocc, pearson_corr