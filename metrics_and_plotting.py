import networkx as nx
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt
from scipy import stats

def test_cluster_similarity(clusters1, clusters2):
    ''' Calculates the proportion of pairs that are correctly paired/unpaired
    clusters1 is a list of sets of node names
    clusters2 is a list of sets of node names
    '''
    nodes = [item for s in clusters1 for item in s]
    assert set(nodes) == set([item for s in clusters2 for item in s])

    # Helper function to check if two nodes are in the same cluster
    def in_same_cluster(node1, node2, clusters):
        for cluster in clusters:
            if node1 in cluster and node2 in cluster:
                return True
        return False
    
    shared_pairs = 0
    both_non_pairs = 0
    c2_not_c1 = 0
    c1_not_c2 = 0
    for node1 in nodes:
        for node2 in nodes:
            if node1 != node2:
                if in_same_cluster(node1, node2, clusters1) and in_same_cluster(node1, node2, clusters2):
                    shared_pairs += 1
                elif not in_same_cluster(node1, node2, clusters1) and not in_same_cluster(node1, node2, clusters2):
                    both_non_pairs += 1
                elif not in_same_cluster(node1, node2, clusters1) and in_same_cluster(node1, node2, clusters2):
                    c2_not_c1 += 1
                else: # in_same_cluster(node1, node2, clusters1) and not in_same_cluster(node1, node2, clusters2):
                    c1_not_c2 += 1

    assert (shared_pairs + both_non_pairs + c1_not_c2 + c2_not_c1) == len(nodes)**2 - len(nodes)

    prop_correct = (shared_pairs + both_non_pairs) / (len(nodes)**2 - len(nodes))
    return shared_pairs, both_non_pairs, c1_not_c2, c2_not_c1, prop_correct


def test_bipartite_cluster_similarity(clusters1, clusters2):
    ''' Calculates the proportion of pairs that are correctly paired/unpaired (clusters1 is original interaction matrix, clusters2 is cooccurr.)
    clusters1 is a list of sets of node names
    clusters2 is a list of sets of node names
    '''
    nodes1 = [item for s in clusters1 for item in s]
    nodes2 = [item for s in clusters2 for item in s]

    # Helper function to check if two nodes are in the same cluster
    def in_same_cluster(node1, node2, clusters):
        for cluster in clusters:
            if node1 in cluster and node2 in cluster:
                return True
        return False
    
    shared_pairs = 0
    both_non_pairs = 0
    c2_not_c1 = 0
    c1_not_c2 = 0
    for node1 in nodes1:
        for node2 in nodes2:
            if node1 != node2:
                if in_same_cluster(node1, node2, clusters1) and in_same_cluster(node1, node2, clusters2):
                    shared_pairs += 1
                elif not in_same_cluster(node1, node2, clusters1) and not in_same_cluster(node1, node2, clusters2):
                    both_non_pairs += 1
                elif not in_same_cluster(node1, node2, clusters1) and in_same_cluster(node1, node2, clusters2):
                    c2_not_c1 += 1
                else: # in_same_cluster(node1, node2, clusters1) and not in_same_cluster(node1, node2, clusters2):
                    c1_not_c2 += 1

    #correctness based on clustered nodes from original interaction matrix
    prop_correct = (shared_pairs + both_non_pairs) / (len(nodes1)**2 - len(nodes1))
    return shared_pairs, both_non_pairs, c1_not_c2, c2_not_c1, prop_correct

def test_bipartite_cluster_similarity_posonly(clusters1, clusters2):
    ''' Calculates the proportion of pairs that are correctly paired/unpaired (clusters1 is original interaction matrix, clusters2 is cooccurr.)
    clusters1 is a list of sets of node names
    clusters2 is a list of sets of node names
    '''
    nodes1 = [item for s in clusters1 for item in s]
    nodes2 = [item for s in clusters2 for item in s]

    # Helper function to check if two nodes are in the same cluster
    def in_same_cluster(node1, node2, clusters):
        for cluster in clusters:
            if node1 in cluster and node2 in cluster:
                return True
        return False
    
    shared_pairs = 0
    both_non_pairs = 0
    c2_not_c1 = 0
    c1_not_c2 = 0
    for node1 in nodes1:
        for node2 in nodes2:
            if node1 != node2:
                if in_same_cluster(node1, node2, clusters1) and in_same_cluster(node1, node2, clusters2):
                    shared_pairs += 1
                elif not in_same_cluster(node1, node2, clusters1) and not in_same_cluster(node1, node2, clusters2):
                    both_non_pairs += 1
                elif not in_same_cluster(node1, node2, clusters1) and in_same_cluster(node1, node2, clusters2):
                    c2_not_c1 += 1
                else: # in_same_cluster(node1, node2, clusters1) and not in_same_cluster(node1, node2, clusters2):
                    c1_not_c2 += 1

    #correctness based on clustered nodes from original interaction matrix
    prop_correct = (shared_pairs) / (shared_pairs+c1_not_c2)
    return shared_pairs, both_non_pairs, c1_not_c2, c2_not_c1, prop_correct

def calculate_centrality_metrics(G, nodes, centrality_metrics = {}):
    """Calculate various centrality metrics for all nodes in the graph.
        nodes = list(G.nodes()) to make sure the nodes are always in the same order
    """
    
    # Degree centrality
    degree_cent = nx.degree_centrality(G)
    degree = np.array([degree_cent[node] for node in nodes])
    if 'degree' not in centrality_metrics:
        centrality_metrics['degree'] = degree
    else:
        centrality_metrics['degree'] = np.append(centrality_metrics['degree'], degree)
    
    # Betweenness centrality
    betweenness_cent = nx.betweenness_centrality(G)
    betweenness = np.array([betweenness_cent[node] for node in nodes])
    if 'betweenness' not in centrality_metrics:
        centrality_metrics['betweenness'] = betweenness
    else:
        centrality_metrics['betweenness'] = np.append(centrality_metrics['betweenness'], betweenness)
    
    # Closeness centrality
    closeness_cent = nx.closeness_centrality(G)
    closeness = np.array([closeness_cent[node] for node in nodes])
    if 'closeness' not in centrality_metrics:
        centrality_metrics['closeness'] = closeness
    else:
        centrality_metrics['closeness'] = np.append(centrality_metrics['closeness'], closeness)
    
    # Eigenvector centrality
    try:
        eigenvector_cent = nx.eigenvector_centrality(G, max_iter=10000, tol=1e-03)
    except nx.PowerIterationFailedConvergence:
        eigenvector_cent = nx.eigenvector_centrality_numpy(G)
    eigenvector = np.array([eigenvector_cent[node] for node in nodes])
    if 'eigenvector' not in centrality_metrics:
        centrality_metrics['eigenvector'] = eigenvector
    else:
        centrality_metrics['eigenvector'] = np.append(centrality_metrics['eigenvector'], eigenvector)
    
    return centrality_metrics

def plot_hub_roc(centrality_metrics, hub_status, title="ROC Curves for Hub Metrics"):
    """Plot ROC curves for different centrality metrics."""
    
    plt.figure(figsize=(10, 8))
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
    
    auc_scores = {}
    
    for i, (metric_name, centrality_values) in enumerate(centrality_metrics.items()):
        if centrality_values is None:
            continue
            
        # Calculate ROC curve
        fpr, tpr, _ = metrics.roc_curve(hub_status, centrality_values)
        roc_auc = metrics.auc(fpr, tpr)
        auc_scores[metric_name] = roc_auc
        
        # Plot ROC curve
        plt.plot(fpr, tpr, color=colors[i % len(colors)], lw=2,
                label=f'{metric_name.title()} (AUC = {roc_auc:.3f})')
    
    # Plot diagonal line (random classifier)
    plt.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--', 
             label='Random (AUC = 0.500)')
    
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', size = 16)
    plt.ylabel('True Positive Rate', size = 16)
    plt.title(title)
    plt.legend(loc="lower right", fontsize=16)
    plt.grid(True, alpha=0.3)
    
    return auc_scores

def calculate_nodf(matrix):
    """
    Calculates the NODF (Nestedness based on Overlap and Decreasing Fill), Almeida-Neto (2008). Can take a matrix with weighted entries, will convert to binary.
    """
    matrix = (matrix != 0).astype(int)
    # Step 1: Sort the matrix by marginal totals to find NODFmax
    # For rows (descending sum)
    row_sums = matrix.sum(axis=1)
    row_sort_indices = np.argsort(row_sums)[::-1]
    sorted_matrix = matrix[row_sort_indices, :]

    # For columns (descending sum)
    col_sums = sorted_matrix.sum(axis=0)
    col_sort_indices = np.argsort(col_sums)[::-1]
    sorted_matrix = sorted_matrix[:, col_sort_indices]

    # Step 2: Calculate nestedness for pairs of rows
    n_rows, n_cols = sorted_matrix.shape
    n_pairs_rows = n_rows * (n_rows - 1) / 2
    n_pairs_cols = n_cols * (n_cols - 1) / 2

    paired_nestedness_rows = []
    for i in range(n_rows):
        for j in range(n_rows):
            if i < j:
                k_i = sorted_matrix[i].sum()
                k_j = sorted_matrix[j].sum()
                # Check decreasing fill property
                if k_i > k_j:
                    # Calculate paired overlap (nij/kj)
                    nij = np.sum(sorted_matrix[i] * sorted_matrix[j])
                    if k_j > 0:
                        paired_nestedness_rows.append(nij / k_j)
                    else:
                        paired_nestedness_rows.append(0)
                else:
                    paired_nestedness_rows.append(0)

    # Step 3: Calculate nestedness for pairs of columns
    paired_nestedness_cols = []
    for k in range(n_cols):
        for l in range(n_cols):
            if k < l:
                k_k = sorted_matrix[:, k].sum()
                k_l = sorted_matrix[:, l].sum()
                # Check decreasing fill property
                if k_k > k_l:
                    # Calculate paired overlap (nkl/kl)
                    nkl = np.sum(sorted_matrix[:, k] * sorted_matrix[:, l])
                    if k_l > 0:
                        paired_nestedness_cols.append(nkl / k_l)
                    else:
                        paired_nestedness_cols.append(0)
                else:
                    paired_nestedness_cols.append(0)

    # Step 4: Combine and normalize the results
    n_paired_sum = np.sum(paired_nestedness_rows) + np.sum(paired_nestedness_cols)
    denominator = n_pairs_rows + n_pairs_cols

    if denominator == 0:
        return 0.0
    
    nodf = n_paired_sum/denominator
    
    return nodf

def calculate_barber_modularity_bipartite(B, partition_row, partition_col,numgroup):
    """
    Calculates the Barber bipartite modularity score for a given biadjacency matrix (lower left quadrant of original interaction matrix), partitions, and group size.
    Assumes evenly sized groups (i.e. numgroup members of type A, numgroup members of type B).
    B: biadjacency matrix (1s/0s) with group A on the rows, group B on columns, 1 indicates interaction between A and B.T
    partition_row: list of lists; internal lists are identified modules, in order that matches partition_col
    numgroup: size of membership in type A or type B (half the total species)

    """
    E = B.sum() # Total number of edges
    if E == 0:
        return 0
        
    q = B.sum(axis=1) # Degrees of nodes in the first set (rows)
    d = B.sum(axis=0) # Degrees of nodes in the second set (columns)
    
    Q_b = 0.0
    for i in range(B.shape[0]): # Iterate over rows (first set)
        for j in range(B.shape[1]): # Iterate over columns (second set)
            for k in range(len(partition_row)):
                if i in partition_row[k] and (j+numgroup) in partition_col[k]:  # Check if in the same community
                    Q_b += (B[i, j] - (q[i] * d[j]) / E)   
    return Q_b / E

def rewire_modular(matrix, modules, nEdges, interactStren):
    '''
    takes nEdges within-module interactions and moves them to a random other edge

    matrix is the GLV interaction matrix (numpy array), modules is a list of sets, and nEdges is how many edges to move
    '''
    moved = 0
    while moved < nEdges:
        # choose an edge to remove
        mod1 = np.random.choice(range(len(modules)))
        node1_1 = np.random.choice(range(len(modules[mod1])))
        node1_2 = np.random.choice(range(len(modules[mod1])))

        # choose an edge to add
        node2_1 = np.random.choice(range(matrix.shape[0]))
        node2_2 = np.random.choice(range(matrix.shape[0]))

        e1 = matrix[modules[mod1][node1_1], modules[mod1][node1_2]]
        e2 = matrix[node2_1, node2_2]
        if node1_1 != node1_2 and node2_1 != node2_2 and e1 != 0 and e2 == 0:
            matrix[modules[mod1][node1_1], modules[mod1][node1_2]] = 0
            matrix[node2_1, node2_2] = np.random.uniform(low = -interactStren, high = interactStren)
            moved += 1

    return matrix


def plot_corr_vs_interact(corrMx, interactMx, cutoff_null_95, title):
    # exclude diagonal elements
    mask = ~np.eye(interactMx.shape[0], dtype=bool)
    interactMx_masked = interactMx[mask]
    coAbd_masked = corrMx[mask]
    
    # Separate zero and non-zero entries
    nonzero_mask = interactMx_masked != 0
    zero_mask = interactMx_masked == 0
    
    # Non-zero entries for scatter plot
    interactMx_nonzero = interactMx_masked[nonzero_mask]
    coAbd_nonzero = coAbd_masked[nonzero_mask]
    
    # Zero entries for histogram
    coAbd_zeros = coAbd_masked[zero_mask]
    
    # linear regression on non-zero data
    slope, intercept, r_value, p_value, std_err = stats.linregress(interactMx_nonzero, coAbd_nonzero)
    r_squared = r_value**2
    
    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8), gridspec_kw={'width_ratios': [3, 1]})
    
    # Main scatter plot (left subplot)
    scatter = ax1.scatter(interactMx_nonzero, coAbd_nonzero, alpha=0.3, s=40)
    
    if type(cutoff_null_95) == float:
            lower_cutoff = -cutoff_null_95
            upper_cutoff = cutoff_null_95
    elif len(cutoff_null_95) == 2:
            lower_cutoff = cutoff_null_95[0]
            upper_cutoff = cutoff_null_95[1]    

    ax1.axhline(y=upper_cutoff, color='black', linestyle='--', alpha=0.7, linewidth=1, label='cutoff for link significance')
    ax1.axhline(y=lower_cutoff, color='black', linestyle='--', alpha=0.7, linewidth=1)
    ax1.axhline(y=0, color='gray', linestyle='-', alpha=0.5, linewidth=2)

    # Plot linear fit
    x_range = np.array([interactMx_nonzero.min(), interactMx_nonzero.max()])
    y_fit = slope * x_range + intercept
    ax1.plot(x_range, y_fit, 'r-', linewidth=2, label=f'Linear fit (${{R^2}}$ = {r_squared:.4f})')
    
    #ax1.text(0.05, 0.95, f'R² = {r_squared:.4f}', transform=ax1.transAxes,
    #        fontsize=14, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Formatting for scatter plot
    ax1.set_xlabel('Pairwise interaction strength \n (Non-zero)', fontsize=16)
    ax1.set_ylabel('Pairwise Pearson correlation of abundances \n (Interaction strength > 0)', fontsize=16)
    ax1.set_title(title, fontsize=14)
    ax1.legend(fontsize=16)
    ax1.grid(True, alpha=0.3)

    # center the plot
    coabd_max = np.max(np.abs(coAbd_masked))
    ax1.set_ylim((-coabd_max, coabd_max))
    
    # Histogram combining interactions and non-interactions (right subplot)
    # Separate positive and negative correlations for non-zero interactions
    positive_corr = coAbd_nonzero[coAbd_nonzero >= 0]
    negative_corr = coAbd_nonzero[coAbd_nonzero < 0]

    # Create bins that span the full range
    all_data = np.concatenate([coAbd_nonzero, coAbd_zeros])
    all_bins = np.linspace(all_data.min(), all_data.max(), 21)

    # Prepare data for stacked histogram
    data_to_plot = []
    labels_to_plot = []
    colors_to_plot = []

    if len(positive_corr) > 0:
        data_to_plot.append(positive_corr)
        labels_to_plot.append(f'Positive interactions \n (n={len(positive_corr)})')
        colors_to_plot.append('blue')

    if len(negative_corr) > 0:
        data_to_plot.append(negative_corr)
        labels_to_plot.append(f'Negative interactions \n  (n={len(negative_corr)})')
        colors_to_plot.append('red')

    if len(coAbd_zeros) > 0:
        data_to_plot.append(coAbd_zeros)
        labels_to_plot.append(f'No interaction \n (n={len(coAbd_zeros)})')
        colors_to_plot.append('orange')

    # Plot stacked histogram
    ax2.hist(data_to_plot, bins=all_bins, orientation='horizontal',
            stacked=True, alpha=0.7, color=colors_to_plot,
            edgecolor='black', linewidth=0.5, label=labels_to_plot)

    # Add significance cutoff lines
    ax2.axhline(y=upper_cutoff, color='black', linestyle='--', alpha=0.7, linewidth=1)
    ax2.axhline(y=lower_cutoff, color='black', linestyle='--', alpha=0.7, linewidth=1)
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=2)

    ax2.set_xlabel('Count (log scale)', fontsize=16)
    ax2.set_ylabel('Pairwise Pearson correlation of abundances', fontsize=16)
    ax2.set_xscale('log')
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3)

    # Match y-axis range with scatter plot
    ax2.set_ylim(ax1.get_ylim())
    
    plt.tight_layout()
    plt.show()
    
    return {"n_nonzero": len(interactMx_nonzero), "n_zero": len(coAbd_zeros), 
            "r_sq": r_squared, "corr_coef (r)": r_value, "p-value": p_value, 
            "slope": slope, "intercept": intercept}


def plot_corr_vs_interact2(corrMx, interactMx, cutoff_null_95, title):
    # exclude diagonal elements
    mask = ~np.eye(interactMx.shape[0], dtype=bool)
    interactMx_masked = interactMx[mask]
    coAbd_masked = corrMx[mask]
    
    # Separate zero and non-zero entries
    nonzero_mask = interactMx_masked != 0
    zero_mask = interactMx_masked == 0
    
    # Non-zero entries for scatter plot
    interactMx_nonzero = interactMx_masked[nonzero_mask]
    coAbd_nonzero = coAbd_masked[nonzero_mask]
    
    # Zero entries for histogram
    coAbd_zeros = coAbd_masked[zero_mask]
    
    # linear regression on non-zero data
    slope, intercept, r_value, p_value, std_err = stats.linregress(interactMx_nonzero, coAbd_nonzero)
    r_squared = r_value**2
    
    # Create subplots with three panels
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6), 
                                        gridspec_kw={'width_ratios': [3, 1, 1]})
    
    # Main scatter plot (left subplot)
    scatter = ax1.scatter(interactMx_nonzero, coAbd_nonzero, alpha=0.3, s=40)
    ax1.axhline(y=cutoff_null_95, color='black', linestyle='--', alpha=0.7, linewidth=1, 
                label='0.95 cutoff for link significance')
    ax1.axhline(y=-cutoff_null_95, color='black', linestyle='--', alpha=0.7, linewidth=1)
    
    # Plot linear fit
    x_range = np.array([interactMx_nonzero.min(), interactMx_nonzero.max()])
    y_fit = slope * x_range + intercept
    ax1.plot(x_range, y_fit, 'r-', linewidth=2, 
             label=f'Linear fit: y = {slope:.3f}x + {intercept:.3f}')
    ax1.text(0.05, 0.95, f'R² = {r_squared:.4f}', transform=ax1.transAxes,
             fontsize=14, verticalalignment='top', 
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Formatting for scatter plot
    ax1.set_xlabel('Pairwise Interaction Strength \n (non-zero)', fontsize=16)
    ax1.set_ylabel('Pairwise Pearson Correlation of Abundance \n (Interaction strength > 0)', fontsize=16)
    ax1.set_title(title, fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Histogram of coAbd where interactMx = 0 (middle subplot)
    ax2.hist(coAbd_zeros, orientation='horizontal', bins=20, alpha=0.7,
                color='orange', edgecolor='black', linewidth=0.5)
    ax2.axhline(y=cutoff_null_95, color='black', linestyle='--', alpha=0.7, linewidth=1)
    ax2.axhline(y=-cutoff_null_95, color='black', linestyle='--', alpha=0.7, linewidth=1)
    ax2.set_xlabel('Count', fontsize=12)
    ax2.set_ylabel('Co-occurrence Pearson Correlation (No interaction)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    # Match y-axis range with scatter plot
    ax2.set_ylim(ax1.get_ylim())
    
    # Third subplot: Histogram of correlations for non-zero interactions, colored by sign
    # Separate positive and negative correlations
    positive_corr = coAbd_nonzero[coAbd_nonzero >= 0]
    negative_corr = coAbd_nonzero[coAbd_nonzero < 0]
    
    # Create bins that span the full range
    all_bins = np.linspace(coAbd_nonzero.min(), coAbd_nonzero.max(), 21)
    
    # Plot histograms
    if len(positive_corr) > 0:
        ax3.hist(positive_corr, bins=all_bins, orientation='horizontal', 
                    alpha=0.7, color='blue', edgecolor='black', linewidth=0.5,
                    label=f'Positive (n={len(positive_corr)})')
    
    if len(negative_corr) > 0:
        ax3.hist(negative_corr, bins=all_bins, orientation='horizontal', 
                    alpha=0.7, color='red', edgecolor='black', linewidth=0.5,
                    label=f'Negative (n={len(negative_corr)})')
    
    # Add significance cutoff lines
    ax3.axhline(y=cutoff_null_95, color='black', linestyle='--', alpha=0.7, linewidth=1)
    ax3.axhline(y=-cutoff_null_95, color='black', linestyle='--', alpha=0.7, linewidth=1)
    ax3.axhline(y=0, color='gray', linestyle='-', alpha=0.5, linewidth=1)
    
    ax3.set_xlabel('Count', fontsize=12)
    ax3.set_ylabel('Pearson Correlation (Non-zero interactions)', fontsize=12)
    ax3.set_title('Correlation Distribution\n(Non-zero interactions)', fontsize=12)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Match y-axis range with scatter plot
    ax3.set_ylim(ax1.get_ylim())
    
    plt.tight_layout()
    plt.show()
    
    # Enhanced return statistics
    n_positive = len(coAbd_nonzero[coAbd_nonzero >= 0])
    n_negative = len(coAbd_nonzero[coAbd_nonzero < 0])
    
    return {
        "n_nonzero": len(interactMx_nonzero), 
        "n_zero": len(coAbd_zeros),
        "n_positive_corr": n_positive,
        "n_negative_corr": n_negative,
        "r_sq": r_squared, 
        "corr_coef (r)": r_value, 
        "p-value": p_value,
        "slope": slope, 
        "intercept": intercept
    }

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

def make_random_mistakes(true_interact_network, cooc_net):
    ''' Takes the interaction network and the co-occurrence network (symmetric adjacency matrices)
    1) computes the confusion matrix (false and true positives)
    2) modifies the interaction network such that there are the same number of false negatives and false positives,
    but they are in random positions. Then outputs that null cooccurrence network.
    '''
    true_bin = true_interact_network != 0
    cooc_bin = cooc_net != 0
    n = true_bin.shape[0]

    iu = np.triu_indices(n, k=1)
    true_upper = true_bin[iu]
    cooc_upper = cooc_bin[iu]

    # get confusion mx
    true_pos = np.sum(true_upper & cooc_upper)
    false_pos = np.sum(~true_upper & cooc_upper)
    false_neg = np.sum(true_upper & ~cooc_upper)
    true_neg = np.sum(~true_upper & ~cooc_upper)

    pos_pair_indices = np.where(true_upper)[0]
    neg_pair_indices = np.where(~true_upper)[0]

    # pick random fps and fns to add
    fn_choice = np.random.choice(pos_pair_indices, size=false_neg, replace=False)
    fp_choice = np.random.choice(neg_pair_indices, size=false_pos, replace=False)

    null_cooc_net = true_bin.astype(int)

    rows_fn, cols_fn = iu[0][fn_choice], iu[1][fn_choice]
    null_cooc_net[rows_fn, cols_fn] = 0
    null_cooc_net[cols_fn, rows_fn] = 0

    rows_fp, cols_fp = iu[0][fp_choice], iu[1][fp_choice]
    null_cooc_net[rows_fp, cols_fp] = 1
    null_cooc_net[cols_fp, rows_fp] = 1

    # check that the null network has the same number of edges as the real cooc_net
    null_cooc_bin = null_cooc_net != 0
    assert np.sum(null_cooc_bin[iu]) == np.sum(cooc_upper)

    # check that the null network has the same number of false positives/negatives as cooc_net
    null_upper = null_cooc_bin[iu]
    null_false_pos = np.sum(~true_upper & null_upper)
    null_false_neg = np.sum(true_upper & ~null_upper)
    assert null_false_pos == false_pos
    assert null_false_neg == false_neg

    confusion_matrix = (true_pos, false_pos, false_neg, true_neg)
    return null_cooc_net, confusion_matrix

