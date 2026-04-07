# this makes a whole bunch of ROC curves for direct edge prediciton when running things with different parameters

# similar to code from ./run_interactDetection.ipynb
# adjusts samplesize to 25, 50, and 75
# adjust numtrials to 50, 100, 250, 500, 1000
# adjust pInteract to 0.01, 0.05, 0.1

import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics

import metaCommunityMx
import community_sim

# fixed parameters 
nSpec = 100
pPositive = 0.5
alpha = 0.001
interactFactor = 1/5
interactStren = alpha * interactFactor
n_runs = 5

# used as a common x-axis grid so that ROC curves from different runs can be averaged together.
base_fpr = np.linspace(0, 1, 101)

# parameters to vary
samplesize_L = [25, 50, 75]
numtrials_L = [50, 100, 250, 500, 1000]
pInteract_L = [0.01, 0.05, 0.1]

for pInteract in pInteract_L:
    fig, axes = plt.subplots(
        len(samplesize_L), len(numtrials_L),
        figsize=(4 * len(numtrials_L), 4 * len(samplesize_L)),
        dpi=100
    )

    for row_idx, samplesize in enumerate(samplesize_L):
        for col_idx, numtrials in enumerate(numtrials_L):
            ax = axes[row_idx, col_idx]
            tprs = []

            for run in range(n_runs):
                # generate a new random matrix and b vector each run
                randMx = metaCommunityMx.getRandomMatrix(
                    maxInteractionStrength=interactStren,
                    pPositive=pPositive,
                    pInteract=pInteract,
                    nSpecies=nSpec,
                    pBidirectional=0
                )
                b = np.ones(nSpec) - np.random.uniform(low=0, high=0.1, size=(nSpec))

                # run GLV simulation and compute co-abundance matrix
                rand_res = community_sim.simulate_coOccurrence(
                    samplesize, numtrials, a=randMx, b=b, alpha=alpha, tmax=100, y0=None
                )
                rand_coAbd = np.round(rand_res[0], 0)

                # pearson correlation as predictor
                pearson_rand = np.corrcoef(rand_coAbd, rowvar=False)
                pearson_rand_abs = np.abs(pearson_rand)

                # binary ground truth: edge exists or not
                binary_rand = (randMx != 0).astype(int)
                mask = ~np.eye(nSpec, dtype=bool)
                binary_flat = binary_rand[mask]
                pearson_flat = pearson_rand_abs[mask]

                # compute and store ROC curve
                fpr, tpr, _ = metrics.roc_curve(binary_flat, pearson_flat)
                ax.plot(fpr, tpr, color='blue', alpha=0.2, lw=1)
                tpr_interp = np.interp(base_fpr, fpr, tpr)
                tpr_interp[0] = 0.0
                tprs.append(tpr_interp)

            # mean +/- std across runs
            tprs_arr = np.array(tprs)
            mean_tpr = tprs_arr.mean(axis=0)
            std_tpr = tprs_arr.std(axis=0)
            mean_auc = metrics.auc(base_fpr, mean_tpr)

            ax.plot(base_fpr, mean_tpr, color='blue', lw=2, label=f'Mean (AUC={mean_auc:.3f})')
            ax.fill_between(
                base_fpr,
                np.maximum(mean_tpr - std_tpr, 0),
                np.minimum(mean_tpr + std_tpr, 1),
                color='blue', alpha=0.25
            )
            ax.plot([0, 1], [0, 1], 'gray', lw=1, linestyle='--')
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1.05)
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=7, loc='lower right')

            # column headers (numtrials = M) on top row only
            if row_idx == 0:
                ax.set_title(f'M={numtrials}', fontsize=9)
            # row labels (S_m = samplesize) on left column only
            if col_idx == 0:
                ax.set_ylabel(f'S_m={samplesize}\nTPR', fontsize=8)
            else:
                ax.set_ylabel('')
            ax.set_xlabel('FPR', fontsize=8)

    fig.suptitle(
        f'ROC Curves - Direct Edge Prediction (p_ij={pInteract})\nRows: S_m, Cols: M, {n_runs} averaged runs each',
        fontsize=13
    )
    plt.tight_layout()
    fname = f'roc_direct_edge_pInteract{pInteract}.png'
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {fname}')