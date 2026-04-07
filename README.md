# Code for "Co-occurrence networks can preserve emergent properties of ecological communities"

DOI: [To be added]

## Key Question

Can co-occurrence patterns (statistical correlations in species abundances) reliably reveal larger properties of ecological interaction networks?

## Approach

We simulate ecological communities with known interaction structures using generalized Lotka-Volterra dynamics, then analyze the co-occurrence patterns that emerge. By comparing simulated interaction networks with derived co-occurrence networks, we test which characteristics of the interaction network are preserved in the co-occurrence network.

## Requirements

- Python
    - **[requirements.txt](requirements.txt)**
- R
    - poweRlaw
    - dplyr
    - gridExtra
    - ggplot2

## General Usage

1. Generate interaction networks with specific structures using `metaCommunityMx.py`
2. Simulate community dynamics using `run_cooccurrence.py`
3. Analyze and compare networks using network tools

See Jupyter notebooks for detailed examples and analyses.

## Main Components

### Simulation Framework
- **[community_sim.py](community_sim.py)**: Generalized Lotka-Volterra population dynamics with deterministic and stochastic simulations
- **[metaCommunityMx.py](metaCommunityMx.py)**: Generation of interaction matrices with controlled network structures (modular, hub-based, nested, scale-free, random)

### Analysis Pipeline
- **[run_cooccurrence.py](run_cooccurrence.py)**: Calculates co-occurrence networks from simulated abundances using Pearson correlations and statistical thresholds
- **[metrics_and_plotting.py](metrics_and_plotting.py)**: Network analysis tools including centrality metrics and visualization

## Paper sections (see paper methods and results)

### Detecting direct interactions
**[run_interactDetection.ipynb](run_interactDetection.ipynb)**

### Hub species and node centrality
- **[run_hub.ipynb](run_hub.ipynb)**

### Detecting modules and modularity
- **[run_modularity_clustering_correctness.ipynb](modularity_clustering_correctness.ipynb)**

### Detecting nestedness and modularity in bipartite networks
- **[run_nested_modular.ipynb](run_nested.ipynb)**

### Differentiating between node degree distributions
- See **[run_powerlaw.sh](run_powerlaw.sh)** for usage.
    - **[process_data_degSeq.py](process_data_degSeq.py)**
    - **[powerLaw.r](powerLaw.r)**

- **[run_degree_sequence_rank.ipynb](degree_sequence_rank.ipynb)**
- **[run_degree_not_binom.ipynb](run_degree_not_binom.ipynb)**

### Supplemental results
#### Role of the environment
- **[run_environ.ipynb](run_environ.ipynb)**

#### Frequency of interaction types
- **[run_frequency_interaction_types.ipynb](frequency_interaction_types.ipynb)**

#### Testing samples needed for power law detection
- **[poweRlaw_test.r](poweRlaw_test.r)**

#### Testing different sets of parameters
- **[loop_parms.py](loop_parms.py)**