# Code for "Pigeon Dies Alone on a Rock and other opinions"

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
3. Analyze and compare networks using metrics tools

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

### Frequency of interaction types
**[frequency_interaction_types.ipynb](frequency_interaction_types.ipynb)**

### Hub species and node centrality
**[run_hub.ipynb](run_hub.ipynb)**

### Detecting modules and modularity
**[modularity_clustering_correctness.ipynb](modularity_clustering_correctness.ipynb)**

### Detecting nestedness
**[run_nested.ipynb](run_nested.ipynb)**

### Differentiating between node degree distributions
**[degree_sequence_rank.ipynb](degree_sequence_rank.ipynb)**
**[process_data_degSeq.py](process_data_degSeq.py)**
**[powerLaw.r](powerLaw.r)**

### Role of the environment, Detecting shifting interactions
**[run_cooccurrence.py](run_cooccurrence.py)**



