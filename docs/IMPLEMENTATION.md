# Implementation details

## Workflow structure

ParetoMarker is organized as three connected computational stages.

### 1. R preprocessing and network refinement

`R/preprocess.R` performs expression preprocessing, optional limma-based target selection, PCA-CMI network refinement, and state-specific mutual-information calculations. Edge capacity for the downstream information-flow analysis is defined as the absolute difference in mutual information between the two phenotype states.

`R/PCA_CMI.R` contains the PCA-CMI functions used by the preprocessing stage.

### 2. Python information-flow module construction

`python/flow.py` constructs bounded source-target paths, forms source-centered induced subnetworks, adds a virtual sink, and applies maximum-flow/minimum-cut calculations. Genes incident to positive-flow edges are retained as candidate genes. The implementation includes reverse residual edges and checks agreement between maximum-flow and minimum-cut values.

For each source-centered module, the Python stage exports a binary candidate-gene mask. These masks are written to `initial_population.csv` and passed directly to MATLAB.

### 3. MATLAB multi-objective optimization and voting

`matlab/run_paretomarker.m` loads the expression matrix, phenotype labels, candidate network, and `initial_population.csv`. The module masks are used directly as the initial evolutionary population. No random-initialization fallback is used when the manuscript configuration is selected.

The optimization evaluates five objectives in manuscript order:

1. ROC AUC;
2. negative AIC-like complexity criterion;
3. phenotype association score;
4. expression entropy;
5. network modularity.

The evolutionary loop uses `ModelLearning`/`MMLP` and `DEgenerator` for MLP-guided candidate updates, followed by alpha/beta dynamic selection. The manuscript configuration uses alpha = 0.4, beta = 3, population size = 8, and 800 candidate-solution evaluations. The final recurrence-based voting rule retains genes present in at least four of the eight final candidate solutions.

## Input consistency checks

The stage interfaces enforce alignment of sample identifiers, gene identifiers, network nodes, and binary phenotype labels. For the manuscript configuration, the Python stage must produce eight non-empty source-centered modules. A module-count mismatch or empty module raises an error instead of changing the initialization strategy.

## Scope

This repository contains the core ParetoMarker implementation used for candidate generation, information-flow module construction, five-objective optimization, dynamic selection, and recurrence-based voting. Downstream cohort validation, statistical summaries, and manuscript figure/table generation are separate analyses and are not required to execute the core algorithm.
