# Stage interfaces

## R inputs

- `expression.csv`: first column `gene`; remaining columns are sample IDs; one row per gene. Values are numeric, with missing values allowed. This is an expression/abundance matrix, not FASTQ/BAM or a raw-sequencing alignment workflow. Set its measurement scale before use.
- `metadata.csv`: `sample_id,state,label`. IDs uniquely match expression columns. `state` identifies the two conditions used for differential MI. `label` is 0/1 for the optimization subset, defined before analysis. Repeated patients require a scientifically justified design; this package does not infer identity from sample-name substrings.
- `network.csv`: `source,target` gene identifiers for a directed prior network.
- `candidate_genes.csv`: optional `gene` column constraining the upstream candidate universe. Set its config field to NULL to use all expressed network genes.
- `sources.csv`, `targets.csv`: `gene` column for disease-related sources and response-related targets. Sources should come from the intended prior biological knowledge. With `target_mode="limma"`, targets instead come from the configured expression contrast.

When using limma, declare the comparison/reference groups and optional patient-block column. With `de_scale="minmax"`, limma operates on min–max scaled values: its coefficient column is not a conventional log2 fold change. Thresholds must be justified on the chosen scale. The example preserves the recovered threshold as a configurable value, not a general recommendation.

The R stage writes `weighted_network.csv` (`source,target,capacity`), `sources.csv`, `targets.csv`, `expression.csv` (sample rows, gene columns), and `labels.csv` (`sample_id,label`). Here `capacity=abs(MI_state_a−MI_state_b)`.

## Python to MATLAB

Python retains the union of nodes on bounded simple paths per source, then constructs the induced directed subgraph as the recovered code did. It adds virtual target-to-sink edges with configurable capacity. Genes incident to positive-flow edges become candidates, excluding the virtual sink. Flow values and minimum-cut edges are algorithm diagnostics; selected genes are derived from positive flow, not merely cut-edge endpoints.

`max_path_nodes` counts nodes; the recovered limit of six nodes corresponds to at most five edges. `sink_capacity` is explicit. `expected_modules` can be set in the flow configuration; the manuscript-oriented example sets it to 8 and requires all eight source-centered modules to be non-empty. Flow solutions may be nonunique, and the implementation uses sorted traversal for deterministic tie handling.

`outputs/flow/matlab_input/` contains:

- `expression.csv`: sample rows, first column `sample_id`, followed by candidate genes.
- `labels.csv`: matching sample IDs and binary labels.
- `network.csv`: directed edges among candidate genes.
- `initial_population.csv`: one row per source-centered module; the first column is `module_source` and the remaining columns are binary masks in candidate-gene space.

MATLAB reads this directory and uses `initial_population.csv` directly as the evolutionary initial population. For the manuscript configuration, `population_size` is 8, so the flow stage must provide eight non-empty source-centered modules. A module-count mismatch or an empty module raises an error rather than silently falling back to random initialization.

## Data scope

The public archive contains synthetic examples only. It contains no real expression matrices, real labels or real per-patient prediction scores. Algorithm fitness on discovery data is not independent validation; upstream selection is not nested by this pipeline.
