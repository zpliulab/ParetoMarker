list(
  expression="examples/input/expression.csv",
  metadata="examples/input/metadata.csv",
  network="examples/input/network.csv",
  candidate_genes="examples/input/candidate_genes.csv",
  sources="examples/input/sources.csv",
  targets="examples/input/targets.csv",
  output_dir="outputs/preprocessing",
  missing_fraction=0.75, knn_k=3, seed=7,
  cmi_threshold=0.5, cmi_order=1,
  state_a="Baseline", state_b="Cycle1Day3",
  optimization_state="Baseline",
  target_mode="supplied",
  de_subset_state="Baseline", de_group_column="label",
  de_reference="0", de_comparison="1", de_block_column=NULL,
  de_scale="minmax", de_fdr=0.01, de_abs_logfc=1.5
)
