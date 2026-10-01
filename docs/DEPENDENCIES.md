# Dependencies

Python uses pandas and its NumPy dependency; install the versions in requirements.txt. The active flow algorithm uses the Python standard library and does not require NetworkX. The preserved legacy scripts reference additional old dependencies; they are reference sources, not the entry points.

In R, install infotheo from CRAN and impute/limma through Bioconductor using an R-compatible Bioconductor version. Example commands in an R session:

```r
install.packages(c("infotheo", "BiocManager"))
BiocManager::install(c("impute", "limma"))
```

The R wrapper records `sessionInfo()` when run so the execution environment can be archived with each analysis.

The MATLAB entry point requires Statistics and Machine Learning Toolbox (`fitcsvm`, `fitPosterior`, `perfcurve`, `kmeans`). Bundled source dependencies and notices are preserved; `run_info.mat` records the MATLAB version and installed products for each run.
