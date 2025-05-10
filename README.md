# MSstats_thermal_profiling

#The original work is in github.com/rokapre/MSstatsThermalProfiler and was moved

#Please use this repository as it contains the most recent work to date.
# Install from GitHub
# install.packages("devtools")
```r
devtools::install_github("CleanValidation/MSstatsThermalProfile")


##Step 1: Install renv: install.packages(renv)
##Step 2: Initialize renv: renv::init()
##Step 3: restore the package versions for the manuscript: renv::restore()
##Step 4: initialize git lfs in terminal: git lfs install
##Step 5: track rda files in terminal: git lfs track "*.rda"
##Step 6: pull files from lfs in terminal: git lfs fetch
##Step 7: update git lfs in terminal: git lfs update

#Please install MSstats
```r
if (!requireNamespace("MSstats", quietly = TRUE)) {
    BiocManager::install("MSstats")
    BiocManager::install("MSstatsConvert")
    BiocManager::install("MSstatsTMT")
}

### Windows users


To install this package from source, you’ll need [Rtools 4.2](https://cran.r-project.org/bin/windows/Rtools/rtools42/) installed and properly configured for R 4.2.
You can test your setup with:

```r
pkgbuild::check_build_tools(debug = TRUE)
