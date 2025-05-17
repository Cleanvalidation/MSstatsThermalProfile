# MSstats Thermal Profiling

The original work is at:  
[rokapre/MSstatsThermalProfiler](https://github.com/rokapre/MSstatsThermalProfiler)  

This repository contains the **most recent updates** and improvements,  
so please use this one instead.

---

## Installation

### 1. Install `devtools` (if not already installed)


if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

### 2. Install MSstatsThermalProfile package from GitHub
devtools::install_github("CleanValidation/MSstatsThermalProfile")

### 3. Install MSstats and related packages if needed
if (!requireNamespace("MSstats", quietly = TRUE)) {
  BiocManager::install("MSstats")
  BiocManager::install("MSstatsConvert")
  BiocManager::install("MSstatsTMT")
}

## --- Using renv for reproducible environments ---

### Step 1: Install renv if needed
 install.packages("renv")

### Step 2: Initialize renv in your project directory
 renv::init()

### Step 3: Restore package versions from the lockfile
 renv::restore()

##--- Git Large File Storage (LFS) setup ---

### Step 4: Initialize git lfs (run this in your terminal)
 git lfs install

### Step 5: Track .rda files with git lfs (terminal)
git lfs track "*.rda"

### Step 6: Fetch .rda files from git lfs (terminal)
 git lfs fetch

### Step 7: Update git lfs (terminal)
 git lfs update

## --- Windows-specific note ---

# To build packages from source on Windows, you need Rtools 4.2:
###  https://cran.r-project.org/bin/windows/Rtools/rtools42/

# Check if your build tools are properly configured:
pkgbuild::check_build_tools(debug = TRUE)
