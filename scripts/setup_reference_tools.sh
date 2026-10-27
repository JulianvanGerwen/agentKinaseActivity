#!/usr/bin/env bash
# Reproduce reference-tool environment (R: JUMPsem, lavaan, KSEAapp; Python: decoupler, gseapy)
set -euo pipefail
W=/workspace/external
export MAMBA_ROOT_PREFIX=$W/mamba
MM=/home/tools/bin/micromamba   # or download: https://micro.mamba.pm/api/micromamba/linux-64/latest
# 1) R env (conda-forge). 'compilers' needed to build CRAN pkgs; r-efatools is NOT on conda-forge -> CRAN.
$MM create -y -p $W/envs/r -c conda-forge r-base=4.4 r-lavaan r-psych r-remotes r-tidyr r-stringr \
  r-dplyr r-data.table r-mass r-matrix r-ggplot2 r-plyr r-gdata r-rcurl r-rcpp r-rcpparmadillo \
  r-tidyverse r-devtools compilers
export PATH=$W/envs/r/bin:$PATH     # REQUIRED so R finds x86_64-conda-linux-gnu-cc
export LC_ALL=C.UTF-8 LANG=C.UTF-8  # REQUIRED: JUMPsem lazydata build fails on BOM csv in the C/POSIX locale
Rscript -e 'options(repos="https://cloud.r-project.org"); install.packages(c("EFAtools","KSEAapp"))'
# 2) JUMPsem from GitHub (official: Wanglab-UTHSC/JUMPsem, commit 52919f5 ~ v1.1.0)
mkdir -p $W/src
[ -d $W/src/JUMPsem ] || git clone https://github.com/Wanglab-UTHSC/JUMPsem.git $W/src/JUMPsem
(cd $W/src/JUMPsem && git checkout 52919f560b9ac038f3641b8794ddacfd3eca40c6)
Rscript -e 'remotes::install_local("'$W'/src/JUMPsem",dependencies=FALSE,upgrade="never",force=TRUE)'
# 3) Python references
/home/plm_env/bin/python -m pip install decoupler gseapy
