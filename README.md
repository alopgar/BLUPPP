[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![Last version](https://img.shields.io/github/tag/alopgar/BLUPPP.svg)](https://img.shields.io/github/tag/alopgar/BLUPPP.svg)

# BLUPPP
Pipeline that prepares genotypes files for BLUPF90

## 0. Software requirements:
For the correct functioning of this scripts, the installation of several software is required:
- **gfortran-9**: https://fortran-lang.org/en/learn/os_setup/install_gfortran/
- **PLINK 1.9**: https://www.cog-genomics.org/plink/
- **R 4.1.3**: https://cran.r-project.org/bin/windows/base/old/4.1.3/

### R dependencies:
- **tidyverse**: `install.packages("tidyverse")`

## 1. Installation:
a) Download all the ./bin files in your installation directory.  
b) Add your installation directory to your `PATH` variable in `~/.bashrc` file.  
c) Change `BINPATH` variable inside BLUPPP_exe.sh to point your BLUPPP installation directory.

## 2. Pipeline execution:
a) Use `BLUPPP_exe.sh -h` for more information.  
b) Download the parameter file model (`BLUPPP_parameters.par`) and fill variables.  
c) Run as: `BLUPPP_exe.sh -i path/to/BLUPPP_parameters.par [additional options]`
