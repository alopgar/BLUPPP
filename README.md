[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![Last version](https://img.shields.io/github/tag/alopgar/BLUPPP.svg)](https://img.shields.io/github/tag/alopgar/BLUPPP.svg)

# BLUPPP
This pipeline processes a file of SNP genotypes and prepares it for its use with [BLUPF90](http://nce.ads.uga.edu/wiki/doku.php) or [ASRGenomics](https://asreml.kb.vsni.co.uk/asrgenomics-download-success/).

### Processing steps included:
#### 1) Add new genotypes:
Two SNP files might be combined if new animals are to be included.
#### 2) Custom processing:
Removal of non-SNP and non-animal ID columns, or other custom operations can be specified in the parameter file.
#### 3) Combination or removal of duplicated animals:
If 2 copies of an animal ID are detected, BLUPPP will combine their SNP information, if the correspondence is higher than 95%, or remove both duplicates if the correspondence is lower than 95%. If 3 or more duplicates are detected, BLUPPP wil warn the user to check these IDs.
#### 4) PLINK custom filters:
An additional script with plink processing (`plink_filt.sh`) is included in BLUPPP, which can be modified according to the user's preferences.
#### 5) Conversion of SNP file to BLUPF90/ASRGenomics formats:
This step includes:
a) The removal/keeping of headers (i.e., SNP names or codes).  
b) The mutation of NA codes to 5.  
c) The collapse or spaced separation of SNP values.  
Following the requirements of BLUPF90/ASRGenomics, respectively.

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
