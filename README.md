# Documentation for MutPred2 repository
## Requirements
* MATLAB R2017b or earlier (not tested on earlier versions)
* PSI-BLAST executable and associated files (bare-bones version included here but full version available [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.18/)
* At least 4 GB RAM
* At least 7 GB disk space (code + [data and model files](http://mutpred.mutdb.org/model_and_data_files.tar.gz))

## Description of files and sub-directories
* *mutpred2.m*: the main MATLAB function that makes predictions (calls the different functions in *all_functions*).
* *startup.m*: contains code to include the all files in the repo in MATLAB path (**edit file before using** - change the path to the location of the MutPred2 repo).
* *excluded_builtin_functions.txt*: all MATLAB functions that are part of the Bioinformatics or Statistics toolboxes that cannot be included under MATLAB's licensing restrictions.
* *all_functions* - directory containing all helper functions called by *mutpred2.m*.
* *blast-2.2.18* - legacy version of the PSI-BLAST executable that is compatible with MutPred2 (*bin* subdirectory) and data files associated with PSI-BLAST (*data* subdirectory).

## Citation
If you use code from this repository, the web-server or the standalone executable in your work, please cite:

**MutPred2: inferring the molecular and phenotypic impact of amino acid variants**

Vikas Pejaver, Jorge Urresti, Jose Lugo-Martinez, Kymberleigh A. Pagel, Guan Ning Lin, Hyun-Jun Nam, Matthew Mort, David N. Cooper, Jonathan Sebat, Lilia M. Iakoucheva, Sean D. Mooney, Predrag Radivojac
*bioRxiv* 134981; doi: [https://doi.org/10.1101/134981](https://doi.org/10.1101/134981)

## Additional notes
### Data and code
* The code here corresponds to MutPred2.0, the version used for all analyses in the [bioRxiv](https://www.biorxiv.org/content/10.1101/134981v1) paper (and the main text of the published paper).
* There is currently no internal organization within the *all_functions* directory, as these needed to be placed in one location when compiling the standalone executable.
* If you are looking for MutPred2's training set, a subset of the training data for MutPred2 is available [here](http://mutpred.mutdb.org/wo_exclusive_hgmd_mp2_training_data.txt).

### License
* We cannot share MATLAB-related files or any of MATLAB's proprietary functionsand have excluded them to the best extent possible. However, please be aware of MATLAB's licensing restrictions when modifying and/or redistributing any code.
* Although MutPred2 requires an unsupported legacy version of PSI-BLAST, please refer to the BLAST suite's license before modifying and/or redistributing.

