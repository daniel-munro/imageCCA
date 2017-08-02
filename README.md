# imageCCA scripts

These scripts were used in the study "Joint analysis of gene expression levels and histological images identifies genes associated with tissue morphology" (under submission).

### `extract_image_from_SVS.ijm`

This is an ImageJ macro that reads images stored in SVS format and extracts a sampled representative region while attempting to avoid empty (background) space. It requires the Bio-Formats plugin for ImageJ.


### `prepare_CCA.R`

This R script contains functions for preparing the image feature and expression data for CCA. This involves lining up observations corresponding to the same sample, removing any features with zero variance, and saving the ordered sample and gene names for analysis of the CCA output.


### `run_CCA.R`

This R script contains functions for running sparse CCA (using the `PMA` package) and extracting relevant output. This output includes the canonical variables, which are the image and expression samples transformed into the CCA component space, and the sparse CCA coefficients, which are the nonzero factors used for the transformation that indicate which image and expression features contributed to each CCA component. At the end of the script is code for the permutation analysis used in the study.


### `CAE/CAE_[dataset].py`

These python scripts contain code for building and training a convolutional autoencoder (CAE) using the `keras` package. It is recommended to run these scripts using GPUs for faster execution.


### `CAE/CAE_output_[dataset].py`

These python scripts contain code for loading the trained CAE network made by `CAE_[dataset].py` and producing image representations to be used for CCA.


### `CAE/CAE_classify_[dataset].py`

These python scripts are similar to `CAE_[dataset].py` but train a CAE that is supervised on the sample classes. It is recommended to run these scripts using GPUs for faster execution.


### `CAE/CAE_classify_output_[dataset].py`

These python scripts are similar to `CAE_output_[dataset].py` but load the trained supervised CAE network made by `CAE_classify_[dataset].py`.