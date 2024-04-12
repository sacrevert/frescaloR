# frescaloR
FRESCALO implemented in R

FRESCALO is a method to correct biological recording data for variation in recorder effort. The method is described in Hill, MO (2012) Local frequency as a key to interpreting species occurrence data when recording effort is not known, Methods in Ecology and Evolution, doi: 10.1111/j.2041-210X.2011.00146.x

The original method is a FORTRAN program. 
This R script takes the same input files as FRESCALO:   

1. a weights file (defining a neighbourhood around a focal region)
2. a species file giving the species observed at each region, and the date of the observation

More information about FRESCALO at http://www.brc.ac.uk/biblio/frescalo-computer-program-analyse-your-biological-records  and https://github.com/BiologicalRecordsCentre/sparta/wiki/frescalo

## Updates Apr. 2024 by Oli Pescott
- Added calculation of standard deviations and fixed some bugs encountered with sparse data.
- Included dataset and weights from [sparta](https://github.com/biologicalRecordsCentre/sparta) used to check all outputs. Matched against original fortran and this independent R coding based on direct translation from fortran: https://github.com/sacrevert/fRescalo 
- `unicorn_TF.rda` includes the results from running the original fortran code on the test data `clusterTestDat.csv` with the weights in `GB_LC_Wts.txt`.
