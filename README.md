# ACM_Results_ExcessEnthalpy
Results and code from the Masters thesis 'Array completion methods for thermodynamic property estimation' by FS Middleton (2023).

## Folders present
* The Code folder contains all code used in the project. 
* The SmallResults folder contains all results for the smallest array completed, a subset of the array at 298.15K. This was the toy problem for this work.
* The Results folder contains all results for the MCM predictions, and some results for the 3ACM predictions.
* The UNIFAC_Dortmund_preds folder contains all UNIFAC (Dortmund) predictions, using the same figures available in Results. 

## Navigation of results
The way to navigate the Results folder is by using the search function from 'Go to file' on GitHub. 

### High level metrics
There are results for all arrays completed that show a high level of the error metrics. They are named according to the following convention:
* 'filling scheme'-'Plot generated'

The filling scheme can be: 
* uni
* tri 
* mix

The plot generated can be: 
* ErrorMetricsPerFracObserved
* ErrorMetricsPerNumCompounds
* ErrorMetricsPerTemp

For example: uniErrorMetricsperTemp


### Metrics per isothermal array
There are also figures comparing the results of each filling scheme for an isothermal array. These are labelled according to the convention:
*  'Temperature'-'r1'-'Plot generated'

The plot generated can be:
* 2Hist-FillSchemes
* 2Hist-FillSchemes+UNIFAC
* CompareFillMethod

For example: 
298.15-ComapreFillMethod 
298.15K-2Hist-FillSchemes

The Excel file called 'TypesOfMixtures' also contains the errors per type of mixture in a type of system, these being type 1 to 5 and describing the shapes of the systems. All results are held in tables in this file.
These results are all for the optimal rank of the array using a given filling scheme and _t_ = 40%.

### All results of interest
The MCM results for each isothermal array and filling scheme are named according to the convention: 
* 'Temperature'-'filling scheme'-'r1'-'Plot generated'
  
The ACM results are named according to the convention:
* 'Number of temperature slices'-'Temperature of the results'-'filling scheme'-fn3='r3'-fn1=-'r1'-'Plot generated'

The names of the plots generated are shortened descriptions:
* ErrorMetrics = the error metrics of the SMSE, wSMSE, AARD, and wAARD for each rank tested
* Parity = Parity plot
* 2Hist = 2-way histogram - errors
* 3Hist = 3-way histogram - errors and ARDs
* Err-3way = error per prediction in a 3-way cubic plot
* UppervsLower = most accurate predictions from the Upper or Lower triangular arrays
* UNIFACvsMCM = most accurate predictions from the MCM or UNIFAC
* UNIFACvs3ACM =  most accurate predictions from the 3ACM or UNIFAC
* MCMvs3ACMvxUNIFAC = most accurate predictions from the MCM, 3ACM, or UNIFAC
* BAC = error metrics for each BAC group including the SMSE, wSMSE, AARD, and wAARD
* Func = error metrics for each combination of functional groups including the SMSE, wSMSE, AARD, and wAARD
* MixturesInterest = predictions for the mixtures of interest found in a set of results
* Mixtures = predictions of all mixtures that were observed in the array. LOOCV results

The BAC and Func figures also differentiate into the 4 different recorded error metrics:
* SMSE
* wSMSE
* AARD
* wAARD


And the Mixtures figures have a number following them. 

Examples of these names include:
* 293.15K-mix-6-Parity
* 293.15K-mix-6-BAC-AARD
* 293.15K-mix-6-Mixtures0
