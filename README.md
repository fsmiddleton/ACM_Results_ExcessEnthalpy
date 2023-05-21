# ACM_Results_ExcessEnthalpy
Results and code from the Masters thesis 'Array completion methods for thermodynamic property estimation' by FS Middleton (2023).

The Code folder contains all code used in the project. 

The SmallResults folder contains all results for the smallest array completed, a subset of the array at 298.15K. This was the toy problem for this work.

The Results folder contains all results for the MCM predictions, and some results for the 3ACM predictions. 

The way to navigate the Results folder is by using the search function from 'Go to file'. 
There are results for all arrays completed named according to the following convention:
* 'filling scheme'-'Plot generated'
The names of the plots generated can be: 
* ErrorMetricsPerFracObserved
* ErrorMetricsPerNumCompounds
* ErrorMetricsPerTemp

The filling schemes used can be: 
* uni
* tri 
* mix 


The MCM results are named according to the convention: 
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
* CompareFillMethod = more accurate predictions from either the uni, mix, or tri filling scheme
