# New version

The code has been retooled and a new version can be found at: https://github.com/CancerRxGene/gdscIC50 


# Estimating IC50 values using a mixed effects model.

_djvMixedIC50_ allows you to estimate IC50 and AUC values using a non-linear mixed model. It estimates all the responses simultaneously and capitalizes on the entire set of responses for each cell line to infer its shape. This significantly reduces the ambiguity of the cell line's shape and thus improves the estimate.

# etails
The package was initially designed for the GDSC compound sensitivity screen produced by the Wellcome Trust Sanger Institute. Their typical assay is a 9 point assay with 2-fold dilutions. However, all screen comprising sets of cell lines exposed to a series of compounds can benefit from this modeling approach.

As a design choice, the IC50 is not directly estimated on the concentration range, but on the dilution steps. In this way, the internal IC50 estimates are comparable in terms of response range. Internally the fold-dilution is constant at 2, other constants can be accommodated by using a custom predictor vector (see also _getXfromConcSeries_).

After the initialization of the data by _initIC50_, and the curve fitting of the data by _fitModel_, the data frames containing the resultant statistics are constructed by _gatherModelStats_. In this function, the internal IC50 estimates are translated to actual (log) concentration values based on the highest test concentration. The individual data matrices holding these statistics can be accessed the named list. For more information on how to format the data for a custom screen, please consult the relevant sections in _initIC50_.

The resultant matrices hold the natural logarithm of the estimated micro Molar IC50 concentration. The area under the curve (AUC) is calculated using the integral of the estimated model parameters for the dose-response curve and normalized by the total area. The numerically integrated version is designated AUCtrap, calculated using the trapezoid rule. The residual of the model fit is summarized by the root mean square error.

See also: Pharmacogenomics, May 2016, Vol. 17, No. 7, Pages 691-700, or http://www.ncbi.nlm.nih.gov/pubmed/27180993

# Examples
gDat      <- initIC50(szFileName='~/myFile.csv')

fmMod1    <- fitModel(gDat)

outStats  <- gatherModelStats(gDat=gDat,fmMod1=fmMod1)

dfIC50    <- outStats$IC50

dfRMSE    <- outStats$RMSE

dfAUC     <- outStats$AUC

dfAUCtrap <- outStats$AUCtrap

The input data looks as follows:

Grouped Data: y ~ x | drug/CL

  x          y     CL drug maxc
  
1 9 0.88079708 MC-CAR    1    2

2 8 0.79139147 MC-CAR    1    2

3 7 0.66075637 MC-CAR    1    2

4 6 0.50000000 MC-CAR    1    2

5 5 0.33924363 MC-CAR    1    2

6 4 0.20860853 MC-CAR    1    2

7 3 0.11920292 MC-CAR    1    2

8 2 0.06496917 MC-CAR    1    2

9 1 0.03444520 MC-CAR    1    2

