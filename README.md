# A continuous amyloid-beta CSF/PET imbalance model to capture Alzheimer's disease heterogeneity 

This code allows to estimate a continuous model for soluble/aggregated amyloid-beta imbalance. 

The function "mincurve1" iteratively fits a hyperbolic regression model between baseline CSF-Aβ42 and Aβ-PET CL data, by minimizing the sum of the Euclidean distance of the experimental points to the fitted line. 

The function "mincurve2" allows to derive two subject-specific measures of:
(i) the relative imbalance between soluble (CSF) and aggregated (PET) Aβ, termed Aβ-aggregation 
(ii) the extent of Aβ-pathology, termed Aβ-severity. 

Aβ-aggregations scores were calculated as the difference between the observed and predicted data point (i.e., the standardized Euclidean distance of the observed data point to the fitted line).  A positive score denotes a higher PET CL value than expected for a given value of CSF-Aβ42 (i.e., more aggregated relative to soluble Aβ), whereas a negative score denotes lower CSF-Aβ42 values than expected for a given value of PET CL (i.e., more soluble relative to aggregated Aβ). In turn, Aβ-severity indicates where along the hyperbolic regression line the data point is located (i.e., the standardized Euclidean distance between the individual predicted value and the median of all predicted values). A positive score reflects more advanced Aβ-pathology and a negative score reflects less advanced Aβ-pathology. Figure-1A shows several hypothetical data points and their respective Aβ-aggregation and Aβ-severity scores.





More details are provided in Mastenbroek, Sala et al., A continuous amyloid-beta CSF/PET imbalance model
to capture Alzheimer's disease heterogeneity (2024). Neurology.


## Set-up

Matlab installation is required to run the functions mincurve1 and mincurve2. These functions have been tested on Matlab 2018b and more recent releases.


## How to use

The user should run functions mincurve1, followed by mincurve2. 

In order to run mincurve1, the user should provide some provisional parameters to initialize the hyperbolic model. We recommend plotting the observed PET against CSF data and then visually estimate some plausible parameters for the hyperbolic model, according to the formula:

$$ CSF =a+ { b*PET \over PET-c} $$

   
## How to acknowledge

The code is provided under GLP-3.0 license. The code was authored by Juan Doming Gispert, BBRC-Foundation Pasqual Maragall (mincurve1; mincurve2) and Arianna Sala, Karolinska Institutet (mincurve2). 

Any publication based on this code should cite: Mastenbroek, Sala et al., A continuous amyloid-beta CSF/PET imbalance model
to capture Alzheimer's disease heterogeneity (2024). Neurology.
