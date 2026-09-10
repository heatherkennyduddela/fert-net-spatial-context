# fert-net-spatial-context



This file contains details about the data processing and analysis of extra-pair fertilization networks for North American barn swallows in Colorado in 2021-2023.



These results are for the preprint titled "Spatial context shapes fertilization networks at both landscape and local scales in wild barn swallows".

Article authors: Heather V. Kenny-Duddela, Drew R. Schield, Zachary M. Laubach, Iris I. Levin, Kayleigh P. Keller, Rebecca J. Safran



The steps described require the following software:



* R and RStudio





## Contents



### 01 paternity assignment

This directory contains subdirectories for each year of paternity data. The scripts assign genetic sires using lcMLkin kinship output, classify fertilizations by type (within-pair, extra-pair) and calculate the total number of unknown (unsampled) sires in each year.



### 02 timing and synchrony

This directory contains scrips for calculating dyad synchrony and Rest-of-Colony Synchrony Index (RoCSI).



### 03 make full data table

This directory contains a script for combining data about distances between nests, timing of breeding, and fertilizations into a table where each row is a dyad. 



### 04 exploratory plots

This directory contains scripts for adding phenotype data to the main table, and calculating summary statistics on the data that is used for modeling. 


### 05 distance intx categorical

This directory contains scripts for doing final filtering steps on the data and testing support for distance as a categorical variable. 



### 06 explore temporal overlap
This directory contains a script for testing the predictive power of different metrics for assessing temporal overlap in breeding activity. 

### 07 between site ferts
This directory contains scripts for categorizing fertilizations as within-pair or extra-pair including unsampled sires. There is also a script for plotting sibling categories for offspring from unsampled sires, calculating summary statistics for the landscape-level model, and fitting a Generalized Linear Mixed-Effect Model to test for an effect of distance on between-structure fertilizations. 

### 08 within site ferts
This directory contains the script for fitting a Generalized Linear Mixed-Effect Model for same-structure fertilizations, along with two helper scripts for calculating confidence intervals for the GLMM and generating custom new data for plotting results. There is also a script that calculates summary statistics for the dyad data used in the model. 

