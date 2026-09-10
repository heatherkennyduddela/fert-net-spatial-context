# fert-net-spatial-context



This file contains details about the data processing and analysis of extra-pair fertilization networks for North American barn swallows in Colorado in 2021-2023.



These results are for the preprint titled "Spatial context shapes fertilization networks at both landscape and local scales in wild barn swallows".

Article authors: Heather V. Kenny-Duddela, Drew R. Schield, Zachary M. Laubach, Iris I. Levin, Kayleigh P. Keller, Rebecca J. Safran



The steps described require the following software:



* R and RStudio





## Contents

* [01\_paternity-assignment](#01_paternity-assignment)
* [02\_timing-and-synchrony](#02_timing-and-synchrony)
* [03\_make-full-data-table](#03_make-full-data-table)
* [04\_exploratory-plots](#04_exploratory-plots)
* [05\_distance-intx-categorical](#05_distance-intx-categorical)
* [06\_explore-temporal-overlap](#06_explore-temporal-overlap)
* [07\_between-site-ferts](#07_between-site-ferts)
* [08\_within-site-ferts](#08_between-site-ferts)



### 01 paternity assignment

This directory contains subdirectories for each year of paternity data. The scripts assign genetic sires using lcMLkin kinship output, classify fertilizations by type (within-pair, extra-pair) and calculate the total number of unknown (unsampled) sires in each year.



### 02 timing and synchrony

