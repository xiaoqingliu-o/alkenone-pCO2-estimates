Alkenone-pCO2 estimates aim to compute atmospheric CO2 levels and associated uncertainties, as published by Zhang et al. (2020), an accepted manuscript at Geochimica et Cosmochimica Acta (https://www.sciencedirect.com/science/article/pii/S0016703720303069). This repository contains MATLAB code to compute pCO2 and an excel template of the input data to help users prepare the data. 

pCO2 and its uncertainty are calculated using a Monte Carlo approach by 10,000 interactions of randomly sampling the input parameters within their 2 standard deviations of the mean (shown in MATLAB code) assuming normal distributions of these parameters, including sea surface temperature, delta13C of C37:2 alkenones, delta13C of carbonates, coccolith mean length, salinity, and pH (shown in the excel template).

The main script is CO2_Estimate.m and it requires the input data of each parameter within the same column as that in the excel template. Please refer to the excel file named as Template-input data. When the input data is ready, import data to the MATLAB Program as ‘Numeric Matrix’, and then rename the xxxx.mat file as data.mat.

Note: delta13C of carbonates are derived from delta13C of planktonic foraminiferal shells.
