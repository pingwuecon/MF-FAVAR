This zip file contains the Matlab code for replicating the results 
in Koop et al (2023). Specifically, it contains two parts:

The file Data_Source.xlsx provides the details of the underlying data sources. 

Most of these data are freely downloadable from the internet. In the case of the CBI and PMI survey data and the Oil Price data used in this paper, these are available with a subscription to the relevant data service. 

In the case of the VAT data, these were accessed via a data sharing agreement with the Office for National Statistics.

In CSV form: There are 12 CSV. Each is for one region. The predictors are in Table 1 (main paper).

VAT.mat: This is the VAT data. Note VAT data is not publicly available.

Quarterly.mat: This is the quarterly data that we use in the main paper (Oil price is empty since the data is from subscription, although similar oil price series are freely available online).

Part II: codes

VBFAVARforecastingALACP_VAT.m: This is the main code to estimate the model. It allows for switching the method to estimate the factors (EMPCA, TW, TP).

This code is free to use for academic purposes only, and with the following attribution: 

Gary Koop, Stuart McIntyre, James Mitchell, Aubrey Poon, and Ping Wu (2023). Incorporating Short Data into Large Mixed Frequency VARs for Regional Nowcasting. 
Journal of the Royal Statistical Society, Series A, Forthcoming.

This code comes without technical support of any kind. It is expected to reproduce the results reported in the paper. Under no circumstances will the author be held responsible for any use (or misuse) of this code in any way.

For questions, please contact Ping Wu: ping.wu@strath.ac.uk