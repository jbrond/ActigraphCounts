# ActigraphCounts
This repository holds the code for generating the ActiGraph physical activity counts from acceleration measured with an alternative device. The implementation was described in the publication:

Generating ActiGraph Counts from Raw Acceleration Recorded by an Alternative Monitor, MSSE 2017 November

https://www.ncbi.nlm.nih.gov/pubmed/28604558

The Open source device used in the study was the Axivity AX3 which can be found here.

https://axivity.com/product/ax3

The matlab code is easy to use. The filter coefficients (A and B) are stored in the agcoefficeints.mat file and the generation of ActiGraph counts is done by invoking the agfilt function.

Example:

counts = agfilt(raw_acceleration_data,sampling_frequency,B,A);

Update: 03/06-2019

R and Python scripts has been added to the repository.
These scripts were devloped thanks to Ruben Brondeel (ruben.brondeel@umontreal.ca)


Jan Christin Brønd
