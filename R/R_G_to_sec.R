###############################################################
# Script with article: 
# (Title here)
# Corresponding author Ruben Brondeel (ruben.brondeel@umontreal.ca)
# 
# This script calculates the Python version of activity counts, 
# The functions are a translation of the Matlab function presented in 
# Brond JC, Andersen LB, Arvidsson D. Generating ActiGraph Counts from 
# Raw Acceleration Recorded by an Alternative Monitor. Med Sci Sports Exerc. 2017.
# 
# R version 3.5.1 (run in RStudio Version 1.1.453)
###############################################################

library(data.table)
library(seewave)
library(signal)

# predefined filter coefficients 
A_coeff = c(1,-4.1637,7.5712,-7.9805,5.385,-2.4636,0.89238,0.06361,-1.3481,2.4734,-2.9257,2.9298,-2.7816,2.4777,-1.6847,0.46483,0.46565,-0.67312,0.4162,-0.13832,0.019852)
B_coeff = c(0.049109,-0.12284,0.14356,-0.11269,0.053804,-0.02023,0.0063778,0.018513,-0.038154,0.048727,-0.052577,0.047847,-0.046015,0.036283,-0.012977,-0.0046262,0.012835,-0.0093762,0.0034485,-0.00080972,-0.00019623)


pptrunc = function(data, max_value){
  
  # Saturate a vector such that no element's absolute value exceeds max_value.
  #   :param data: a vector of any dimension containing numerical data
  #   :param max_value: a float value of the absolute value to not exceed
  #   :return: the saturated vector
  
  outd = ifelse(data > max_value, max_value, data)
  return(ifelse(outd < - max_value, - max_value, outd))
}

trunc = function(data, min_value){
  
  # Truncate a vector such that any value lower than min_value is set to 0.
  #   :param data: a vector of any dimension containing numerical data
  #   :param min_value: a float value the elements of data should not fall below
  #   :return: the truncated vector
  
  return(ifelse(data < min_value, 0, data))
}

runsum = function(data, len, threshold){
  
  # Compute the running sum of values in a vector exceeding some threshold within a range of indices.
  # Divides the data into len(data)/length chunks and sums the values in excess of the threshold for each chunk.
  #   :param data: a 1D numerical vector to calculate the sum of
  #   :param len: the length of each chunk to compute a sum along, as a positive integer
  #   :param threshold: a numerical value used to find values exceeding some threshold
  #   :return: a vector of length len(data)/length containing the excess value sum for each chunk of data
  
  N = length(data)
  cnt = ceiling(N/len)
  rs = rep(0, cnt)
  
  for(n in 1:cnt){
    for(p in (1+len*(n-1)):(len*n)){
      if(p < N & data[p] >= threshold){
        rs[n] = rs[n] + data[p] - threshold
      }
    }  
  }
  return(rs)
}

# actual function to calculate counts from an array of raw accelerometer data
counts = function(data, filesf, B=B_coeff, A=A_coeff){
  
  # Get activity counts for a set of accelerometer observations.
  # First resamples the data frequency to 30Hz, then applies a Butterworth filter to the signal, then filters by the
  # coefficient matrices, saturates and truncates the result, and applies a running sum to get the final counts.
  #   :param data: the vertical axis of accelerometer readings, as a vector
  #   :param filesf: the number of observations per second in the file
  #   :param a: coefficient matrix for filtering the signal
  #   :param b: coefficient matrix for filtering the signal
  #   :return: a vector containing the final counts
  
  deadband = 0.068
  sf = 30
  peakThreshold = 2.13
  adcResolution = 0.0164
  integN = 10
  gain = 0.965
  
  if(filesf>sf) {
    datares = resamp(data, filesf, sf, 'matrix')
  }
  
  datab = bwfilter(datares, f = sf, n=4, from = 0.01, to = 7, bandpass = TRUE)
  
  B = B * gain
  
  fx8up = filter(B, A, datab)
  
  fx8 = pptrunc(fx8up[seq(1, length(fx8up), 3)], peakThreshold) #downsampling is replaced by slicing with step parameter
  
  out = runsum(floor(trunc(abs(fx8), deadband)/adcResolution), integN, 0)
  return(out)
}



main = function(file, folderInn, folderOut, filesf){
  # main function to run the previous functions over a list of data files
  # :param file: file name (.csv)
  # :param folderInn: directory of input (raw accelerometer data in g-units)
  # :param folderOut: directory of outcome (counts per second)
  # :param filesf: sampling frequency of data collection
  # :return none (instead writes a .Rdata file)
  
  fileInn = paste0(folderInn, file)

  # file name of activity count file (.Rdata) 
  fileR = gsub('.csv', '.Rdata', file)
  fileOut = paste0(folderOut, fileR)

  # read raw data
  dt = fread(fileInn, showProgress=FALSE)
  
  # calculate counts
  cx = counts(dt[, x], filesf)
  cy = counts(dt[, y], filesf)
  cz = counts(dt[, z], filesf)
  print('     Counts calculated')
  
  # combine counts data table
  counts = data.table(x = cx,  y = cy, z = cz)

  # save 'counts' object in a .Rdata file
  save(list = 'counts', file = fileOut)
  print('     Finished')
  
}


###############################################
# execute
###############################################

path = 'your_path/'
folderInn = paste0(path, "raw_accelerometer_data/")
folderOut = paste0(path, "activity_count_sec_R/")
files = dir(folderInn)

for(file in files){
  main(file, folderInn, folderOut, filesf = 50)
}




