# coding=utf-8

'''
Script with article: 
Converting raw accelerometer data to activity counts 
using open source code in Matlab, Python, and R ? a comparison to ActiLife activity counts
Corresponding author Ruben Brondeel (ruben.brondeel@umontreal.ca)

This script calculates the Python version of activity counts, 
The functions are a translation of the Matlab function presented in 
Br?nd JC, Andersen LB, Arvidsson D. Generating ActiGraph Counts from 
Raw Acceleration Recorded by an Alternative Monitor. Med Sci Sports Exerc. 2017.

Python (3.6); run in Eclipse (Oxygen.3a Release (4.7.3a))

'''

% read in libraries needed for the functions
import math, os
import numpy as np
from scipy import signal
import pandas as pd
import resampy

##predefined filter coefficients, as found by Jan Brond
A_coeff = np.array(
    [1, -4.1637, 7.5712,-7.9805, 5.385, -2.4636, 0.89238, 0.06361, -1.3481, 2.4734, -2.9257, 2.9298, -2.7816, 2.4777,
     -1.6847, 0.46483, 0.46565, -0.67312, 0.4162, -0.13832, 0.019852])
B_coeff = np.array(
    [0.049109, -0.12284, 0.14356, -0.11269, 0.053804, -0.02023, 0.0063778, 0.018513, -0.038154, 0.048727, -0.052577,
     0.047847, -0.046015, 0.036283, -0.012977, -0.0046262, 0.012835, -0.0093762, 0.0034485, -0.00080972, -0.00019623])

def pptrunc(data, max_value):
    '''
    Saturate a vector such that no element's absolute value exceeds max_abs_value.
    Current name: absolute_saturate().
      :param data: a vector of any dimension containing numerical data
      :param max_value: a float value of the absolute value to not exceed
      :return: the saturated vector
    '''
    outd = np.where(data > max_value, max_value, data)
    return np.where(outd < -max_value, -max_value, outd)

def trunc(data, min_value):
  
    '''
    Truncate a vector such that any value lower than min_value is set to 0.
    Current name zero_truncate().
    :param data: a vector of any dimension containing numerical data
    :param min_value: a float value the elements of data should not fall below
    :return: the truncated vector
    '''

    return np.where(data < min_value, 0, data)

def runsum(data, length, threshold):
    '''
    Compute the running sum of values in a vector exceeding some threshold within a range of indices.
    Divides the data into len(data)/length chunks and sums the values in excess of the threshold for each chunk.
    Current name run_sum().
    :param data: a 1D numerical vector to calculate the sum of
    :param len: the length of each chunk to compute a sum along, as a positive integer
    :param threshold: a numerical value used to find values exceeding some threshold
    :return: a vector of length len(data)/length containing the excess value sum for each chunk of data
    '''
    
    N = len(data)
    cnt = int(math.ceil(N/length))

    rs = np.zeros(cnt)

    for n in range(cnt):
        for p in range(length*n, length*(n+1)):
            if p<N and data[p]>=threshold:
                rs[n] = rs[n] + data[p] - threshold

    return rs

def counts(data, filesf, B=B_coeff, A=A_coeff):
    '''
    Get activity counts for a set of accelerometer observations.
    First resamples the data frequency to 30Hz, then applies a Butterworth filter to the signal, then filters by the
    coefficient matrices, saturates and truncates the result, and applies a running sum to get the final counts.
    Current name get_actigraph_counts()
    :param data: the vertical axis of accelerometer readings, as a vector
    :param filesf: the number of observations per second in the file
    :param a: coefficient matrix for filtering the signal, as found by Jan Brond
    :param b: coefficient matrix for filtering the signal, as found by Jan Brond
    :return: a vector containing the final counts
    '''
    
    deadband = 0.068
    sf = 30
    peakThreshold = 2.13
    adcResolution = 0.0164
    integN = 10
    gain = 0.965

    if filesf>sf:
        data = resampy.resample(np.asarray(data), filesf, sf)

    B2, A2 = signal.butter(4, np.array([0.01, 7])/(sf/2), btype='bandpass')
    dataf = signal.filtfilt(B2, A2, data)

    B = B * gain

    #NB: no need for a loop here as we only have one axis in array
    fx8up = signal.lfilter(B, A, dataf)

    fx8 = pptrunc(fx8up[::3], peakThreshold) #downsampling is replaced by slicing with step parameter

    return runsum(np.floor(trunc(np.abs(fx8), deadband)/adcResolution), integN, 0)


def main(file, folderInn, folderOut, filesf):
    
    '''
    Creates activity counts per second from raw acceleromter data (g-units) 
    This function:
      - reads in data into a pandas dataFrame
      - Calculates activity counts per axis
      - combines the axis in a pandas dataFrame
    :param file: file name of both input and output file
    :param folderInn: directory with input files, containing raw accelerometer data
    :param folderOut: directory with out files, containing activity counts.
    :param filesf: sampling frequency of raw accelerometer data
    :return: none (writes .csv file instead)
    '''
    
    # read raw accelerometer data
    fileInn = folderInn + file
    dt = pd.read_table(fileInn, delimiter=',', header=0)
    
    # calculate counts per axis
    c1_1s = counts(dt.x, filesf)
    c2_1s = counts(dt.y, filesf)
    c3_1s = counts(dt.z, filesf)
    
    # combine counts in pandas dataFrame
    c_1s = pd.DataFrame(data = {'axis1' : c1_1s, 'axis2' : c2_1s, 'axis3' : c3_1s})
    c_1s = c_1s.astype(int)

    # write to output folder
    fileOut = folderOut + file
    c_1s.to_csv(fileOut, sep=',', index = False) 
    

#-------------------------------------------------------------------------------
# Execute in loop
#-------------------------------------------------------------------------------

path = 'your_path/'
folderInn = path + "raw_accelerometer_data/"
folderOut = path + "count_sec_python/" 
files = os.listdir(folderInn)
      
# Loop over .csv files
[main(file, folderInn, folderOut, filesf = 50) for file in files]


