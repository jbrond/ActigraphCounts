A Java port of Jan Brond's Matlab code to produce Actigraph counts from raw accelerometer records.

Written by Timo Rantalainen 2019 (tjrantal at gmail dot com).

This Java port is released to the public domain with a CC0-BY Creative Commons Attribution license. Butterworth filter coefficient calculations were ported from http://www-users.cs.york.ac.uk/~fisher/mkfilter/ and filtering and zero-lag filtering were ported from Octave source code https://github.com/greenm01/forc/blob/master/sandbox/filter.m and https://searchcode.com/codesearch/view/9521190/. The licenses that apply to those codes apply to this code as well where appropriate.

Depends on Apache Math 3.6.1 library

Use gradle to build (gradle build in this folder). Gradle will download the Apache Math library from maven repository. Once the jar is built, test with Matlab (src/matlabTest/TestActigraph.m) to confirm you get the same results from the Matlab implementation and this Java port. Note that the results are not exactly the same if full java analysis is applied because the java implements resampling with linear interpolation, whereas Brond's matlab code uses matlab resample command.
