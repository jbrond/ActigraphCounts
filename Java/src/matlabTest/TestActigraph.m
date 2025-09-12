%%Matlab file to test whether Java port and Brond ActigraphCounts matlab implementation give the same result
%%Spoiler, the results are not exactly the same because the java
%%implementation does not initialise filter states similarly to matlab

fclose all;
close all;
clear all;
clc;

addpath('functions');
addpath('../../../Matlab');	%Brond Matlab implementation
load('../../../Matlab/agcoefficients.mat');	%Load the filter coefficients

%writetable(array2table(num2cell([B;A])),'coeffs.csv');

if exist('OCTAVE_VERSION', 'builtin') 
	pkg load signal	%Filter functions (butter) are in the signal package
	javaaddpath('../../build/libs/BrondActigraphCounts-1.0.jar')
else
	javaclasspath('../../build/libs/BrondActigraphCounts-1.0.jar'); %Add the java implementation into classpath
end
%Read a sample file recorded with a smartphone worn in front pocket. The file contains data from roughly a 25 min jog.
dataFolder = '../res/';
accFile = getFilesAndFolders([dataFolder]);
data = readLog([dataFolder accFile(1).name]);
sRate = round(1/median(diff(data.data(:,1)./1000))); %Time stamps are milliseconds
acc = data.data(:,2:4)./9.81;	%Acceleration in g
resultant = sqrt(sum(acc.^2,2));

%Full java analysis, resampling different from agfilt (agfilt uses
%resample, Java uses linear interpolation).
javaAGC = javaObject('timo.jyu.ActigraphCounts',resultant,sRate);
countsFJ = javaMethod('getCounts',javaAGC);

%Anti-alising filter
b = javaMethod('getB',javaAGC);
a = javaMethod('getA',javaAGC);
aFiltered = filtfilt(b,a,resultant);

if size(aFiltered,2) > size(aFiltered,1)
	aFiltered = aFiltered';
end
countsMat = agfilt(aFiltered,100,B,A);


%Replicate Matlab analysis
deadband = 0.068;
peakThreshold = 2.13;
adcResolution = 0.0164;

agResampledData = resample(aFiltered,30,sRate);
javaAGC = javaObject('timo.jyu.ActigraphCounts',30);
agFilt =javaMethod('getAGFiltered',javaAGC,agResampledData);
down_10 =javaMethod('down10',javaAGC,agFilt);
trunc8bit =javaMethod('pptruncDeadband8bit',javaAGC,down_10,peakThreshold,deadband,adcResolution);
countsJ =javaMethod('getCounts',javaAGC,trunc8bit);

difference = countsMat-countsFJ;
disp(sprintf('Full java counts diffMax %f diffMin %f meanDiff %f',max(difference),min(difference),mean(difference)));


difference = countsMat-countsJ;
disp(sprintf('Matched java counts diffMax %f diffMin %f meanDiff %f',max(difference),min(difference),mean(difference)));


t = ([1:length(resultant)]-1)./sRate;
indices = 1:sRate:length(t);

fh = figure;
ah =  []; cnt = 0;
subplot(2,1,1)
plot(t,resultant,'k','linewidth',3);
title('Resultant acceleration');
xlabel('Time [s]');
ylabel('[g]');
cnt = cnt+1; ah(cnt) = gca();
subplot(2,1,2)
plot(t(indices),countsMat,'k','linewidth',3);
hold on;
plot(t(indices(1:length(countsJ))),countsJ,'r','linewidth',3);
plot(t(indices(1:length(countsJ))),countsFJ,'b','linewidth',3);
xlabel('Time [s]');
ylabel('one second count sums [count]');
cnt = cnt+1; ah(cnt) = gca();
linkaxes(ah,'x');
waitfor(fh);
