function counts = agfilt(data,filesf,B,A)

deadband = 0.068;

sf = 30;
peakThreshold = 2.13;
adcResolution = 0.0164;
integN = 10;
gain = 0.965;

if (filesf>sf)
    dataf = resample(data,sf,filesf);
else
    %Aliasing Filter
    [B2,A2] = butter(4,[0.01 7]./(sf/2));
    dataf = filtfilt(B2,A2,data);
end

S = size(dataf);

B = B * gain;

for n=1:S(2)
    
    fx8up = filter(B,A,dataf(:,n));
    fx810 = downsample(fx8up,3);
	
    fx8 = pptrunc(downsample(fx8up,3),peakThreshold);
    
    counts(:,n) = runsum(floor(trunc(abs(fx8),deadband)./adcResolution),integN,0);
end
    