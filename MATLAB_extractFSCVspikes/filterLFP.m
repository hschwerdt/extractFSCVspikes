function filteredData=filterLFP(data,samplingRate,filterBand)
%forward/backward filter, butterworth
%data is the data to be sampled
%samplingRate is the sampling rate
%filterBand is the frequency band [LP HP] in Hz
filter_band=filterBand;
order = 4;
buttertype = '';
freqlim=filter_band;
result=0;
if freqlim(1)==0
    buttertype='low';
    freqlim=freqlim(2);
end
if length(freqlim)>1
%Below if statement added sometime after 2021 may be, error if supply [0
%100] so added conditional above this...
if freqlim(2)==inf
    buttertype='high';
    freqlim=freqlim(1);
end
end
fwdflag = true;
revflag = true; %Always to backward/forward filter
samplePeriod=1/samplingRate;
% 'butter' requires freqs spec'd with 1.0 corresponding to half the sample rate. 
%freqlim must be 0.0 < freqlim < 1.0, with 1.0 corresponding to half the sample rate.
%freqlim(1) cannot be zero for default argument of bandpass
freqlim = freqlim*samplePeriod*2;

if isempty(buttertype)
    [z, p, k] = butter(4, freqlim);
else
    [z, p, k] = butter(4, freqlim, buttertype);
end
[sos,g]=zp2sos(z,p,k);
h2=dfilt.df2sos(sos,g);
if revflag
    result = filter(h2, data(end:-1:1));%Reverse filter
    if ~fwdflag
        result = result(end:-1:1);
    end
end
if fwdflag
    if revflag
        result = filter(h2, result(end:-1:1));%Forward filter by reversing reversed data above
    else
        result = filter(h2, data(:));
    end
end


filteredData=result;
end