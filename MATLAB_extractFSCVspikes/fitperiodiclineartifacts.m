function [artpeakids,fractionbad] = fitperiodiclineartifacts(ts, samples, thresh, maxblanks)
%2021-2023 Helen largely replaced Dan code (dg_fitperiodicartifacts4)
%similar to fitFSCVart, but customized for detecting 60/120 hz line noise
%Here the Input samples should be  the absolute value of the highpass filtered data 
% (300 hz cut off used and works well)> thresh is 5-8 times the mean of that filtered
%absolute value signal
%%Get midpoints of artifacts, fit 120hz peak between 60 hz and also create
%60 hz peaks around 15 or 30 hz peaks (e.g. if no peaks because of fscv
%artifact)
%Explicit nbefore and nafter of 3 and 10 samples now for interp
%Maxblanks should be 8 (short) for interval between detected xingids so
%that don't overestimate number of allpeakids in getperiodic peaks. For
%actual interpolating, want to overestimate for smaller ones so make it 16,
%but not too high.
%INPUTS
% ts: timestamps in seconds of all samples, same size as <samples>.
% samples: array of sampled values, one per timestamp.
% thresh: positive-going threshold to trigger interpolation.
%OUTPUTS
% artpeakids: the index into <samples> of the center of each identified
% artifact
% fractionbad: fraction of bad harmonicpeakid intervals (i.e., intervals >
% 60 hz period)
%   
%NOTES
%   First part to approximate spike detection based on threshold crossing
%   of an high-pass filtered signal (absolute value), very simialr to
%   fitFSCVart.m but modified based on different frequency content and
%   parameters of the 60/120 hz noise
%   The sample scanperiod is estimated based only on the first 10,000 samples.

artpeakids=[];
fnoise_0=60;    %frequency of fundamental noise 60 hz
fnoise_1=120;   %harmonic of noise 120 hz
sampleperiod = median(diff(ts(1:1e4))); %   The sample period is estimated based only on the first 10,000 samples., rate = 1/sampleperiod
intervalnoise=ceil(1/fnoise_0*1/sampleperiod);%60 hz noise interval
fractionbad=0;
maxwidth = 10;  %Max artifact width in samples = 10/30e3 = .3 ms
periodicitytol = 0.05; % max fractional deviation from the excpected period, don't check for 60 hz noise %MAYBE MAKE 5 %%%%% used in getperiodicpeaks, need tight

samples_orig=samples;   %Store origianl samples
samples=abs(samples);   %Get absolute value for xing crossing

% Find threshold crossings ("xings"), scanperiod, and start of first scanperiod.
% <xingidx> points to the first sample that is above <thresh> in each
% putative FSCV artifact.
isxing = [ false
    reshape( samples(2:end) >= thresh & samples(1:end-1) < thresh, ...
    [], 1 ) ]; %Only gets first threshold crossing because of logical ids
xingidx = reshape(find(isxing), [], 1 );
peakidx=zeros(size(xingidx));

badslope=[];
% Get peak ids from xing ids (i.e. find local max after xing)
for idx2 = 1:length(xingidx)
    slope = samples(xingidx(idx2)) - samples(xingidx(idx2) - 1);
    if slope < 0.02 * thresh
        % Bogus xing, skip:
        badslope=[badslope; xingidx(idx2)];
        continue
    end

    %Instead of looking for true art start, look for local peak
    if xingidx(idx2)+maxwidth>length(samples)
        %skip near end
        continue
    else
        arttrace=samples(xingidx(idx2):xingidx(idx2)+maxwidth);
        %Multiple ways to record artifacts, center or peak.
        % Center means looking for midpoint between positive going slope 
        % and negative going slope. 
        % Here we use peak, the maximum within max window past the positive
        % crossing xingidx. this works better than center. Similar to
        % fitFSCVart.m
        [peakmax,peakid]=max(arttrace);
        artpeak=peakid+xingidx(idx2)-1;
        peakidx(idx2)=artpeak;
    end
end
% Delete the entries for the bogus xings:
peakidx(peakidx == 0) = [];

%Now go through peakidx's (midpoints of artifacts) to check periodicity,    
% add artifacts in between stable harmonics,
%Get all 120 Hz peaks fit between those 60 harmonics found
%Max blank periods defines window of 60 hz or 120 hz intervals to skip filling of artifact
% ids as there may be periods in recording where source of 60 hz noise turned off or 
%larger srouces of interference that we don't want to unnecessarily interpolate over, 
%these periods are checked for in subseq programs that look for nonperiodic
% movement artifacts or pump related artifacts
sessiondur=length(samples)*sampleperiod;
fractfoundline=length(peakidx)/(sessiondur*120);    %Expected is sessiondur*120
if fractfoundline>.05
    %At least 10 percent of expected artifacts found based on just threshold
    %crossing, then continue
    %Max blank periods = 32 or 32*1/120 ~ 0.25 s, if 64 creates over fitting in session 125 for exmaple, 
    % 16 for interpolating with grand average of artifacts, and use 8 short 
    % for identifying artifact channels to avoid finding artifactual artifacts
    %Get timestamp ids for periodic 60/120 hz artifacts
    harmonicpeakids=getperiodicpeaks(sampleperiod,peakidx,periodicitytol,fnoise_0,fnoise_1,maxblanks);
    artpeakids=harmonicpeakids;

    %Find how many of the harmonic ids identified in algo match original
    %ids with +/-1 sample tolerance
    ismemberharmonics=find(ismember(harmonicpeakids,[peakidx-1; peakidx; peakidx+1;]));  %Find what harmonicpeakids contain original xingidx ids (with +/-1 sample offset for some slight differences in identifying peaks vs periodicity
    origxingsinharmonicids=[]; %original peakidx values also in harmonicpeakids
    if ~isempty(ismemberharmonics)
        origxingsinharmonicids=harmonicpeakids(ismemberharmonics);
    end
    %Find how many of the original peak ids identified based on xing crossing match 
    % newly identified harmonic ids with +/-1 sample tolerance
    ismemberoriginalxings=ismember(peakidx,[harmonicpeakids-1; harmonicpeakids; harmonicpeakids+1;]); %Logical ids of harmonic peaks identified in original peakidx values
    harmonicsinorigxings=[];
    harmonicsnotinorigxings=[];
    if ~isempty(ismemberoriginalxings)
        harmonicsinorigxings=peakidx(ismemberoriginalxings);
        %original peak ids that were not identified as harmonics in
        %getperiodicpeaks
        harmonicsnotinorigxings=peakidx(~ismemberoriginalxings);
        %IDs of peakidx not found in harmonicpeakids, to compare if a lot 
        % of offset artifacts may be should retain and include with the 
        % harmonics as these might represeent the 180 Hz artifacts seen eg session 83 and 83b
    end
    %fraction of original peak ids that were not identified as harmonics
    fractionnotacctinharm=length(harmonicsnotinorigxings)/length(peakidx);  
    disp(['fraction of peakidx thresh crossings originally identified and not ' ...
        'accounted for in harmonic ids = ' num2str(fractionnotacctinharm,5)]);
    if fractionnotacctinharm>.03
        disp('Greater than 3%,');
    end
    disp(['include original peakidx with the harmonic ids to get allpeakids']);      

    artpeakids=sort([harmonicsnotinorigxings; harmonicpeakids]);

    %percentage of bad harmonicpeakid intervals (or shifts in 60 hz noise)    
    fractionbad=length(find(diff(harmonicpeakids)>intervalnoise))/length(harmonicpeakids);
    if fractionbad/fractfoundline>.25
        disp(['fraction of bad intervals detected in getperiodicpeaks is close ' ...
            'to the initial xing crossings detected.' ...
            ' Most likely overfitting harmonic ids into non-harmonic noise.' ...
            ' Making artpeakids empty for this channel'])
        artpeakids=[];
    end
    if fractionbad>.001
        %more than 0.1% then may want to check manually
        warning(['fraction of out of sync 120 hz peak ids = ' num2str(fractionbad)])
    end

else
    disp('insufficient number of artifactual peaks found returning no artifacts from fitperiodicartifacts')
end


end
