function [periodicpeaks,peaksnotcounted]=getperiodicpeaks(sampleperiod,peakidx,periodicitytol,fa,ft,maxblankperiods)
%10/2023 cleaned up
%Function to create a template of timestamp ids of when artifacts should occur 
% according to the fundamental and harmonic frequencies of measured
% artifact ids (input peakidx), ignoring non-harmonic artifacts and long
% delays between artifacts (based on maxblankperiods)
%peakidx is original artifact peaks detected based on threshold crossing in
%fitFSCVart
%Fa is frequency of artifact to add artifacts as needed into
%Ft is the higher harmonic or same frequency of primary artifact, which
%helps to define timepoints of periodic artifacts if these are also
%threshold crossing. Fa = Ft when looking at 10 Hz FSCV artifacts only
%sampleperiod is sampling period of original signal based on sampling rate (e.g. 1/30k)
%periodicitytol = tolerance of periodicity (usually 5%) (not perfect
%because of low sampling rate of FSCV system in comparison to EPhys and
%potentially some inconsistencies in timestamps)
%Now go through peak (midpoints of artifacts, check periodicity, 
% add artifacts in between stable harmonics)
sampperiod=round(1/sampleperiod*1/fa);%500 samples for 60 Hz at 30e3 sample rate
target_period=round(1/sampleperiod*1/ft);%For 10 hz fscv artifact ft (target freq) = fa
if sampperiod/target_period<1 || round(sampperiod/target_period)~=sampperiod/target_period
    %If targetperiod not equal to or less than sampperiod than cannot
    %interpolate between sampperiod. If sampperiod cannot be divided by
    %targetperiod then also not a harmonic
    error('not a periodic target frequency');
end

allpeakids=zeros(length(peakidx)*2,1);  %Initiate allpeakids variable to contain harmonic peaks
countpeaks=1;   %Harmonic peaks counted
maxblankflag=0; %Flag if already had an interval above maxblankperiods (prolonged period of no evident artifact)
lastxingidx=peakidx(1); %Set previous threshold crossing idx as first peakidx to initiate
peaksnotcounted=zeros(length(peakidx),1);   %Record peakidx's not fitting definition of a consistent periodic
nc=1;
for idx2 = 2:length(peakidx)    
    %loop through each identified peakidx from fitFSCVart and
    %retrospectively add harmonic peaks as needed for artifacts that don't
    %cross previously applied threshold (in fitFSCVart) and to remove
    %peaks that don't fit the criterion for being a periodic signal
    if countpeaks==1
        interval = peakidx(idx2) - lastxingidx;   %Likely using bad lastxingidx to begin with, 
        % %will need to do this a couple times until a good one pops up, 
        % %then we can actually start looking assuming harmonic intervals or skip bad ones
        remainder = rem(interval, sampperiod);%Should be close to zero if period matches 1/fa
        if (remainder/(sampperiod)>periodicitytol && ...
             remainder/(sampperiod)<(1-periodicitytol))
            %Is not a harmonic within tolerance store 
            lastxingidx=peakidx(idx2); %Use current peak id as next lastxingidx to check next run
            peaksnotcounted(nc)=peakidx(idx2);
            nc=nc+1;
            continue
        end
        %Otherwise Is harmonic, loops out below to store in
        %allpeakids and increment countpeaks, but no intervals
        %added in betweeen because this is just about the
        %lastxingidx, not currentpeak. Next conditional adds.
        allpeakids(countpeaks)=lastxingidx;%Add lastxingidx since it is a harmonic
        countpeaks=countpeaks+1;
    end
    if countpeaks>1
        %This happens the first time when above conditional ran also.
        %Check if harmonic based on last stored id
         interval = peakidx(idx2) - allpeakids(countpeaks-1); 
         %Interval will never be a harmonic, in theory, if at some point 
         % there is delay shift in noise timing, hence prospective check added below
         remainder = rem(interval, target_period);%Accounts for lower/higher harmonics (ft)
         if (remainder/(sampperiod)<=periodicitytol || ...
                 remainder/(sampperiod)>(1-periodicitytol))
             %Tolerance based on fundamental fa frequency (sampperiod)
            %Check it is a harmonic/periodic 
            numperiods=floor(interval/(target_period));%How many 1/ft intervals to be added
            if remainder/sampperiod>periodicitytol
                %If the number of remainder is greater than it should be,
                %include an additional number to be added
                numperiods=numperiods+1;
            end
            if numperiods>maxblankperiods
                %Too long of a gap between artifacts based on maxblankperiods threshold
                % may be no artifacts anymore in the signal?
                %temporarily store ids to check in next loop
                % If last peakidx was also after a maxblankperiod then we
                %are just storing ids of long quiet intervals that we dont' want to interpolate over
                % In this case, remove the last peakidx stored in
                % allpeakids, replace with current peakidx
                if maxblankflag==1
                    %Already flagged for quiet period (no artifacts) in
                    %previous loop
                    allpeakids(countpeaks-1)=peakidx(idx2);%Add current peak idx to replace last one that was found after long interval (Ie remove last id)
                    lastxingidx=allpeakids(countpeaks-1);%Store current peak as the last good for next loop
                    peaksnotcounted(nc)=peakidx(idx2);
                    nc=nc+1;
                    %Dont increment countpeaks, just keep replacing old id
                    %until over multiple maxblank periods intervals
                    continue
                else
                    %First time, do not fill retrospectively with periodic intervals
                    % set maxblankflag=1 to check next run and if this
                    % peakidx still produces a long interval prospectively
                    % also then remove it from the allpeakids in the
                    % conditional above, otherwise proceed with filling
                    % below
                    allpeakids(countpeaks)=peakidx(idx2);%Add current peak idx 
                    countpeaks=countpeaks+1;    %Increment counts for allpeakids
                    lastxingidx=allpeakids(countpeaks-1);%Store current peak as the last good for next loop
                    peaksnotcounted(nc)=peakidx(idx2);%Store in variable as does not match definition of a consistent periodic artifact
                    nc=nc+1;
                    maxblankflag=1;                %FLAG too long of a blank periodic interval            
                    continue
                end
            end
            if numperiods>=2
                %if numperiods==1 then interval already
                %Add peaks up until current peak peakidx(idx2), but not
                %current peak since that's added by default at end of loop
                for ii=1:numperiods-1
                    allpeakids(countpeaks)=allpeakids(countpeaks-1)+target_period;  %Fill higher harmonic (e.g., 120 hz) peak timestamps predicted
                    countpeaks=countpeaks+1;          
                end
            end
         else             
            %Not a harmonic based on previous idx and current index 
            % --> prospective check ahead next peakidx interval
            %  what if signal was delayed/shifted off of non-periodic noise glitch-type signal
            if idx2<length(peakidx)
                %Check prospective peakidx relative to current peakidx
                interval_next = peakidx(idx2+1) - peakidx(idx2);
                remainder_next = rem(interval_next, target_period);%Accounts for lower/higher harmonics          
                if remainder_next/sampperiod>periodicitytol && ...
                    remainder_next/sampperiod<(1-periodicitytol)
                    %Curent idx minus previous idx not a harmonic and next idx
                    %- current indx also not a harmonic, skip this peakidx
                    peaksnotcounted(nc)=peakidx(idx2);
                    nc=nc+1;
                    continue
                end
                %If harmonic next prospective index then peak will be added
                %to allpeakids below, so will be checked out next loop for
                %filling peak intervals as needed, but no retrospective
                %fill this loop
            end
         end
        maxblankflag=0; %Reset flag for long blank interval
        allpeakids(countpeaks)=peakidx(idx2);%Add current peak idx finally
        countpeaks=countpeaks+1;    
        lastxingidx=allpeakids(countpeaks-1);%Store current peak as the last good for next loop
    end

end
allpeakids=unique(allpeakids);
allpeakids(allpeakids==0)=[];
periodicpeaks=allpeakids;%output return
peaksnotcounted=unique(peaksnotcounted);    %for record keeping, original peakidxs that are not harmonics identified