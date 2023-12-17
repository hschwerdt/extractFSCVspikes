function allpeakids = fitFSCVart(ts, samples)
%10/2023 cleaned up from original code dg_fitFSCVart3, by Helen and Dan
%Assuming 10 Hz FSCV scan frequency, find artifact center indices based on
%signal periodicity
%Accounts for when baseline DC shift or prolonged low frequency changes
% that would place the signal above the threshold for a long time and prevent
% detection based on harmonics. Instead look for peaks during such baseline shifts using
%highpass filtered signal for 60 hz harmonic noise
%INPUTS
% ts: timestamps in seconds of all samples, same size as <samples>.
% samples: array of sampled values, one per timestamp.
%OUTPUTS
% allpeakids: the ids of <samples> of the "center" of each identified FSCV
%   artifact. Center based on peak of band-passed waveform
%NOTES
%   Positive-going threshold crossings are used as a rough first 
% approximation to the FSCV scan times.  If the slope of the two samples at and
% immediately before the threshold crossing is less than 2% of <thresh> per
% sample, then the threshold crossing is assumed not to actually be an FSCV
% artifact and is removed from the list.  (Note that this criterion also
% removes crossings where the baseline is so close to threshold that there
% is less than 300 microseconds between the start of the scan and the
% threshold crossing.) For the remaining threshold crossings, the peak is 
% taken for the peak artifact idx. These peak ids are then fed into
% getperiodicpeaks.m to get the actual periodic artifacts and remove
% non-periodic artifact peaks identified from initial rough detection
% within this code
%   A simple linear interpolator is provided in dg_rmFSCVart.
%   The sample scanperiod is estimated based only on the first 10,000 samples.

maxblanks=20;   %Max blank periods/artifacts to interpolate sometimes if low freq bias 
% then the apparent freq shifts and will be filtered out with 10-100 bandpass below, 
% but easier to work with this than to make wider bandwidth which would allow glitches in
sampleperiod=mean(diff(ts(1:1e4)));%This is more accurate scanperiod estimated based only on the first 10,000 samples.
samplerate=round(1/sampleperiod);
% These constants might have to become arguments:
fitwindow = floor(3e-3 / sampleperiod); % for piecewise linear fit
minscansamp = round(.001/sampleperiod); % min number suprathreshold samples
periodicitytoltight=0.05;   %max fractional deviation from the excpected period
waveformperiod=round(8.5e-3/sampleperiod);  %theoretical fscv waveform period in samples
waveformtol=0.5;   %Low band pass filter applied below + filtering affect of tissue widens ~ 3.3 ms from 8.5 ms theoretical
maxwidth=round(waveformperiod*waveformtol+waveformperiod);

%Remove DC offsets to help detect artifacts based on threshold crossing
origsamples=samples;    %store original samples
%Band pass filter high-pass to remove DC bias and low pass to avoid 
% high-frequency phase distortions (pg 65 10/15/2022 notes)
lbpsamples=filterLFP(abs(origsamples),samplerate,[10 100]);    
%120 Hz is theroetical primary freq cut off for fscv period of 8.5 ms, but 
% RC from tissue makes it slower in appareance, also want to make sure to remove 
% higher frequency omponent glitches and artifacts from 60 hz
thresh= std(lbpsamples,'omitnan')*1.75; %Create threshold for artifact detection
samples=lbpsamples;



%Find positive going threshold crossings from bandpass filtered data
isxing = [ false
    reshape( samples(2:end) >= thresh & samples(1:end-1) < thresh, ...
    [], 1 ) ];
xingidx = reshape(find(isxing), [], 1 );    %Threshold crossing ids
isxingdown = [ false
reshape( samples(1:end-1) >= thresh & samples(2:end) < thresh, ...
[], 1 ) ];  %Not used
xingidxdown = reshape(find(isxingdown), [], 1 );    %Not used

%Realign to real "peak" points of the xingidx artifacts detected above
peakidx=zeros(size(xingidx));   %initiate variable for "peak" idx of artifacts
bogusxings=[];  %Record "bogus" xings based on low slope
toobrief=[];    %Record artifacts whose waveforms are too brief to be fscv artifact
for idx2 = 1:length(xingidx)
    slope = samples(xingidx(idx2)) - samples(xingidx(idx2) - 1);
    if slope < 0.005 * thresh
        % Bogus xing, skip:
        bogusxings=[bogusxings xingidx(idx2)];
        continue
    end
    % Check to make sure there really is something that could possibly be
    % an FSCV artifact:
    if sum(samples(xingidx(idx2) + (0:(2*fitwindow))) > thresh) ...
            < minscansamp
        % Ignore this one, it's too brief.
        toobrief=[toobrief xingidx(idx2)];
        continue
    end

    %Instead of looking for true art start, look for local peak
    arttrace=samples(xingidx(idx2):xingidx(idx2)+maxwidth); 
    %Multiple ways to record artifacts, center or peak.
    % Center means looking for midpoint between positive going slope 
    % and negative going slope. 
    % Here we use peak, the maximum within max window past the positive
    % crossing xingidx
    [peakmax,peakid]=max(arttrace);
    artpeak=peakid+xingidx(idx2)-1; %true peak idx of artifact
    peakidx(idx2)=artpeak;

end
% Delete the entries for the bogus xings:
peakidx(peakidx==0)=[];



%Get all 10 Hz peaks fit between those found in case miss any
%Why we set max blank periods to 10? May be FSCV turned off or some other
%larger srouces of interference that we don't want to artificially
%generating interpolated signals over, these larger srouces of interference
%are checked for in subseq algorithms that look for movement artifacts or
%pump related artifacts
[allpeakids,nonperiodicpeaks]=getperiodicpeaks(sampleperiod,peakidx,periodicitytoltight,10,10,maxblanks);


%%Test figures

%Setup test figure
ax=testAxes(8);%Create test figure with 8 subplots

        %Test Figure 1 filtered siganl adn threshold
        curax=ax{1};
        title(curax,'Band pass filtered (10 - 100) of signal magnitude')
        hold(curax,'on')
        ylabel(curax,'Voltage (V)')
        xlabel(curax,'Time (s)')
        sampleIDs=samplerate*125*60:samplerate*125*60+samplerate*1;
        %%Window 30 min + 1 s clean sess 83
        %sampleIDs=samplerate*237.25*60:samplerate*237.25*60+samplerate*1; 
        %sampleIDs=samplerate*12.6*60:samplerate*12.6*60+samplerate*1; %example LFP, sess83
        %sampleIDs=samplerate*104.84*60:samplerate*104.84*60+samplerate*1; %example FSCV glitch, sess83 nonperiodic peaks
        %sampleIDs=samplerate*160.56*60:samplerate*160.56*60+samplerate*1; %example FSCV glitch, sess83 nonperiodic peaks
        % sampleIDs=samplerate*558.5:samplerate*558.5+samplerate*1; %sess 127 example large glitch
                          %  sampleIDs=samplerate*4904.2:samplerate*4904.2+samplerate*1;      %Sess 127 csc 12 spike example    
                                            sampleIDs=samplerate*4903.97:samplerate*4903.97+samplerate*1;      %Sess 127 csc 12 spike example          

        plot(curax,ts(sampleIDs)-ts(1),origsamples(sampleIDs));
        yyaxis(curax,'right');
        plot(curax,ts(sampleIDs)-ts(1),lbpsamples(sampleIDs));
        plot(curax,[ts(sampleIDs(1))-ts(1) ts(sampleIDs(end))-ts(1)]...
            ,[thresh thresh],'g--');
        xlim(curax,[ts(sampleIDs(1))-ts(1) ts(sampleIDs(end))-ts(1)])

        %Test Figure 2 initial peak ids
        curax=ax{2};
        title(curax,'Initial pass peaks based on threshold crossing')
        hold(curax,'on')
        ylabel(curax,'Voltage (V)')
        xlabel(curax,'Time (s)')
        %sampleIDs=samplerate*237.25*60:samplerate*237.25*60+samplerate*1; 
        %sampleIDs=samplerate*12.6*60:samplerate*12.6*60+samplerate*1; 
        % sampleIDs=samplerate*558.5:samplerate*558.5+samplerate*1; %sess 127 example large glitch
        plot(curax,ts(sampleIDs)-ts(1),origsamples(sampleIDs));
        samppeakIDs=sampleIDs((ismember(sampleIDs,peakidx)));
        plot(curax,ts(samppeakIDs)-ts(1),origsamples(samppeakIDs),'ko')
        yyaxis(curax,'right');
        plot(curax,ts(sampleIDs)-ts(1),lbpsamples(sampleIDs));
                plot(curax,[ts(sampleIDs(1))-ts(1) ts(sampleIDs(end))-ts(1)]...
            ,[thresh thresh],'g--');
        plot(curax,ts(samppeakIDs)-ts(1),lbpsamples(samppeakIDs),'ko')
        xlim(curax,[ts(sampleIDs(1))-ts(1) ts(sampleIDs(end))-ts(1)])

        %Test Figure 3 showing FSCV glitch, nonperiodicpeak (sess 83,
        %104.84 min)
        curax=ax{3};
        title(curax,'Initial pass peaks based on threshold crossing')
        hold(curax,'on')
        ylabel(curax,'Voltage (V)')
        xlabel(curax,'Time (s)')
      %  sampleIDs=samplerate*104.84*60:samplerate*104.84*60+samplerate*1; %example FSCV glitch, sess83 nonperiodic peaks
       %  sampleIDs=samplerate*558.5:samplerate*558.5+samplerate*1; %sess 127 example large glitch
        plot(curax,ts(sampleIDs)-ts(1),origsamples(sampleIDs));
        samppeakIDs=sampleIDs((ismember(sampleIDs,peakidx)));
        plot(curax,ts(samppeakIDs)-ts(1),origsamples(samppeakIDs),'ko')
        yyaxis(curax,'right');
        plot(curax,ts(sampleIDs)-ts(1),lbpsamples(sampleIDs));
        plot(curax,ts(samppeakIDs)-ts(1),lbpsamples(samppeakIDs),'ko')
        xlim(curax,[ts(sampleIDs(1))-ts(1) ts(sampleIDs(end))-ts(1)])
        samppeakIDsCorrected=sampleIDs((ismember(sampleIDs,allpeakids)));
        plot(curax,ts(samppeakIDs)-ts(1),lbpsamples(samppeakIDsCorrected),'ro')
        yyaxis(gca,'left')

        %Test Figure 4 FFT
        %get 1 min period of signal
     %   sampleIDs=samplerate*100*60:samplerate*100*60+samplerate*1*60;
       %  sampleIDs=samplerate*558.5:samplerate*558.5+samplerate*1; %sess 127 example large glitch
       sampleIDs=samplerate*125*60:samplerate*125*60+samplerate*60;
                           sampleIDs=samplerate*4904.2:samplerate*4904.2+samplerate*60;      %Sess 127 csc 12 spike example          
                sampleIDs=samplerate*4903.97:samplerate*4903.97+samplerate*60;      %Sess 127 csc 12 spike example          

        %%Window 30 min + 1 s clean sess 83 (1 minute for FFT)
        samplesnippet=abs(origsamples(sampleIDs));
        samplesFFT = reshape(fft(samplesnippet), [], 1);  %get fft
        samplesFFTshift=fftshift(samplesFFT);       %shift spectra to normal frequency scale (ie. mirroring 0)
        df=samplerate./length(samplesnippet);
        freq = (-samplerate/2:df:samplerate/2-df)' + mod(length(samplesnippet),2)*df/2;
        curax=ax{4};
        title(curax,'FFT (over one min)')
        hold(curax,'on')
        plot(curax,freq,20*log10(abs(samplesFFTshift)))
        xlabel(curax,'Frequency (Hz)')
        ylabel(curax,'Magnitude (dB)')
        xlim(curax,[0 300])
        %legend({'original signal' ,'FSCV harmonics subtracted'});
        ylim(curax,[-50 50]); 

        %Test Figure 5 bandpass filter freq response
        curax=ax{5};
        title(curax,'Filter Response')
         hold(curax,'on')
         plot(curax,freq,20*log10(abs(samplesFFTshift)))
         ylim(curax,[-20 50]); 
         ylabel(curax,'Magnitude (dB)')
         yyaxis(curax,'right')
         freqlim = [10 100]/samplerate*2;
        [zb, pb, kb] = butter(8, freqlim);
        [bb,ab] = zp2tf(zb,pb,kb);
        [hb,wb] = freqs(bb,ab,4096);
        [sos,g]=zp2sos(zb,pb,kb);
        [h,w] = freqz(sos,5000) ;
        plot(curax,w*samplerate/(2*pi), ...
            20*log10(abs(h*g)),'k--');    %Need to multiply gain factor g
         hold(curax,'on')
        axis(curax,[0 300 -50 10])
        xlabel(curax,'Frequency (Hz)')
        ylabel(curax,'Attenuation (dB)')
        yyaxis(gca,'left')

        %Test Figure 6 filterede
        samplesnippetLBP=lbpsamples(sampleIDs);
        samplesFFTLBP = reshape(fft(samplesnippetLBP), [], 1);  %get fft
        samplesFFTshiftLBP=fftshift(samplesFFTLBP);       %shift spectra to normal frequency scale (ie. mirroring 0)
        df=samplerate./length(samplesnippetLBP);
        curax=ax{6};
        title(curax,'Filtered Signal')
        hold(curax,'on')
        plot(curax,freq,20*log10(abs(samplesFFTshift)))   %Original
        plot(curax,freq,20*log10(abs(samplesFFTshiftLBP)))         
        ylabel(curax,'Magnitude (dB)')                
        xlim(curax,[0 300]);
        ylim(curax,[-50 50]);         
        yyaxis(curax,'right')
        plot(curax,w*samplerate/(2*pi), ...
            20*log10(abs(h*g)),'k--');    %Need to multiply gain factor g
        hold(curax,'on')
        ylim(curax,[-50 10])
        xlabel(curax,'Frequency (Hz)')
        ylabel(curax,'Attenuation (dB)')
        legend(curax,'Original','Filtered','Filter Response') 

        
%Save Test Figure
pdfSave(['fitFSCVartsession83example_samplesmean_normal'], [8 11.5], gcf);
pdfSave(['fitFSCVartsession127example_samplesmean_spikeexample'], [8 11.5], gcf);
