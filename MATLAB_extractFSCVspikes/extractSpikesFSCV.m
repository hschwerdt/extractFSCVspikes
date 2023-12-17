function extractSpikesFSCV(subject,sessnum,varargin)
%10/09/2023 HNS revise formating and make sure code all exists for manusciprt
%configdir needs to be set manually for config files, otherwise checks current directory
%FSCV artifacts detected using fitFSCVart applying periodic template
% Interpolation of artifacts done in rmFSCVart
%made for entire csc file recording, also output interpolated data as
%temporary storage unless checked otherwise varargin
%use chronicXXchconfig.m to include
%ncsnoartifacts={'cl1', 's4','s3'};  %to denote what channels don't have
%FSCV %artifacts
%Whole file
%Also check for 60/120 Hz noise and interpolate these after fscv
%interpolation, not assuming what channels have this, explicit check

%%Variables that should be changed according to experimental config
samplerate=30e3;           %sample rate for neuralynx data, 
railvalue = 1e-3;   % Neuralynx clips at 1 mV
%default on putamen pc in lab, get data directory from config file
configdir='C:\Users\putamen\Documents\MATLAB\fscv\analysis\config\';
if ~isfolder(configdir)
    %separate config directory doesn't exist, use current workspace
    warning('default configdir does not exist, using current workspace')
    configdir=[pwd filesep]; 
end
saveNameLabel='fscv_multi_';    %Initial label of save name, followed by number defined by trial
FSCVfreq=10;    %FSCV sampling rate
forwardwin=7e-3*samplerate; %Forward window from FSCV artifact to interpolate = 7 ms (see 10/2022 notes on this)
backwardwin=5e-3*samplerate;    %Backward window from FSCV artifact to interpolate = 5 ms
forwardwinHF=5; %Forward window (in samples) ~ 0.166 ms for 60/120 hz artifacts
backwardwinHF=5;

%Setup test figure
ax=testAxes(8);%Create test figure with 8 subplots

%%Initialize main variables and load experimental settings
logdoc=[];  %logging errors and details to be exported
flagNLX=0;
subjectname='patra';        %Database identifier
if strcmp(subject,'cleo')
    subjectname='cleo';
end
sessid='';  %session number for experiment data
if ~isnumeric(sessnum)
    sessid=sessnum;
else
    sessid=num2str(sessnum);
end

ncschannels={};     %CSC channels (neuralynx channel #s)
csc_map={};     %Map variable to assign channel #s to specific site names from map file below
paths={};
event_codes={}; %Event codes from map (not used)
%Load experimental settings
switch subjectname
    case 'patra'
        try
            run([configdir ['chronic' sessid 'chconfig.m']]);        %get 'paths'
            run([configdir 'patra_map_bipolar.m']);      %get csc_map
        catch
            %in current directory?
            disp('config files not in configdir directory, trying current directory')
            run([pwd filesep ['chronic' sessid 'chconfig.m']]);        %get 'paths'
            run([pwd filesep 'patra_map_bipolar.m']);      %get csc_map
        end
    case 'cleo'
        run([configdir ['cleo_chronic' sessid '.m']]);        %get 'paths'
        run([configdir ['cleo_map_bipolar' sessid '.m']]);      %get csc_map & event_codes
end

intervalsplit=[];       %default full file for fixed interval default
mergeflag=1;            %default merge flag (i.e. not task separated)
alignname='fixedintervals'; %default full file for fixed interval default
    
%%Check for input arguments
argnum=1;
cleolickflag=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'nomerge'
            %don't merge entire recording into 1 file
            mergeflag=0;
        case 'cleolick37'
            cleolickflag=1;     %for chronic 25 where lick rec
        case 'interval'
            argnum=argnum+1;
            intervalsplit=varargin{argnum};
        case 'nlx'
            %convert to nlx
            flagNLX=1;
           % argnum=argnum+1;
           % cscadd=varargin{argnum};
        case 'sites'
            %explicitly indicate sites to convert, not from
            %chronicXXchconfig file
            %should be in format {'sitename1','sitename2',...}
            argnum=argnum+1;
            ncschannels=varargin{argnum};
    end
    argnum=argnum+1;
end
if isempty(intervalsplit)
    %if empty argument, then means entire file
    mergeflag=1;
end
preOffset=0;
durationTrial=0;

%Get csc channel numbers based on user defined ncschannels and csc_map
if cleolickflag==0
    [ephyschs,bipolarchs,otherchs,sumchs]=getcscids(ncschannels,csc_map);
else
    [ephyschs,bipolarchs,otherchs,sumchs]=getcscids(ncschannels,csc_map,'cleolick37');
end
selectcsc=[ephyschs];
noartcsc=[];
if exist('ncsnoartifacts')
    %channels listing no artifacts from chronicXXchconfig.m
    %These should not be included in signal used to get global artifact template
    noartcsc=getcscids(ncsnoartifacts,csc_map);
end

%Define savefolder within original fscv directory from paths{1} defined in
%chronicXXchconfig
if ~isdir(paths{1})
    %fscv folder may not exist since only lfps
    mkdir(paths{1});
end
dirfscv=dir([paths{1} '*_cvtotxt']);       %get cvtotxt files in path
filesfscv={dirfscv.name};
pathfscv=paths{1};
sep=findstr(filesep,pathfscv);      %indexes of separator '\' or '/' in pathname
savefolder=fullfile(pathfscv,'matlab',filesep);
if ~isdir(savefolder)
    status = mkdir(savefolder);
end

%Load ephys samples from directory paths{2} defined in chronicXXchconfig
pathnlx=paths{2};
sep2=findstr(filesep,pathnlx);
[samples, nlxFileNames, TS,header,nlxjunk]=ephys_getCSCwheader(selectcsc,pathnlx);
%TS here is vector nlx timestamps in seconds
%if TS is negative (session 134b, issue with frame conversion on the first
%~530 samples)
if any(TS<0)
    warning('Negative TSs present, attempting to fix based on non-negative TS')
    xid=find(TS>0);
    id1=xid(1);
    TSnormalstart=TS(id1);  %First normal timestamp
    mediansampleperiod=median(TS(id1+2:end)-TS(id1+1:end-1));
    medianrate=1/mediansampleperiod;
    duration=mediansampleperiod*(id1-1);
    fillTS(id1-1:-1:1)=TS(id1)-mediansampleperiod:-mediansampleperiod:TS(id1)-duration;
    TS(1:id1-1)=fillTS;
end

numChs=length(samples);

%Create save folder
subfolder='';
savepath=fullfile(pathfscv,'matlab','spikes',subfolder,filesep);
if ~isdir(savepath)
    status = mkdir(savepath);
end
temppath=fullfile(pathfscv,'matlab','spikes2','interp',subfolder,filesep);
%storage for interpolated raw data
if ~isdir(temppath)
    status = mkdir(temppath);
end
initialNum=100;
logdoc=[];

%Get header information for writing CSC
ADBitVoltstr={};
ADBitVolts=[];
for k = 1:length(header)
    if regexp(header{k}, '^\s*-ADBitVolts\s+')
        ADBitVoltstr = regexprep(header{k}, ...
            '^\s*-ADBitVolts\s+', '');
        ADBitVolts = str2num(ADBitVoltstr);                           
    end
end

durationrecording=(TS(end)-TS(1))/60;
writelog=['recorded ' num2str(durationrecording,5) ' minutes' newline ...
    'in session ' num2str(sessnum) newline...
    'writing interpolated spike files for csc ' ...
    num2str(selectcsc)];

%Create logdoc file
if isempty(logdoc)
    logdoc=writelog;
else
    logdoc=[logdoc sprintf('\n') writelog];
end

origTS=TS;  %Store original timestamps
origSamples={};

saveName=[saveNameLabel num2str(initialNum)];     
%get start IDs of sample
currentSampleIndex=1;   %first index of samples
lastSampleIdx=length(origTS);  %to end of recording     
nextSampleIdx=lastSampleIdx;    %get to end of recording   
if nextSampleIdx>length(samples{1})
    nextSampleIdx=length(samples{1});
end
for EphysCh=1:numChs
    origSamples{EphysCh}=samples{EphysCh}(currentSampleIndex:nextSampleIdx);
end
TS=TS(currentSampleIndex:nextSampleIdx);    %Timestamps defined in window (if needed)
samplesMod=origSamples; %Samples to be modified
    
%%Spike Extraction by interpolating FSCV artifacts 
%First get threshold by averaging all samples all chs with artifacts

%First find & remove channels that are open circuit (disconnected) based on 
% how often the signal hits the rail 1mv (1e-3)
opench=[];
for iachs=1:length(selectcsc)
    railnums=find(abs(origSamples{iachs})>=railvalue-0.01e-3);%find number of samples at rail close to 1e-3 V
    if length(railnums)/length(origSamples{iachs})>.4   
        %open channel can also be closer to 25%, this is an underestimate. 
        %some open channels will carry fscv artifacts so stil useful for
        %artifact removal
        opench=selectcsc(iachs);
        writelog=['open channel: 40% or more of csc ' num2str(selectcsc(iachs)) ' is hitting the rail, do not use, not interpolated or used for furhter processing'];
        disp(writelog);
        %record in logdoc file
        if isempty(logdoc)
            logdoc=writelog;
        else
            logdoc=[logdoc newline writelog];
        end
    end
end

%Collect ch ids of only artifact containing neural signals & not disconnected chs
artchs=find(~ismember(selectcsc,noartcsc) & ~ismember(selectcsc,opench));
sumsamples=samplesMod{artchs(1)};
%Get sum of all art chs, with polarity
for iachs=2:length(artchs)
    sumsamples=sumsamples+samplesMod{artchs(iachs)};
end
samplesmean=sumsamples./length(artchs);
        %TEST figure, Step 1 - Find channels with artifacts, and get mean
        curax=ax{1};
        title(curax,'Channels with FSCV Artifacts')
        hold(curax,'on')
        ylabel(curax,'Voltage (V)')
        xlabel(curax,'Time (s)')
        sampleIDs=samplerate*125*60:samplerate*125*60+samplerate*1;    %Window 30 min + 1 s
        sampleIDs=samplerate*237.25*60:samplerate*237.25*60+samplerate*1; %sess 83
                sampleIDs=samplerate*100*60:samplerate*100*60+samplerate*1;
                sampleIDs=samplerate*100.5*60:samplerate*100.5*60+samplerate*1; 
        sampleIDs=samplerate*558.5:samplerate*558.5+samplerate*1; %sess 127 example large glitch
        sampleIDs=find(TS>=9.2584349015e4 & TS<=9.259977444199999e+04); %sess 127 example trial 190 big reward
                sampleIDs=samplerate*4903.97:samplerate*4903.97+samplerate*1;      %Sess 127 csc 12 spike example          
        plot(curax,TS(sampleIDs)-TS(1),samplesMod{artchs(1)}(sampleIDs));    
        chstr={['csc' num2str(selectcsc(artchs(1)))]};
        for ix=2:length(artchs)
            plot(curax,TS(sampleIDs)-TS(1),samplesMod{artchs(ix)}(sampleIDs));
            chstr{ix}=['csc' num2str(selectcsc(artchs(ix)))];
        end
        legend(curax,chstr);
        xlim(curax,[TS(sampleIDs(1))-TS(1) TS(sampleIDs(end))-TS(1)])
        %TEST figure 2
        curax=ax{2};
        title(curax,'Average Signal')
        hold(curax,'on')
        ylabel(curax,'Voltage (V)')
        xlabel(curax,'Time (s)')
        plot(curax,TS(sampleIDs)-TS(1),samplesmean(sampleIDs));
        xlim(curax,[TS(sampleIDs(1))-TS(1) TS(sampleIDs(end))-TS(1)])

%Some samples will have NaN values if frames have excessively variable duration
%as found in ephys_getCSCwheader
%e.g. SESSION 126, nan signals found in less than 1 s period, but need
%to remove otherwise cannot filter without creating all nan signal
if any(isnan(samplesmean))
    numnan=length(find(isnan(samplesmean)));
        writelog=[num2str(numnan) '# samples in the signals are nan, most likely frames with excessively variable duration as found in ephys_getCSCheader'];
        warning(writelog);
        %record in logdoc file
        if isempty(logdoc)
            logdoc=writelog;
        else
            logdoc=[logdoc newline writelog];
        end
    samplesmean(isnan(samplesmean))=min(samplesmean); %Set signal to rail if is nan, any ways it would be disaqualified as an artifact whether LFP or spike analysis
end

%Find FSCV artifact peaks based on known FSCV parameters applied
allpeakids = fitFSCVart(TS, samplesmean); %Thres is made in function to detect artifact peaks
disp(['found ' num2str(length(allpeakids)) ' artifacts']);   
sessiondur=length(samplesmean)/samplerate;
fractfoundfscv=length(allpeakids)/(sessiondur*FSCVfreq);  
%Expected fraction of artifacts based on known 10 Hz frequency of FSCV artifact
disp(['fraction of expected 10 hz artifacts found: ' num2str(fractfoundfscv,4)])
if fractfoundfscv<.01
    %Less than 1% of recording, don't interpolate
    writelog=['No fscv interpolation occured as less than 1% of expected'];
    disp(writelog);
    %record in logdoc file
    if isempty(logdoc)
        logdoc=writelog;
    else
        logdoc=[logdoc newline writelog];
    end
else
    %> 1% of recording has FSCV artifacts
                %TEST figure, before interpolation, Peaks, example channel
                curax=ax{3};
                title(curax,'FSCV Peaks')
                hold(curax,'on')
                ylabel(curax,'Voltage (V)')
                xlabel(curax,'Time (s)')
                plot(curax,TS(sampleIDs)-TS(1),samplesMod{(selectcsc==12)}(sampleIDs));
                samppeakIDs=sampleIDs((ismember(sampleIDs,allpeakids)));
                plot(curax,TS(samppeakIDs)-TS(1),samplesMod{(selectcsc==12)}(samppeakIDs),'ko')
                xlim(curax,[TS(sampleIDs(1))-TS(1) TS(sampleIDs(end))-TS(1)])
    for iachs=1:length(artchs)
        %interpolate artifacts away with fixed 12 ms window for all chs, 5
        %ms backward and 7 ms forward
        samplesinterp=rmFSCVart(samplesMod{artchs(iachs)},allpeakids,backwardwin,forwardwin); %See notes on how this window was determined from 10/2022
        samplesMod{artchs(iachs)}=samplesinterp;
    end
                %TEST figure, after interpolation, example channel
                example_csc=12;
                curax=ax{4};
                title(curax,'After FSCV Interpolation')
                hold(curax,'on')
                ylabel(curax,'Voltage (V)')
                xlabel(curax,'Time (s)')
                plot(curax,TS(sampleIDs)-TS(1),origSamples{(selectcsc==12)}(sampleIDs));
                plot(curax,TS(sampleIDs)-TS(1),samplesMod{(selectcsc==12)}(sampleIDs),'k');
                legend(curax,{['Original Signal, csc' num2str(example_csc)],'After FSCV Interpolation'})
                xlim(curax,[TS(sampleIDs(1))-TS(1) TS(sampleIDs(end))-TS(1)])


    samplesinterp=rmFSCVart(samplesmean,allpeakids,backwardwin,forwardwin); 
    %Get mean signal without fscv artifacts, to isolate 60Hz+harmonics
    %could be for next section of code 

    samplesmean=samplesinterp;  %Replace samplesmean with interpolated mean
end

%%
%Remove 60/120 hz harmonic noise based on high pass filtered signal
%First, identify channels that contain 60 hz noise since the
%ones with fscv artifacts don't necessiarly contain 60 hz noise
art60chs=[];    %Channels with 60 hz noise
count=1;
filtsamplesabs={};
filtsamplessum=[];  %To get average of all channels to find consistent 60 hz noise
%Include noartchs (don't have FSCV artifacts, but could have 60 hz) 
artchs2=find(~ismember(selectcsc,opench));%all potential chs with 60 hz noise
disp('searching for channels with 60 hz artifacts');
for iachs=1:length(artchs2)
    tstart=tic;%timer function
    %look at all channels to find those that have 60 hz artifacts
    %after the fscv artifacts were interpolated above
    samplestemp=samplesMod{artchs2(iachs)};
    if any(isnan(samplesMod{artchs2(iachs)}))
        numnan=length(find(isnan(samplesMod{artchs2(iachs)})));
        warning([num2str(numnan) '# samples in the signals csc' ...
            num2str(selectcsc(artchs2(iachs)))   ' are nan, remove nans by setting as' ...
            ' min of signal before filtering for artifact detection'])
        samplestemp(isnan(samplestemp))=min(samplestemp); %Set signal to rail if is nan, any ways it would be disaqualified as an artifact whether LFP or spike analysis
    end
    %hp filter to enhance spike artifact detection 
    % (60/120 hz noise look very similar to neural spikes) 
    filtsamples=filterLFP(samplestemp,samplerate,[300 inf]); 
    %Get absolute value (for single-polarity threshold detection)
    filtsamplesabstemp=abs(filtsamples);
    %Define a low threshold for this first pass of detection 
    % in individual channesls (i.e. csc2, p1 has smaller artifacts)
    thres2=mean(filtsamplesabstemp, 'omitnan')*5;
    %Determine if harmonic noise signal by trying to fit harmonic peaks
    %in between peaks detected
    %Max blanks = 8 for initial check 
    [artpeakids,fractionbad] = fitperiodiclineartifacts(TS, filtsamplesabstemp, thres2,8);
    fractfoundline=length(artpeakids)/(sessiondur*120);           
    if fractfoundline<.05 || fractionbad>0.20
        %Less than 5% of recording, don't use
        %OR, fraction out of sync > 20% meaning not really a 60 or 120 hz
        %periodic signal
        writelog=['60 hz artifacts is ' num2str(fractfoundline*100,2) ...
            '% of expected in CSC ' num2str(selectcsc(artchs2(iachs))) ...
            ' Percent out of sync is ' num2str(fractionbad*100,2) ...
            ' - either less than 5% artifacts expected or out of sync ' ...
            '> 20% so not including'];
        disp(writelog);
        %record in logdoc file
        if isempty(logdoc)
            logdoc=writelog;
        else
            logdoc=[logdoc newline writelog];
        end
    else
        %Channel contains legit 60 hz noise
        % include in art60chs
        % Add to cumulative sum of channel signals with 60 hz noise
        writelog=['60 hz artifacts is ' num2str(fractfoundline*100,2) ...
            '% of expected in CSC ' num2str(selectcsc(artchs2(iachs))) ...
            ' - adding to process, add filtered to mean artifacts'];
        disp(writelog);
        %record in logdoc file
        if isempty(logdoc)
            logdoc=writelog;
        else
            logdoc=[logdoc newline writelog];
        end

        filtsamplesabs{count}=filtsamplesabstemp;

        % Add filtered signal to cumulative sum of channel signals with 60 hz noise
        if count==1
            filtsamplessum=filtsamplesabstemp;  %Using high pass filtered data
        else
            %Get mean from additional channels high pass filtered
            filtsamplessum=filtsamplessum+filtsamplesabstemp;
        end
        count=count+1;
        % include this channel in art60chs
        art60chs=[art60chs artchs2(iachs)];
    end
    tend=toc(tstart);   %timer end of loop
    disp(['Time for iding 60 hz in CSC ' num2str(selectcsc(artchs2(iachs))) ...
        ' : ' num2str(tend,3) ' seconds'])
end

%If only 4 channels or less with 60 hz noise, may need to check summed
%signal, are there sufficient recordings to get an accurate noise
%measurement
if length(art60chs)<5 && ~isempty(art60chs)
    warning('Less than 5 channels for identifying commons 60 hz related or other HF transients')
    warning('Check signals before proceeding');
    keyboard;       %Give control to user to debug
end

%60 hz noise exists, remove these artifacts
if ~isempty(art60chs)
    %If 60 hz artifacts exists on any of the channels
    disp('Identifying 60 hz artifacts on prospective channels');
    samplesmean=filtsamplessum./length(art60chs);
    %hp filter for spike artifact detection 
    filtsamplesmean=filterLFP(samplesmean,samplerate,[300 inf]);
    %Absolute value to constrain to positive crossings only
    filtsamplesmeanabs=abs(filtsamplesmean);
    %Higher threshold defined for actual detection and fitting of
    %harmonics to be more precise than initial detection pass
    thres2=mean(filtsamplesmeanabs,'omitnan')*8;
    %Get all harmonic peaks based on periodicity of expected 60/120 hz noise 
    %and fit these in between original peaks detected, keep original
    %peaks also as these might be transient glitches to remove
    %Maxblanks=16 for actual interp, can be bigger now since we already identified those channels with artifacts
    artpeakids = fitperiodiclineartifacts(TS, filtsamplesmeanabs, thres2,16); 

                %Test Figure 60 hz samples hpf mean, thres peaks
                curax=ax{5};
                title(curax,'Mean 2nd-order HPF magnitude After FSCV Interpolation')
                hold(curax,'on')
                ylabel(curax,'Voltage (V)')
                xlabel(curax,'Time (s)')
                %sampleIDs=samplerate*125*60:samplerate*125*60+samplerate*1;    %Window 30 min + 1 s
                sampleIDs=samplerate*237.25*60:samplerate*237.25*60+samplerate*.5;      %sess 83
                sampleIDs=samplerate*558.75:samplerate*558.75+samplerate*.5;      %Sess 127 example glitch
                    sampleIDs=samplerate*4904.2:samplerate*4904.2+samplerate*.5;      %Sess 127 csc 12 spike example          
                %filtsamples=filterLFP(samplesMod{3},samplerate,[300 inf]); 
                %Get absolute value (for single-polarity threshold detection)
                %filtsamplesabstemp2=abs(filtsamples);
                %plot(curax,TS(sampleIDs)-TS(1),filtsamplesabstemp2(sampleIDs));
                plot(curax,TS(sampleIDs)-TS(1),filtsamplesmeanabs(sampleIDs));
                plot(curax,[TS(sampleIDs(1))-TS(1) TS(sampleIDs(end))-TS(1)]...
                    ,[thres2 thres2],'g--');
                samppeakIDs=sampleIDs((ismember(sampleIDs,artpeakids)));
                plot(curax,TS(samppeakIDs)-TS(1),filtsamplesmeanabs(samppeakIDs),'ko')                    
                legend(curax,{'Mean HPF signal','Threshold'})
                xlim(curax,[TS(sampleIDs(1))-TS(1) TS(sampleIDs(end))-TS(1)])

    fractfoundline=length(artpeakids)/(sessiondur*120);
    writelog=['fraction of expected 60/120 hz artifacts found: ' num2str(fractfoundline,4)];
    disp(writelog);
    %record in logdoc file
    if isempty(logdoc)
        logdoc=writelog;
    else
        logdoc=[logdoc newline writelog];
    end
    if fractfoundline<.01
        %Less than 1% of recording, don't interpolate
        disp('No interpolation occured as less than 1% of expected')
    else
        %Interpolate 60/120 hz artifacts away
        samplesInterp=samplesMod;   %initiate variable
        writelog=['Interpolating ' num2str(selectcsc(art60chs))];
        disp(writelog);
        %record in logdoc file
        if isempty(logdoc)
            logdoc=writelog;
        else
            logdoc=[logdoc newline writelog];
        end
        for iachs=1:length(art60chs)
            %interpolate artifacts away with fixed 11 ms window for all chs, 
            % 5 samples backward and 5 samples forward from center
            % peak
            %interpolate 5 before and 5 samples after (333 microseconds)
            samplestempint=rmFSCVart(samplesMod{art60chs(iachs)},artpeakids,backwardwinHF,forwardwinHF); 
            samplesInterp{art60chs(iachs)}=samplestempint;
        end       
                %%Test Figure 6 and 7 after interpolating 60 hz      
                cla(ax{7})
                cla(ax{6})
                curax=ax{6};
                title(curax,'After 60 Hz Interpolation')
                hold(curax,'on')
                ylabel(curax,'Voltage (V)')
                xlabel(curax,'Time (s)')
                sampleIDs=samplerate*237.25*60:samplerate*237.25*60+samplerate*.5;   %sess 83
                sampleIDs=samplerate*558.75:samplerate*558.75+samplerate*.5;   %Sess 127  
                plot(curax,TS(sampleIDs)-TS(1),samplesMod{(selectcsc==example_csc)}(sampleIDs));
                plot(curax,TS(sampleIDs)-TS(1),samplesInterp{(selectcsc==example_csc)}(sampleIDs),'k');
                legend(curax,{['Original Signal, csc' num2str(example_csc)],'After 60 Hz Interpolation'},'location','northeast')
                xlim(curax,[TS(sampleIDs(1))-TS(1) TS(sampleIDs(end))-TS(1)])
                %Test Figure spikes example after interpolating 60 hz                                
                curax=ax{7};
                title(curax,'Filtered (Forward only Butter) at 250 hz')
                hold(curax,'on')
                ylabel(curax,'Voltage (V)')
                xlabel(curax,'Time (s)')
                filtsamples=filterForwardOnly(samplesInterp{(selectcsc==example_csc)},samplerate,[250 inf]);    %Forward only to emulate plexon
                filtsamples=filterForwardOnly(samplesMod{(selectcsc==example_csc)},samplerate,[250 inf]);    %Forward only to emulate plexon, sess 127
               % xspikeids=find(filtsamples<-30e-6);%find high amp spikes
                plot(curax,TS(sampleIDs)-TS(1),filtsamples(sampleIDs));
                xlim(curax,[TS(sampleIDs(1))-TS(1) TS(sampleIDs(end))-TS(1)])
               

        samplesMod=samplesInterp;%Store 60 hz removed signals
    end
else
    disp('no significant 60 hz artifacts detected, no interpolation')
end

                %%Test Figure 6 and 7 showing close up of spike waveforms after interpolating and hpf   
                cla(ax{7})
                cla(ax{6})
                curax=ax{6};
                title(curax,'After 60 Hz Interpolation')
                hold(curax,'on')
                ylabel(curax,'Voltage (V)')
                xlabel(curax,'Time (s)')
                example_csc=12;
                sampleIDs=samplerate*4903.97:samplerate*4903.97+samplerate*.2;      %Sess 127 csc 12 spike example          
                plot(curax,TS(sampleIDs)-TS(1),origSamples{(selectcsc==example_csc)}(sampleIDs));
                plot(curax,TS(sampleIDs)-TS(1),samplesMod{(selectcsc==example_csc)}(sampleIDs),'k');
                legend(curax,{['Original Signal, csc' num2str(example_csc)],'After 60 Hz Interpolation'},'location','northeast')
                xlim(curax,[TS(sampleIDs(1))-TS(1) TS(sampleIDs(end))-TS(1)])
                ylim(curax,[-3e-4 1e-4])
                %Test Figure spikes example after interpolating 60 hz                                
                curax=ax{7};
                title(curax,'Filtered (Forward only Butter) at 250 hz')
                hold(curax,'on')
                ylabel(curax,'Voltage (V)')
                xlabel(curax,'Time (s)')
                filtsamples=filterForwardOnly(samplesMod{(selectcsc==example_csc)},samplerate,[250 inf]);    %Forward only to emulate plexon, sess 127
                plot(curax,TS(sampleIDs)-TS(1),filtsamples(sampleIDs));
                xlim(curax,[TS(sampleIDs(1))-TS(1) TS(sampleIDs(end))-TS(1)])

%%Save interpolated data as Neuralynx NCS CSC files
%save in nlx format compatible with dg / lfp lib
dg_Nlx2Mat_Timestamps=round(TS.*1e6);     %convert to Microseconds * MAKE SURE TO ROUND OTHERWISE GET ROUNDOFF ERROR
dg_Nlx2Mat_Samples={};
dg_Nlx2Mat_SamplesUnits='microVolts'; %convert to microvolts for Plexon offline sorter

for nlxch=1:numChs
    %save each nlx channel as separate CSC neuralynx file
    %convert to microvolts for Plexon offline sorter
    dg_Nlx2Mat_Samples=samplesMod{nlxch}.*1e6;   
    savecsc=['csc' num2str(selectcsc(nlxch)) '_' num2str(initialNum)];
    if ~flagNLX
        save([temppath savecsc],'dg_Nlx2Mat_Samples','dg_Nlx2Mat_SamplesUnits','dg_Nlx2Mat_Timestamps','-v7.3');
    else
        %convert to NLX using dg_writeCSC
        %first convert to frames of 512 samples per frame,
        %timestamp for last sample
        %save in nlx format compatible with dg / lfp lib
        dg_Nlx2Mat_Timestamps=round(TS.*1e6);     %in microseconds for nlx
        %convert to ad units for nlx, see header
        samplesSave={};                    
        samplesSave{nlxch} = samplesMod{nlxch}./ADBitVolts; 
        dg_Nlx2Mat_Samples=round(samplesSave{nlxch});   
        dg_Nlx2Mat_SamplesUnits='AD'; 
        nlx_TS=dg_Nlx2Mat_Timestamps(1:512:end);    %get TS every 512 timestamps to put into frames
        samplesnlx=reshape(dg_Nlx2Mat_Samples,512,length(dg_Nlx2Mat_Samples)/512);   %reshape samples into 512 length frames
        disp(['writing ' savecsc '.ncs to ' temppath]);
        dg_writeCSC([temppath savecsc '.ncs'], nlx_TS, samplesnlx, header,'nlxjunk',nlxjunk);   %need to write junk that came with origianl nlx file in order to read into plexon
    end                
end
         
savelogname=['log_extractSpikesFSCV'];
save([temppath, savelogname],'logdoc'); %save logdoc of latencies

%%Copy events files to savepath if not already in new folder
%useful for offline sorter to have this file in the new path
%get events file from path{2} from chronicXXchconfig.m
direvents=dir([paths{2} '*.nev']);
fileevents=[paths{2} direvents.name];       %get events file name, this gets re-exported to new folder also
if ~exist([temppath 'events.nev'],'file')
    status=copyfile(fileevents, temppath); %Copy events .nev file to new directory  
    if status
        disp(['copied ' fileevents ' to ' temppath]);
    elseif ~status
        disp('failed to copy events .nev file')
    end
    status2=copyfile([[fileevents(1:end-4)] '.mat'], temppath); %Copy converted events .mat file to new directory  
    if status2
        disp(['copied ' [[fileevents(1:end-4)] '.mat'] ' to ' temppath]);
    elseif ~status
        disp('failed to copy events .mat file')
    end
end

%Save Test Figure
%pdfSave([temppath 'extractSpikesfig3session127_trial190-3'], [8 11.5], gcf);
pdfSave([temppath 'extractSpikesfig3session127_4'], [15 11.5], gcf);

end