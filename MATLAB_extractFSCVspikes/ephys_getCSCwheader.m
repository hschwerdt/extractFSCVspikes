function [samples, fileNames, TS,header,nlxjunk]=ephys_getCSCwheader(chnums, sessiondir)
%08/2021 include header so that write back with eheader in
%extractSpikesFSCV
%Requires dg_readCSC and lfp_reconcileFrames
%also output nlxjunk from dg_readcsc necessary for plexon read if exporting
%back to nlx format
%convert signals to dg style mat
%INPUTS
% 
% sessiondir: absolute or relative path to directory containing CSC files.
% 
%OUTPUTS are only to .mat files.
%OPTIONS
% 'srcdir', srcdir - the raw CSC files are in <srcdir>. Default value is
%   <sessiondir>.

srcdir = sessiondir;

%search directory for default naming convention of cheetah files based on
%selected channel input in chnums
currDir=dir(srcdir);
postFormat='.ncs';  preFormat='csc';
nameDefault=[preFormat num2str(chnums(1)) postFormat];
%defaultNameFound=contains({currDir.name},nameDefault,'IgnoreCase',true);
%Contains does not work in matlab 2013
defaultNameFound=strfind({currDir.name},nameDefault);     %case-sensitive
defaultNameFound=regexpi({currDir.name},nameDefault);    %Case insensitive 03/2023
defaultNameFound=find(~cellfun(@isempty,defaultNameFound));

%if sum(defaultNameFound)<1
if isempty(defaultNameFound)
    %not found try other types of names eg. downsampled format in mat based
    %on 'down' separator typically found in these names to identiy "format"
    %before/after csc name
    nameTry=['csc' num2str(chnums(1))];
    %find file index of files in directory with nameTry name
    %defaultNameFoundTry=contains({currDir.name},nameTry,'IgnoreCase',true);
    defaultNameFoundTry=strfind({currDir.name},nameTry);     %case-sensitive
    otherNameID=find(~cellfun(@isempty,defaultNameFoundTry));

    %get file index identified iwth nameTry name
    %otherNameID=find(defaultNameFoundTry>0);
    %for first identified Name, get name index in name where digit occurs
    %(should be channel #)
    digitFound=regexp(currDir(otherNameID(1)).name,'\d');
    %get name before channel #
    preFormat=currDir(otherNameID(1)).name(1:digitFound(1)-1);
    %get name after # (that is not part of channel # and not part of other
    %#s in name eg downsamplingrate)
    postFormatID=digitFound(1);
    %if channel identifying # greater than 1 digit, need to get last digit
    if length(digitFound)>1
        if digitFound(2)==digitFound(1)+1
            postFormatID=digitFound(2);
        end
    end
    postFormat=currDir(otherNameID(1)).name(postFormatID+1:end);
   
end

%get list of files in directory to load
fileNames={};
for refidx=1:length(chnums)
    getName=[preFormat num2str(chnums(refidx)) postFormat];
    %check named file exists in directory
    %nameFound=contains({currDir.name},getName,'IgnoreCase',true);
    nameFound=strfind({currDir.name},getName);     %case-sensitive
        nameFound=regexpi({currDir.name},getName);     %case-insensitive 03/2023

    nameFound=find(~cellfun(@isempty,nameFound));

    %if sum(nameFound)<1
    if isempty(nameFound)
        error(['file ''' getName ''' does not exist']);
    end
    %get actual name if case differs
   % actualName=currDir(nameFound>0).name;
    actualName=currDir(nameFound).name;
    fileNames{refidx}=actualName;    
end

timestamps={};
lfp_Samples={};
lfp_TimeStamps={};
header={};
nlxjunk={};
samples={};     %non-frame samples
TS=[];          %non-frame TS
for refidx = 1:length(chnums)
    display(['reading ' fileNames{refidx} ' | ' num2str(refidx) '/' num2str(length(chnums))]);
    CSCfilespec = chnums(refidx);
    if isempty(CSCfilespec)
        continue
    end
   
    clear dg_Nlx2Mat_Samples dg_Nlx2Mat_Timestamps;
    pathnameCSC=[srcdir fileNames{refidx}];
    read_mode=[];
    [pathup, infilenameup, ext] = fileparts(fileNames{refidx});
    switch lower(ext)
        case {'.ncs' '.dat'}
            read_mode = 'nlx';
        case '.mat'
            read_mode = 'mat';
        otherwise
            error('lfp_read2:badPreset', ...
                'Preset CSC files must be .NCS or .MAT' );
    end
    
    filenum=refidx;
    switch read_mode
    case 'nlx'
        [timestamps{filenum}, lfp_Samples{filenum}, header,nlxjunk] ...
            = dg_readCSC(pathnameCSC);
        if length(timestamps{filenum}) ~= ...
                size(lfp_Samples{filenum},2)
            error('lfp_read2:badNlx', ...
                'The file %s contains different numbers of frames and timestamps', ...
                src );
        end
        lfp_SamplesUnits{filenum} = 'AD';
        if ~exist('UnitlessFileRegexp') || ...
                isempty(regexpi(pathnameCSC, UnitlessFileRegexp))
            for k = 1:length(header)
                if regexp(header{k}, '^\s*-ADBitVolts\s+')
                    ADBitVoltstr = regexprep(header{k}, ...
                        '^\s*-ADBitVolts\s+', '');
                    ADBitVolts = str2num(ADBitVoltstr);
                    if isempty(ADBitVolts)
                        warning('lfp_read2:badADBitVolts', ...
                            'Could not convert number from:\n%s', ...
                            header{k} );
                    else
                        lfp_Samples{filenum} = ADBitVolts ...
                            * lfp_Samples{filenum};
                        lfp_SamplesUnits{filenum} = 'V';
                    end
                end
            end
        end
    case 'mat'
        load('-MAT', pathnameCSC);
        timestamps{filenum}= double(dg_Nlx2Mat_Timestamps);
        clear dg_Nlx2Mat_Timestamps;
        lfp_Samples{filenum} = double(dg_Nlx2Mat_Samples);
        clear dg_Nlx2Mat_Samples;
        if exist('dg_Nlx2Mat_SamplesUnits') == 1
            lfp_SamplesUnits{filenum} = dg_Nlx2Mat_SamplesUnits;
        end
    end
samples{refidx}= reshape(lfp_Samples{refidx}, 1, []);     
    
    
end
    
%get info and consolidate all channel data
lfp_TimeStamps = lfp_reconcileFrames(timestamps,'path',srcdir); 
%lfp_TimeStamps and timestamps{1} should be same if one channel only
if isempty(lfp_TimeStamps)
    error('lfp_read2:frames2', ...
        'Frame timestamps are irreconcilable! ');
end
% We assume that the first two frames were recorded continuously, and
% calculate the sample period from them:
lfp_SamplesPerFrame = size(lfp_Samples{1}, 1);
%Convert lfp_TimeStamps from usec to seconds:
lfp_TimeStamps = 1.0e-6 * lfp_TimeStamps;       %THIS IS WHERE TIMESTAMPS GET CHANGED
%If convert back to microseconds make sure to use round to get back to
%integer and remove truncation error, because 0.1 is infinite binary 
% Find the sample period, ignoring timestamp intervals that
% indicate a break in recording:
framedurs = lfp_TimeStamps(2:end)-lfp_TimeStamps(1:end-1);
minfd = min(framedurs);
medianfd = median(framedurs);
mediansample = medianfd/lfp_SamplesPerFrame;
%framedurs = framedurs(framedurs < 2*minfd);        %helen took out 10/2018
badframes=find(framedurs<medianfd - 50*mediansample);
%threshold anything lower than average frame duration - 50*average sig
%duration
if std(framedurs(framedurs < minfd + 1.01*mediansample)) > 1e-6
warning('lfp_read2:frames1', ...
    'CSC frames have sub-sample jitter.');
end
if std(framedurs) > mediansample || medianfd > minfd + mediansample
warning('lfp_read2:frames2', ...
    'CSC frames are of excessively variable duration.');
end
lfp_SamplePeriod = median(framedurs) / lfp_SamplesPerFrame;
lfp_SampleRate = 1/lfp_SamplePeriod;
%get vector of timestamps relative to non-frame samples relative to
%absolute time stamp from neuralynx system
allrawTS = repmat(lfp_TimeStamps, lfp_SamplesPerFrame, 1) + ...
    repmat((0 : lfp_SamplesPerFrame - 1)' * lfp_SamplePeriod, 1, length(lfp_TimeStamps));
%take frame timestamp and then add from 0 - 511 times the sample period so
%you get individual timestamps for each sample, but resets to individual
%time stamps at top of individual columns
allsubTS = reshape(allrawTS, 1, []);

if ~isempty(badframes)
    warning([num2str(length(badframes)) ' bad frames to nan ']);
    
for badid = 1:length(badframes)
for refidx = 1:length(chnums)
    %nan all bad frames signals in both preceding and current frame
    samples{refidx}(badframes(badid)*lfp_SamplesPerFrame + ...
        (-lfp_SamplesPerFrame-1:lfp_SamplesPerFrame-1)) =nan;
end
%tsnan=allsubTS(badframes(badid)*lfp_SamplesPerFrame + (0:lfp_SamplesPerFrame-1));
%nan time points of bad frame preceding bad ts & current
allsubTS(badframes(badid)*lfp_SamplesPerFrame +...
    (-lfp_SamplesPerFrame-1:lfp_SamplesPerFrame-1))=nan;
end
lengthnan=length(find(isnan(allsubTS)));
warning([num2str(round(lengthnan/lfp_SampleRate*1000)) ' ms nanned']);    

end

TS=allsubTS;

end

