function dg_writeCSC(filename, TS, Samples, Hdr, varargin)
%dg_writeCSC(filename, TS, Samples, Hdr)
% INPUTS
%   TS:  Timestamps in clock ticks (units as recorded).
%   Samples:  Guess what!
%   Hdr:  Neuralynx format header (<16kB); may be omitted or empty.
% OPTIONS
% 'nlxjunk', nlxjunk - <nlxjunk> is a struct as returned from 'dg_readCSC',
%   and it is used to populate all the other CSC fields that for our usual
%   purposes are junk.  This option is only available in Windows.
% NOTES
%   The dg_write* series of functions for creating Neuralynx format files
% only works reliably on Windows, and this one doesn't even try.  The
% dg_save* series can be used to save .mat files in lfp_lib-compatible
% format.

%$Rev: 278 $
%$Date: 2021-10-04 18:41:51 -0400 (Mon, 04 Oct 2021) $
%$Author: dgibson $

if nargin < 4 || isempty(Hdr)
    Hdr = {'######## Neuralynx Data File Header '
        sprintf('## File Name: %s ', filename)
        sprintf('## Time Opened: %s ', datestr(now))
        '## written by dg_writeCSC '
        ' '
        };
end

nlxjunk = [];
argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'nlxjunk'
            argnum = argnum + 1;
            nlxjunk = varargin{argnum};
        otherwise
            error('dg_writeCSC:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

if ~ispc
    error('dg_writeCSC:notPC', ...
        'Neuralynx format files can only be created on Windows machines.');
end
if isempty(nlxjunk)
    Mat2NlxCSC_411(filename, 0, 1, 1, length(TS), [1 0 0 0 1 1], TS, ...
        Samples, Hdr);
else
    % In 'Mat2NlxCSC_411.m', "Help for file Mat2NlxCSC" it is written:
    % NumFields = 6;
    % Field List
    %     1. Timestamps
    %     2. Channel Numbers
    %     3. Sample Frequency
    %     4. Number of Valid Samples
    %     5. Samples
    %     6. Header
    % ...
    %         Function( Filename, AppendFile, ExtractMode, ModeArray,
    %         NumRecs, FieldSelection, Timestamps, ChannelNumbers,
    %         SampleFrequency, NumberValidSamples, Samples, ...
    %         Header );
    
    Mat2NlxCSC_411(filename, 0, 1, 1, ...
        length(TS), [1 1 1 1 1 1], TS, nlxjunk.ChannelNumbers, ...
        nlxjunk.SampleFrequencies, nlxjunk.NumberValidSamples, Samples, ...
        Hdr);
end


