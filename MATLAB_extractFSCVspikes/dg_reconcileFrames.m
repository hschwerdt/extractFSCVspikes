function [canonicalstamps, activefnums, includedfile] = ...
    dg_reconcileFrames(timestamps, varargin)
% When Cheetah recording stops, it sometimes happens that some channels
% record one frame more than the other channels.  Rodent data may contain
% many starts and stops, so there can be frame discrepancies in the middles
% of files.  We make the files uniform by deleting the extra frames
% wherever they occur.
%INPUT
% timestamps - a cell vector containing the list of timestamps from each of
%   the files that are to be reconciled, with one file's timestamps in each
%   cell.  May also contain empty cells (e.g. from empty files), which are
%   ignored.
%OUTPUTS
% canonicalstamps - the list of timestamps common to all files.  If the
%   vast majority of files have a common set of timestamps but a small
%   minority substantially disagree, then the small minority are skipped
%   unless the 'allfiles' option is specified.  "Vast majority" and
%   "substantially disagree" are both defined in terms of <threshfrac> (see
%   code).
% activefnums - the equivalent of lfp_ActiveFilenums for the current
%   value of <timestamps>, i.e. a list of indices of non-empty cells in
%   <timestamps>.
% includedfile - logical vector of same size as <timestamps>, which is true
%   for elements of <timestamps> that were included in the computation of
%   <canonicalstamps> and false for empty or skipped elements.
%OPTIONS
% 'allfiles' - strictly requires <canonicalstamps> to contain timestamps
%   that are common to literally ALL of the elements of <timestamps>.


%$Rev: 270 $
%$Date: 2020-07-11 14:29:51 -0400 (Sat, 11 Jul 2020) $
%$Author: dgibson $

threshfrac = 0.9; % definition of "vast majority"
argnum = 0;
allfilesflag = false;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'allfiles'
            allfilesflag = true;
        otherwise
            error('dg_reconcileFrames:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

if isempty(timestamps)
    error('dg_reconcileFrames:timestamps', ...
        '<timestamps> is empty.');
end
if ~iscell(timestamps) || ~isvector(timestamps)
    error('dg_reconcileFrames:timestamps2', ...
        '<timestamps> must be a cell vector.');
end

% Find <activefnums>, 
gotstuff = false(size(timestamps));
for k = 1:length(timestamps)
    gotstuff(k) = ~isempty(timestamps{k});
end
activefnums = find(gotstuff);
if length(activefnums) < 1
    canonicalstamps = [];
    return
end

if length(activefnums) == 1
    canonicalstamps = timestamps{activefnums(1)};
    return
end
% Use activefnums(1) by default to seed the reconciliation process:
[canonicalstamps, includedfile] = computecanonical(activefnums, ...
    timestamps, 1, threshfrac, allfilesflag);
if sum(includedfile)/length(activefnums) < threshfrac
    % Too many files are excluded.  This might simply indicate that the
    % first file is bad, so try using a different <seedidx>.  But we want a
    % <seedidx> for a file that has a median number of timestamps, not
    % another defective file that is way too short.
    timestamplengths = cellfun(@length, timestamps(activefnums));
    mednumframes = median(timestamplengths);
    foundone = false;
    for seedidx = 2:length(activefnums)
        if length(timestamps{activefnums(seedidx)}) >= mednumframes
            foundone = true;
            break
        end
    end
    if foundone
        % see if it works better than seedidx = 1:
        [canonicalstamps, includedfile] = computecanonical(activefnums, ...
            timestamps, seedidx, threshfrac, allfilesflag);
    end
    if sum(includedfile)/length(activefnums) < threshfrac
        % The problem is more complicated, and calls for manual examination
        % of the sets of timestamps:
        error( 'dg_reconcileFrames:failed', ...
            'Failed to reconcile frames; %d%% of files are marked bad:\n%s', ...
            round(100 - 100 * sum(includedfile)/length(includedfile)), ...
            dg_canonicalSeries(find(~includedfile)) );
    end
end
end

function [canonicalstamps, includedfile] = computecanonical(activefnums, ...
    timestamps, seedidx, threshfrac, allfilesflag)
% Computes canonical timestamps starting with <seedidx> as the initial
% list. <seedidx> is an index into <activefnums>.
includedfile = false(size(timestamps));
includedfile(activefnums) = true;
canonicalstamps = timestamps{activefnums(seedidx)};
for filenum = activefnums(1:end)
    if ~isempty(timestamps{filenum})
        newcanonicalstamps = intersect(canonicalstamps, ...
            timestamps{filenum});
        if length(newcanonicalstamps) == length(canonicalstamps)
            % There are no mismatches, and no need to update
            % <newcanonicalstamps>.
        else
            if ~allfilesflag && ...
                    ( length(newcanonicalstamps)/length(canonicalstamps) ...
                    < threshfrac )
                % More than 10% of the old <canonicalstamps> fail to match
                % this file; mark this as a bad file and skip it:
                includedfile(filenum) = false;
            else
                % accept the reduction in size of <canonicalstamps>:
                canonicalstamps = newcanonicalstamps;
            end
        end
    end
end
end
