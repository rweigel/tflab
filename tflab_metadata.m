function S = tflab_metadata(S)

if iscell(S)
    for s = 1:length(S)
        S{s} = tflab_metadata(S{s});
    end
    return
end

meta = struct();

% Needed for calculations

meta.timedelta = 1;  % Measurement cadence.
meta.timeunit  = ''; % Unit for timedelta
meta.frequnit  = ''; % Unit for frequency
meta.freqsf    = 1;  % Multiply frequencies by this scale factor (sf) to get frequnit
meta.zunit     = ''; % Unit for Z

% Needed for all plots.

meta.instr = 'In'; % Cell array or string. 
% Or ,{'$B_x$', ...} (1 cell element per column in In)

meta.outstr = 'Out'; % Cell array or string. 
% Or, {'$E_x$', ...} (1 cell element per column in Out)

% Timestart is a time string of the form
%   'yyyy-mm-ddTHH:MM:SS.FFF'.
%   Example: 
%     timestart = '2001-01-01T00:00:00.000';
% or an integer. 
meta.timestart = 1;
if ~iscell(S.In)
    % TODO: If S.In is cell array, so should time{start/stop}
    meta.timestop = size(S.In,1);
end

meta.inunit    = '';
meta.outunit   = '';

meta.chainid   = '';   % Usualy over-arching project name
meta.stationid = '';   % Usually an abbreviation, e.g., VAQ58

if exist('S','var')
    if isfield(S,'Metadata')
        S.Metadata = updatedefaults(meta,S.Metadata);
    else
        S.Metadata = meta;
    end
else
    S = meta;
end

if ischar(S.Metadata.timestart) && ischar(S.Metadata.timeunit)
    if strcmp(S.Metadata.timeunit,'s')
        startdn = datenum(S.Metadata.timestart,'yyyy-mm-ddTHH:MM:SS.FFF');
        stopdn  = startdn + S.Metadata.timedelta*size(S.In,1)/86400;
        S.Metadata.timestop = datestr(stopdn,'yyyy-mm-ddTHH:MM:SS.FFF');
    end
end
