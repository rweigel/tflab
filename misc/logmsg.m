function logmsg(fmt, varargin)
%LOGMSG
%
%  LOGMSG(dbstack, FORMAT, ...) Creates a console log message with name of
%  calling file, line number and link, and indentation depending on call
%  stack depth.
%
%  The first argument must always be dbstack.
%
%  FORMAT is modified and the arguments after it are passed to FPRINTF.
%
%  See also LOGMSG_TEST.

% TODO: Allow logmsg(fid, dbstack, fmt, ....) to write to file (and omit
% hyperline).

stack = dbstack;

if length(stack) > 1
    stack = stack(2:end); % Remove entry for logmsg.m
    % Indent one space per call stack entry.
    indent = repmat(' ',1,length(stack)-1);
    str = sprintf('%s (line %d): ', stack(1).file, stack(1).line);
    link = sprintf(...
            '%s%s (<a href="matlab: matlab.desktop.editor.openAndGoToLine(''%s'', %d);">line %d</a>)',...
            indent,...
            stack(1).file,...
            which(stack(1).file),...
            stack(1).line,...
            stack(1).line);
    fmtx = fmt(end);
    % Replace newlines not at end of string with a newline
    % and then indentation.
    fmt = strrep(fmt(1:end-1), '\n',...
                ['\n',repmat(' ', 1, length(str)+length(indent))]);
    fprintf([link, ': ', fmt, fmtx], varargin{:});
else
    % Called from command line.
    fprintf(fmt, varargin{:});
end
