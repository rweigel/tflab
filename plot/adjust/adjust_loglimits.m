function adjust_loglimits(direction, force)

% Not implemented.

if nargin < 2
    force = 0; 
    % For use if algorthim for detecting offsetted x10^{N} does not catch
    % all cases.
end

if nargin == 0
    adjust_exponent('x');
    adjust_exponent('y');
    adjust_exponent('z');
    return;
end

assert(any(strcmp(direction,{'x','y','z'})), 'dir must be x, y, or z');

if strcmp(get(gca, [direction,'Scale']), 'log')
    labels = get(gca, [direction,'TickLabel']);
    if isempty(labels)
        return;
    end
end