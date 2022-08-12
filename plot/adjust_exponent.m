function adjust_exponent(direction, force, listen)
%ADJUST_EXPONENT - Relabel axes number with exponents
%
%  ADJUST_EXPONENT() Relabels x, y, and z labels
%
%  ADJUST_EXPONENT(s) Relabels s labels only, where s = 'x',
%  'y', or 'z'.
%
%  On log axes, relabels
%     10^{-1} to 0.1
%     10^{0} to 1
%     10^{1} to 10
%     10^{2} to 100
%
%  On linear axes, removes an offsetted x10^{N} if it appears near last
%  axis label and appends $\cdot 10^{N}$ to the last axis label.
%
%  Call adjust exponent after limits are set because MATLAB does not update
%  labels when the limits change if the labels have been modified.

if nargin < 3
    listen = 1;
end

if nargin < 2
    force = 0; 
    % For use if algorthim for detecting offsetted x10^{N} does not work.
end

if nargin == 0
    adjust_exponent('x');
    adjust_exponent('y');
    %adjust_exponent('z'); % Not tested
    return;
end

assert(any(strcmp(direction,{'x','y'})), 'dir must be x or y');

drawnow;

% Relabel 10^{-1} to 0.1, 10^{0} to 1, 10^{1} to 10, and 10^{2} to 100.
if strcmp(get(gca, [direction,'Scale']), 'log')
    labels = get(gca, [direction,'TickLabel']);
    if isempty(labels)
        return;
    end
    found = 0;
    for i = 1:length(labels)
        if strcmp(labels{i}(1),'$')
            found = 1;
        end
        if strcmp(labels{i},'10^{-1}')
            labels{i} = '0.1';
        end
        if strcmp(labels{i},'$10^{-1}$')
            labels{i} = '$0.1$';
        end
        if strcmp(labels{i},'10^{0}')
            labels{i} = '1';
        end
        if strcmp(labels{i},'$10^{0}$')
            labels{i} = '$1$';
        end
        if strcmp(labels{i},'10^{1}')
            labels{i} = '10';
        end
        if strcmp(labels{i},'$10^{1}$')
            labels{i} = '$10$';
        end
        if strcmp(labels{i},'10^{2}')
            labels{i} = '100';
        end
        if strcmp(labels{i},'$10^{2}$')
            labels{i} = '$100$';
        end
    end
    if found
        % Only change if found = 1. If range of y is < 10, no exponents are
        % used to label each tick. If there was an offset x10^{N} shown,
        % doing the following would drop it.
        set(gca,[direction,'TickLabel'], labels);
    end
end

% Remove the offsetted x10^{N} notation that appears on last axis label
% and add it to the last (top) label.
if strcmp(get(gca, [direction,'Scale']),'linear')
    labels = get(gca, [direction,'TickLabel']);
    ticks = get(gca, [direction,'Tick']);
    if isempty(labels)
        return;
    end
    ax = get(gca,[direction,'Axis']);
    if length(ax) > 1
        % If yyaxis used, ax will have two elements.
        % This will only adjust left axis.
        % TODO: Loop over both.
        ax = ax(1);
    end
    if isprop(ax,'Exponent') && ax.Exponent ~= 0
        % Newer versions of MATLAB (when?)
        force = 1;
        ed = ax.Exponent;
    else
        ed = NaN;
    end
    if force || ticks(end) >= 1000 || ticks(end) <= 0.001 % Check 1
        % There does not seem to be a direct way of determining if the
        % offset notation is used (or what it is) in older versions of
        % MATLAB, so Check 1 and Check 2 are used then.

        if ~iscell(labels)
            for i = 1:length(ticks)
                labelsc{i} = labels(i,:);
            end
            labels = labelsc;
        end
        if contains(labels{end},'$')
            return
        end
        if isnan(ed)
            r = abs(ticks(end)/str2double(labels{end}));
            if force == 0
                if (r < 1.1 && r > 0.9) % Check 2.
                    % E.g., ticks(end) = 2000 and labels{end} = '2';
                    return;
                end
            end
        
            % Exponent digit
            ed = floor(log10(r));
            if ed == 0
                return
            end
            if abs(r-10^ed) > eps
                % We have exact value already
                warning('Relabeling failed')
                fprintf('Original: %.6e, New: %.6e\n',r,10^ed);
                return;
            end
        end
        for i = 1:length(ticks)-1
            labels_new{i} = sprintf('%s', labels{i});
        end
        labels_new{i+1} = sprintf('%s$\\cdot 10^{%d}$', labels{i+1}, ed);
        set(gca, [direction,'TickLabel'], labels_new);
        if direction == 'x'
            set(get(gca,'XLabel'),'HorizontalAlignment','center')
        end
        if direction == 'y'
            set(get(gca,'YLabel'),'VerticalAlignment','top')
        end
    end
end

drawnow

if listen
    % On zoom, compute default tick labels.
    % Based on
    % https://blogs.mathworks.com/loren/2015/12/14/axes-limits-scream-louder-i-cant-hear-you/

    % Ideally we would determine if this listener has already been attached
    % to gca and then return if it has. Then we would not need to pass the
    % 'listen' argument. To do this, would need to loop through the
    % AutoListeners__ cell array (undocumented). See commented out code below
    % and https://stackoverflow.com/questions/48345311/find-and-delete-listeners
    addlistener(gca, [upper(direction),'Lim'], 'PostSet', @(a,b)reset(a,b));
    %tmp = gca;
    %tmp.AutoListeners__
end

function reset(~,~)
    set(gca,[upper(direction), 'TickLabelMode'],'auto');
    adjust_exponent(direction, force, 0);
end
end