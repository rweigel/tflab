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

debug = 1;

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

pause(0.1)
drawnow;

% Relabel 10^{-1} to 0.1, 10^{0} to 1, 10^{1} to 10, and 10^{2} to 100.
if strcmp(get(gca, [direction,'Scale']), 'log')
    if debug
        fprintf('Log scale\n');
    end
    labels = get(gca, [direction,'TickLabel']);
    if debug
        fprintf('Current %s labels:\n',direction);
        labels
    end
    if isempty(labels)
        if debug
            fprintf('No labels. Returning.\n');
        end
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
    if found == 0 && debug
        fprintf('No labels found with exponential notation. No modifications made.\n');
    end
    if found
        % Only change if found = 1. If range of y is < 10, no exponents are
        % used to label each tick. If there was an offset x10^{N} shown,
        % doing the following would drop it.
        if debug
            fprintf('Setting modified labels:\n');
            labels
        end
        set(gca,[direction,'TickLabel'], labels);
    end
end

% Remove the offsetted x10^{N} notation that appears on last axis label
% and add it to the last (top) label.
if strcmp(get(gca, [direction,'Scale']),'linear')
    if debug
        fprintf('Linear scale for %s\n',direction);
    end
    labels = get(gca, [direction,'TickLabel']);
    if debug
        labels
    end
    if isempty(labels)
        fprintf('No labels. Returning.\n');
        return;
    end
    ticks = get(gca, [direction,'Tick']);
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
        % MATLAB, so Check 1 and Check 2 are used.

        if ~iscell(labels)
            for i = 1:length(ticks)
                labelsc{i} = labels(i,:);
            end
            labels = labelsc;
        end
        if contains(labels{end},'$')
            %return
        end
        if isnan(ed)
            r = abs(ticks(end)/str2double(labels{end}));
            if isnan(r) % Catch 0/0 case
                return;
            end
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
        if debug
            ticks
        end
        for i = 1:length(ticks)-1
            labels_new{i,1} = sprintf('%s', labels{i});
        end
        if debug
            labels_new
        end
        labels_new{i+1,1} = sprintf('%s$\\cdot 10^{%d}$', labels{i+1}, ed);
        if debug
            labels_new
        end
        %fprintf('Setting labels_new in 5 s\n');
        %pause(5)
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
    if debug
        fprintf('Setting Listener for %sLim change.\n',direction);
    end
    ax = gca();
    if isprop(ax.YAxis,'LimitsChangedFcn')
        ax.YAxis.LimitsChangedFcn = @(src,evt)reset(src,evt);
    else
        % For before 2021a when LimitsChangedFcn introduced.
        % When "Restore View" is clicked, ticks will be wrong due to bug:
        % https://www.mathworks.com/matlabcentral/discussions/highlights/134586-new-in-r2021a-limitschangedfcn
        % Uncomment to use.
        %addlistener(gca, [upper(direction),'Lim'], 'PostSet', @(src,evt)reset(src,evt));
    end
end

function reset(obj,src,evt)
    if debug
        fprintf('Reset called. Setting TickLabelMode to auto and deleting listener for %sLim change.\n',direction);
    end
    set(gca,[upper(direction), 'TickLabelMode'],'auto');
    delete(obj);
    %fprintf('Calling adjust_exponent in 5 s\n');
    adjust_exponent(direction, force, 0);
end
end