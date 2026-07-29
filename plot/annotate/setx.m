function setx(opts,last,frequnit,fontsize)
if nargin < 4
    fontsize = 12;
end
% TODO: Support other units.
assert(strcmp(frequnit,'Hz') || strcmp(frequnit,''),...
    sprintf('frequnit must be ''Hz'' or '''', not ''%s''',frequnit));

periodunit = '';
if ~isempty(frequnit)
    periodunit = 's';
end

if isfield(opts,'vs_period') && opts.vs_period
    if isfield(opts,'period_range') && ~isempty(opts.period_range)
        set(gca,'XLim',opts.period_range);
    end
    if ~isempty(periodunit)
        periodunit = sprintf(' [%s]', periodunit);
        period_lines(last, fontsize);
    end
    xlabel(sprintf('$T$%s',periodunit));
else
    if isfield(opts,'frequency_range') && ~isempty(opts.frequency_range)
        set(gca,'XLim',opts.frequency_range);
    end
    if ~isempty(frequnit)
        frequnit = sprintf(' [1/%s]', frequnit);
    end
    xlabel(sprintf('$f$%s',frequnit));
end

% Force xtick labels to be at 1, 10, etc. if log scale.

is_log = get(gca,'XScale');
is_log = strcmp(is_log,'log');
if is_log
    set(gca,'XScale','log');
    xl = get(gca,'XLim');
    xticks = 10.^(floor(log10(xl(1))):ceil(log10(xl(2))));
    xticks = xticks(xticks >= xl(1) & xticks <= xl(2));
    set(gca,'XTick',xticks);
end

if last == 0
    set(gca,'XTickLabel',[]);
    set(gca,'XLabel',[]);
    return
end