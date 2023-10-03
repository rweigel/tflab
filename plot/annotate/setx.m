function setx(opts,last,frequnit,freqsf)

if nargin < 4
    freqsf = 1;
end

% TODO: Support other units.
assert(strcmp(frequnit,'Hz') || strcmp(frequnit,''),...
    sprintf('frequnit must be ''Hz'' or '''', not ''%s''',frequnit));
periodunit = 's';

if opts.vs_period 
    if ~isempty(opts.period_range)
        set(gca,'XLim',opts.period_range);
    end
    if ~isempty(periodunit)
        periodunit = sprintf(' [%s]', periodunit);
        period_lines(freqsf);
    end
    xlabel(sprintf('$T$%s',periodunit));
else
    if ~isempty(opts.frequency_range)
        set(gca,'XLim',opts.frequency_range);
    end    
    if ~isempty(frequnit)
        frequnit = sprintf(' [1/%s]', frequnit);
    end
    xlabel(sprintf('$f$%s',frequnit));
end

if last == 0
    set(gca,'XTickLabel',[]);
    set(gca,'XLabel',[]);
    return
end

end