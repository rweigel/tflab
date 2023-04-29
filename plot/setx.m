function setx(opts,last,timeunit)

if opts.vs_period 
    if ~isempty(opts.period_range)
        set(gca,'XLim',opts.period_range);
    end
    if ~isempty(timeunit)
        timeunit = sprintf(' [%s]', timeunit);
        period_lines();
    end
    xlabel(sprintf('$T$%s',timeunit));
else
    if ~isempty(opts.frequency_range)
        set(gca,'XLim',opts.frequency_range);
    end    
    if ~isempty(timeunit)
        timeunit = sprintf(' [1/%s]', timeunit);
    end
    xlabel(sprintf('$f$%s',timeunit));
end

if last == 0
    set(gca,'XTickLabel',[]);
    set(gca,'XLabel',[]);
    return
end

end