function axes_handles = spectrum_plot(S,popts)
%SPECTRUM_PLOT

% Default options
opts = struct();
    opts.type = 'raw';
    opts.title = '';
    opts.vs_period = 0;
    if iscell(S)
        for s = 1:length(S)
            lens(s) = size(S{s}.In,1);
        end
        opts.period_range = [1, 2*max(lens)];
    else        
        opts.period_range = [1, 2*size(S.In,1)];
    end
    opts.frequency_range = [0, 0.5];
    opts.filename = 'spectrum';
    opts.savefmt = {}; % One or more extensions allowed by export_fig
                       % e.g. {'pdf'} or {'svg','png'}.
    
% Use default options if options not given
if nargin > 1
    fns = fieldnames(popts);
    for i = 1:length(fns)
        opts.(fns{i}) = popts.(fns{i});
    end
end

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];
figprep();

% Title string
if ~isempty(opts.title)
    ts = opts.title;
else
    ts = '';
    if strcmp(opts.type,'raw')
        ts = '';
    end
    if strcmp(opts.type,'error')
        ts = '';
    end
    if strcmp(opts.type,'error-smoothed')
        ts = 'Smoothed (in FD)';
    end
    if strcmp(opts.type,'zeropadded')
        if isfield(S,'Zeropad')
            ts = sprintf('Padded with %d zeros',S.Options.td.zeropad);
        end
    end    
    if strcmp(opts.type,'windowed')
        if isfield(S,'Window')
            ts = [S.Options.td.window.functionstr, '-windowed (in TD)'];
        end
    end
    if strcmp(opts.type,'smoothed')
        ts = 'Smoothed (in FD)';
    end
end

% Line options
lnopts = {'marker','.','markersize',10,'linewidth',2};
% Legend options
lgopts = {'Location','NorthWest','Orientation','Horizontal'};

if strcmp(opts.type,'error')
    if iscell(S)
        for s = 1:length(S)
            fe{s} = S{s}.Metrics.PSD.Raw.fe;
            y1{s} = S{s}.Metrics.PSD.Raw.Error;
            y2a{s} = S{s}.Metrics.SN.Raw;
            y2b{s} = S{s}.Metrics.Coherence.Raw;
        end
    else
        fe = S.Metrics.PSD.Raw.fe;
        y1 = S.Metrics.PSD.Raw.Error;
        y2a = S.Metrics.SN.Raw;
        y2b = S.Metrics.Coherence.Raw;
    end
end
if strcmp(opts.type,'error-smoothed')
    if iscell(S)
        if ~isfield(S{1}.Metrics.PSD,'Smoothed')
            logmsg('First PSD was not smoothed. Not plotting smoothed PSDs.\n');
            return
        end                
        for s = 1:length(S)
            fe{s} = S{s}.Metrics.PSD.Smoothed.fe;
            y1{s} = S{s}.Metrics.PSD.Smoothed.Error;
            y2a{s} = S{s}.Metrics.SN.Smoothed;
            y2b{s} = S{s}.Metrics.Coherence.Smoothed;
        end
    else
        if ~isfield(S.Metrics.PSD,'Smoothed')
            logmsg('PSD was not smoothed. Not plotting smoothed PSDs.\n');
            return
        end        
        fe = S.Metrics.PSD.Smoothed.fe;
        y1 = S.Metrics.PSD.Smoothed.Error;
        y2a = S.Metrics.SN.Smoothed;
        y2b = S.Metrics.Coherence.Smoothed;
    end
end
if strcmp(opts.type,'raw')
    if iscell(S)
        for s = 1:length(S)
            fe{s} = S{s}.Metrics.PSD.Raw.fe;
            y1{s} = S{s}.Metrics.PSD.Raw.In;
            y2{s} = S{s}.Metrics.PSD.Raw.Out;
        end
    else
        fe = S.Metrics.PSD.Raw.fe;
        y1 = S.Metrics.PSD.Raw.In;
        y2 = S.Metrics.PSD.Raw.Out;
    end
end
if strcmp(opts.type,'raw-phase')
    if iscell(S)
        for s = 1:length(S)
            fe{s} = S{s}.Metrics.DFT.Raw.fe;
            y1{s} = (180/pi)*angle(S{s}.Metrics.DFT.Raw.In);
            y2{s} = (180/pi)*angle(S{s}.Metrics.DFT.Raw.Out);
        end
    else
        fe = S.Metrics.PSD.Raw.fe;
        y1 = (180/pi)*angle(complex(S.Metrics.DFT.Raw.In));
        y2 = (180/pi)*angle(complex(S.Metrics.DFT.Raw.Out));
    end
end

if any(strcmp(opts.type,{'windowed','prewhitened','zeropadded','smoothed'}))
    if strcmp(opts.type,'windowed')
        el1 = 'Window';
        el2 = 'Raw';
    elseif strcmp(opts.type,'prewhitened')
        el1 = 'Prewhiten';
        el2 = 'Raw';
    elseif strcmp(opts.type,'zeropadded')
        el1 = 'Zeropad';
        el2 = 'Raw';
    else
        el1 = 'Metrics';    
        el2 = 'Smoothed';
    end
    if iscell(S)
        cnt = 0;
        for s = 1:length(S)
            if ~isfield(S{s},el1) || ~isfield(S{s}.(el1).PSD,el2)
                logmsg(sprintf('Time series %d was not %s. Plotting raw PSD.\n',s,opts.type));
                fe{s} = S{s}.Metrics.PSD.Raw.fe;
                y1{s} = S{s}.Metrics.PSD.Raw.In;
                y2{s} = S{s}.Metrics.PSD.Raw.Out;
            else
                cnt = cnt + 1;
                fe{s} = S{s}.(el1).PSD.(el2).fe;
                y1{s} = S{s}.(el1).PSD.(el2).In;
                y2{s} = S{s}.(el1).PSD.(el2).Out;
            end
        end
        if cnt == 0
            return
        end
    else
        if ~isfield(S,el1) || ~isfield(S.(el1).PSD,el2)
            logmsg(sprintf('Time series was not %s. Not plotting %s PSDs.\n',opts.type,opts.type));
            return
        end
        fe = S.(el1).PSD.(el2).fe;
        y1 = S.(el1).PSD.(el2).In;
        y2 = S.(el1).PSD.(el2).Out;
    end
end

if startsWith(opts.type,'error')
    y1 = ftrim(fe,y1);
    y2a = ftrim(fe,y2a);
    y2b = ftrim(fe,y2b);
else
    y1 = ftrim(fe,y1);
    y2 = ftrim(fe,y2);
end

if iscell(S)
    for s = 1:length(S)
        if opts.vs_period
            x{s} = S{s}.Options.info.timedelta./fe{s};
        else
            x{s} = fe{s};
        end
    end
else
    if opts.vs_period
        x = S.Options.info.timedelta./fe;
    else
        x = fe;
    end
end

if any(strcmp(opts.type,{'raw','raw-phase','windowed','prewhitened','zeropadded','smoothed'}))

    if iscell(S)
        j = 1; % Only plots first component (need to generalize).
        for s = 1:length(S)
            if isempty(x{s})
                continue;
            end
            info = S{s}.Options.info;
            inunit = '';
            if ~isempty(info.inunit)
                inunit = sprintf(' [%s]', info.inunit);
            end
            lg1{s} = S{s}.Options.description;
            if iscell(info.instr)
                %lg1{s} = sprintf('%s(:,%d) PSD%s%s',info.instr{j},j,inunit);
            else
                %lg1{s} = sprintf('%s PSD%s',info.instr,inunit);
            end
            outunit = '';
            if ~isempty(info.outunit)
                outunit = sprintf(' [%s]', info.outunit);
            end
            lg2{s} = S{s}.Options.description;
            if iscell(info.outstr)
                %lg2{s} = sprintf('%s(:,%d) PSD%s',info.outstr{j},j,outunit);
            else
                %lg2{s} = sprintf('%s PSD%s',info.outstr,outunit);
            end
        end
        lg1 = lg1(~cellfun('isempty',lg1));
        lg2 = lg2(~cellfun('isempty',lg2));
    else
        info = S.Options.info;
        inunit = '';
        if ~isempty(info.inunit)
            inunit = sprintf(' [%s]', info.inunit);
        end
        lab = 'PSD';
        if endsWith(opts.type,'phase')
            lab = 'Phase';
        end
        for j = 1:size(y1,2)
            if iscell(info.instr)
                lg1{j} = sprintf('%s(:,%d) %s%s',info.instr{j},j,lab,inunit);
            else
                lg1{j} = sprintf('%s %s%s',info.instr,lab,inunit);
            end
        end
        outunit = '';
        if ~isempty(info.outunit)
            outunit = sprintf(' [%s]', info.outunit);
        end
        for j = 1:size(y2,2)
            if iscell(info.outstr)
                lg2{j} = sprintf('%s(:,%d) %s%s',info.outstr{j},j,lab,outunit);
            else
                lg2{j} = sprintf('%s %s%s',info.outstr,lab,outunit);
            end
        end
    end
    
    subplot('Position',PositionTop);
        plot_(x,y1,lnopts)
        if endsWith(opts.type,'phase')
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
        else
            set(gca,'YScale','log');
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if ~iscell(S)
            lg1 = plotnoise(lg1,'InNoisePSD', inunit);
        end
        if iscell(S)
            ylabel(sprintf('In PSD %s', inunit));
        end
        title(ts,'FontWeight','Normal');
        grid on;box on;
        legend(lg1,lgopts{:});
        setx(opts,0);
        adjust_ylim('both');
        adjust_yticks();
        adjust_exponent('y');
    subplot('Position',PositionBottom);
        plot_(x,y2,lnopts)
        set(gca,'YScale','log');
        if endsWith(opts.type,'phase')
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
        else
            set(gca,'YScale','log');
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if ~iscell(S)
            lg2 = plotnoise(lg2,'OutNoisePSD', outunit);
        end        
        if iscell(S)
            ylabel(sprintf('Out PSD %s', outunit));
        end
        grid on;box on;
        legend(lg2,lgopts{:});
        setx(opts, 1);
        adjust_ylim('both');
        adjust_yticks();
        adjust_exponent('y');
else
    subplot('Position',PositionTop)
        if iscell(S)
            j = 1; % Only plots first component (need to generalize).
            for s = 1:length(S)
                info = S{s}.Options.info;
                outunit = '';
                if ~isempty(info.inunit)
                    outunit = sprintf(' [%s]', info.outunit);
                end
                ls{s} = S{s}.Options.description;
                if ~iscell(info.instr) && size(S{s}.Metrics.PSD.Raw.Error,2) > 1
                    %ls{s} = sprintf('Error(:,%d) PSD%s%s',j,outunit);
                else
                    %ls{s} = sprintf('Error PSD%s',outunit);
                end
            end            
        else
            info = S.Options.info;
            outunit = '';
            if ~isempty(info.inunit)
                outunit = sprintf(' [%s]', info.outunit);
            end
            for j = 1:size(S.Metrics.PSD.Raw.Error,2)
                if ~iscell(info.instr) && size(S.Metrics.PSD.Raw.Error,2) > 1
                    ls{j} = sprintf('Error(:,%d) PSD%s%s',j,outunit);
                else
                    ls{j} = sprintf('Error PSD%s',outunit);
                end
            end
        end
        plot_(x,y1,lnopts)
        set(gca,'YScale','log');
        if opts.vs_period
            set(gca,'YScale','log');
        end
        title(ts,'FontWeight','Normal');
        grid on;box on;
        legend(ls,lgopts{:});
        if iscell(S)
            if size(S{s}.Metrics.PSD.Raw.Error,2) > 1
                ylabel('PSD of Error(:,1)');
            else
                ylabel('PSD of Error');                
            end
        end
        setx(opts,0);
        adjust_ylim();
        adjust_yticks();
        adjust_exponent('y');
    subplot('Position',PositionBottom)
        if iscell(S)
            q = 1;
            for k = 1:2
                for s = 1:length(S)
                    ls2{q} = S{s}.Options.description;
                    q = q + 1;
                end
            end
        end
        yyaxis left
        plot_(x,y2a,lnopts);
        set(gca,'YScale','log');
        if opts.vs_period
            set(gca,'YScale','log');
        end
        if iscell(S)
            ylabel('[PSD Out]/[PSD Error]');
        else
            ls2{1} = '[PSD Out]/[PSD Error]';
        end
        adjust_ylim();
        %adjust_yticks();
        adjust_exponent('y');
        
        yyaxis right
        set(gca,'YScale','linear');
        plot_(x,y2b,lnopts);
        if opts.vs_period
            set(gca,'YScale','log');
        end
        set(gca,'YLim',[-0.01,1.19]);
        if iscell(S)
            ylabel('Coherence');
        else
            ls2{2} = 'Coherence';
        end
        legend(ls2,lgopts{:});
        grid on;box on;
        setx(opts,1);

end

function plot_(x,y,lnopts)
    if iscell(x) && iscell(y)
        hold on;
        for c = 1:length(x)
            if ~isempty(y)
                plot(x{c},y{c},lnopts{:});
            end
        end
        hold off;
    else
        plot(x,y,lnopts{:});
    end
end

for i = 1:length(opts.savefmt)
    figsave([opts.filename,'.',opts.savefmt{i}]);
end

function ls1 = plotnoise(ls1,comp,unit)

    if isfield(S,comp) && strcmp(opts.type,'raw')
        hold on;
        y = ftrim(S.fe,S.(comp));
        loglog(x,y,lnopts{:});
        jl = size(S.(comp),2);
        if strcmp(comp,'InNoisePSD')
            compstrs = info.instr;
        else
            compstrs = info.outstr;
        end
        for j = 1:jl
            if iscell(info.instr)
                ls1{j+jl} = sprintf('%s(:,%d) Noise PSD%s',compstrs{j},j,unit);
            else
                ls1{j+jl} = sprintf('%s Noise PSD%s',compstrs,unit);
            end
        end
    end
    
end

function setx(opts,last)
    if opts.vs_period && ~isempty(opts.period_range)
        set(gca,'XLim',opts.period_range);
    end
    if ~opts.vs_period && ~isempty(opts.frequency_range)
        set(gca,'XLim',opts.frequency_range);
    end
    if last == 0
        set(gca,'XTickLabel',[]);
        return
    end
    timeunit = '';
    if opts.vs_period 
        if ~isempty(info.timeunit)
            timeunit = sprintf(' [%s]', info.timeunit);
        end
        xlabel(sprintf('$T$%s',timeunit));
    else
        if ~isempty(info.timeunit)
            timeunit = sprintf(' [1/%s]', info.timeunit);
        end
        xlabel(sprintf('$f$%s',timeunit));
    end
end

end

function y = ftrim(fe,y)
    if iscell(y)
        % For 1-D case, f = 0 and f = 0.5 are always computed
        % using only one frequency. As a result, in-sample error will
        % be near zero if Z is computed using Z(f) = E(f)/B(f) where
        % E(f) and B(f) are FFTs of time series with no windowing,
        % zero padding, or other processing steps that modify the 
        % Fourier coefficients. These frequencies are set to NaN
        % so that large values do not make the scale span many 
        % order of magnitudes.
        for s = 1:length(y)
            if fe{s}(1) == 0
                y{s}(1,:)  = NaN;
            end
            if fe{s}(end) == 0.5
                y{s}(end,:)  = NaN;
            end
        end
    else
        % See comment above for motivation for NaNs.
        if fe(1) == 0
            y(1,:) = NaN;
        end
        if fe(end) == 0.5
            y(end,:) = NaN;
        end    
    end
end