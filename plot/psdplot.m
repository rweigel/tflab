function psdplot(S,popts)
%PSDPLOT

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

% Default options
opts = struct();
    opts.title = '';
    opts.type = 'raw';
    opts.print = 0;
    opts.printname = 'psd';
    opts.printdir = '';
    opts.printfmt = {'pdf'};

    opts.vs_period = 1;
    if iscell(S)
        for s = 1:length(S)
            lens(s) = size(S{s}.In,1);
        end
        opts.period_range = [1, 2*max(lens)];
    else        
        opts.period_range = [1, 2*size(S.In,1)];
    end
    opts.frequency_range = [0, 0.5];
    
% Use default options if options not given
if nargin > 1
    fns = fieldnames(popts);
    for i = 1:length(fns)
        opts.(fns{i}) = popts.(fns{i});
    end
end

% Line options
lnopts = {'marker','.','markersize',20,'linestyle','none'};

% Legend options
lgopts = {'Location','NorthWest','Orientation','Horizontal'};

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];
figprep();

if iscell(S)
    % TODO: Check all same.
    timeunit = S{1}.Options.info.timeunit;
else
    timeunit = S.Options.info.timeunit;
end

if strcmp(opts.type,'error')
    if iscell(S)
        for s = 1:length(S)
            fe{s} = S{s}.Metrics.PSD.Raw.fe;
            y1{s} = sqrt(S{s}.Metrics.PSD.Raw.Error);
            y2{s} = S{s}.Metrics.SN.Raw;
        end
    else
        fe = S.Metrics.PSD.Raw.fe;
        y1 = sqrt(S.Metrics.PSD.Raw.Error);
        y2 = S.Metrics.SN.Raw;
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
            y1{s} = sqrt(S{s}.Metrics.PSD.Smoothed.Error);
            y2{s} = S{s}.Metrics.SN.Smoothed;
        end
    else
        if ~isfield(S.Metrics.PSD,'Smoothed')
            logmsg('PSD was not smoothed. Not plotting smoothed PSDs.\n');
            return
        end        
        fe = S.Metrics.PSD.Smoothed.fe;
        y1 = sqrt(S.Metrics.PSD.Smoothed.Error);
        y2 = S.Metrics.SN.Smoothed;
    end
end
if strcmp(opts.type,'raw')
    if iscell(S)
        for s = 1:length(S)
            fe{s} = S{s}.Metrics.PSD.Raw.fe;
            y1{s} = sqrt(S{s}.Metrics.PSD.Raw.In);
            y2{s} = sqrt(S{s}.Metrics.PSD.Raw.Out);
        end
    else
        fe = S.Metrics.PSD.Raw.fe;
        y1 = sqrt(S.Metrics.PSD.Raw.In);
        y2 = sqrt(S.Metrics.PSD.Raw.Out);
    end
end
if strcmp(opts.type,'smoothed')
    if iscell(S)
        for s = 1:length(S)
            fe{s} = S{s}.Metrics.PSD.Smoothed.fe;
            y1{s} = sqrt(S{s}.Metrics.PSD.Smoothed.In);
            y2{s} = sqrt(S{s}.Metrics.PSD.Smoothed.Out);
        end
    else
        fe = S.Metrics.PSD.Smoothed.fe;
        y1 = sqrt(S.Metrics.PSD.Smoothed.In);
        y2 = sqrt(S.Metrics.PSD.Smoothed.Out);
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
                logmsg(sprintf('Time series %d was not %s. Plotting raw Fourier amplitudes.\n',s,opts.type));
                fe{s} = S{s}.Metrics.PSD.Raw.fe;
                y1{s} = sqrt(S{s}.Metrics.PSD.Raw.In);
                y2{s} = sqrt(S{s}.Metrics.PSD.Raw.Out);
            else
                cnt = cnt + 1;
                fe{s} = S{s}.(el1).PSD.(el2).fe;
                y1{s} = sqrt(S{s}.(el1).PSD.(el2).In);
                y2{s} = sqrt(S{s}.(el1).PSD.(el2).Out);
            end
        end
        if cnt == 0
            return
        end
    else
        if ~isfield(S,el1) || ~isfield(S.(el1).PSD,el2)
            logmsg(sprintf('Time series was not %s. Not plotting %s Fourier amplitudes.\n',opts.type,opts.type));
            return
        end
        fe = S.(el1).PSD.(el2).fe;
        y1 = sqrt(S.(el1).PSD.(el2).In);
        y2 = sqrt(S.(el1).PSD.(el2).Out);
    end
end

if startsWith(opts.type,'error')
    y1 = ftrim(fe,y1);
    y2 = ftrim(fe,y2);
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

if any(strcmp(opts.type,...
      {'raw','raw-phase','windowed','prewhitened','zeropadded','smoothed'}))
    
    [lg1, lg2] = legend_(S);
    
    subplot('Position',PositionTop);
        plot_(x,y1,lnopts)
        if endsWith(opts.type,'phase')
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
            ylabel('$[^\circ]$');            
        else
            ylabel('$|$Fourier Amplitude$|$');
            set(gca,'YScale','log');
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if iscell(S)
            ylabel(sprintf('In %s', S{1}.Options.info.inunit));
        else
            lg1 = plotnoise(lg1,'InNoisePSD', S.Options.info.inunit);
            tflab_title(S,opts,'psd');
        end
        grid on;box on;
        legend(lg1,lgopts{:});
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(opts, 0, timeunit);

    subplot('Position',PositionBottom);
        plot_(x,y2,lnopts)
        set(gca,'YScale','log');
        if endsWith(opts.type,'phase')
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
            ylabel('$[^\circ]$');
        else
            ylabel('$|$Fourier Amplitude$|$');
            set(gca,'YScale','log');
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if iscell(S)
            ylabel(sprintf('Out %s', S{1}.Options.info.inunit));
        else
            lg2 = plotnoise(lg2,'OutNoisePSD', S.Options.info.inunit);
        end
        grid on;box on;
        legend(lg2,lgopts{:});
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(opts, 1, timeunit);        
else
    subplot('Position',PositionTop)
        [~,lg2] = legend_(S);
        plot_(x,y1,lnopts)
        set(gca,'YScale','log');
        if opts.vs_period
            set(gca,'XScale','log');
        end
        grid on;box on;
        ylabel('$|$Fourier Amp of Error$|$')
        legend(lg2,lgopts{:});
        tflab_title(S,opts,'psd');
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(opts,0,timeunit);
        
    subplot('Position',PositionBottom)
        plot_(x,y2,lnopts);
        set(gca,'YScale','log');
        if opts.vs_period
            set(gca,'XScale','log');
            set(gca,'YScale','log');
        end
        ylabel('SNR');
        grid on;
        legend(lg2,lgopts{:});
        adjust_ylim();
        adjust_exponent('y');
        setx(opts,1,timeunit);

end

if opts.print
    for i = 1:length(opts.printfmt)
        fname = sprintf('%s_%s.%s', opts.printname, opts.type, opts.printfmt{i});
        figsave(fullfile(opts.printdir, fname));
    end
end

function [lg1, lg2] = legend_(S)

    if iscell(S)
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
            lg2{s} = S{s}.Options.description;
        end
        lg1 = lg1(~cellfun('isempty',lg1));
        lg2 = lg2(~cellfun('isempty',lg2));
    else        
        info = S.Options.info;
        inunit = '';
        if ~isempty(info.inunit)
            inunit = sprintf(' [%s]', info.inunit);
        end
        for j = 1:size(S.In,2)
            if iscell(info.instr)
                lg1{j} = sprintf('%s %s',info.instr{j},inunit);
            else
                lg1{j} = sprintf('%s(:,%d) %s',info.instr,j,inunit);
            end
        end
        outunit = '';
        if ~isempty(info.outunit)
            outunit = sprintf(' [%s]', info.outunit);
        end
        for j = 1:size(S.Out,2)
            if iscell(info.outstr)
                lg2{j} = sprintf('%s %s',info.outstr{j},outunit);
            else
                lg2{j} = sprintf('%s(:,%d) %s',info.outstr,j,outunit);
            end
        end
    end
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
                ls1{j+jl} = sprintf('%s(:,%d) Noise Amplitudes%s',compstrs{j},j,unit);
            else
                ls1{j+jl} = sprintf('%s Noise Amplitudes%s',compstrs,unit);
            end
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

end

