function psdplot(S,popts)
%PSDPLOT

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

if nargin < 2
    popts = struct();
end

opts = tflabplot_options(S, popts, 'raw', 'psdplot');

if iscell(S) && length(S) == 1
    S = S{1};
end

if iscell(S)
    for s = 1:length(S)
        S{s} = defaultinfo(S{s});
    end
    % TODO: Check all same.
    timeunit = S{1}.Options.info.timeunit;
    timedelta = S{1}.Options.info.timedelta;
else
    S = defaultinfo(S);
    timeunit = S.Options.info.timeunit;
    timedelta = S.Options.info.timedelta;
end

if strcmp(opts.type,'error')
    if iscell(S)
        for s = 1:length(S)
            fe{s} = S{s}.Metrics.PSD.Raw.fe;
            y1{s} = sqrt(S{s}.Metrics.PSD.Raw.Error);
            y2{s} = S{s}.Metrics.SN.Raw;
            y2{s} = (180/pi)*angle(S{s}.Metrics.DFT.Raw.Error);            
        end
    else
        fe = S.Metrics.PSD.Raw.fe;
        y1 = sqrt(S.Metrics.PSD.Raw.Error);
        y2 = (180/pi)*angle(S.Metrics.DFT.Raw.Error);            
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
            y2{s} = (180/pi)*angle(S{s}.Metrics.DFT.Smoothed.Error);            
        end
    else
        if ~isfield(S.Metrics.PSD,'Smoothed')
            logmsg('PSD was not smoothed. Not plotting smoothed PSDs.\n');
            return
        end        
        fe = S.Metrics.PSD.Smoothed.fe;
        y1 = sqrt(S.Metrics.PSD.Smoothed.Error);
        y2 = (180/pi)*angle(S.Metrics.DFT.Raw.Error);            
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

figprep();

if ~any(strcmp(opts.type,{'error','error-smoothed'}))
    
    [lg1, lg2] = legend_(S);
    [yl1, yl2] = ylabel_(S);
    
    subplot('Position',opts.PositionTop);
        plot_(x,y1,opts)
        grid on;box on;
        if endsWith(opts.type,'phase')
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
            ylabel('$\phi$ $[^\circ]$');            
        else
            set(gca,'YScale','log');
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if ~isempty(lg1)
            legend(lg1,opts.legend{:});
        end
        ylabel(yl1);
        tflab_title(S, opts, 'psd');
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(opts, 0, timeunit);

    subplot('Position',opts.PositionBottom);
        plot_(x,y2,opts)
        grid on;box on;
        set(gca,'YScale','log');
        if endsWith(opts.type,'phase')
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
            ylabel('$\phi$ $[^\circ]$');
        else
            set(gca,'YScale','log');
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if ~isempty(lg2)
            legend(lg2,opts.legend{:});
        end
        ylabel(yl2);
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(opts, 1, timeunit);        
else
    [~,lg2] = legend_(S);
    subplot('Position',opts.PositionTop);
        plot_(x,y1,opts);
        grid on;box on;
        set(gca,'YScale','log');
        if opts.vs_period
            set(gca,'XScale','log');
        end
        grid on;box on;
        outunit = '';
        if iscell(S)
            outunit = S{1}.Options.info.outunit;
            legend(lg2,opts.legend{:});
        else
            outunit = S.Options.info.outunit;
            tflab_title(S,opts,'psd');
        end
        if ~isempty(outunit)
            outunit = sprintf(' [%s]',outunit);
        end
        tflab_title(S, opts, 'psd');        
        ylabel(['$|\widetilde{\mbox{Error}}|$', outunit]);
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(opts,0,timeunit);
        
    subplot('Position',opts.PositionBottom);
        plot_(x,y2,opts);
        grid on;box on;        
        if opts.vs_period
            set(gca,'XScale','log');
        end
        set(gca,'YScale','linear');
        set(gca,'YTick',-180:45:180);
        ylabel('$\phi$ of Error')
        if iscell(S)
            legend(lg2,opts.legend{:});
        end
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
            info = S{s}.Options.info;
            lg1{s} = S{s}.Options.description;
            lg2{s} = S{s}.Options.description;
        end
    else
        info = S.Options.info;
        lg1 = '';
        lg2 = '';
        if size(S.In) > 1
            if length(info.instr) == 1
                lg1 = labelstr_(info.instr);
            else
                for j = 1:length(info.instr)
                    lg1{j} = labelstr_(info.instr{j});
                end
            end
            if length(info.outstr) == 1
                lg2 = labelstr_(info.outstr);
            else
                for j = 1:length(info.outstr)
                    lg2{j} = labelstr_(info.outstr{j});
                end
            end
        end
    end    
end

function [yl1, yl2] = ylabel_(S)

    if iscell(S)
        info = S{1}.Options.info;
        yl1 = [labelstr_(info.instr),unitstr_(info.inunit)];
        yl2 = [labelstr_(info.outstr),unitstr_(info.outunit)];
    else
        info = S.Options.info;        
        if size(S.In) > 1
            yl1 = unitstr_(info.inunit);
            yl2 = unitstr_(info.outunit);
        else
            yl1 = [labelstr_(info.instr),unitstr_(info.inunit)];
            yl2 = [labelstr_(info.outstr),unitstr_(info.outunit)];
        end
    end
end

function s = labelstr_(labelstr)
    if contains(labelstr,'$')
        s = ['$|\widetilde{',replace(labelstr,'$',''),'}|$'];
    else
        s = ['$|\widetilde{\mbox{',labelstr,'}}|$'];
    end
end

function s = unitstr_(unit, padstart)
    prepad = ' ';
    if nargin < 2 || padstart == 0
        prepad = '';
    end
    s = '';
    if ~isempty(unit)
        s = sprintf('%s[%s]', prepad, unit);
    end
end

function plot_(x,y,opts)
    if iscell(x) && iscell(y)
        hold on;
        for c = 1:length(x)
            if ~isempty(y)
                plot(x{c},y{c},opts.line{:});
            end
        end
        hold off;
    else
        plot(x,y,opts.line{:});
    end
end

function ls1 = plotnoise(ls1,comp,unit)

    if isfield(S,comp) && strcmp(opts.type,'raw')
        hold on;
        y = ftrim(S.fe,S.(comp));
        loglog(x,y,opts.line{:});
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

