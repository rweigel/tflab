function ax = tsplot(S,popts)
%TSPLOT - Plot timeseries in output of TRANSFERFNFD.
%
%  TSPLOT(S), where S is the output of TRANSFERFNFD.
%  TSPLOT(S, opts)

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

if nargin < 2
    popts = struct();
end

opts = tflabplot_options(S, popts, 'raw', 'tsplot');

if iscell(S) && length(S) == 1
    S = S{1};
end

if iscell(S) && length(S) == 1
    S = S{1};
end

if iscell(S)
    info = S{1}.Options.info;
    t = S{1}.Time;
    if length(S) > 1
        for j = 2:length(S)
            if ~all(size(S{1}.Time) == size(S{j}.Time))
                error('Time arrays must have same size.\n');
            end
            d = S{1}.Time == S{j}.Time;
            if ~all(d(:))
                error('Time arrays must have same content.\n');
            end
        end
    end
else
    info = S.Options.info;
    t = S.Time;
end

if ~isempty(info.timestart)
    try
        fmt = 'yyyy-mm-ddTHH:MM:SS.FFF';
        to = datenum(info.timestart,fmt);
    catch
        warning(['Could not parse Options.info.timestart. Format must be ',fmt]);
        info.timestart = '';
    end
end

if ~isempty(info.timestart)
    dt = info.timedelta;
    t = [0:length(t)-1]';
    if strcmp(info.timeunit,'ms')
        ppd = 86400000/info.timedelta;
    elseif strcmp(info.timeunit,'s')
        ppd = 86400/info.timedelta;
    elseif strcmp(info.timeunit,'m')
        ppd = 1440/info.timedelta;
    else
        warning('Options.td.timeunit = %s not recognized. Must be ms, s, or m\n',...
            info.timeunit);
        info.timestart = '';
    end
    t = to + t/ppd;
end

figprep();

if ~iscell(S) && any(strcmp(opts.type,{'raw','windowed','prewhitened'}))

    if strcmp(opts.type,'raw')
        In = S.In;
        Out = S.Out;
    end
    
    if strcmp(opts.type,'windowed')
        if ~isfield(S,'Window')
            logmsg('Data were not windowed. Not plotting windowed timeseries.\n');
            return
        end
        In = S.Window.In;
        Out = S.Window.Out;
    end
    if strcmp(opts.type,'prewhitened')
        if ~isfield(S,'Prewhiten')
            logmsg('Data were not prewhitened. Not plotting prewhitened timeseries.\n');
            return
        end
        In = S.Prewhiten.In;
        Out = S.Prewhiten.Out;
    end
    
    ax(1) = subplot('Position',opts.PositionTop);
        plot(t,In);
        grid on;grid minor;box on;
        lg = legend_('In');
        if ~isempty(info.inunit)
           ylabel(sprintf('[%s]', info.inunit));
        end
        if isfield(S,'InNoise') && ~strcmp(opts.type,'windowed')
            hold on;
            plot(t,S.InNoise);
        end
        tflab_title(S,opts,'ts');
        [~, lo] = legend(lg,opts.legend{:});
        adjust_legend_lines(lo);
        adjust_ylim();
        adjust_exponent('y');
        setx(0,info,[t(1),t(end)]);
        
    ax(2) = subplot('Position',opts.PositionBottom);
        plot(t,Out);
        grid on;grid minor;box on;    
        lg = legend_('Out');
        if ~isempty(info.outunit)
           ylabel(sprintf('[%s]', info.outunit));
        end
        if isfield(S,'OutNoise') && ~strcmp(opts.type,'windowed')
            hold on;
            plot(t,S.OutNoise);
        end
        [~, lo] = legend(lg,opts.legend{:});
        adjust_legend_lines(lo);
        adjust_ylim();
        adjust_exponent('y');
        setx(1,info,[t(1),t(end)]);  
    
    if opts.print
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s_%s.%s',opts.printname, opts.type, opts.printfmt{i});
            figsave(fullfile(opts.printdir, fname), opts);
        end
    end
end

if ~iscell(S) && strcmp(opts.type,'error')
    
    for j = 1:size(S.Out,2)
        if j > 1
            figure();
            figprep();
        end
        outunit = '';
        if ~isempty(info.outunit)
            outunit = sprintf(' [%s]', info.outunit);
        end
        if iscell(info.outstr)
            outstrerr = info.outstr{j};
            outstr = sprintf('%s%s',info.outstr{j},outunit);
        else
            outstr = info.outstr;
            outstrerr = info.outstr;
            if size(S.Out,2) > 1
                outstr = sprintf('%s(:,%d)%s',outstr,j,outunit);
                outstrerr = sprintf('%s(:,%d)',outstr,j);
            end
        end
        metrics = sprintf('PE/CC/MSE = %.3f/%.3f/%.3f',...
                        S.Metrics.PE(j),...
                        S.Metrics.CC(j),...
                        S.Metrics.MSE(j));
        desc = S.Options.description;
        lg1{1} = 'Measured';
        lg1{2} = 'Predicted';
        lg2 = sprintf('%s Error; %s',outstrerr,metrics);

        ax(1,j) = subplot('Position',opts.PositionTop);
            plot(t,S.Out(:,j));
            grid on;grid minor;box on;hold on;
            plot(t,S.Metrics.Predicted(:,j));
            [~, lo] = legend(lg1(:),opts.legend{:});
            ylabel(outstr);
            adjust_legend_lines(lo);
            adjust_ylim();
            adjust_exponent('y');
            setx(0,info,[t(1),t(end)]);
            tflab_title(S,opts,'ts');
        ax(2,j) = subplot('Position',opts.PositionBottom);
            plot(t,S.Metrics.Predicted(:,j)-S.Out(:,j));
            grid on;grid minor;box on;
            ylabel(outunit);            
            %adjust_ylim();
            adjust_exponent('y')            
            [~, lo] = legend(lg2,opts.legend{:});
            adjust_legend_lines(lo);
            setx(1,info,[t(1),t(end)]);

        if opts.print
            for i = 1:length(opts.printfmt)
                fname = sprintf('%s_%s.%s',opts.printname,opts.type, opts.printfmt{i});
                figsave(fullfile(opts.printdir, fname));
            end
        end
    end

end

c = {'k','r','g','b'}; % TODO: Define more colors

% Compare
if iscell(S)
    if ~strcmp(opts.type,'error')
        error('tsplot for cell array input must be ''error''');
    end

    subplot('Position',opts.PositionTop);
        plot(t,S{1}.Out(:,1),c{1});
        grid on;grid minor;box on;hold on;
        ylabel(sprintf('%s [%s]',...
                    S{1}.Options.info.outstr{1},S{1}.Options.info.outunit));
        lg0 = sprintf('Observed at %s\n', S{1}.Options.info.stationid);

        for j = 1:length(S)
            plot(t,S{j}.Metrics.Predicted(:,1),c{j+1});
            lg{j} = sprintf('Predicted %s\n',...
                        S{j}.Options.description);
        end
        if ~isempty(opts.title)
            title(opts.title);
        end
        [~, lo] = legend({lg0,lg{:}},...
                        'Location','NorthEast','Orientation','Vertical');
        adjust_legend_lines(lo);
        adjust_ylim();
        adjust_exponent('y');
        setx(0,info,[t(1),t(end)]);
    subplot('Position',opts.PositionBottom);
        plot(t,S{1}.Out(:,2),c{1});
        grid on;grid minor;box on;hold on;
        ylabel(sprintf('%s [%s]',...
                    S{1}.Options.info.outstr{2},S{1}.Options.info.outunit));
        lg0 = sprintf('Observed at %s\n', S{1}.Options.info.stationid);

        for j = 1:length(S)
            for j = 1:length(S)
                plot(t,S{j}.Metrics.Predicted(:,2),c{j+1});
                lg{j} = sprintf('Predicted %s\n',...
                            S{j}.Options.description);
            end
        end
        [lh, lo] = legend({lg0,lg{:}},...
                        'Location','NorthEast','Orientation','Vertical');
        adjust_legend_lines(lo);
        adjust_ylim();
        adjust_exponent('y');
        setx(1,info,[t(1),t(end)]);
    
    if opts.print
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s_compare_%s.%s',...
                        opts.printname, opts.type, opts.printfmt{i});
            figsave(fullfile(opts.printdir, fname));
        end
    end
end


function ls = legend_(comp)

    if strcmp(comp,'In')
        compstrs = info.instr;
    else
        compstrs = info.outstr;
    end

    for j = 1:size(S.(comp),2)
        if ~iscell(compstrs) && size(S.(comp),2) > 1
            ls{j} = sprintf('%s(:,%d)\n',compstrs,j);
        else
            if iscell(compstrs)
                ls{j} = sprintf('%s\n', compstrs{j});
            else
                ls{j} = sprintf('%s\n', compstrs);
            end
        end
    end
    
    return;
    if isfield(S,[comp,'Noise'])
        jl = j;
        for j = 1:size(S.([comp,'Noise']),2)
            if iscell(compstrs)        
                ls{j+jl} = sprintf('%s(:,%d) Noise%s\n',compstrs{j},j,unitstr);    
            else
                ls{j+jl} = sprintf('%s Noise%s\n',compstrs,unitstr);
            end
        end
    end
end


function setx(last,info,tl)
    if ~isempty(info.timestart)
        % Set tick positions
        set(gca,'XLim',tl);
        adjust_xlim();
        datetick('x','keeplimits');
        xt = get(gca,'XTick');
        io = 1;
        if xt(1) < tl(1) % First label is not visible.
            io = 2; 
        end
        xtl = cellstr(get(gca,'XTickLabel'));
        xtl{io} = sprintf('$$\\begin{array}{c}\\mbox{%s} \\\\ %s\\end{array}$$',...
                            xtl{io},datestr(tl(1),'yyyy'));
        set(gca,'XTickLabel',xtl, 'TickLabelInterpreter', 'latex');
    end
    if ~last
        % Hide tick labels
        set(gca,'XTickLabel',[]);
    end
    if last && isempty(info.timestart)
        % Use default label
        if isempty(info.timeunit)
            xlabel('t');
        else
            xlabel(sprintf('%s since start', info.timeunit));
        end
        adjust_exponent('x');
    end
        
end

end % function

