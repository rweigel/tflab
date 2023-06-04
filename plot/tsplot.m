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

opts = tflabplot_options(S, popts, 'original', 'tsplot');

if iscell(S) && length(S) == 1
    S = S{1};
end

S = defaultinfo(S);
if iscell(S)
    info = S{1}.Options.info;
    timeunit  = S{1}.Options.info.timeunit;
    timedelta = S{1}.Options.info.timedelta;

    t = timeVector_(S{1});
    timeunits = {};
    timestarts = {};
    for s = 1:length(S)
        timeunits{s}  = S{s}.Options.info.timeunit;
        timestarts{s} = S{s}.Options.info.timestart;
    end
    % TODO: Allow different timeunits and timestarts.
    if length(unique(timeunits)) > 1
        error('Time units must all be the same');
    end
    if length(unique(timestarts)) > 1
        error('Time starts must all be the same');
    end
else
    info = S.Options.info;
    timeunit  = S.Options.info.timeunit;
    timedelta = S.Options.info.timedelta;
    t = timeVector_(S);
end

figprep();

if ~iscell(S) && ~strcmp(opts.type,'error')

    if strcmp(opts.type,'original')
        In = S.In;
        Out = S.Out;
    elseif strcmp(opts.type,'final')
        In = S.In_.Final;
        Out = S.Out_.Final;
    else    
        typeuc = [upper(opts.type(1)),opts.type(2:end)];
        if ~isfield(S.In_,typeuc)
            error('Data were not %s. No plot created.',opts.type);
        end
        In = S.In_.(typeuc);
        Out = S.Out_.(typeuc);
    end
    
    [yl1, yl2] = ylabel_(S);
    [lg1, lg2] = legend_(S);
    
    ax(1) = subplot('Position',opts.PositionTop);
        plot(t,In);
        colororder_(ax(1), S.In);
        grid on;grid minor;box on;
        titlestr(S,opts,'ts');
        ylabel(yl1);
        if ~isempty(lg1)
            [~, lo] = legend(lg1,opts.legend{:});
            adjust_legend_lines(lo);
        end
        if isfield(S,'InNoise') && ~strcmp(opts.type,'windowed')
            %hold on;
            %plot(t,S.InNoise);
        end
        %adjust_ylim();
        adjust_exponent('y');
        setx_(0,info,[t(1),t(end)]);
        
    ax(2) = subplot('Position',opts.PositionBottom);
        plot(t,Out);
        colororder_(ax(2), S.Out);
        grid on;grid minor;box on;    
        ylabel(yl2);
        if ~isempty(lg2)
            [~, lo] = legend(lg2,opts.legend{:});
            adjust_legend_lines(lo);
        end
        if isfield(S,'OutNoise') && ~strcmp(opts.type,'windowed')
            %hold on;
            %plot(t,S.OutNoise);
        end
        %adjust_ylim();
        adjust_exponent('y');
        setx_(1,info,[t(1),t(end)]);  
    
    if opts.print
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s_%s.%s',...
                        opts.printname, opts.type, opts.printfmt{i});
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
            plot(t,S.Out_.Predicted(:,j));
            [~, lo] = legend(lg1(:),opts.legend{:});
            ylabel(outstr);
            adjust_legend_lines(lo);
            adjust_ylim();
            adjust_exponent('y');
            setx_(0,info,[t(1),t(end)]);
            titlestr(S,opts,'ts');
        ax(2,j) = subplot('Position',opts.PositionBottom);
            plot(t,S.Out_.Predicted(:,j)-S.Out(:,j));
            grid on;grid minor;box on;
            ylabel(outunit);            
            %adjust_ylim();
            adjust_exponent('y')            
            [~, lo] = legend(lg2,opts.legend{:});
            adjust_legend_lines(lo);
            setx_(1,info,[t(1),t(end)]);

        if opts.print
            for i = 1:length(opts.printfmt)
                fname = sprintf('%s_%s.%s',...
                            opts.printname,opts.type, opts.printfmt{i});
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
            plot(t,S{j}.Out_.Predicted(:,1),c{j+1});
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
        setx_(0,info,[t(1),t(end)]);
    subplot('Position',opts.PositionBottom);
        plot(t,S{1}.Out(:,2),c{1});
        grid on;grid minor;box on;hold on;
        ylabel(sprintf('%s [%s]',...
                    S{1}.Options.info.outstr{2},S{1}.Options.info.outunit));
        lg0 = sprintf('Observed at %s\n', S{1}.Options.info.stationid);

        for j = 1:length(S)
            for j = 1:length(S)
                plot(t,S{j}.Out_.Predicted(:,2),c{j+1});
                lg{j} = sprintf('Predicted %s\n',...
                            S{j}.Options.description);
            end
        end
        [lh, lo] = legend({lg0,lg{:}},...
                        'Location','NorthEast','Orientation','Vertical');
        adjust_legend_lines(lo);
        adjust_ylim();
        adjust_exponent('y');
        setx_(1,info,[t(1),t(end)]);
    
    if opts.print
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s_compare_%s.%s',...
                        opts.printname, opts.type, opts.printfmt{i});
            figsave(fullfile(opts.printdir, fname));
        end
    end
end

end % function

function setx_(last,info,tl)

    if ~isempty(info.timestart)

        % Set tick positions
        set(gca,'XLim',tl);
        adjust_xlim();

        if ischar(info.timestart)
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
    end
    
    if last == 0
        % Hide tick labels
        set(gca,'XTickLabel',[]);
    end

    if last == 1 && ~ischar(info.timestart)
        % Use default label
        if isempty(info.timeunit)
            xlabel('t');
        else
            xlabel(sprintf('%s since start', info.timeunit));
        end
        adjust_exponent('x');
    end
        
end

function [lg1, lg2] = legend_(S)

    info = S.Options.info;
    lg1 = '';
    lg2 = '';
    if size(S.In) > 1
        for j = 1:length(info.instr)
            lg1{j} = info.instr{j};
        end
    end
    if size(S.Out) > 1
        for j = 1:length(info.outstr)
            lg2{j} = info.outstr{j};
        end
    end    
end

function [yl1, yl2] = ylabel_(S)    
    info = S.Options.info;
    if size(S.In,2) > 1
        yl1 = unitstr(info.inunit);
    else
        yl1 = [info.instr,' ',unitstr(info.inunit)];
    end
    if size(S.Out,2) > 1
        yl2 = unitstr(info.outunit);
    else
        yl2 = [info.outstr,' ',unitstr(info.outunit)];
    end
end

function t = timeVector_(S)

    info = S.Options.info;
    nt = size(S.In,1);
    if ischar(info.timestart)
        try
            fmt = 'yyyy-mm-ddTHH:MM:SS.FFF';
            to = datenum(info.timestart,fmt);
        catch
            warning(['Could not parse Options.info.timestart. Format must be ',fmt]);
            info.timestart = '';
        end
        dt = info.timedelta;
        % Determine ppd (point per day)
        if strcmp(info.timeunit,'ms') || startsWith(info.timeunit,'millis')
            ppd = 86400000/dt;
        elseif  startsWith(info.timeunit,'s')
            ppd = 86400/dt;
        elseif strcmp(info.timeunit,'m') || startsWith(info.timeunit,'min')
            ppd = 1440/dt;
        else
            error(['Options.td.timeunit = %s not recognized. ',...
                   'Must be "millseconds", "seconds", or "minutes".'],...
                    info.timeunit);
        end
        t = to + (0:nt-1)'/ppd;
    else
        t = (info.timestart:info.timedelta:nt)';
    end
end

