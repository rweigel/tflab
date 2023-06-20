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

% Apply default metadata for fields not specified in S.Metadata.
S = tflab_metadata(S);

% Apply default plot options for fields in popts not specified.
popts = tflabplot_options(S, popts, 'original', 'tsplot');

if iscell(S) && length(S) == 1
    S = S{1};
end

if iscell(S)
    info = S{1}.Metadata;
    timeunit  = S{1}.Metadata.timeunit;
    timedelta = S{1}.Metadata.timedelta;

    t = timeVector_(S{1}, size(S{1}.In,1));
    timeunits = {};
    timestarts = {};
    for s = 1:length(S)
        timeunits{s}  = S{s}.Metadata.timeunit;
        timestarts{s} = S{s}.Metadata.timestart;
    end
    % TODO: Allow different timeunits and timestarts.
    if length(unique(timeunits)) > 1
        timeunits
        error('Time units must all be the same');
    end
    if length(unique(timestarts)) > 1
        error('Time starts must all be the same');
    end
else
    info = S.Metadata;
    timeunit  = S.Metadata.timeunit;
    timedelta = S.Metadata.timedelta;
    if strcmp(popts.type,'error')
        y1{1} = S.Out;
        y1{2} = S.Out_.Predicted;
        t1 = timeVector_(S, size(y1{1},1));

        y2 = S.Out_.Error;
        t2 = timeVector_(S, size(y2,1));

        lg1 = {'Measured', 'Predicted'};
    else
        if strcmp(popts.type,'original')
            y1 = S.In;
            y2 = S.Out;
        elseif strcmp(popts.type,'final')
            y1 = S.In_.Final;
            y2 = S.Out_.Final;
        else    
            typeuc = [upper(popts.type(1)),popts.type(2:end)];
            if ~isfield(S.In_,typeuc)
                error('Data were not %s. No plot created.',popts.type);
            end
            y1 = S.In_.(typeuc);
            y2 = S.Out_.(typeuc);
        end
        t1 = timeVector_(S, size(y1,1));
        t2 = timeVector_(S, size(y2,1));
        [yl1, yl2] = ylabel_(S);
        [lg1, lg2] = legend_(popts);
    end
end

figprep();

if ~iscell(S) && ~strcmp(popts.type,'error')
    
    ax(1) = subplot('Position',popts.PositionTop);
        plot(t1,y1);
        colororder_(ax(1), S.In);
        grid on;grid minor;box on;
        titlestr(S,popts,'ts');
        ylabel(yl1);
        if ~isempty(lg1)
            [~, lo] = legend(lg1,popts.legend{:});
            adjust_legend_lines(lo);
        end
        %adjust_ylim();
        adjust_exponent('y');
        setx_(0,info,[t1(1),t1(end)]);
        
    ax(2) = subplot('Position',popts.PositionBottom);
        plot(t1,y2);
        colororder_(ax(2), S.Out);
        grid on;grid minor;box on;    
        ylabel(yl2);
        if ~isempty(lg2)
            [~, lo] = legend(lg2,popts.legend{:});
            adjust_legend_lines(lo);
        end
        if isfield(S,'OutNoise') && ~strcmp(popts.type,'windowed')
            %hold on;
            %plot(t,S.OutNoise);
        end
        %adjust_ylim();
        adjust_exponent('y');
        setx_(1,info,[t2(1),t2(end)]);  
    
    if popts.print
        for i = 1:length(popts.printfmt)
            fname = sprintf('%s_%s.%s',...
                        popts.printname, popts.type, popts.printfmt{i});
            figsave(fullfile(popts.printdir, fname), popts);
        end
    end
end

if ~iscell(S) && strcmp(popts.type,'error')
    
    for j = 1:size(S.Out,2)
        if j > 1
            figure();
            figprep();
        end
        outunit = '';
        if ~isempty(info.outunit)
            outunit = sprintf(' [%s]', info.outunit);
        end
        if iscell(popts.outstr)
            outstrerr = popts.outstr{j};
            outstr = sprintf('%s%s',popts.outstr{j},outunit);
        else
            outstr = popts.outstr;
            outstrerr = popts.outstr;
            if size(S.Out,2) > 1
                outstr = sprintf('%s(:,%d)%s',outstr,j,outunit);
                outstrerr = sprintf('%s(:,%d)',outstr,j);
            end
        end
        metrics = sprintf('PE/CC/MSE = %.3f/%.3f/%.3f',...
                        S.Metrics.PE(j),...
                        S.Metrics.CC(j),...
                        S.Metrics.MSE(j));

        lg2 = sprintf('%s Error; %s',outstrerr,metrics);

        ax(1,j) = subplot('Position',popts.PositionTop);
            plot(t1,y1{1}(:,j));
            grid on;grid minor;box on;hold on;
            plot(t1,y1{2}(:,j));
            [~, lo] = legend(lg1(:),popts.legend{:});
            ylabel(outstr);
            adjust_legend_lines(lo);
            adjust_ylim();
            adjust_exponent('y');
            setx_(0,info,[t1(1),t1(end)]);
            titlestr(S,popts,'ts');
        ax(2,j) = subplot('Position',popts.PositionBottom);
            plot(t2,y2(:,j));
            grid on;grid minor;box on;
            ylabel(outunit);            
            %adjust_ylim();
            adjust_exponent('y')            
            [~, lo] = legend(lg2,popts.legend{:});
            adjust_legend_lines(lo);
            setx_(1,info,[t2(1),t2(end)]);

        if popts.print
            for i = 1:length(popts.printfmt)
                fname = sprintf('%s_%s.%s',...
                            popts.printname,popts.type, popts.printfmt{i});
                figsave(fullfile(popts.printdir, fname));
            end
        end
    end

end

c = {'k','r','g','b'}; % TODO: Define more colors

% Compare
if iscell(S)
    if ~strcmp(popts.type,'error')
        error('tsplot for cell array input must be ''error''');
    end

    subplot('Position',popts.PositionTop);
        plot(t,S{1}.Out(:,1),c{1});
        grid on;grid minor;box on;hold on;
        ylabel(sprintf('%s [%s]',...
                    popts.outstr{1},S{1}.Metadata.outunit));
        lg0 = sprintf('Observed at %s\n', S{1}.Metadata.stationid);

        for j = 1:length(S)
            plot(t,S{j}.Out_.Predicted(:,1),c{j+1});
            lg{j} = sprintf('Predicted %s\n',...
                        S{j}.Options.description);
        end
        if ~isempty(popts.title)
            title(popts.title);
        end
        [~, lo] = legend({lg0,lg{:}},...
                        'Location','NorthEast','Orientation','Vertical');
        adjust_legend_lines(lo);
        adjust_ylim();
        adjust_exponent('y');
        setx_(0,info,[t(1),t(end)]);
    subplot('Position',popts.PositionBottom);
        plot(t,S{1}.Out(:,2),c{1});
        grid on;grid minor;box on;hold on;
        ylabel(sprintf('%s [%s]',...
                    popts.outstr{2},S{1}.Metadata.outunit));
        lg0 = sprintf('Observed at %s\n', S{1}.Metadata.stationid);

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
    
    if popts.print
        for i = 1:length(popts.printfmt)
            fname = sprintf('%s_compare_%s.%s',...
                        popts.printname, popts.type, popts.printfmt{i});
            figsave(fullfile(popts.printdir, fname));
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

function [lg1, lg2] = legend_(opts)

    lg1 = '';
    lg2 = '';
    if length(opts.instr) > 1
        lg1 = opts.instr;
    end
    if length(opts.outstr) > 1
        lg2 = opts.outstr;
    end    
end

function [yl1, yl2] = ylabel_(S)    
    meta = S.Metadata;
    if size(S.In,2) > 1
        yl1 = unitstr(meta.inunit);
    else
        yl1 = [meta.instr,' ',unitstr(meta.inunit)];
    end
    if size(S.Out,2) > 1
        yl2 = unitstr(meta.outunit);
    else
        yl2 = [meta.outstr,' ',unitstr(meta.outunit)];
    end
end

function t = timeVector_(S, nt)

    meta = S.Metadata;
    if ischar(meta.timestart)
        try
            fmt = 'yyyy-mm-ddTHH:MM:SS.FFF';
            to = datenum(meta.timestart,fmt);
        catch
            warning(['Could not parse Metadata.timestart. Format must be ',fmt]);
            meta.timestart = '';
        end
        dt = meta.timedelta;
        % Determine ppd (points per day)
        if strcmp(meta.timeunit,'ms') || startsWith(meta.timeunit,'millis')
            ppd = 86400000/dt;
        elseif  startsWith(meta.timeunit,'s')
            ppd = 86400/dt;
        elseif strcmp(meta.timeunit,'m') || startsWith(meta.timeunit,'min')
            ppd = 1440/dt;
        else
            error(['Metadata.timeunit = %s not recognized. ',...
                   'Must be "millseconds", "seconds", or "minutes".'],...
                    meta.timeunit);
        end
        t = to + (0:nt-1)'/ppd;
    else
        t = (meta.timestart:meta.timedelta:nt)';
    end
end

