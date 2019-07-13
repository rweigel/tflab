function fh = timeseries_plot(S,pt)
%TIMESERIES_PLOT

opts = S.Options;
t = S.Time;

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

fh = figure();
figprep();

if strcmp(pt,'raw')
    if size(S.In,2) > 1,
        s1='s';
    else
        s1='';
    end
    ts = sprintf('Raw Input%s (top) and Output (bottom)',s1);
    
    subplot('Position',PositionTop);
        plot(t,S.In);
        grid on;box on;
        for j = 1:size(S.In,2)
            if iscell(opts.info.instr)
                instr = opts.info.instr{j};
                ls{j} = sprintf('%s [%s]\n',...
                            instr,...
                            opts.info.inunit);
            else
                instr = opts.info.instr;
                ls{j} = sprintf('%s(:,%d) [%s]\n',...
                            instr,...
                            1:size(S.In,2),...
                            opts.info.inunit);
            end
        end
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        adjust_ylim();
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')
    subplot('Position',PositionBottom);
        plot(t,S.Out);
        grid on;box on;    
        if iscell(opts.info.outstr)
            ls = sprintf('%s [%s]\n',...
                    opts.info.outstr{1},...
                    opts.info.outunit);
        else
            ls = sprintf('%s(:,1) [%s]\n',...
                    opts.info.outstr,...
                    opts.info.outunit);
        end
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        adjust_ylim();
        if isfield(opts,'td') && ischar(opts.td.start)
            datetick('x');
        else
            xlabel(sprintf('%s [%s]',...
                            opts.info.timestr,...
                            opts.info.timeunit));        
        end
elseif strcmp(pt,'windowed')
    if size(S.In,2) > 1,s1='s';else,s1='';end  
    if strcmp(pt,'windowed')
        ts = sprintf('%s-windowed Input%s (top) and Output (bottom)',...
                     opts.td.window.functionstr,...
                     s1);
    end
    
    subplot('Position',PositionTop);
        plot(t,S.Window.In);
        grid on;box on;
        for j = 1:size(S.Window.In,2)
            if iscell(opts.info.instr)
                instr = opts.info.instr{j};
                ls{j} = sprintf('%s [%s]\n',...
                                instr,...
                                opts.info.inunit);
            else
                instr = opts.info.instr;
                ls{j} = sprintf('%s(:,%d) [%s]\n',...
                                instr,1:size(S.In,2),...
                                opts.info.inunit);
            end
        end
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        adjust_ylim();
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')
    subplot('Position',PositionBottom);
        plot(t,S.Window.Out);
        grid on;box on;    
        if iscell(opts.info.outstr)
            ls = sprintf('%s [%s]\n',...
                        opts.outstr{1},...
                        opts.info.outunit);
        else
            ls = sprintf('%s(:,1) [%s]\n',...
                    opts.info.outstr,...
                    opts.info.outunit);
        end
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        adjust_ylim();        
        xlabel(sprintf('%s [%s]',...
                        opts.info.timestr,...
                        opts.info.timeunit));
        if ischar(opts.td.start)
            datetick('x');
        end
elseif strcmp(pt,'error')
    ts = sprintf('Method: %s',opts.description);
    subplot('Position',PositionTop);
        plot(t,S.Out);
        grid on;box on;    
        if iscell(opts.info.outstr)        
            ls = sprintf('%s [%s]\n',...
                        opts.info.outstr{1},...
                        opts.info.outunit);
        else
            ls = sprintf('%s(:,1) [%s]\n',...
                        opts.info.outstr,...
                        opts.info.outunit);    
        end
        title(['Raw output (top) predicted raw output (bottom)',ts],...
                'FontWeight','Normal');
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        adjust_ylim();
        set(gca,'XTickLabel',[]);
    subplot('Position',PositionBottom);
        plot(t,S.Out-S.Predicted);
        grid on;box on;
        if iscell(opts.info.outstr)        
            ls = sprintf('%s - %s [%s]',...
                        opts.info.outstr{1},...
                        opts.info.outstr{1},...
                        opts.info.outunit);
        else
            ls = sprintf('%s(:,1) meas. - %s(:,1) pred. [%s]',...
                        opts.info.outstr,...
                        opts.info.outstr,...
                        opts.info.outunit);    
        end
        ls = sprintf('%s PE/CC/MSE = %.3f/%.3f/%.3f',...
                        ls,...
                        S.Metrics.PE,...
                        S.Metrics.CC,...
                        S.Metrics.MSE);
                 
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        adjust_ylim();
        xlabel(sprintf('%s [%s]',...
                    opts.info.timestr,...
                    opts.info.timeunit));
        if ischar(opts.td.start)
            datetick('x');
        end
end

try
    if opts.transferfnFD.plot.timeseries(2)
        % Print png
    end
    if opts.transferfnFD.plot.timeseries(3)
        % Print pdf
    end
catch
end
