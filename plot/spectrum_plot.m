function fh = spectrum_plot(S,pt)
%SPECTRUM_PLOT

opts = S.Options;
t = S.Time;

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

fh = figure();

if strcmp(pt,'raw')
    if size(S.In,2) > 1,s1='s';else,s1='';end  
    ts = sprintf('Smoothed PSD of Raw Input%s (top) and Output (bottom)',s1);
    
    subplot('Position',PositionTop);
        loglog(S.fe,S.PSD.In);
        grid on;box on;
        for j = 1:size(S.PSD.In,2)
            if iscell(opts.info.instr)
                instr = opts.info.instr{j};
                ls{j} = sprintf('%s [%s]\n',...
                        instr,opts.info.inunit);
            else
                instr = opts.info.instr;
                ls{j} = sprintf('%s(:,%d) [%s]\n',...
                            instr,1:size(S.In,2),opts.info.inunit);
            end
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')
    subplot('Position',PositionBottom);
        loglog(S.fe,S.PSD.Out);
        grid on;box on;
        if iscell(opts.info.outstr)
            ls = sprintf('%s [%s]\n',...
                    opts.info.outstr{1},opts.info.outunit);
        else
            ls = sprintf('%s(:,1) [%s]\n',...
                    opts.info.outstr,opts.info.outunit);
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        xlabel(sprintf('Frequency [1/%s]',opts.info.timeunit));
elseif strcmp(pt,'windowed')
    if size(S.In,2) > 1,s1='s';else,s1='';end  
    ts = sprintf('Smoothed PSD of %s-windowed Input%s (top) and Output (bottom)',...
                    opts.td.window.functionstr,s1);
    subplot('Position',PositionTop);
        loglog(S.fe,S.Window.PSD.In);
        grid on;box on;
        for j = 1:size(S.In,2)
            if iscell(opts.info.instr)
                instr = opts.info.instr{j};
                ls{j} = sprintf('%s [%s]\n',...
                        instr,opts.info.inunit);
            else
                instr = opts.info.instr;
                ls{j} = sprintf('%s(:,%d) [%s]\n',...
                            instr,1:size(S.In,2),opts.info.inunit);
            end
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')
    subplot('Position',PositionBottom);
        loglog(S.fe,S.Window.PSD.Out);
        grid on;box on;    
        if iscell(opts.info.outstr)
            ls = sprintf('%s [%s]\n',...
                    opts.outstr{1},opts.info.outunit);
        else
            ls = sprintf('%s(:,1) [%s]\n',...
                    opts.info.outstr,opts.info.outunit);
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        xlabel(sprintf('Frequency [1/%s]',opts.info.timeunit));
elseif strcmp(pt,'error')
    ts = sprintf('PSD(error) (top) and PSD(output)./PSD(error) (bottom)');
    subplot('Position',PositionTop);
        loglog(S.fe,S.PSD.Error);
        grid on;box on;    
        if iscell(opts.info.outstr)
            ls = sprintf('PSD(%s meas.) - PSD(%s pred.)',...
                        opts.info.outstr{1},opts.info.outstr{1});
        else
            ls = sprintf('PSD(%s meas.) - PSD(%s pred.)',...
                        opts.info.outstr,opts.info.outstr);
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')        
    subplot('Position',PositionBottom);
        loglog(S.fe,S.PSD.Out./S.PSD.Error);
        grid on;box on;
        if iscell(opts.info.outstr)
            ls = sprintf('PSD(%s output)./PSD(%s predicted)',...
                        opts.info.outstr{1},opts.info.outstr{1});
        else
            ls = sprintf('PSD(%s output)./PSD(%s predicted)',...
                        opts.info.outstr,opts.info.outstr);
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        xlabel(sprintf('Frequency [1/%s]',opts.info.timeunit));
end
    
if opts.transferfnFD.plot.spectrum(2)
    % Print png
end
if opts.transferfnFD.plot.spectrum(3)
    % Print pdf
end
