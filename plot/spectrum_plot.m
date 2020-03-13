function fh = spectrum_plot(S,popts)
%SPECTRUM_PLOT

opts = struct();
    opts.type = 'raw';
    opts.title = '';
    opts.period = 1;
    
% Use default options if options not given
if nargin > 1
    fns = fieldnames(popts);
    for i = 1:length(fns)
        if isfield(opts,fns{i})
           opts.(fns{i}) = popts.(fns{i});
        end
    end
end

info = S.Options.info;
t = S.Time;

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

figprep();

if strcmp(opts.type,'raw')
    if size(S.In,2) > 1,s1='s';else,s1='';end  
    if isempty(opts.title)
        ts = sprintf('Smoothed PSD of Raw Input%s (top) and Output (bottom)',s1);
    else
        ts = opts.title;
    end
    subplot('Position',PositionTop);
        if opts.period
            loglog(1./S.fe,S.PSD.In);
        else
            loglog(S.fe,S.PSD.In);
        end
        grid on;box on;
        for j = 1:size(S.PSD.In,2)
            if iscell(info.instr)
                ls{j} = sprintf('%s [%s]\n',...
                        info.instr{j},info.inunit);
            else
                ls{j} = sprintf('%s(:,%d) [%s]\n',...
                            info.instr,j,info.inunit);
            end
        end
        legend(ls,'Location','NorthWest','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')
        adjust_exponent('y');
    subplot('Position',PositionBottom);
        if opts.period
            loglog(1./S.fe,S.PSD.Out);
        else
            loglog(S.fe,S.PSD.Out);
        end
        grid on;box on;
        for j = 1:size(S.PSD.In,2)        
            if iscell(info.outstr)
                ls{j} = sprintf('%s [%s]\n',...
                            info.outstr{j},info.outunit);
            else
                ls{j} = sprintf('%s(:,1) [%s]\n',...
                           info.outstr,info.outunit);
            end
        end
        legend(ls,'Location','NorthWest','Orientation','Horizontal');
        if opts.period
            xlabel(sprintf('Period [%s]',info.timeunit));
        else
            xlabel(sprintf('Frequency [1/%s]',info.timeunit));
        end
        adjust_exponent('y');
elseif strcmp(opts.type,'windowed')
    if size(S.In,2) > 1,s1='s';else,s1='';end  
    if size(S.In,2) > 1,s1='s';else,s1='';end  
    if isempty(opts.title)
        ts = sprintf('Smoothed PSD of %s-windowed Input%s (top) and Output (bottom)',...
                        S.Options.td.window.functionstr,s1);
    else
        ts = opts.title;
    end
    subplot('Position',PositionTop);
        if opts.period
            loglog(1./S.fe,S.Window.PSD.In);
        else
            loglog(S.fe,S.Window.PSD.In);
        end
        grid on;box on;
        for j = 1:size(S.In,2)
            if iscell(info.instr)
                instr = info.instr{j};
                ls{j} = sprintf('%s [%s]\n',...
                        instr,info.inunit);
            else
                instr = info.instr;
                ls{j} = sprintf('%s(:,%d) [%s]\n',...
                            instr,j,info.inunit);
            end
        end
        legend(ls,'Location','NorthWest','Orientation','Horizontal');
        adjust_exponent('y');
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')
    subplot('Position',PositionBottom);
        if opts.period
            loglog(1./S.fe,S.Window.PSD.Out);
        else
            loglog(S.fe,S.Window.PSD.Out);
        end
        grid on;box on;    
        for j = 1:size(S.Out,2)
            if iscell(info.outstr)
                ls{j} = sprintf('%s [%s]\n',...
                        info.outstr{j},info.outunit);
            else
                ls{j} = sprintf('%s(:,%d) [%s]\n',...
                            info.outstr,j,info.outunit);
            end
        end
        legend(ls,'Location','NorthWest','Orientation','Horizontal');
        if opts.period
            xlabel(sprintf('Period [%s]',info.timeunit));
        else
            xlabel(sprintf('Frequency [1/%s]',info.timeunit));
        end
elseif strcmp(opts.type,'error')
    ts = sprintf('PSD(error) (top) and PSD(output)./PSD(error) (bottom)');
    subplot('Position',PositionTop);
        loglog(S.fe,S.PSD.Error);
        grid on;box on;    
        if iscell(info.outstr)
            ls = sprintf('PSD(%s meas.) - PSD(%s pred.)',...
                        info.outstr{1},info.outstr{1});
        else
            ls = sprintf('PSD(%s meas.) - PSD(%s pred.)',...
                        info.outstr,info.outstr);
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        title(ts,'FontWeight','Normal')        
    subplot('Position',PositionBottom);
        loglog(S.fe,S.PSD.Out./S.PSD.Error);
        grid on;box on;
        if iscell(info.outstr)
            ls = sprintf('PSD(%s output)./PSD(%s predicted)',...
                        info.outstr{1},info.outstr{1});
        else
            ls = sprintf('PSD(%s output)./PSD(%s predicted)',...
                        info.outstr,info.outstr);
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        xlabel(sprintf('Frequency [1/%s]',info.timeunit));
end
