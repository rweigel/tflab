function transferfnZ_plot(S1,S2,popts)

% Default options
opts = struct();
    opts.title = '';
    opts.period = 1;
    opts.unwrap = 1;
    opts.magnitude = 0;
    opts.period_range = [];

% Use default options if options not given    
if nargin > 2
    fns = fieldnames(popts);
    for i = 1:length(fns)
        if isfield(opts,fns{i})
           opts.(fns{i}) = popts.(fns{i});
        end
    end
end
    
figprep();

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

if size(S1.Z,2) > 1
    Zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};
    Phistrs = {'\phi_{xx}','\phi_{xy}','\phi_{yx}','\phi_{yy}'};
else
    Zstrs = {'Z'};
    Phistrs = {'Z'};
end

if nargin < 2
    subplot('Position', PositionTop);
        if opts.period
            loglog(1./S1.fe,abs(S1.Z),...
                    'marker','.','markersize',10,'linewidth',2);
        else
            loglog(S1.fe,abs(S1.Z),...
                    'marker','.','markersize',10,'linewidth',2);
        end
        grid on;box on;hold on;
        unitstr = sprintf('[%s/(%s)]',...
                            S1.Options.info.inunit,...
                            S1.Options.info.outunit);
        title(ts,'FontWeight','Normal');
        if size(S1.Z,2) == 1
            legend(sprintf('$|Z|$ %s', unitstr),'Location','NorthEast');
        else
            for j = 1:size(S1.Z,2)
                ls{j} = sprintf('$|%s|$ %s',Zstrs{j},unitstr);
            end
            legend(ls,'Location','NorthEast');
        end
        set(gca,'XTickLabel',[]);
        adjust_exponent('y');
    subplot('Position', PositionBottom);
        if opts.unwrap
            ang = (180/pi)*unwrap(atan2(imag(S1.Z),real(S1.Z)));
        else
            ang = (180/pi)*(atan2(imag(S1.Z),real(S1.Z)));    
        end
        if opts.period
            semilogx(1./S1.fe,ang,...
                        'linewidth',2,'marker','.','markersize',10);
        else
            semilogx(S1.fe, (180/pi)*(atan2(imag(S1.Z),real(S1.Z))),...
                        'linewidth',2,'marker','.','markersize',10);
        end
        grid on;box on;hold on;
        if ~opts.unwrap
            set(gca,'YLim',[-185,185])
            set(gca,'YTick',[-180:60:180]);
        end
        if opts.period
            xlabel(sprintf('$T$ [%s]', S1.Options.info.timeunit));
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
        else
            xlabel(sprintf('$f$ [1/%s]', S1.Options.info.timeunit));
        end
        if size(S1.Z,2) == 1
            legend(sprintf('$|\\phi|$ %s', unitstr),'Location','NorthEast');
        else
            for j = 1:size(S1.Z,2)
                %ls{j} = sprintf('$%s$ [$^\\circ$]',Phistrs{j});
                ls{j} = sprintf('$%s$',Phistrs{j});
            end
            legend(ls,'Location','NorthEast');
        end
        ylabel('[$^\circ$]');
        %legend('$\phi$ [degrees]','Location','NorthEast');
        adjust_exponent('x');
end

if nargin > 1
    %assert(all(size(S1.Z) == size(S2.Z)), 'Required: size(S1.Z) == size(S2.Z)');
    % Assumes units are the same. TODO: Allow different unit labels.
    for j = 1:size(S1.Z,2)
        figure();
        figprep();
        if opts.period
            x1 = 1./S1.fe;
            x2 = 1./S2.fe;
        else
            x1 = S1.fe;
            x1 = S2.fe;
        end
        subplot('Position', PositionTop);        
            if opts.magnitude
                y1 = abs(S1.Z(:,j));
                y2 = abs(S2.Z(:,j));
                ls = {sprintf('$|%s|$ %s',Zstrs{j},S1.Options.description),...
                      sprintf('$|%s|$ %s',Zstrs{j},S2.Options.description)};
                loglog(x1, y1);
                grid on;box on;hold on;
                loglog(x2, y2);
            else
                y1 = real(S1.Z(:,j));
                y2 = real(S2.Z(:,j));
                ls = {sprintf('$\\Re(%s)$ %s',Zstrs{j},S1.Options.description),...
                      sprintf('$\\Re(%s)$ %s',Zstrs{j},S2.Options.description)};
                semilogx(x1, y1);
                grid on;box on;hold on;
                semilogx(x2, y2);
            end
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
            adjust_exponent('y');
            ylabel(sprintf('[%s/(%s)]',...
                        S1.Options.info.inunit,...
                        S1.Options.info.outunit));
            legend(ls,'Location','NorthEast');
            set(gca,'XTickLabel',[]);
            title(opts.title);
        subplot('Position', PositionBottom);
            if opts.magnitude
                y1 = (180/pi)*unwrap(atan2(imag(S1.Z(:,j)),real(S1.Z(:,j))));
                y2 = (180/pi)*unwrap(atan2(imag(S2.Z(:,j)),real(S2.Z(:,j))));
                ylabel('[degrees]');
                ls = {sprintf('$|%s|$ $\\,$ %s',Phistrs{j},S1.Options.description),...
                      sprintf('$|%s|$ $\\,$ %s',Phistrs{j},S2.Options.description)},...
                semilogx(x1, y1);
                grid on;box on;hold on;
                semilogx(x2, y2);
            else
                y1 = imag(S1.Z(:,j));
                y2 = imag(S2.Z(:,j));
                ylabel(sprintf('[%s/(%s)]',...
                       S1.Options.info.inunit,...
                       S1.Options.info.outunit));                
                ls = {sprintf('$\\Im(%s)$ %s',Zstrs{j},S1.Options.description),...
                      sprintf('$\\Im(%s)$ %s',Zstrs{j},S2.Options.description)};
                semilogx(x1, y1);
                grid on;box on;hold on;
                semilogx(x2, y2);
                ylabel(sprintf('[%s/(%s)]',...
                        S1.Options.info.inunit,...
                        S1.Options.info.outunit));
            end
            if opts.magnitude && ~opts.unwrap
                set(gca,'YLim',[-185,185])
                set(gca,'YTick',[-180:60:180]);
            end
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
            if ~opts.magnitude
                adjust_exponent('y');
            end
            if opts.period
                xlabel(sprintf('$T$ [%s]', S1.Options.info.timeunit));
                if ~isempty(opts.period_range)
                    set(gca,'XLim',opts.period_range);
                end
            else
                xlabel(sprintf('$f$ [1/%s]', S1.Options.info.timeunit));
            end
            legend(ls,'Location','NorthEast');
            if isfield(popts,'filename') && ~isempty(popts.filename)
                figsave(1,sprintf('%s-%d.pdf',popts.filename,j));
            end
    end
end