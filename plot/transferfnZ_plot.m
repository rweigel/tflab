function transferfnZ_plot(S,popts)

% Default options
opts = struct();
    opts.title = '';
    opts.period = 1;
    opts.unwrap = 1;
    opts.magnitude = 0;
    opts.period_range = [];

% Use default options if options not given    
if nargin > 1
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

if iscell(S) && length(S) == 1
    S = S{1};
end

if (iscell(S) && size(S{1}.Z,2) > 1) || size(S.Z,2) > 1
    Zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};
    Phistrs = {'\phi_{xx}','\phi_{xy}','\phi_{yx}','\phi_{yy}'};
else
    Zstrs = {'Z'};
    Phistrs = {'Z'};
end

if isstruct(S)
    subplot('Position', PositionTop);
        if opts.period
            x = 1./S.fe;
        else
            x = S.fe;
        end
        if opts.magnitude
            y = abs(S.Z);
            loglog(x, y,'marker','.','markersize',10,'linewidth',2);
            grid on;
            for j = 1:size(S.Z,2)
                ls{j} = sprintf('$|%s|$',Zstrs{j});
            end
        else
            y = real(S.Z);
            semilogx(x,y,'marker','.','markersize',10,'linewidth',2);
            grid on;
            for j = 1:size(S.Z,2)
                ls{j} = sprintf('$\\Re(%s)$',Zstrs{j});
            end
        end
        unitstr = sprintf('[%s/(%s)]',...
                            S.Options.info.inunit,...
                            S.Options.info.outunit);
        ylabel(unitstr);
        title(S.Options.description,'FontWeight','Normal');
        if opts.period
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        adjust_exponent();
    subplot('Position', PositionBottom);
        if opts.magnitude
            if opts.unwrap
                y = (180/pi)*unwrap(atan2(imag(S.Z),real(S.Z)));
            else
                y = (180/pi)*(atan2(imag(S.Z),real(S.Z)));    
            end
            semilogx(x, y,'linewidth',2,'marker','.','markersize',10);
            grid on;
            ylabel('[$^\circ$]');
            for j = 1:size(S.Z,2)
                ls{j} = sprintf('$%s$',Phistrs{j});
            end
        else
            y = imag(S.Z);
            semilogx(x,y,'linewidth',2,'marker','.','markersize',10);
            grid on;
            ylabel(unitstr);            
            for j = 1:size(S.Z,2)
                ls{j} = sprintf('$\\Im(%s)$',Zstrs{j});
            end
        end

        if ~opts.unwrap
            set(gca,'YLim',[-185,185])
            set(gca,'YTick',[-180:60:180]);
        end
        if opts.period
            xlabel(sprintf('$T$ [%s]', S.Options.info.timeunit));
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
        else
            xlabel(sprintf('$f$ [1/%s]', S.Options.info.timeunit));
        end
        legend(ls,'Location','Best','Orientation','Horizontal');
        adjust_exponent();
        if isfield(popts,'filename') && ~isempty(popts.filename)
            [fpath,fname,fext] = fileparts(popts.filename);
            figsave(1,sprintf('%s/%s%s',fpath,fname,fext));
        end
        % If S.Segment, create 4-panel plot of Z with all segment Z's.
else
    %TODO:
    % assert(all(size(S1.Z) == size(S2.Z)), 'Required: size(S1.Z) == size(S2.Z)');
    % Assumes units are the same. Check this.
    for j = 1:size(S{1}.Z,2)
        figure();
        figprep();
        subplot('Position', PositionTop);
            for s = 1:length(S)
                if opts.period
                    x = 1./S{s}.fe;
                else
                    x = S{s}.fe;
                end
                if opts.magnitude
                    y = abs(S{s}.Z(:,j));
                    ls{s} = sprintf('$%s$ %s',Zstrs{j},S{s}.Options.description);
                    loglog(x, y);
                else
                    y = real(S{s}.Z(:,j));
                    ls{s} = sprintf('$\\Re(%s)$ %s',Zstrs{j},S{s}.Options.description);
                    semilogx(x, y);
                end
                if s == 1
                    grid on;box on;hold on;
                end
            end
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
            adjust_exponent('y');
            % Assumes all units are the same
            ylabel(sprintf('[%s/(%s)]',...
                        S{1}.Options.info.inunit,...
                        S{1}.Options.info.outunit));
            legend(ls,'Location','NorthEast');
            set(gca,'XTickLabel',[]);
            title(opts.title);
        subplot('Position', PositionBottom);
            for s = 1:length(S)
                if opts.period
                    x = 1./S{s}.fe;
                else
                    x = S{s}.fe;
                end
                if opts.magnitude
                    y = (180/pi)*unwrap(atan2(imag(S{s}.Z(:,j)),real(S{s}.Z(:,j))));
                    yl = '[degrees]';
                    ls{s} = sprintf('$%s$ $\\,$ %s',Phistrs{j},S{s}.Options.description);
                    semilogx(x, y);
                else
                    y = imag(S{s}.Z(:,j));
                    yl = sprintf('[%s/(%s)]',...
                                 S{1}.Options.info.inunit,...
                                 S{1}.Options.info.outunit);
                    ls{s} = sprintf('$\\Im(%s)$ %s',Zstrs{j},S{s}.Options.description);
                    semilogx(x, y);
                end
                if s == 1
                    grid on;box on;hold on;
                end
            end
            ylabel(yl);
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
                xlabel(sprintf('$T$ [%s]', S{1}.Options.info.timeunit));
                if ~isempty(opts.period_range)
                    set(gca,'XLim',opts.period_range);
                end
            else
                xlabel(sprintf('$f$ [1/%s]', S{1}.Options.info.timeunit));
            end
            legend(ls,'Location','NorthEast');
        if isfield(popts,'filename') && ~isempty(popts.filename)
            [fpath,fname,fext] = fileparts(popts.filename);
            figsave(1,sprintf('%s/%s-%s%s',fpath,fname,Zstrs{j},fext));
        end
    end
end