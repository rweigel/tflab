function transferfnZ_plot(S,popts)

% Default options
opts = struct();
    opts.title = '';
    opts.period = 1;
    opts.unwrap = 1;
    opts.plottype = 1; % 1 = Z,phi, 2 = rho,phi; 3 = Re,Im
    opts.period_range = [];
    opts.savefmt = {}; % Any extension allowed by export_fig

% Use default options if options not given    
if nargin > 1
    fns = fieldnames(popts);
    for i = 1:length(fns)
        opts.(fns{i}) = popts.(fns{i});
    end
end

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

if iscell(S) && length(S) == 1
    S = S{1};
end

if (iscell(S) && size(S{1}.Z,2) > 1) || size(S.Z,2) > 1
    Zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};
    Rhostrs = {'\rho^a_{xx}','\rho^a_{xy}','\rho^a_{yx}','\rho^a_{yy}'};
    Phistrs = {'\phi_{xx}','\phi_{xy}','\phi_{yx}','\phi_{yy}'};
else
    Zstrs = {'Z'};
    Phistrs = {'Z'};
end

if isstruct(S)
    figure();
    figprep();
    S = defaultmeta(S);
    if opts.period
        x = S.Options.info.timedelta./S.fe;
    else
        x = S.Options.info.timedelta;
    end
    subplot('Position', PositionTop);
        switch opts.plottype
            case 1
                y = abs(S.Z);
                for j = 1:size(S.Z,2)
                    ls{j} = sprintf('$%s$',Zstrs{j});
                end
                loglog(x, y, 'linewidth',2,'marker','.','markersize',10);
            case 2
                % See Egbert, Booker, and Schultz 1992 pg 15,116.
                % rho_ij = (mu_o/omega)*|E_i|^2/|B_j|^2
                % rho_ij = (mu_o/omega)*|Z_ij|^2
                % mu_o = 4*pi*e-7 N/A^2
                % rho_ij = |Z_ij|^2/(5*f)
                % when f in Hz, Z in (mV/km)/nT
                y = abs(S.Z).^2./(5*S.fe/S.Options.info.timedelta);
                for j = 1:size(S.Z,2)
                    ls{j} = sprintf('$%s$',Rhostrs{j});
                end
                loglog(x, y, 'linewidth',2,'marker','.','markersize',10);
            case 3
                y = real(Z);
                ls{s} = sprintf('Re$(%s)$ %s',Zstrs,S.Options.description);
                semilogx(x, y, 'linewidth',2,'marker','.','markersize',10);
        end
        grid on;
        unitstr = sprintf('[(%s)/%s]',S.Options.info.outunit,S.Options.info.inunit);
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
        if opts.plottype ~= 3
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
        rep = '\{|\}';
        pre = [opts.filename,'-',regexprep(Zstrs{j},rep,'')];
        switch opts.plottype
            case 1
                fname = [pre,'_Magnitude_Phase'];
            case 2
                fname = [pre,'_Rhoa_Phase'];            
            case 3
                fname = [pre,'_Real_Imaginary'];            
        end
        for i = 1:length(opts.savefmt)
            figsave([fname,'.',opts.savefmt{i}]);
        end
else
    %TODO:
    % assert(all(size(S1.Z) == size(S2.Z)), 'Required: size(S1.Z) == size(S2.Z)');
    % Assumes units are the same. Check this.
    for j = 1:size(S{1}.Z,2)
        figure();
        figprep();
        for s = 1:length(S)
            S{s} = defaultmeta(S{s});
        end        
        subplot('Position', PositionTop);
            for s = 1:length(S)
                if opts.period
                    x = S{s}.Options.info.timedelta./S{s}.fe;
                else
                    x = S{s}.fe/S{s}.Options.info.timedelta;
                end
                switch opts.plottype
                    case 1
                        y = abs(S{s}.Z(:,j));
                        ls{s} = sprintf('$%s$ %s %s',Zstrs{j},...
                            S{s}.Options.info.stationid, S{s}.Options.description);
                        loglog(x, y, 'linewidth',2,'marker','.','markersize',10);
                    case 2
                        % See Egbert, Booker, and Schultz 1992 pg 15,116.
                        % rho_ij = (mu_o/omega)*|E_i|^2/|B_j|^2
                        % rho_ij = (mu_o/omega)*|Z_ij|^2
                        % mu_o = 4*pi*e-7 N/A^2
                        % rho_ij = |Z_ij|^2/(5*f)
                        % when f in Hz, Z in (mV/km)/nT
                        y = (abs(S{s}.Z(:,j))).^2./(5*S{s}.fe/S{s}.Options.info.timedelta);
                        ls{s} = sprintf('$%s$ %s',Rhostrs{j},S{s}.Options.description);
                        loglog(x, y, 'linewidth',2,'marker','.','markersize',10);
                    case 3
                        y = real(S{s}.Z(:,j));
                        ls{s} = sprintf('Re$(%s)$ %s',Zstrs{j},S{s}.Options.description);
                        semilogx(x, y, 'linewidth',2,'marker','.','markersize',10);                        
                end
                if s == 1
                    grid on;box on;hold on;
                end
            end
            if opts.plottype == 2
                yl = '$\Omega\cdot$m';
            else
                yl = sprintf('[(%s)/%s]',...
                    S{1}.Options.info.outunit,...
                    S{1}.Options.info.inunit);                
            end
            ylabel(yl);
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
            adjust_exponent();
            legend(ls,'Location','NorthEast');
            set(gca,'XTickLabel',[]);
            title(opts.title,'FontWeight','Normal');
        subplot('Position', PositionBottom);
            for s = 1:length(S)
                if opts.period
                    x = S{s}.Options.info.timedelta./S{s}.fe;
                else
                    x = S{s}.fe/S{s}.Options.info.timedelta;
                end
                if opts.plottype == 3
                    y = imag(S{s}.Z(:,j));
                    yl = sprintf('[(%s)/%s]',...
                                 S{1}.Options.info.outunit,...
                                 S{1}.Options.info.inunit);
                    ls{s} = sprintf('Im$(%s)$ %s',Zstrs{j},S{s}.Options.description);
                    semilogx(x, y, 'linewidth',2,'marker','.','markersize',10);
                else
                    y = (180/pi)*unwrap(atan2(imag(S{s}.Z(:,j)),real(S{s}.Z(:,j))));
                    yl = '[degrees]';
                    ls{s} = sprintf('$%s$ %s %s',Phistrs{j},...
                        S{s}.Options.info.stationid,S{s}.Options.description);
                    semilogx(x, y, 'linewidth',2,'marker','.','markersize',10);
                end
                if s == 1
                    grid on;box on;hold on;
                end
            end
            ylabel(yl);
            if opts.plottype ~= 3 && ~opts.unwrap
                set(gca,'YLim',[-185,185])
                set(gca,'YTick',[-180:60:180]);
            end
            if opts.period
                xlabel(sprintf('$T$ [%s]', S{1}.Options.info.timeunit));
                if ~isempty(opts.period_range)
                    set(gca,'XLim',opts.period_range);
                end
            else
                if ~isempty(opts.period_range)
                    set(gca,'XLim',[1/opts.period_range(2),1/opts.period_range(1)]);
                end
                xlabel(sprintf('$f$ [1/%s]', S{1}.Options.info.timeunit));
            end
            legend(ls,'Location','NorthEast');
            adjust_exponent();
        rep = '\{|\}';
        pre = [opts.filename,'-',regexprep(Zstrs{j},rep,'')];
        switch opts.plottype
            case 1
                fname = [pre,'_Magnitude_Phase'];
            case 2
                fname = [pre,'_Rhoa_Phase'];            
            case 3
                fname = [pre,'_Real_Imaginary'];            
        end
        for i = 1:length(opts.savefmt)
            figsave([fname,'.',opts.savefmt{i}]);
        end
        end
    end
end


function S = setsubfield(S,varargin)

    if length(varargin) == 1
        error('At least three arguments needed.');
    end
    if length(varargin) == 2
        if ~isfield(S,varargin{1})
            S.(varargin{1}) = varargin{2};
        end
        return
    end
    if ~isfield(S,varargin{1})
        S.(varargin{1}) = struct(varargin{2},varargin{3});        
    end

    S.(varargin{1}) = setsubfield(S.(varargin{1}),varargin{2:end});
    
end

function S = defaultmeta(S)
    S = setsubfield(S,'Options','description','');
    S = setsubfield(S,'Options','info','timedelta',1);
    S = setsubfield(S,'Options','info','timedelta',1);
    S = setsubfield(S,'Options','info','outunit','?');
    S = setsubfield(S,'Options','info','inunit','?');
    S = setsubfield(S,'Options','info','timeunit','?');
end