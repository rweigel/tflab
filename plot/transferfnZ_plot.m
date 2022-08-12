function transferfnZ_plot(S,popts)

% Default options
opts = struct();
    opts.title = '';
    opts.filename = 'transferfnZ';
    opts.vs_period = 0;
    opts.unwrap = 0;
    opts.plottype = 1; % 1 = Z,phi, 2 = rho,phi; 3 = Re,Im
    if isstruct(S)
        opts.period_range = [1, 2*size(S.In,1)];
        opts.frequency_range = [0, 0.5];
    else
        if opts.vs_period
            Nt = size(S{1}.In,1);
            for i = 2:length(S)
                Nt = max(Nt,size(S{i}.In,1));
            end
            opts.period_range = [1, 2*Nt];
        else
            opts.frequency_range = [0, 0.5];
        end
    end
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
% Line options
lnopts = {'marker','.','markersize',10,'linewidth',2};

figprep();

if iscell(S) && length(S) == 1
    S = S{1};
end

if (iscell(S) && size(S{1}.Z,2) > 1) || (isstruct(S) && size(S.Z,2) > 1)
    Zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};
    Rhostrs = {'\rho^a_{xx}','\rho^a_{xy}','\rho^a_{yx}','\rho^a_{yy}'};
    Phistrs = {'\phi_{xx}','\phi_{xy}','\phi_{yx}','\phi_{yy}'};
else
    Zstrs = {'Z'};
    Phistrs = {'\phi'};
end

% Single transfer function
if isstruct(S)
    ts = opts.title;
    if isempty(ts)
        sta = '';
        if isfield(S.Options.info,'stationid')
            sta = S.Options.info.stationid;
        end
        ts = sprintf('Site: %s',sta);
    end
    figprep();
    S = defaultmeta(S);
    if opts.vs_period
        x = S.Options.info.timedelta./S.fe;
        xi = S.Options.info.timedelta./S.fi;
    else
        x = S.Options.info.timedelta*S.fe;
        xi = S.Options.info.timedelta*S.fi;
    end
    subplot('Position', PositionTop);
        switch opts.plottype
            case 1
                y = abs(S.Z);
                yi = abs(S.Zi);
                for j = 1:size(S.Z,2)
                    ls{j} = sprintf('$%s$',Zstrs{j});
                end
            case 2
                % See Egbert, Booker, and Schultz 1992 pg 15,116.
                % rho_ij = (mu_o/omega)*|E_i|^2/|B_j|^2
                % rho_ij = (mu_o/omega)*|Z_ij|^2
                % mu_o = 4*pi*e-7 N/A^2
                % rho_ij = |Z_ij|^2/(5*f)
                % when f in Hz, Z in (mV/km)/nT
                y = abs(S.Z).^2./(5*S.fe/S.Options.info.timedelta);
                yi = abs(S.Zi).^2./(5*S.fi/S.Options.info.timedelta);
                for j = 1:size(S.Z,2)
                    ls{j} = sprintf('$%s$',Rhostrs{j});
                end
            case 3
                y = real(S.Z);
                yi = real(S.Zi);
                for j = 1:size(S.Z,2)
                    ls{j} = sprintf('Re$(%s)$ Interpolated',Zstrs{j});
                end
                jx = size(S.Z,2);
                for j = 1:size(S.Z,2)
                    ls{jx+j} = sprintf('Re$(%s)$',Zstrs{j});
                end
        end
        plot(xi, yi, '.','markersize',20);
        hold on;grid on;
        plot(x, y, '.','markersize',15);
        if opts.plottype < 3
            set(gca,'YScale','log');
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        unitstr = '';
        if ~isempty(S.Options.info.outunit)
            unitstr = sprintf('[(%s)/%s]',S.Options.info.outunit,S.Options.info.inunit);
            ylabel(unitstr);
        end
        if opts.vs_period
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
        end
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        adjust_ylim('both');
        adjust_yticks(1e-4);
        %period_lines(max_T);
        adjust_exponent('y');
    subplot('Position', PositionBottom);
        if opts.plottype ~= 3
            if opts.unwrap
                y = (180/pi)*unwrap(atan2(imag(S.Z),real(S.Z)));
                yi = (180/pi)*unwrap(atan2(imag(S.Z),real(S.Z)));
            else
                y = (180/pi)*(atan2(imag(S.Zi),real(S.Zi)));    
                yi = (180/pi)*unwrap(atan2(imag(S.Zi),real(S.Zi)));
            end
            for j = 1:size(S.Z,2)
                ls{j} = sprintf('Re$(%s)$',Phistrs{j});
            end
            ylabel('[$^\circ$]');
        else
            y = imag(S.Z);
            yi = imag(S.Zi);
            ylabel(unitstr);
            for j = 1:size(S.Z,2)
                ls{j} = sprintf('Im$(%s)$ Interpolated',Zstrs{j});
            end
            jx = size(S.Z,2);
            for j = 1:size(S.Z,2)
                ls{jx+j} = sprintf('Im$(%s)$',Zstrs{j});
            end
        end
        plot(xi, yi, '.','markersize',20);
        hold on;grid on;
        plot(x, y, '.','markersize',15);
        if opts.plottype < 3
            set(gca,'YScale','log');
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end

        if ~opts.unwrap && opts.plottype ~= 3
            set(gca,'YLim',[-185,185])
            set(gca,'YTick',[-180:60:180]);
        end
        unitstr = '';
        if opts.vs_period
            if ~isempty(S.Options.info.timeunit)
                unitstr = sprintf(' [%s]', S.Options.info.timeunit);
            end
            xlabel(sprintf('$T$%s', unitstr));
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
        else
            if ~isempty(S.Options.info.timeunit)
                unitstr = sprintf(' [1/%s]', S.Options.info.timeunit);
            end
            xlabel(sprintf('$f$%s',unitstr));
        end
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        %period_lines(max_T);
        adjust_exponent('y');

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
        if j > 1
            figure();
            figprep();
        end
        set(gcf,'DefaultLegendAutoUpdate','off')
        for s = 1:length(S)
            S{s} = defaultmeta(S{s});
        end                
        subplot('Position', PositionTop);
            max_T = 0;
            for s = 1:length(S)
                ms = max(30-8*s,1);
                if size(S{s}.Z,1) == 1
                    ms = 30;
                    lineopts = {'.','markersize',ms};
                else
                    lineopts = {'linewidth',max(4-s,1),'marker','.','markersize',ms};
                end
                if opts.vs_period
                    x = S{s}.Options.info.timedelta./S{s}.fe;
                    max_T = max(max_T,max(x(x < Inf)));
                else
                    x = S{s}.fe/S{s}.Options.info.timedelta;
                end
                yl = '';                
                switch opts.plottype
                    case 1
                        yl = sprintf('$%s$',Zstrs{j});
                        y = abs(S{s}.Z(:,j));
                        ls{s} = sprintf('%s %s',...
                            S{s}.Options.info.stationid, S{s}.Options.description);
                    case 2
                        % See Egbert, Booker, and Schultz 1992 pg 15,116.
                        % rho_ij = (mu_o/omega)*|E_i|^2/|B_j|^2
                        % rho_ij = (mu_o/omega)*|Z_ij|^2
                        % mu_o = 4*pi*e-7 N/A^2
                        % rho_ij = |Z_ij|^2/(5*f)
                        % when f in Hz, Z in (mV/km)/nT
                        y = (abs(S{s}.Z(:,j))).^2./(5*S{s}.fe/S{s}.Options.info.timedelta);
                        ls{s} = sprintf('$%s$ %s',Rhostrs{j},S{s}.Options.description);
                        yl = '$\Omega\cdot$m';
                    case 3
                        y = real(S{s}.Z(:,j));
                        if ~isempty(S{1}.Options.info.outunit)
                            yl = sprintf('[(%s)/%s]',...
                                         S{1}.Options.info.outunit,...
                                         S{1}.Options.info.inunit);
                        end
                        ls{s} = sprintf('%s',S{s}.Options.description);
                        yl = sprintf('%sRe$(%s)$ %s',yl,Zstrs{j});
                end
                if opts.vs_period
                    h(s) = semilogx(x, y, lineopts{:});
                else
                    h(s) = plot(x, y, lineopts{:});
                end
                if s == 1
                    grid on;box on;hold on;
                end
                %if isfield(S{s},'ZVAR') && ~isfield(S{s},'ZCL')
                    %errorbars(x,y,S{s}.ZVAR(:,j))
                %end
                if isfield(S{s},'ZCL')
                    yl = y-squeeze(abs(S{s}.ZCL.Bootstrap.x_95(:,1,j)));
                    yu = -y+squeeze(abs(S{s}.ZCL.Bootstrap.x_95(:,2,j)));
                    errorbars(x,y,yl,yu);
                end
            end
            if opts.vs_period && ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
            ylabel(yl);            
            legend(h,ls,'Location','NorthEast','Orientation','Horizontal');
            set(gca,'XTickLabel',[]);
            adjust_ylim('both');
            adjust_yticks(1e-4);
            %period_lines(max_T);
            adjust_exponent('y');
        subplot('Position', PositionBottom);
            for s = 1:length(S)
                ms = max(30-8*s,1);
                if size(S{s}.Z,1) == 1
                    ms = 30;
                    lineopts = {'.','markersize',ms};
                else
                    lineopts = {'linewidth',max(4-s,1),'marker','.','markersize',ms};
                end                
                if opts.vs_period
                    x = S{s}.Options.info.timedelta./S{s}.fe;
                else
                    x = S{s}.fe/S{s}.Options.info.timedelta;
                end
                yl = '';
                if opts.plottype == 3
                    y = imag(S{s}.Z(:,j));
                    yl = '';
                    if ~isempty(S{1}.Options.info.outunit)
                        yl = sprintf('[(%s)/%s]',...
                                     S{1}.Options.info.outunit,...
                                     S{1}.Options.info.inunit);
                    end
                    ls{s} = sprintf('%s',S{s}.Options.description);
                    yl = sprintf('%sIm$(%s)$ %s',yl,Zstrs{j});
                else
                    if opts.unwrap
                        y = (180/pi)*unwrap(atan2(imag(S{s}.Z(:,j)),real(S{s}.Z(:,j))));
                    else
                        y = (180/pi)*atan2(imag(S{s}.Z(:,j)),real(S{s}.Z(:,j)));
                    end
                    yl = sprintf('$%s$ [degrees]',Phistrs{j});
                    ls{s} = sprintf('%s %s',...
                        S{s}.Options.info.stationid,S{s}.Options.description);
                end
                if opts.vs_period
                    semilogx(x, y,lineopts{:});
                else
                    plot(x, y,lineopts{:});
                end

                if s == 1
                    grid on;box on;hold on;
                end
            end
            ylabel(yl);
            if opts.plottype ~= 3 && ~opts.unwrap
                set(gca,'YLim',[-185,185]);
                set(gca,'YTick',[-180:60:180]);
            end
            unitstr = '';
            if opts.vs_period
                if ~isempty(S{1}.Options.info.timeunit)
                    unitstr = sprintf(' [%s]', S{1}.Options.info.timeunit);
                end
                xlabel(sprintf('$T$%s',unitstr));
                if ~isempty(opts.period_range)
                    set(gca,'XLim',opts.period_range);
                end
            else
                if opts.vs_period && ~isempty(opts.period_range)
                    set(gca,'XLim',[1/opts.period_range(2),1/opts.period_range(1)]);
                end
                if ~isempty(S{1}.Options.info.timeunit)
                    unitstr = sprintf(' [1/%s]', S{1}.Options.info.timeunit);
                end
                xlabel(sprintf('$f$%s',unitstr));
            end
            legend(ls,'Location','NorthEast','Orientation','Horizontal');
            adjust_ylim('upper');
            adjust_yticks(1e-4);
            %period_lines(max_T);
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

end % function


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