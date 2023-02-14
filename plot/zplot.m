function [ax1,ax2] = zplot(S,popts)
%ZPLOT

% Default options
opts = struct();
    opts.title = '';
    opts.type = 1; % 1 = Z,phi, 2 = rho,phi; 3 = Re,Im
    opts.print = 0;
    opts.printname = 'zplot';
    opts.printdir = '';
    opts.printfmt = {'pdf'};
    if strcmp(opts.printname,'zplot')
        switch opts.type
            case 1
                opts.printname = [opts.printname,'_magnitude_phase'];
            case 2
                opts.printname = [opts.printname,'_rhoa_ohase'];            
            case 3
                opts.printname = [opts.printname,'_real_imaginary'];            
        end
    end
    opts.vs_period = 1;
    opts.unwrap = 0;
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
    Rhostrs = {'\rho^a'};
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
    interp_Z_exists = 0;
    if isfield(S,'Zi')
        interp_Z_exists = 1;
    end
    
    figprep();
    S = defaultmeta(S);
    if opts.vs_period
        x = S.Options.info.timedelta./S.fe;
        if interp_Z_exists
            xi = S.Options.info.timedelta./S.fi;
        end
    else
        x = S.Options.info.timedelta*S.fe;
        if interp_Z_exists
            xi = S.Options.info.timedelta*S.fi;
        end
    end
    ax1 = subplot('Position', PositionTop);
        switch opts.type
            case 1
                y = abs(S.Z);
                if interp_Z_exists
                    yi = abs(S.Zi);
                end
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
                if interp_Z_exists
                    yi = abs(S.Zi).^2./(5*S.fi/S.Options.info.timedelta);
                end
                for j = 1:size(S.Z,2)
                    ls{j} = sprintf('$%s$',Rhostrs{j});
                end
            case 3
                y = real(S.Z);
                for j = 1:size(S.Z,2)
                    ls{j} = sprintf('Re$(%s)$ Estimated',Zstrs{j});
                end
                if interp_Z_exists
                    yi = real(S.Zi);
                    jx = size(S.Z,2);
                    for j = 1:size(S.Z,2)
                        ls{jx+j} = sprintf('Re$(%s)$ Interpolated',Zstrs{j});
                    end
                end
        end
        if opts.vs_period && ~isempty(opts.period_range)
            idx = find(x <= opts.period_range(1) | x >= opts.period_range(2));
            y(idx) = NaN;
            if interp_Z_exists
                idx = find(xi <= opts.period_range(1) | xi >= opts.period_range(2));
                yi(idx) = NaN;
            end
        end
        plot(x, y, '.','markersize',20);
        hold on;grid on;
        if interp_Z_exists
            plot(xi, yi, '.','markersize',15);
        end
        if opts.type < 3
            set(gca,'YScale','log');
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        unitstr = '';
        if ~isempty(S.Options.info.outunit)
            unitstr = sprintf('[(%s)/%s]',...
                        S.Options.info.outunit,S.Options.info.inunit);
        end
        ylabel(unitstr);
        
        if opts.vs_period
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
        end
        legend(ls,'Location','NorthEast','Orientation','Horizontal');
        set(gca,'XTickLabel',[]);
        adjust_ylim('both');
        adjust_yticks(1e-4);
        adjust_exponent('y');
    ax2 = subplot('Position', PositionBottom);
        if opts.type ~= 3
            if opts.unwrap
                y = (180/pi)*unwrap(atan2(imag(S.Z),real(S.Z)));
                if interp_Z_exists
                    yi = (180/pi)*unwrap(atan2(imag(S.Zi),real(S.Zi)));
                end
            else
                y = (180/pi)*(atan2(imag(S.Z),real(S.Z)));
                if interp_Z_exists
                    yi = (180/pi)*(atan2(imag(S.Zi),real(S.Zi)));
                end
            end
            for j = 1:size(S.Z,2)
                ls{j} = sprintf('Re$(%s)$',Phistrs{j});
            end
            ylabel('[$^\circ$]');
        else
            y = imag(S.Z);
            ylabel(unitstr);
            for j = 1:size(S.Z,2)
                ls{j} = sprintf('Im$(%s)$ Estimated',Zstrs{j});
            end
            if interp_Z_exists
                yi = imag(S.Zi);
                jx = size(S.Z,2);
                for j = 1:size(S.Z,2)
                    ls{jx+j} = sprintf('Im$(%s)$ Interpolated',Zstrs{j});
                end
            end
        end
        if opts.vs_period && ~isempty(opts.period_range)
            idx = find(x <= opts.period_range(1) | x >= opts.period_range(2));
            y(idx) = NaN;
            if interp_Z_exists
                idx = find(xi <= opts.period_range(1) | xi >= opts.period_range(2));
                yi(idx) = NaN;
            end
        end
        plot(x, y, '.','markersize',20);
        hold on;grid on;
        if interp_Z_exists
            plot(xi, yi, '.','markersize',15);
        end
        if ~opts.unwrap && opts.type ~= 3
            set(gca,'YLim',[-185,185])
            set(gca,'YTick',[-180:60:180]);
        end
        unitstr = '';
        if opts.vs_period
            set(gca,'XScale','log');
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
        adjust_exponent();

        if opts.print
            for i = 1:length(opts.printfmt)
                fname = sprintf('%s.%s',opts.printname, opts.printfmt{i});
                figsave(fullfile(opts.printdir, fname), opts);
            end
        end        
else
    %TODO:
    % assert(all(size(S1.Z) == size(S2.Z)), 'Required: size(S1.Z) == size(S2.Z)');
    % Assumes units are the same. Check this.
    for j = 1:size(S{1}.Z,2)
        figprep();
        if j > 1
            figure();
        end
        for s = 1:length(S)
            S{s} = defaultmeta(S{s});
        end                
        ax1(j) = subplot('Position', PositionTop);
            for s = 1:length(S)
                ls{s} = sprintf('%s',S{s}.Options.description);
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
                yebl = [];
                yebu = [];
                switch opts.type
                    case 1
                        y = abs(S{s}.Z(:,j));
                        if isfield(S{s},'ZCL')
                            %yebl =  y-squeeze(S{s}.ZCL.Magnitude.Normal.x_1sigma(:,j,1));
                            %yebu = -y+squeeze(S{s}.ZCL.Magnitude.Normal.x_1sigma(:,j,2));
                            yebl =  y-squeeze(S{s}.ZCL.Magnitude.Bootstrap.x_1sigma(:,j,1));
                            yebu = -y+squeeze(S{s}.ZCL.Magnitude.Bootstrap.x_1sigma(:,j,2));
                        end
                    case 2
                        % See Egbert, Booker, and Schultz 1992 pg 15,116.
                        % rho_ij = (mu_o/omega)*|E_i|^2/|B_j|^2
                        % rho_ij = (mu_o/omega)*|Z_ij|^2
                        % mu_o = 4*pi*e-7 N/A^2
                        % rho_ij = |Z_ij|^2/(5*f)
                        % when f in Hz, Z in (mV/km)/nT
                        y = (abs(S{s}.Z(:,j))).^2./(5*S{s}.fe/S{s}.Options.info.timedelta);
                    case 3
                        y = real(S{s}.Z(:,j));
                        if isfield(S{s},'ZCL')
                            yebl =  y-squeeze(real(S{s}.ZCL.Z.Normal.x_1sigma(:,j,1)));
                            yebu = -y+squeeze(real(S{s}.ZCL.Z.Normal.x_1sigma(:,j,2)));
                        end
                end
                if opts.vs_period && ~isempty(opts.period_range)
                    idx = find(x <= opts.period_range(1) | x >= opts.period_range(2));
                    y(idx) = NaN;
                end                
                if opts.vs_period
                    h(s) = semilogx(x, y, lineopts{:});
                else
                    h(s) = plot(x, y, lineopts{:});
                end
                if s == 1
                    grid on;box on;hold on;
                end
                if ~isempty(yebl)
                    errorbars(x,y,yebl,yebu);
                    %return
                    %b = squeeze(S{s}.Segment.Z(:,j,:));
                    %for k = 2:length(x)
                    %    plot(x(k),abs(b(k,:)),'k.','MarkerSize',1);
                    %end
                end
            end

            yunitstr = '';
            if ~isempty(S{1}.Options.info.outunit) ...
                    && ~isempty(S{1}.Options.info.inunit)
                yunitstr = sprintf(' [(%s)/%s]',...
                                S{1}.Options.info.outunit,...
                                S{1}.Options.info.inunit);
            end
            if opts.type == 1
                yl = sprintf('$|%s|$%s',Zstrs{j},yunitstr);
            end
            if opts.type == 2
                yl = sprintf('$%s% [$\Omega\cdot$m]',Rhostrs{j});
            end
            if opts.type == 3
                yl = sprintf('Re$(%s)$ %s',Zstrs{j},yunitstr);
            end
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
            if opts.type < 3
                %set(gca,'YScale','log');
            end
            ylabel(yl);
            legend(h,ls,'Location','NorthEast','Orientation','Horizontal');
            set(gca,'XTickLabel',[]);
            adjust_ylim('both');
            adjust_yticks(1e-4);
            if ~isempty(S{1}.Options.info.timeunit) && opts.vs_period
                period_lines();
            end
            adjust_exponent('y');
        ax2(j) = subplot('Position', PositionBottom);
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
                yebl = [];
                yebu = [];
                if opts.type == 3
                    y = imag(S{s}.Z(:,j));
                    yl = sprintf('Im$(%s)$',Zstrs{j});
                    if ~isempty(S{1}.Options.info.outunit) ...
                            && ~isempty(S{1}.Options.info.inunit)
                        yl = sprintf('%s [(%s)/%s]',yl,...
                                     S{1}.Options.info.outunit,...
                                     S{1}.Options.info.inunit);
                    end
                    ls{s} = sprintf('%s',S{s}.Options.description);
                    if isfield(S{s},'ZCL')
                        yebl = squeeze(imag(S{s}.ZCL.Z.Normal.x_1sigma(:,j,1)));
                        yebu = squeeze(imag(S{s}.ZCL.Z.Normal.x_1sigma(:,j,2)));
                    end
                else
                    if opts.unwrap
                        y = (180/pi)*unwrap(atan2(imag(S{s}.Z(:,j)),real(S{s}.Z(:,j))));
                    else
                        y = (180/pi)*atan2(imag(S{s}.Z(:,j)),real(S{s}.Z(:,j)));
                    end
                    yl = sprintf('$%s$ [$^\\circ]$',Phistrs{j});
                    ls{s} = sprintf('%s',S{s}.Options.description);
                end
                if opts.vs_period && ~isempty(opts.period_range)
                    idx = find(x <= opts.period_range(1) | x >= opts.period_range(2));
                    y(idx) = NaN;
                end
                if opts.vs_period
                    semilogx(x, y, lineopts{:});
                else
                    plot(x, y, lineopts{:});
                end
                if s == 1
                    grid on;box on;hold on;
                end
                if ~isempty(yebl)
                    errorbars(x,y,yebl,yebu);
                end                
            end
            ylabel(yl);
            if opts.type ~= 3 && ~opts.unwrap
                set(gca,'YLim',[-185,185]);
                set(gca,'YTick',[-180:60:180]);
            end
            if opts.vs_period
                unit = '';
                if S{1}.Options.info.timeunit
                    unit = sprintf(' [%s]',S{1}.Options.info.timeunit);
                end
                xlabel(sprintf('$T$%s', unit));
                if ~isempty(opts.period_range)
                    set(gca,'XLim',opts.period_range);
                end
            else
                if opts.vs_period && ~isempty(opts.period_range)
                    set(gca,'XLim',[1/opts.period_range(2),1/opts.period_range(1)]);
                end
                unit = '';
                if S{1}.Options.info.timeunit
                    unit = sprintf(' [1/%s]',S{1}.Options.info.timeunit);
                end
                xlabel(sprintf('$f$%s', unit));
            end
            legend(ls,'Location','NorthEast','Orientation','Horizontal');
            adjust_ylim('upper');
            adjust_yticks(1e-4);
            adjust_exponent();
            adjust_ylim('upper');
            if ~isempty(S{1}.Options.info.timeunit) && opts.vs_period
                period_lines();
            end
        if opts.print
            comp = regexprep(Zstrs{j},'\{|\}','');
            for i = 1:length(opts.printfmt)
                fname = sprintf('%s-%s.%s',opts.printname, comp, opts.printfmt{i});
                figsave(fullfile(opts.printdir, fname), opts);
            end
        end    
    end % j
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
    S = setsubfield(S,'Options','info','outunit','');
    S = setsubfield(S,'Options','info','inunit','');
    S = setsubfield(S,'Options','info','timeunit','');
end