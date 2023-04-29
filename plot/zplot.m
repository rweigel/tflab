function [ax1,ax2] = zplot(S,popts)
%ZPLOT

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

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
                opts.printname = [opts.printname,'_rhoa_phase'];            
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

% Line options
lnopts = {'marker','.','markersize',10,'linestyle','none'};

% Legend options
lgopts = {'Location','NorthWest','Orientation','Horizontal'};

if iscell(S) && length(S) == 1
    S = S{1};
end

if iscell(S)
    for s = 1:length(S)
        S{s} = defaultinfo(S{s});
    end
    % TODO: Check all same.
    timeunit = S{1}.Options.info.timeunit;
    timedelta = S{1}.Options.info.timedelta;
else
    S = defaultinfo(S);
    timeunit = S.Options.info.timeunit;
    timedelta = S.Options.info.timedelta;
end

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

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

    interp_Z_exists = 0;
    if isfield(S,'Zi')
        interp_Z_exists = 1;
    end
    
    if opts.vs_period
        x = timedelta./S.fe;
        if interp_Z_exists
            xi = timedelta./S.fi;
        end
    else
        x = timedelta*S.fe;
        if interp_Z_exists
            xi = timedelta*S.fi;
        end
    end
    
    figprep();
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
                dt = S.Options.info.timedelta;
                % TODO: Assumes Z in (mV/km)/nT and dt in seconds.
                y = z2rho(S.fe/dt, S.Z);
                if interp_Z_exists
                    yi = z2rho(S.fi/dt, S.Zi);
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
            idx = x <= opts.period_range(1) | x >= opts.period_range(2);
            y(idx) = NaN;
            if interp_Z_exists
                idx = xi <= opts.period_range(1) | xi >= opts.period_range(2);
                yi(idx) = NaN;
            end
        end
        
        plot(x, y, lnopts{:}, 'markersize', 20);
        hold on;grid on;box on;
        if interp_Z_exists
            plot(xi, yi, lnopts{:}, 'markersize', 15);
        end
        if length(ls) > 1
            ylabel(unitstr(S.Options.info));
            legend(ls,lgopts{:});
        else
            ylabel(sprintf('%s %s',ls{1},unitstr(S.Options.info)));
        end
        tflab_title(S,opts,'z');
        
        if opts.type < 3
            set(gca,'YScale','log');
            adjust_ylim('upper');
        else
            adjust_ylim('both');
        end
        adjust_yticks(1e-4);
        adjust_exponent('y');
        setx(opts,0,timeunit);

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
                ls{j} = sprintf('$%s$',Phistrs{j});
            end
            yl = '[$^\circ$]';
        else
            y = imag(S.Z);
            yl = unitstr(S.Options.info);
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
            idx = x <= opts.period_range(1) | x >= opts.period_range(2);
            y(idx) = NaN;
            if interp_Z_exists
                idx = xi <= opts.period_range(1) | xi >= opts.period_range(2);
                yi(idx) = NaN;
            end
        end
        hold on;grid on;box on;
        plot(x, y, lnopts{:}, 'markersize', 20);
        if interp_Z_exists
            plot(xi, yi, lnopts{:}, 'markersize', 15);
        end
        if length(ls) > 1
            ylabel(yl);
            legend(ls,lgopts{:});
        else
            ylabel(sprintf('%s %s',ls{1},yl));
        end
        if ~opts.unwrap && opts.type ~= 3
            set(gca,'YLim',[-185,185])
            set(gca,'YTick',-180:60:180);
        end        
        if opts.vs_period
            set(gca,'XScale','log');
        end
        
        legend(ls,lgopts{:});
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(opts,1,timeunit);

    if opts.print
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s.%s',opts.printname, opts.printfmt{i});
            figsave(fullfile(opts.printdir, fname), opts);
        end
    end
end % if isstruct(S)

% Multiple transfer functions
if iscell(S)
    
    % TODO: Assumes units are the same for all Zs.
    %       Check this before plotting.
    for j = 1:size(S{1}.Z,2)
        figprep();
        if j > 1
            figure();
            figprep();
        end
        ax1(j) = subplot('Position', PositionTop);
            for s = 1:length(S)
                if opts.vs_period
                    x = timedelta./S{s}.fe;
                else
                    x = S{s}.fe/timedelta;
                end
                
                ls{s} = sprintf('%s',S{s}.Options.description);
        
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
                        dt = S{s}.Options.info.timedelta;
                        % TODO: Assumes Z in (mV/km)/nT and dt in seconds.
                        y = z2rho(S{s}.fe/dt, S{s}.Z(:,j));
                    case 3
                        y = real(S{s}.Z(:,j));
                        if isfield(S{s},'ZCL')
                            yebl =  y-squeeze(real(S{s}.ZCL.Z.Normal.x_1sigma(:,j,1)));
                            yebu = -y+squeeze(real(S{s}.ZCL.Z.Normal.x_1sigma(:,j,2)));
                        end
                end
                if opts.vs_period && ~isempty(opts.period_range)
                    idx = x <= opts.period_range(1) | x >= opts.period_range(2);
                    y(idx) = NaN;
                end
                
                lnopts = lineopts(size(S{s}.Z,1), s);
                if opts.vs_period
                    h(s) = semilogx(x, y, lnopts{:});
                else
                    h(s) = plot(x, y, lnopts{:});
                end
                if s == 1
                    grid on;hold on;box on;
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

            yunitstr = unitstr(S{1}.Options.info);
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
            %tflab_title(S,opts,'z');        
            adjust_ylim('both');
            adjust_yticks(1e-4);
            if ~isempty(timeunit) && opts.vs_period
                period_lines();
            end
            adjust_exponent('y');
            setx(opts,0,timeunit);

        ax2(j) = subplot('Position', PositionBottom);
            for s = 1:length(S)
                yebl = [];
                yebu = [];
                if opts.type == 3
                    y = imag(S{s}.Z(:,j));
                    yl = sprintf('Im$(%s)$ %s',Zstrs{j},unitstr(S{1}.Options.info));
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
                
                if opts.vs_period
                    x = timedelta./S{s}.fe;
                else
                    x = S{s}.fe/timedelta;
                end
                if opts.vs_period && ~isempty(opts.period_range)
                    idx = x <= opts.period_range(1) | x >= opts.period_range(2);
                    y(idx) = NaN;
                end
                
                lnopts = lineopts(size(S{s}.Z,1), s);
                if opts.vs_period
                    semilogx(x, y, lnopts{:});
                else
                    plot(x, y, lnopts{:});
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
            legend(ls,'Location','NorthEast','Orientation','Horizontal');
            adjust_ylim('upper');
            adjust_yticks(1e-4);
            adjust_exponent();
            adjust_ylim('upper');
            setx(opts,1,timeunit);

        if opts.print
            comp = regexprep(Zstrs{j},'\{|\}','');
            for i = 1:length(opts.printfmt)
                fname = sprintf('%s-%s.%s',opts.printname, comp, opts.printfmt{i});
                figsave(fullfile(opts.printdir, fname), opts);
            end
        end    
    end % j
end % if iscell(S)

end % function

function lnopts = lineopts(Nz, s)
    ms = max(30-8*s,1);
    if Nz == 1
        ms = 30;
        lnopts = {'.','markersize', ms};
    else
        lnopts = {'linewidth',max(4-s,1),'marker','.','markersize',ms};
    end
end

function str = unitstr(info)

    str = '';
    if isempty(info.inunit) || isempty(info.outunit)
        return
    end
    inunit = sprintf('%s',info.inunit);
    if contains(inunit,'/')
        inunit = sprintf('(%s)',inunit);
    end
    outunit = sprintf('%s',info.outunit);
    if contains(outunit,'/')
        outunit = sprintf('(%s)',outunit);
    end
    str = sprintf('[%s/%s]',outunit,inunit);
end

