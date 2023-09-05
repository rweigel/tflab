function [ax1,ax2] = zplot(S,popts)
%ZPLOT
%
%   ZPLOT(S)
%   ZPLOT(S, popts)
%
%   popts.type is one of, 1, 2, or 3
%
%   1 => Z,phase
%   2 => rho,phase (assumes Z in (mV/km)/nT and frequencies in Hz).
%   3 => Real,Imaginary

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

if nargin < 2
    popts = struct();
end

if iscell(S) && length(S) == 1
    S = S{1};
end

S = tflab_metadata(S);

ptype = 1; % 1 = Z,phi, 2 = rho,phi; 3 = Re,Im
if iscell(S)
    % TODO: Check all same.
    timeunit  = S{1}.Metadata.timeunit;
    timedelta = S{1}.Metadata.timedelta;
    opts = tflabplot_options(S, popts, ptype, 'zplot');
    Zstrs   = opts.zstrs;
    Rhostrs = opts.rhostrs;
    Phistrs = opts.phistrs;
else
    timeunit  = S.Metadata.timeunit;
    timedelta = S.Metadata.timedelta;
    opts = tflabplot_options(S, popts, ptype, 'zplot');
    Zstrs   = opts.zstrs;
    Rhostrs = opts.rhostrs;
    Phistrs = opts.phistrs;
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
    ax(1) = subplot('Position', opts.PositionTop);
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
                % TODO: Assumes Z in (mV/km)/nT and dt in seconds.
                y = z2rho(S.fe/timedelta, S.Z);
                if interp_Z_exists
                    yi = z2rho(S.fi/timedelta, S.Zi);
                end
                for j = 1:size(S.Z,2)
                    ls{j} = sprintf('$%s$',Rhostrs{j});
                end
            case 3
                y = real(S.Z);
                for j = 1:size(S.Z,2)
                    ls{j} = sprintf('Re$(%s)$',Zstrs{j});
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
                idx = xi < opts.period_range(1) | xi > opts.period_range(2);
                yi(idx) = NaN;
            end
        end
        
        plot(x, y, opts.line{:}, 'markersize', 20);
        colororder_(ax(1),y);
        hold on;grid on;box on;
        if interp_Z_exists
            plot(xi, yi, opts.line{:}, 'markersize', 15);
        end
        if length(ls) > 1
            ylabel(unitstr_(S.Metadata));
            legend(ls,opts.legend{:});
        else
            ylabel(sprintf('%s %s',ls{1},unitstr_(S.Metadata)));
        end
        titlestr(S,opts,'z');            
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if opts.type < 3
            set(gca,'YScale','log');            
            adjust_ylim('upper');
        else
            adjust_ylim('both');
        end
        adjust_yticks(1e-4);
        adjust_exponent('y');
        setx(opts,0,timeunit);

    ax(2) = subplot('Position', opts.PositionBottom);
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
            yl = unitstr_(S.Metadata);
            for j = 1:size(S.Z,2)
                ls{j} = sprintf('Im$(%s)$',Zstrs{j});
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
        plot(x, y, opts.line{:}, 'markersize', 20);
        hold on;grid on;box on;
        colororder_(ax(2),y);
        if interp_Z_exists
            plot(xi, yi, opts.line{:}, 'markersize', 15);
        end
        if length(ls) > 1
            ylabel(yl);
            legend(ls,opts.legend{:});
        else
            ylabel(sprintf('%s %s',ls{1},yl));
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if ~opts.unwrap && opts.type ~= 3
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
        end
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
        % Creates new figure for each column.
        figprep();
        if j > 1
            figure();
            figprep();
        end
        ax1(j) = subplot('Position', opts.PositionTop);
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
                        % TODO: Assumes Z in (mV/km)/nT
                        y = z2rho(S{s}.fe/timedelta, S{s}.Z(:,j));
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
                
                opts.line = lineopts_(size(S{s}.Z,1), s);
                h(s) = plot(x, y, opts.line{:});
                if opts.type < 3
                    set(gca,'YScale','log');
                end
                if opts.vs_period
                    set(gca,'XScale','log');
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
            yunitstr_ = unitstr_(S{1}.Metadata);
            if opts.type == 1
                yl = sprintf('$|%s|$%s',Zstrs{j},yunitstr_);
                adjust_ylim('upper');
            end
            if opts.type == 2
                yl = sprintf('$%s$ [$\\Omega\\cdot$m]',Rhostrs{j});
                adjust_ylim('upper');            
            end
            if opts.type == 3
                yl = sprintf('Re$(%s)$ %s',Zstrs{j},yunitstr_);
                adjust_ylim('both');            
            end
            if ~isempty(opts.period_range)
                set(gca,'XLim',opts.period_range);
            end
            ylabel(yl);
            legend(h,ls,opts.legend{:});
            adjust_yticks(1e-4);
            if ~isempty(timeunit) && opts.vs_period
                period_lines();
            end
            adjust_exponent('y');
            setx(opts,0,timeunit);
        ax2(j) = subplot('Position', opts.PositionBottom);
            for s = 1:length(S)
                yebl = [];
                yebu = [];
                if opts.type == 3
                    y = imag(S{s}.Z(:,j));
                    yl = sprintf('Im$(%s)$ %s',Zstrs{j},unitstr_(S{1}.Metadata));
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
                
                opts.line = lineopts_(size(S{s}.Z,1), s);
                if opts.vs_period
                    semilogx(x, y, opts.line{:});
                else
                    plot(x, y, opts.line{:});
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
                set(gca,'YScale','linear');
                set(gca,'YLim',[-180,180]);
                set(gca,'YTick',[-180:45:180]);
                adjust_ylim();
            else
                adjust_ylim('upper');
            end
            legend(ls,opts.legend{:});
            adjust_yticks(1e-4);
            adjust_exponent();
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

function line = lineopts_(Nz, s)
    ms = max(30-8*s,1);
    if Nz == 1
        ms = 30;
        line = {'.','markersize', ms};
    else
        line = {'linewidth',max(4-s,1),'marker','.','markersize',ms};
    end
end

function str = unitstr_(meta)

    str = '';
    if isempty(meta.inunit) || isempty(meta.outunit)
        return
    end
    inunit = sprintf('%s',meta.inunit);
    if contains(inunit,'/')
        inunit = sprintf('(%s)',inunit);
    end
    outunit = sprintf('%s',meta.outunit);
    if contains(outunit,'/')
        outunit = sprintf('(%s)',outunit);
    end
    str = sprintf('[%s/%s]',outunit,inunit);
end

