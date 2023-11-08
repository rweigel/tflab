function [ax1,ax2] = zplot(S,popts,comp)
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
if nargin < 3
    comp = 1;
end
if iscell(S) && length(S) == 1
    S = S{1};
end

S = tflab_metadata(S);

if iscell(S)
    % TODO: Check all same units and sfs or allow different.
    frequnit = S{1}.Metadata.frequnit;
    opts = tflabplot_options(S, popts, 'zplot');
    Zstrs   = opts.zstrs;
    Rhostrs = opts.rhostrs;
    Phistrs = opts.phistrs;
else
    frequnit = S.Metadata.frequnit;
    opts = tflabplot_options(S, popts, 'zplot');
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
        x = 1./(S.fe*S.Metadata.freqsf);
        if interp_Z_exists
            xi = 1./(S.fi*S.Metadata.freqsf);
        end
    else
        x = S.fe*S.Metadata.freqsf;
        if interp_Z_exists
            xi = S.fi*S.Metadata.freqsf;
        end
    end
    
    figprep();
    ax1 = subplot('Position', opts.PositionTop);
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
                % TODO: Assumes Z in (mV/km)/nT and f in Hz.
                y = z2rho(x, S.Z);
                if interp_Z_exists
                    yi = z2rho(x, S.Zi);
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
        colororder_(ax1,y);
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
        setx(opts,0,frequnit);

    ax2 = subplot('Position', opts.PositionBottom);
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
        colororder_(ax2,y);
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
        setx(opts,1,frequnit);

    if opts.print
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s.%s',opts.printname, opts.printfmt{i});
            figsave(fullfile(opts.printdir, fname), opts);
        end
    end
end % if isstruct(S)

% Multiple transfer functions
if iscell(S)

    if nargin < 3
        for comp = 1:size(S{1}.Z,2)
            if comp > 1
                figure;
            end
            [ax1(comp),ax2(comp)] = zplot(S,popts,comp);
        end
        return;
    end
    
    figprep();

    ax1 = subplot('Position', opts.PositionTop);
        [x,y,yebu,yebl] = xycell(S,opts,comp,'top');
        for s = 1:length(S)
            ls{s} = sprintf('%s',S{s}.Options.description);
            opts.line = lineopts_(size(S{s}.Z,1), s);
            h(s) = plot(x{s}, y{s}, opts.line{:});
            if s == 1
                grid on;hold on;box on;
            end
            if ~isempty(yebl{s})
                errorbars(x{s},y{s},yebl{s},yebu{s});
                %return
                %b = squeeze(S{s}.Segment.Z(:,j,:));
                %for k = 2:length(x)
                %    plot(x(k),abs(b(k,:)),'k.','MarkerSize',1);
                %end
            end
        end
        if opts.type < 3
            set(gca,'YScale','log');
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        yunitstr_ = unitstr_(S{1}.Metadata);
        if opts.type == 1
            yl = sprintf('$|%s|$%s',Zstrs{comp},yunitstr_);
        end
        if opts.type == 2
            yl = sprintf('$%s$ [$\\Omega\\cdot$m]',Rhostrs{comp});
        end
        if opts.type == 3
            yl = sprintf('Re$(%s)$ %s',Zstrs{comp},yunitstr_);
        end
        if ~isempty(opts.period_range)
            set(gca,'XLim',opts.period_range);
        end
        ylabel(yl);
        legend(h,ls,opts.legend{:});
        adjust_yticks(1e-4);
        adjust_exponent('y');
        if opts.type == 3
            adjust_ylim('both');
        else
            adjust_ylim('upper');
        end
        setx(opts,0,frequnit);

    ax2 = subplot('Position', opts.PositionBottom);
        [x,y,yebu,yebl] = xycell(S,opts,comp,'bottom');
        for s = 1:length(S)
            ls{s} = sprintf('%s',S{s}.Options.description);

            opts.line = lineopts_(size(S{s}.Z,1), s);
            if opts.vs_period
                semilogx(x{s}, y{s}, opts.line{:});
            else
                plot(x{s}, y{s}, opts.line{:});
            end

            if s == 1
                grid on;box on;hold on;
            end

            if ~isempty(yebl{s})
                errorbars(x{s},y{s},yebl{s},yebu{s});
            end

        end

        yl = sprintf('$%s$ [$^\\circ]$',Phistrs{comp});
        if opts.type == 3
            yl = sprintf('Im$(%s)$ %s',Zstrs{comp},unitstr_(S{1}.Metadata));
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
        setx(opts,1,frequnit);

    if opts.print
        comp = regexprep(Zstrs{j},'\{|\}','');
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s-%s.%s',opts.printname, comp, opts.printfmt{i});
            figsave(fullfile(opts.printdir, fname), opts);
        end
    end
    
end % if iscell(S)

end % function

function [x,y,yebu,yebl] = xycell(S,opts,comp,panel)

    % TODO: Assumes units are the same for all Zs.
    %       Check this before plotting.

    for s = 1:length(S)
        if opts.vs_period
            x{s} = 1./(S{s}.fe*S{s}.Metadata.freqsf);
        else
            x{s} = S{s}.fe*S{s}.Metadata.freqsf;
        end
        
        yebu{s} = [];
        yebl{s} = [];
        if strcmp(panel,'top')
            switch opts.type
                case 1
                    y{s} = abs(S{s}.Z(:,comp));
                    if isfield(S{s},'ZCL')
                        %yebl =  y-squeeze(S{s}.ZCL.Magnitude.Normal.x_1sigma(:,j,1));
                        %yebu = -y+squeeze(S{s}.ZCL.Magnitude.Normal.x_1sigma(:,j,2));
                        yebl{s} =  y-squeeze(S{s}.ZCL.Magnitude.Bootstrap.x_1sigma(:,j,1));
                        yebu{s} = -y+squeeze(S{s}.ZCL.Magnitude.Bootstrap.x_1sigma(:,j,2));
                    end
                case 2
                    % TODO: Check that f in Hz and Z in (mV/km)/nT
                    f = S{s}.fe*S{s}.Metadata.freqsf;
                    y{s} = z2rho(f, S{s}.Z(:,comp));
                case 3
                    y{s} = real(S{s}.Z(:,comp));
                    if isfield(S{s},'ZCL')
                        yebl{s} =  y-squeeze(real(S{s}.ZCL.Z.Normal.x_1sigma(:,comp,1)));
                        yebu{s} = -y+squeeze(real(S{s}.ZCL.Z.Normal.x_1sigma(:,comp,2)));
                    end
            end
            if opts.vs_period && ~isempty(opts.period_range)
                idx = x{s} <= opts.period_range(1) | x{s} >= opts.period_range(2);
                y{s}(idx) = NaN;
            end
        end
        
        if strcmp(panel,'bottom')
            if opts.type == 3
                y{s} = imag(S{s}.Z(:,comp));
                if isfield(S{s},'ZCL')
                    yebl{s} = squeeze(imag(S{s}.ZCL.Z.Normal.x_1sigma(:,comp,1)));
                    yebu{s} = squeeze(imag(S{s}.ZCL.Z.Normal.x_1sigma(:,comp,2)));
                end
            else
                if opts.unwrap
                    y{s} = (180/pi)*unwrap(atan2(imag(S{s}.Z(:,comp)),real(S{s}.Z(:,comp))));
                else
                    y{s} = (180/pi)*atan2(imag(S{s}.Z(:,comp)),real(S{s}.Z(:,comp)));
                end
            end
            if opts.vs_period && ~isempty(opts.period_range)
                idx = x{s} <= opts.period_range(1) | x{s} >= opts.period_range(2);
                y{s}(idx) = NaN;
            end            
        end

    end
end

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

