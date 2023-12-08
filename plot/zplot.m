function [ax1,ax2] = zplot(S,popts,comps)
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

if isstruct(S)
    S = {S};
end
if nargin < 2
    popts = struct();
end
if nargin < 3
    comps = 1:size(S{1}.Z,2);
end
comps = sort(comps);

S = tflab_metadata(S);

% TODO: Check all same units and sfs or allow different.
frequnit = S{1}.Metadata.frequnit;
opts = tflabplot_options(S, popts, 'zplot');
Zstrs   = opts.zstrs;
Rhostrs = opts.rhostrs;
Phistrs = opts.phistrs;

% Single transfer function
if length(S) == 1

    interp_Z_exists = 0;
    if isfield(S,'Zi')
        interp_Z_exists = 1;
    end

    if opts.vs_period
        if interp_Z_exists
            xi = 1./(S{1}.fi*S.Metadata.freqsf);
        end
    else
        if interp_Z_exists
            xi = S{1}.fi*S{1}.Metadata.freqsf;
        end
    end
    
    figprep();
    ax1 = subplot('Position', opts.PositionTop);
        switch opts.type
            case 1
                if interp_Z_exists
                    yi = abs(S.Zi);
                end
                for j = 1:length(comps)
                    ls{j} = sprintf('$|%s|$',Zstrs{comps(j)});
                end
            case 2
                if interp_Z_exists
                    yi = z2rho(x, S.Zi);
                end
                for j = 1:length(comps)
                    ls{j} = sprintf('$%s$',Rhostrs{comps(j)});
                end
            case 3
                for j = 1:length(comps)
                    ls{j} = sprintf('Re$(%s)$',Zstrs{comps(j)});
                end
                if interp_Z_exists
                    yi = real(S.Zi);
                    jx = size(S.Z,2);
                    for j = 1:length(comps)
                        ls{jx+j} = sprintf('Re$(%s)$ Interpolated',Zstrs{comps(j)});
                    end
                end
        end

        for j = 1:length(comps)
            [x,y(:,j),dyebu(:,j),dyebl(:,j)] = xyvals_(S,opts,comps(j),'top');
        end

        if interp_Z_exists
            yi = y(:,comps);
        end
        opts.line = lineopts_(size(S{1}.Z,1), 1);

        h1 = plot(x, y, opts.line{:});
        hold on;grid on;box on;
        titlestr(S,opts,'z');
        colororder_(ax1,y);
        
        if interp_Z_exists
            plot(xi, yi, opts.line{:}, 'markersize', 15);
        end

        if length(ls) > 1
            ylabel(unitstr_(S{1}.Metadata));
            legend(ls,opts.legend{:});
        else
            ylabel(sprintf('%s %s',ls{1},unitstr_(S{1}.Metadata)));
        end
        if opts.vs_period
            set(gca,'XScale','log');
        end
        if opts.type < 3
            set(gca,'YScale','log');            
            adjust_ylim('upper');
        else
            adjust_ylim('both');
        end

        if ~isempty(dyebl)
            colors = get(h1,'Color');
            if ~iscell(colors)
                colors = {colors};
            end
            for j = 1:length(comps)
                errorbars(x,y(:,j),dyebl(:,j),dyebu(:,j),...
                         'y','Color',colors{j},'LineWidth',2);
            end
        end
        adjust_yticks(1e-4);
        adjust_exponent('y');
        setx(opts,0,frequnit);
    ax2 = subplot('Position', opts.PositionBottom);
        if opts.type ~= 3
            yl = '[$^\circ$]';
            if opts.unwrap
                if interp_Z_exists
                    yi = (180/pi)*unwrap(atan2(imag(S.Zi),real(S.Zi)));
                end
            else
                if interp_Z_exists
                    yi = (180/pi)*(atan2(imag(S.Zi),real(S.Zi)));
                end
            end
            for j = 1:length(comps)
                ls{j} = sprintf('$%s$',Phistrs{comps(j)});
            end
        else
            yl = unitstr_(S{1}.Metadata);
            for j = 1:length(comps)
                ls{j} = sprintf('Im$(%s)$',Zstrs{comps(j)});
            end
            if interp_Z_exists
                yi = imag(S.Zi);
                jx = size(S.Z,2);
                for j = 1:length(comps)
                    ls{jx+j} = sprintf('Im$(%s)$ Interpolated',Zstrs{comps(j)});
                end
            end
        end
        for j = 1:length(comps)
            [x(:,j),y(:,j)] = xyvals_(S,opts,comps(j),'bottom');
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
if length(S) > 1

    if length(comps) == 1
        comp = comps;
    else
        for comp = 1:length(comps)
            if comp > 1
                figure;
            end
            [ax1(comp),ax2(comp)] = zplot(S,popts,comp);
        end
        return;
    end
    
    figprep();

    ax1 = subplot('Position', opts.PositionTop);

        [x,y,dyebu,dyebl] = xyvals_(S,opts,comp,'top');

        for s = 1:length(S)
            ls{s} = sprintf('%s',S{s}.Options.description);
            opts.line = lineopts_(size(S{s}.Z,1), s);            
            h(s) = plot(x{s}, y{s}, opts.line{:});
            if s == 1
                grid on;hold on;box on;
            end
        end

        legend(h,ls,opts.legend{:});

        for s = 1:length(S)
            if ~isempty(dyebl{s})
                errorbars(x{s},y{s},dyebl{s},dyebu{s},...
                    'y','Color',get(h(s),'Color'),'LineWidth',2);
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
        adjust_yticks(1e-4);
        adjust_exponent('y');
        if opts.type == 3
            adjust_ylim('both');
        else
            adjust_ylim('upper');
        end
        setx(opts,0,frequnit);

    ax2 = subplot('Position', opts.PositionBottom);

        [x,y,dyebu,dyebl] = xyvals_(S,opts,comp,'bottom');
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
        end

        yl = sprintf('$%s$ [$^\\circ]$',Phistrs{comp});
        if opts.type == 3
            yl = sprintf('Im$(%s)$ %s',Zstrs{comp},unitstr_(S{1}.Metadata));
        end
        ylabel(yl);

        if opts.type ~= 3 && ~opts.unwrap
            set(gca,'YScale','linear');
            set(gca,'YLim',[-180,180]);
            set(gca,'YTick',-180:45:180);
            adjust_ylim();
        else
            adjust_ylim('upper');
        end
        legend(ls,opts.legend{:});
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(opts,1,frequnit);

    if opts.print
        comp = regexprep(Zstrs{comp},'\{|\}','');
        for i = 1:length(opts.printfmt)
            fname = sprintf('%s-%s.%s',opts.printname, comp, opts.printfmt{i});
            figsave(fullfile(opts.printdir, fname));
        end
    end
    
end % if iscell(S)

end % function

function [dyebu,dyebl] = compute_errorbars_(S,opts,panel,comp)

    dyebu = zeros(size(S.Z(:,comp)));
    dyebl = zeros(size(S.Z(:,comp)));

    if isfield(S,'Regression') && isfield(S.Regression,'Bootstrap')
        if opts.type == 1 && strcmp(panel,'top')
            dyebl = abs(S.Z(:,comp)) - S.Regression.Bootstrap.ZCL95l(:,comp);            
            dyebu = S.Regression.Bootstrap.ZCL95u(:,comp) - abs(S.Z(:,comp));
            %dyebl = abs(S.Z(:,comp)) - 2*sqrt(S.Regression.Bootstrap.ZVAR(:,comp));
            %dyebu = 2*sqrt(S.Regression.Bootstrap.ZVAR(:,comp)) - abs(S.Z(:,comp));
        end
        return;
    end

    if ~isfield(S,'ZVAR')
        return;
    end

    % https://ds.iris.edu/spud/resources/pdf/SPUD-XML-change-log.pdf
    % (cached in runs/data/EarthScope/SPUD-XML-change-log.pdf)
    % "To compute the standard error of a real or imaginary part, e.g.,
    % for plotting, one would divide the variance by 2, then take square root
    % (i.e. Std = sqrt(Var/2)). To plot the error bars for a real or imaginary part,
    % one would then multiply that number by 2 ..."

    % http://ds.iris.edu/spudservice/data/15014569
    % Here I multiply by sqrt(2) to get an approxmate visual
    % match to error bars at http://ds.iris.edu/spud/emtf/15014571.
    % The documentation does not indicate if the error bars
    % are 1 or 2 Standard Errors. Also, it may be that the plots
    % have not been updated since 2018 when the
    % reprocessing of the files was performed as mentioned
    % in https://ds.iris.edu/spud/resources/pdf/SPUD-XML-change-log.pdf
    %
    % TODO: Find out what is plotted at http://ds.iris.edu/spud/emtf/15014571
    % so that I can state the error bar length in SEs.

    % See also
    %   https://ds.iris.edu/files/products/emtf/definitions/
    %   http://ds.iris.edu/files/products/emtf/definitions/variance.html
    %   https://library.seg.org/doi/epdf/10.1190/geo2018-0679.1
    
    if strcmp(panel,'top')
        dZ = sqrt(2)*sqrt(S.ZVAR(:,comp));
        if opts.type == 1
            dyebl = dZ;
            dyebu = dZ;
        end
        if opts.type == 2
            % ρ = |Z|^2/(5*f) => Δρ = 2|Z|dZ/(5*f)
            dyebl = 2*abs(S.Z(:,comp)).*dZ./(5*S.fe);
            dyebu = yebl;
        end
        if opts.type == 3
            dyebl = dZ;
            dyebu = dZ;
        end            
    end
    if strcmp(panel,'bottom')
        if opts.type == 3
            dZ = sqrt(2)*sqrt(S.ZVAR(:,comp)/2);
            dyebl = dZ;
            dyebu = dZ;
        end            
    end
end

function [x,y,dyebu,dyebl] = xyvals_(S,opts,comp,panel)

    % TODO: Assumes units are the same for all Zs.
    %       Check this before plotting.

    for s = 1:length(S)
        
        if strcmp(panel,'top')
            if opts.type == 1
                y{s} = abs(S{s}.Z(:,comp));
            end
            if opts.type == 2
                f = S{s}.fe*S{s}.Metadata.freqsf;
                y{s} = z2rho(f, S{s}.Z(:,comp));
            end
            if opts.type == 3
                y{s} = real(S{s}.Z(:,comp));
            end
        end
        
        if strcmp(panel,'bottom')
            if opts.type == 3
                y{s} = imag(S{s}.Z(:,comp));
            else
                if opts.unwrap
                    y{s} = (180/pi)*unwrap(atan2(imag(S{s}.Z(:,comp)),real(S{s}.Z(:,comp))));
                else
                    y{s} = (180/pi)*atan2(imag(S{s}.Z(:,comp)),real(S{s}.Z(:,comp)));
                end
            end
        end

        [dyebu{s},dyebl{s}] = compute_errorbars_(S{s},opts,panel,comp);            

        if opts.vs_period
            x{s} = 1./(S{s}.fe*S{s}.Metadata.freqsf);
        else
            x{s} = S{s}.fe*S{s}.Metadata.freqsf;
        end
        
        if opts.vs_period && ~isempty(opts.period_range)
            idx = x{s} <= opts.period_range(1) | x{s} >= opts.period_range(2);
            y{s}(idx) = NaN;
        end
    end
    if length(S) == 1
        x = x{s};
        y = y{s};
        dyebl = dyebl{s};
        dyebu = dyebu{s};
    end
end

function line = lineopts_(Nz, s);
    if Nz == 1 % Single point
        ms = 30;
        line = {'.','markersize', ms};
    else
        %ms = max(30-8*s,1);
        ms = 30;
        line = {'.','markersize',ms};
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

