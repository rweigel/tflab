function ax = dftplot(S,popts,comp)
%DFTPLOT   Plot DFT derived quantities
%
%   DFTPLOT(S)
%   DFTPLOT(S, popts)
%   DFTPLOT(S, popts, component)
%
%   For non-prediction error related plots, popts.type has form 
%   'a-b-c', where
%
%   a is one of: original, final, detrended, windowed, prewhitened, zeropadded
%   b is one of: raw, averaged
%   c is one of: magnitudes,phases
%
%   For prediction error related plots, popts.type has the form
%   'error-b-c', where
%
%   b is one of: raw, averaged
%   c is one of: magphase, realimag

assert(isstruct(S) || iscell(S), ...
    'S must be a tflab struct or cell array of tflab structs');

if nargin < 2
    popts = struct();
    comp = 1; % Used if iscell(S);
end
if nargin < 3
    comp = 1;
end

popts = tflabplot_options(S, popts, 'original', 'dftplot');


tparts = split(popts.type,'-');
tparts1 = {'error',...
           'original','final','detrended','windowed','prewhitened','zeropadded'};
if ~any(strcmp(tparts{1},tparts1))
    list = join(tparts1,', ');
    error('popts.ptype must start with one of: %s.',list{1});
end

if ~iscell(S)
    S = {S};
end

S = defaultinfo(S);

if size(S{1}.In,2) ~= size(S{1}.Out,2)
    error('Case of size(S.In,2) ~= size(S.Out,2) not handled');
end

% TODO: Check all same timeunit. If not convert to same.
timeunit = S{1}.Options.info.timeunit;
timedelta = S{1}.Options.info.timedelta;

        
for s = 1:length(S)
    opts = S{s}.Options;
    if strcmp(tparts{1},'error')
        % error-{raw,averaged}-{magphase,realimag}
        Error = S{s}.Out_.Predicted - S{s}.Out;

        if strcmp(tparts{2},'averaged')
            [segse, f, fe{s}] = dftsegments(Error, opts);
            w = dftweights(f, [], [], opts);
            dfte = dftaverage(segse, w);
        else
            [dfte,fe{s}] = fftu(Error(:,comp));
        end

        if strcmp(tparts{3},'magphase')
            y1{s} = abs(dfte);
            y2{s} = (180/pi)*angle(y1{s});
        else
            y1{s} = real(dfte);
            y2{s} = imag(dfte);
        end
    else
        % {original,detrended,windowed,prewhitened,zeropadded,final}-{raw,averaged}-{magnitudes,phases}
        if any(strcmp(tparts{1},{'detrended','windowed','prewhitened','zeropadded','final'}))
            % Capatilize first letter
            tparts1uc = [upper(tparts{1}(1)),tparts{1}(2:end)];
            if ~isfield(S{s}.In_,tparts1uc)
                error('Time series was not %s.',tparts{1});
            end
            [segsIn, f, fc] = dftsegments(S{s}.In_.(tparts1uc), opts);
            segsOut = dftsegments(S{s}.Out_.(tparts1uc), opts);
        else
            % original
            [segsIn, f, fc] = dftsegments(S{s}.In, opts);
            segsOut = dftsegments(S{s}.Out, opts);
        end
        
        if strcmp(tparts{2},'averaged')
            w = dftweights(f, segsOut, segsIn, popts);
            DFTIn{s} = dftaverage(segsIn(:,comp), w);
            DFTOut{s} = dftaverage(segsOut(:,comp), w);
            fe{s} = fc;
        else
            tmp = cat(1,f{:});
            fe{s,1} = tmp(:);
            tmp = cat(1,segsIn{:});
            DFTIn{s} = tmp(:,comp);
            tmp = cat(1,segsOut{:});
            DFTOut{s} = tmp(:,comp);
        end
        if strcmp(tparts{3},'phases')
            y1{s} = (180/pi)*angle(DFTIn{s});
            y2{s} = (180/pi)*angle(DFTOut{s});
        else
            y1{s} = abs(DFTIn{s});
            y2{s} = abs(DFTOut{s});
        end
    end
    if popts.vs_period
        x{s} = (timedelta)./fe{s};
    else
        x{s} = fe{s}/(timedelta);
    end
end

figprep();
if strcmp(tparts{1},'error')

    [~,lg2] = legend_(S);
    [yl1, yl2] = ylabelerror_(S,tparts{3});

    ax(1) = subplot('Position',popts.PositionTop);
        plot_(x,y1,popts);
        colororder_(ax(1),y1);
        grid on;box on;
        set(gca,'YScale','log');
        if popts.vs_period
            set(gca,'XScale','log');
        end
        if ~isempty(lg2)
            legend(lg2,popts.legend{:});
        else
            titlestr(S,popts,'psd');
        end
        ylabel(yl1);
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(popts,0,timeunit);

    ax(2) = subplot('Position',popts.PositionBottom);
        plot_(x,y2,popts);
        colororder_(ax(2),y2);
        grid on;box on;        
        if popts.vs_period
            set(gca,'XScale','log');
        end
        ylabel(yl2);
        if ~isempty(lg2)
            legend(lg2,popts.legend{:});
        end
        if strcmp(tparts{1},'magphase')
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
        end
        adjust_ylim();
        adjust_exponent('y');
        setx(popts,1,timeunit);
end

if ~strcmp(tparts{1},'error')
    
    [lg1, lg2] = legend_(S);
    [yl1, yl2] = ylabel_(S,tparts{3});
    
    ax(1) = subplot('Position',popts.PositionTop);
        plot_(x,y1,popts);
        colororder_(ax(1),y1);
        grid on;box on;
        if strcmp(tparts{3},'phases')
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
        else
            set(gca,'YScale','log');
        end
        if popts.vs_period
            set(gca,'XScale','log');
        end
        if ~isempty(lg1)
            legend(lg1,popts.legend{:});
        end
        ylabel(yl1);
        titlestr(S, popts, 'dft');
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(popts, 0, timeunit);

    ax(2) = subplot('Position',popts.PositionBottom);
        plot_(x,y2,popts);
        colororder_(ax(2),y2);
        grid on;box on;
        set(gca,'YScale','log');
        if strcmp(tparts{3},'phases')
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
        else
            set(gca,'YScale','log');
        end
        if popts.vs_period
            set(gca,'XScale','log');
        end
        if ~isempty(lg2)
            legend(lg2,popts.legend{:});
        end
        ylabel(yl2);
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(popts, 1, timeunit);
end

if popts.print
    compstr = '';
    if length(S) > 1 && size(S{1}.In,2) > 1
        compstr = sprintf('_%s',comp);
    end
    for i = 1:length(popts.printfmt)
        fname = sprintf('%s_%s%s.%s',...
                popts.printname, popts.type, compstr, popts.printfmt{i});
        figsave(fullfile(popts.printdir, fname));
    end
end

if length(S) > 1 && comp < size(S{1}.In,2)
    figure();
    psdplot(S,popts,comp+1);
end

end

function plot_(x,y,popts)
    if iscell(x) && iscell(y)
        hold on;
        for s = 1:length(x)
        	plot(x{s},y{s},popts.line{:});
        end
        hold off;
    else
        plot(x,y,popts.line{:});
    end
end

function [yl1, yl2] = ylabel_(S, what)
    info = S{1}.Options.info;
    if strcmp(what,'magnitudes')
        inunit = unitstr(info.inunit);
        outunit = unitstr(info.outunit);
    else
        inunit = '$[^\circ]$';
        outunit = '$[^\circ]$';
    end
    
    if length(S) > 1
        yl1 = [labelstr_(info.instr), ' ',inunit];
        yl2 = [labelstr_(info.outstr),' ',outunit];
    else
        if size(S{1}.In,2) > 1
            yl1 = unitstr(info.inunit);
        else
            yl1 = [labelstr_(info.instr),' ',inunit];
        end
        if size(S{1}.Out,2) > 1
            yl2 = unitstr(info.outunit);
        else
            yl2 = [labelstr_(info.outstr),' ',outunit];
        end
    end
    if strcmp(what,'phases')
        yl1 = ['$\angle$',' ',yl1];
        yl2 = ['$\angle$',' ',yl2];
    end

end

function s = labelstr_(labelstr, prefix, mag)
    if nargin < 2, prefix = ''; end
    if nargin < 3, mag = 1; end
    if iscell(labelstr)
        labelstr = labelstr{comp};
    end
    if contains(labelstr,'$')
        s = replace(labelstr,'$','');
    else
        s = ['\mbox{',labelstr,'}'];
    end
    if mag
        s = ['$|','\widetilde{',prefix, ' ', s, '}|$'];
    else
        s = ['$\angle$ of $','\widetilde{', prefix, ' ', s, '}\;[^\circ]$'];
    end
end

function [yl1, yl2] = ylabelerror_(S,what)
    info = S{1}.Options.info;
    yl1 = [labelstr_(info.outstr,'\Delta',1), ' ',unitstr(info.outunit)];
    yl2 = labelstr_(info.outstr,'\Delta',0);
end

function [lg1, lg2] = legend_(S)

    if length(S) > 1
        for s = 1:length(S)
            info = S{s}.Options.info;
            lg1{s} = S{s}.Options.description;
            lg2{s} = S{s}.Options.description;
        end
    else
        lg1 = '';
        lg2 = '';
        if size(S{1}.In) > 1
            for j = 1:length(info.instr)
                lg1{j} = labelstr_(info.instr{j});
            end
        end
        if size(S{1}.Out) > 1
            for j = 1:length(info.outstr)
                lg2{j} = labelstr_(info.outstr{j});
            end
        end
    end    
end

function ls1 = plotnoise(ls1,comp,unit)

    if isfield(S,comp) && strcmp(popts.type,'raw')
        hold on;
        y = ftrim_(S.fe,S.(comp));
        loglog(x,y,popts.line{:});
        jl = size(S.(comp),2);
        if strcmp(comp,'InNoisePSD')
            compstrs = info.instr;
        else
            compstrs = info.outstr;
        end
        for j = 1:jl
            if iscell(info.instr)
                ls1{j+jl} = sprintf('%s(:,%d) Noise Amplitudes%s',compstrs{j},j,unit);
            else
                ls1{j+jl} = sprintf('%s Noise Amplitudes%s',compstrs,unit);
            end
        end
    end
    
end
