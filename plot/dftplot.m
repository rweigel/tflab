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
%   a is one of: original, final, detrended, windowed, whitened, zeropadded
%   b is one of: raw, averaged
%   c is one of: magnitudes,phases,reals,imaginaries
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
if ~iscell(S)
    S = {S};
end

% Apply default metadata for fields not specified in S.Metadata.
S = tflab_metadata(S);

popts = tflabplot_options(S, popts, 'dftplot');
argcheck_(S, popts)

% TODO: Check all same. If not convert to same.
frequnit = S{1}.Metadata.frequnit;

tparts = split(popts.type,'-');
for s = 1:length(S)
    opts = S{s}.Options;
    if strcmp(tparts{1},'error')
        % error-{raw,averaged}-{magphase,realimag}
        if strcmp(tparts{2},'averaged')
            dftError = dftaverage(S{s}.DFT.Out_.Error);
            fe{s} = S{s}.DFT.fe;            
        else
            tmp = cat(1,S{s}.DFT.f{:});
            fe{s,1} = tmp(:);
            dftError = cat(1,S{s}.DFT.Out_.Error{:});
        end
        dftError = dftError(:,comp);
        if strcmp(tparts{3},'magphase')
            sf = (size(S{s}.Out_.Error,1)-1)/2;
            y1{s} = abs(dftError)/sf;
            y2{s} = (180/pi)*angle(dftError);
            if fe{s}(end) == 0.5
                y1{s}(end) = y1{s}(end)/2;
            end
        else
            y1{s} = real(dftError);
            y2{s} = imag(dftError);
        end
    else
        % {original,detrended,windowed,whitened,zeropadded,final}-{raw,averaged}-{magnitudes,phases}
        if any(strcmp(tparts{1},{'detrended','windowed','whitened','zeropadded','final'}))
            % Capatilize first letter
            tparts1uc = [upper(tparts{1}(1)),tparts{1}(2:end)];
            if ~isfield(S{s}.DFT,'In_') ||  ~isfield(S{s}.DFT.In_,tparts1uc)
                error('Time series was not %s.',tparts{1});
            end
            segsIn = S{s}.DFT.In_.(tparts1uc);
            segsOut = S{s}.DFT.Out_.(tparts1uc);
            f = S{s}.DFT.f_;
            fe{s} = S{s}.DFT.fe_;
        else
            % original
            segsIn = S{s}.DFT.In;
            segsOut = S{s}.DFT.Out;
            f = S{s}.DFT.f;
            fe{s} = S{s}.DFT.fe;
        end

        if strcmp(tparts{2},'averaged')
            [wIn,wOut] = dftweights(f, segsOut, segsIn, opts);
            DFTIn{s}  = dftaverage(segsIn, wIn);
            DFTOut{s} = dftaverage(segsOut, wOut);
        else
            tmp = cat(1,f{:});
            fe{s,1} = tmp(:);
            DFTIn{s} = cat(1,segsIn{:});
            DFTOut{s} = cat(1,segsOut{:});
        end
        
        if length(S) > 1
            % Compare mode. One plot per component.
            DFTIn{s} = DFTIn{s}(:,comp);
            DFTOut{s} = DFTOut{s}(:,comp);
        end

        if strcmp(tparts{3},'phases')
            y1{s} = (180/pi)*angle(DFTIn{s});
            y2{s} = (180/pi)*angle(DFTOut{s});
        end
        if strcmp(tparts{3},'reals')
            y1{s} = real(DFTIn{s});
            y2{s} = real(DFTOut{s});
        end
        if strcmp(tparts{3},'imaginaries')
            y1{s} = imag(DFTIn{s});
            y2{s} = imag(DFTOut{s});
        end
        if strcmp(tparts{3},'magnitudes')
            sf = (size(S{s}.In,1)-1)/2;
            y1{s} = abs(DFTIn{s})/sf;
            y2{s} = abs(DFTOut{s})/sf;
            if fe{s}(end) == 0.5
                y1{s}(end) = y1{s}(end)/2;
                y2{s}(end) = y2{s}(end)/2;
            end
        end
    end
    if popts.vs_period
        x{s} = 1./(fe{s}*S{s}.Metadata.freqsf);
    else
        x{s} = fe{s}*S{s}.Metadata.freqsf;
    end
    if length(S) == 1
        x = x{1};
        y1 = y1{1};
        y2 = y2{1};
    end
end

figprep();
if strcmp(tparts{1},'error')

    lg = '';
    if length(S) > 1
        [lg,~] = legend_(S,tparts{3},popts);
    end
    [yl1, yl2] = ylabelerror_(S,tparts{3},popts,comp);

    ax(1) = subplot('Position',popts.PositionTop);
        plot_(x,y1,popts);
        colororder_(ax(1),y1);
        grid on;box on;
        set(gca,'YScale','log');
        if popts.vs_period
            set(gca,'XScale','log');
        end
        if ~isempty(lg)
            legend(lg,popts.legend{:});
        else
            titlestr(S{1},popts,'dft');
        end
        ylabel(yl1);
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(popts,0,frequnit);

    ax(2) = subplot('Position',popts.PositionBottom);
        plot_(x,y2,popts);
        colororder_(ax(2),y2);
        grid on;box on;        
        if popts.vs_period
            set(gca,'XScale','log');
        end
        ylabel(yl2);
        if ~isempty(lg)
            legend(lg,popts.legend{:});
        end
        if strcmp(tparts{1},'magphase')
            set(gca,'YScale','linear');
            set(gca,'YTick',-180:45:180);
        end
        adjust_ylim();
        adjust_exponent('y');
        setx(popts,1,frequnit);
end

if ~strcmp(tparts{1},'error')

    [lg1, lg2] = legend_(S,tparts{3},popts);
    [yl1, yl2] = ylabel_(S,tparts{3},popts,comp);
    
    ax(1) = subplot('Position',popts.PositionTop);
        plot_(x,y1,popts);
        colororder_(ax(1),y1);
        grid on;box on;
        setscales_(tparts{3},popts.vs_period)
        if ~isempty(lg1)
            legend(lg1,popts.legend{:});
        else
            titlestr(S{1},popts,'dft');
        end
        ylabel(yl1);
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(popts,0,frequnit);

    ax(2) = subplot('Position',popts.PositionBottom);
        plot_(x,y2,popts);
        colororder_(ax(2),y2);
        grid on;box on;
        setscales_(tparts{3},popts.vs_period)
        if ~isempty(lg2)
            legend(lg2,popts.legend{:});
        end
        ylabel(yl2);
        adjust_ylim('upper');
        adjust_yticks(1e-4);
        adjust_exponent();
        setx(popts,1,frequnit);
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

if strcmp(tparts{1},'error')
    if comp < size(S{1}.In,2)
        figure();
        dftplot(S,popts,comp+1);
    end
else
    if length(S) > 1 && comp < size(S{1}.In,2)
        figure();
        dftplot(S,popts,comp+1);
    end
end

end % function

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

function setscales_(tparts3, vs_period)
    if strcmp(tparts3,'magnitudes')
        set(gca,'YScale','log');
    end
    if strcmp(tparts3,'phases')
        set(gca,'YScale','linear');
        set(gca,'YTick',-180:45:180);
    end
    if vs_period
        set(gca,'XScale','log');
    end
end

function [lg1, lg2] = legend_(S,what,popts)

    if length(S) > 1
        % Comparing two tfs. Each plot will be for one component.
        for s = 1:length(S)
            lg1{s} = S{s}.Options.description;
            lg2{s} = lg1{s};
        end
    else
        if strcmp(what,'magnitudes')
            what = 'magnitude';
        end
        if strcmp(what,'phases')
            what = 'angle';
        end
        % Single tf. All components on same plot.
        lg1 = '';
        lg2 = '';
        if size(S{1}.In) > 1
            % One legend entry per per component.
            for j = 1:length(popts.instr)
                lg1{j} = labelstr_(popts.instr{j},'',what);
            end
        end
        if size(S{1}.Out) > 1
            for j = 1:length(popts.outstr)
                lg2{j} = labelstr_(popts.outstr{j},'',what);
            end
        end
    end    
end

function [yl1, yl2] = ylabel_(S,what,popts,comp)

    meta = S{1}.Metadata;
    
    if length(S) > 1
        if strcmp(what,'magnitudes')
            yl1 = [labelstr_(popts.instr{comp}), ' ',meta.inunit];
            yl2 = [labelstr_(popts.outstr{comp}),' ',meta.outunit];
        end
        if strcmp(what,'phases')
            yl1 = labelstr_(popts.instr{comp},'','angle');
            yl2 = labelstr_(popts.outstr{comp},'','angle');
        end
        if strcmp(what,'reals')
            yl1 = [labelstr_(popts.instr{comp},'','real'), ' ',meta.inunit];
            yl2 = [labelstr_(popts.outstr{comp},'','real'),' ',meta.outunit];
        end
        if strcmp(what,'imaginaries')
            yl1 = [labelstr_(popts.instr{comp}),'','imaginary', ' ',meta.inunit];
            yl2 = [labelstr_(popts.outstr{comp}),'','imaginary',' ',meta.outunit];
        end
    else
        yl1 = single_(S{1}.In, what, comp, popts.instr, meta.inunit);
        yl2 = single_(S{1}.Out, what, comp, popts.outstr, meta.outunit);
    end
    
    function yl = single_(X, what, comp, instr, inunit)
        if size(X,2) > 1
            if any(strcmp(what,{'magnitudes','reals','imaginaries'}))
                yl = unitstr(inunit);
            end
            if strcmp(what,'phases')
                yl = '$\angle$ $[^\circ]$';
            end
        else
            if any(strcmp(what,{'magnitudes','reals','imaginaries'}))
                yl = [labelstr_(instr{comp},' ',what),' ',inunit];
            end
            if strcmp(what,'phases')
                yl = labelstr_(instr{comp},'','angle');
            end
        end
    end
end

function [yl1, yl2] = ylabelerror_(S,what,popts,comp)
    meta = S{1}.Metadata;

    if strcmp(what,'magphase')
        l1 = labelstr_(popts.outstr{comp},'\Delta','magnitude');
        yl1 = [l1,' ',unitstr(meta.outunit)];
        yl2 = labelstr_(popts.outstr{comp},'\Delta','angle');
    else
        % Real and Imaginary
        l1 = labelstr_(popts.outstr{comp},'\Delta','real');
        yl1 = [l1,' ',unitstr(meta.outunit)];
        l2 = labelstr_(popts.outstr{comp},'\Delta','imaginary');
        yl2 = [l2,' ',unitstr(meta.outunit)];
    end
end

function s = labelstr_(labelstr, prefix, ltype)

    if nargin < 2, prefix = ''; end
    if nargin < 3, ltype = 'magnitude'; end
    if contains(labelstr,'$')
        s = replace(labelstr,'$','');
    else
        s = ['\mbox{',labelstr,'}'];
    end
    s = ['\widetilde{',prefix,' ',s,'}'];
    if startsWith(ltype,'magnitude')
        s = ['$|',s,'|$'];
    end
    if startsWith(ltype,'angle')
        s = ['$\angle$ of $',s,'$'];
    end
    if startsWith(ltype,'real')
        s = ['Re of $',s,'$'];
    end
    if startsWith(ltype,'imag')
        s = ['Im of $',s,'$'];
    end
end

function argcheck_(S,popts)

    if size(S{1}.In,2) ~= size(S{1}.Out,2)
        error('Case of size(S.In,2) ~= size(S.Out,2) not handled');
    end

    tparts = split(popts.type,'-');
    tparts1 = {'error',...
               'original','final','detrended','windowed','whitened','zeropadded'};
    if ~any(strcmp(tparts{1},tparts1))
        list = join(tparts1,', ');
        error('popts.type must start with one of: %s (not "%s")',list{1},tparts{1});
    end
    if strcmp(tparts{1},'error')
        tparts3 = {'magphase','realimag'};
        if ~any(strcmp(tparts{3},tparts3))
            list = join(tparts3,', ');
            error('pots.type must end with one of %s (not "%s")',list{1},tparts{3});
        end
    else
        tparts3 = {'magnitudes','phases','reals','imaginaries'};
        if ~any(strcmp(tparts{3},tparts3))
            list = join(tparts3,', ');
            error('popts.type must end with one of %s (not "%s")',list{1},tparts{3});
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
            compstrs = popts.instr;
        else
            compstrs = popts.outstr;
        end
        for j = 1:jl
            if iscell(meta.instr)
                ls1{j+jl} = sprintf('%s(:,%d) Noise Amplitudes%s',...
                                        compstrs{j},j,unit);
            else
                ls1{j+jl} = sprintf('%s Noise Amplitudes%s',...
                                        compstrs,unit);
            end
        end
    end
    
end
