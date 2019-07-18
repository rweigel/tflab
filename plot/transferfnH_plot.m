function transferfnH_plot(S1,S2,xl)

figure;
figprep();

if nargin < 2
    S2 = [];
end
if nargin == 2 && ~isstruct(S2)
    xl = S2;
end

Hstrs = {'H_{xx}','H_{xy}','H_{yx}','H_{yy}'};

if isstruct(S1) && ~isstruct(S2)
    H1 = fftshift(S1.H);
    [a,b] = findss(H1);
    tH1 = fftshift(S1.tH);
    
    plot(tH1(a:b), H1(a:b));
    grid on;box on;hold on;

    unitstr = sprintf('[%s/%s]', S1.Options.info.inunit, S1.Options.info.outunit);
    title(sprintf(S1.Options.description),'FontWeight','Normal');
    legend(sprintf('$H$ %s',unitstr), 'Location', 'NorthEast');
    xlabel(sprintf('$t$ [%s]', S1.Options.info.timeunit));
    if nargin > 2
        set(gca(),'XLim',xl);
    end
end

if nargin > 1 && isstruct(S2)
    % Assumes units are the same. TODO: Allow different unit labels.

    H1 = fftshift(S1.H);
    H2 = fftshift(S2.H);
    [a1,b1] = findss(H1);
    [a2,b2] = findss(H2);

    tH1 = fftshift(S1.tH);
    tH2 = fftshift(S2.tH);

    if b1-a1 > 1
        plot(tH1(a1:b1),H1(a1:b1));
    else
        stem(tH1(a1:b1),H1(a1:b1));    
    end
    grid on;box on;hold on;

    if b1-a1 > 1    
        plot(tH2(a2:b2),H2(a2:b2));
    else
        stem(tH2(a2:b2),H2(a2:b2));
    end

    ylabel(sprintf('[%s/%s]', S1.Options.info.inunit, S1.Options.info.outunit));
    legend(['$H$ $\,$ ', S1.Options.description],...
           ['$H$ $\,$ ', S2.Options.description],...
           'Location','NorthEast');
    xlabel(sprintf('$t$ [%s]', S1.Options.info.timeunit));
    if nargin > 2
        set(gca(),'XLim',xl);
    end
end

end

function [a,b] = findss(x, w, t)
    % Find limits above/below which x is "steady state".
    if nargin < 2
        w = 10;
    end
    if nargin < 3
        t = 0.05;
    end
    if length(x) < w
        a = 1;
        b = length(x);
        return;
    end

    % Compute std in non-overlaping windows of length w.
    x = x(1:w*floor(length(x)/w));
    xr = reshape(x,w,round(length(x)/w));
    xs = std(xr);

    % ss start is start of window where this condition true.
    a = w*find(xs/max(xs) > t,1,'first')-w+1;

    % ss end is end of window where this condition is true.
    b = w*find(xs/max(xs) > t,1,'last');
    
    if isempty(a) || isempty(b)
        a = 1;
        b = length(x);
    end
end
