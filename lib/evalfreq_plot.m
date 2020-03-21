function evalfreq_plot(N,varargin)
%EVALFREQ_PLOT
%
%   evalfreq_plot(...) takes the same inputs as EVALFREQ() and creates a
%   plot showing the frequency bands and centers.
%
%   See also EVALFREQ, EVALFREQ_LOG.

logscale = 0;
period = 1;

if nargin > 1
    [fe,Ic,Ne] = evalfreq(N,varargin{:});
else
    [fe,Ic,Ne] = evalfreq(N);
end

% Unique DFT frequencies
[~,f] = fftfreq(N);

if isempty(varargin) || (length(varargin) > 1 && strcmp('logarithmic',varargin{2}))
    logscale = 1;
end

figure;clf;grid on;box on;hold on;
for i = 1:length(fe)
    Np = (2*Ne(i)+1);
    fs = f(Ic(i)-Ne(i):Ic(i)+Ne(i));
    fa{i} = fs;
    if period
        x = 1./fs;
        xe = 1./fe;
    else
        x = fs;
        xe = fe;
    end
    if logscale
        % Shift duplicates a small amount vertically.
        if i > 1 && Ne(i) == Ne(i-1)
            Np = Np + i/10;
        end
        y = Np*ones(size(fs));
        plot(x,y,'go');
        plot(x,y,'k-');
        plot(xe(i),Np,'k.','MarkerSize',20);
    else
        y = i*ones(size(fs));        
        plot(x,y,'go');
        plot(x,y,'k-');        
        plot(xe(i),i,'k.','MarkerSize',20);
    end
    if i == 1
        % Hack for legend.
        plot(nan,nan,'r.');
    end
end
d = setdiff(f,cat(2,fa{:}));
if period
    plot(1./d,length(d)*ones(size(d)),'r.');
    legend('DFT frequencies','Band','Center Frequency',...
           'Un-used DFT frequencies','Location','NorthEast');
    xlabel('period');
else
    plot(d,length(d)*ones(size(d)),'r.');
    legend('DFT frequencies','Band','Center Frequency',...
           'Un-used DFT frequencies','Location','SouthEast');    
    xlabel('frequency');
end

if logscale
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    ylabel(sprintf('# of DFT points associated with band\n(overlapping bands shifted vertically)'));
else
    set(gca,'YTickLabel',[])
    set(gca,'YTick',[])
end