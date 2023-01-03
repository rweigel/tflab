clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100;

%f = 25.25/N; % frequency that is not DFT freq for N or 2*N
%f = 25.5/N;  % frequency that is DFT freq for 2*N but not N
f= 25/N; % Exact DFT frequency

A = 1;
B = 0;
L = 2*N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = (0:N-1);
x = A*cos(2*pi*t*f) + B*sin(2*pi*t*f);
Xu = fft(x); % Unpadded
x = [x,zeros(1,L-N)];

X = fft(x);

fc = (0:L-1)/L;

more off
% Exact DFT for signal with length L
for l = 1:L
    Xc(l) = 0;
    for k = 1:N
        W = ((k-1)/N-(l-1)/L);
        G = (1/N)*exp(pi*1j*(N-1)*W)*sin(pi*N*W)/sin(pi*W);
        if W == 0
            % TODO: If L/N is integer, could use mod to check when W = 0.
            G = 1;
        end
        Xc(l) = Xc(l) + Xu(k)*G;
    end
end

fe = (0:L-1)/L;

PositionTop = [0.1300 0.5100 0.7750 0.38];
PositionBottom = [0.1300 0.1000 0.7750 0.38];

xrange = [0,1];
fname1 = 'zeropad-mag-phase.pdf';
fname2 = 'zeropad-real-imag.pdf';
for z = 1:2
    if z == 2
        fname1 = 'zeropad-mag-phase-zoom.pdf';
        fname2 = 'zeropad-real-imag-zoom.pdf';
        xrange = [0.2,0.3];
    end
    
    figure();figprep();
    subplot('Position', PositionTop);
        semilogy(f,(N/2)*sqrt(A^2 + B^2),'o');
        hold on;grid on;
        semilogy(fe,abs(X),'r*');
        semilogy(fc,abs(Xc),'k.')
        semilogy(1-f,(N/2)*sqrt(A^2 + B^2),'o');
        set(gca,'XLim',xrange);
        set(gca,'XTickLabel',[]);
        box on;
        title(sprintf('$N=%d; f=%.3f; X(f)=\\sum_{t=1}^N x(t)e^{-2\\pi j f(t-1)}$\n\n $x(t) = %d\\sin[2\\pi f (t-1)] + %d\\cos[2\\pi f(t-1)]; t=1,...,N/2;$   $x(t)=0; t=N/2+1,...,N$',2*N,f,A,B))
        ylabel('$|X(f)|$')
        legend({'CFT[$x(t)$]','FFT[$x(t)$]', 'DFT[$x(t)$]'});
    subplot('Position', PositionBottom);
        plot(f,(180/pi)*atan2(-B,A),'o');
        hold on;grid on;
        plot(fe,(180/pi)*angle(X),'r*');
        plot(1-f,(180/pi)*atan2(-B,A),'o');
        set(gca,'YTick',[-180:45:180]);
        set(gca,'XLim',xrange);
        box on;
        ylabel('Phase')
        xlabel('$\displaystyle f = \frac{k-1}{N},\quad k = 1, ...,N$')
        legend({'CFT[$x(t)$]','FFT[$x(t)$]', 'DFT[$x(t)$]'});

     figsave(fullfile(scriptdir(), fname1),'-p10');
        
    figure();figprep();
    subplot('Position', PositionTop);
        plot(f,(N/2)*A,'o');
        hold on;grid on;
        plot(fe,real(X),'r*');
        plot(1-f,(N/2)*A,'o');
        set(gca,'XTickLabel',[]);
        set(gca,'XLim',xrange);
        box on
        ylabel('Re$[X(f)]$')
        legend({'CFT[$x(t)$]','FFT[$x(t)$]', 'DFT[$x(t)$]'});
    subplot('Position', PositionBottom);
        plot(f,-(N/2)*B,'o');
        hold on;grid on;
        plot(fe,imag(X),'r*');
        plot(1-f,-(N/2)*B,'o');
        box on;
        set(gca,'XLim',xrange);
        ylabel('Im$[X(f)]$')
        xlabel('$f$')
        legend({'CFT[$x(t)$]','FFT[$x(t)$]', 'DFT[$x(t)$]'});

     figsave(fullfile(scriptdir(), fname2),'-p10');

end