clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100;

t = (0:N-1);
f = 25.5/N;
% f= 25/N; % Exact DFT frequency

% Note how asymmetry around f occurs in PSD when one of A or B is not zero.
A = 1;
B = 1;
% B = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = A*cos(2*pi*t*f) + B*sin(2*pi*t*f);

% Compute DFT using using slow sum method instead of FFT
for k = 1:N
    Xn = 0;
    for n = 1:N
       Xn = Xn + x(n)*exp(-1j*2*pi*(k-1)*(n-1)/N);
    end
    X(k) = Xn;
end

Ng = 10001;
fg = (0:Ng-1)/Ng;

% Exact DFT for signal with arbitrary frequency f
% For derivation, see slide 14 of
% https://www.eecs.umich.edu/courses/eecs452/Lec/L11SlidesF14.pdf
% or copy named leakage_equation.pdf.
for k = 1:Ng
    Fn = f - (k-1)/Ng;
    Fp = f + (k-1)/Ng;
    Lp(k) = exp(-1j*pi*Fp*(N-1))*sin(pi*Fp*N)/sin(pi*Fp);    
    Ln(k) = exp(1j*pi*Fn*(N-1))*sin(pi*Fn*N)/sin(pi*Fn);    
end
Lx = A*(Ln+Lp)/2 + B*(Ln-Lp)/(2*1j);

clear Lp Ln
% Manual calculation. Should give same result as exact equation used
% above, which is derived by using the formula for a partial sum of a
% geometric series.
for k = 1:Ng
    Ln(k) = 0;
    Lp(k) = 0;
    Fn = f - (k-1)/Ng;
    Fp = f + (k-1)/Ng;
    for n = 1:N
        Ln(k) = Ln(k) + exp(+2*pi*1j*(n-1)*Fn);
        Lp(k) = Lp(k) + exp(-2*pi*1j*(n-1)*Fp);
    end
end
L = A*(Ln+Lp)/2 + B*(Ln-Lp)/(2*1j);

fe = (0:N-1)/N;

PositionTop = [0.1300 0.5100 0.7750 0.38];
PositionBottom = [0.1300 0.1000 0.7750 0.38];

xrange = [0,1];
fname1 = 'leakage-mag-phase.pdf';
fname2 = 'leakage-real-imag.pdf';
for z = 1:2
    if z == 2
        fname1 = 'leakage-mag-phase-zoom.pdf';
        fname2 = 'leakage-real-imag-zoom.pdf';
        xrange = [0.2,0.3];
    end
    
    figure();figprep();
    subplot('Position', PositionTop);
        semilogy(f,(N/2)*sqrt(A^2 + B^2),'o');
        hold on;grid on;
        semilogy(fe,abs(X),'r*');
        %semilogy(fg,abs(L),'k-')
        semilogy(fg,abs(Lx),'b-')
        semilogy(f,(N/2)*sqrt(A^2 + B^2),'o');
        set(gca,'XLim',xrange);
        set(gca,'XTickLabel',[]);
        box on;
        title(sprintf('$N=%d; f=%.3f; X(f)=\\sum_{t=1}^N x(t)e^{-2\\pi j f(t-1)}$\n\n $x(t) = %d\\sin[2\\pi f (t-1)] + %d\\cos[2\\pi f(t-1)]; t=1,...,N$',N,f,A,B))
        ylabel('$|X(f)|$')
        %legend({'CFT[$x(t)$]','FFT[$x(t)$]', 'DFT[$x(t)$] using sum', 'DFT[$x(t)$] w/o using sum'});
        legend({'CFT[$x(t)$]','FFT[$x(t)$]', 'DFT[$x(t)$]'});
    subplot('Position', PositionBottom);
        plot(f,(180/pi)*atan2(-B,A),'o');
        hold on;grid on;
        plot(fe,(180/pi)*angle(X),'r*');
        %plot(fg,(180/pi)*angle(L),'k-')
        plot(fg,(180/pi)*angle(Lx),'b-')
        plot(1-f,(180/pi)*atan2(-B,A),'o');
        set(gca,'YTick',[-180:45:180]);
        set(gca,'XLim',xrange);
        box on;
        ylabel('Phase')
        xlabel('$\displaystyle f = \frac{k-1}{N},\quad k = 1, ...,N$')
        %legend({'CFT[$x(t)$]','FFT[$x(t)$]', 'DFT[$x(t)$] using sum', 'DFT[$x(t)$] w/o using sum'});
        legend({'CFT[$x(t)$]','FFT[$x(t)$]', 'DFT[$x(t)$]'});
    
    figsave(fullfile(scriptdir(), fname1),'-p10');
 
    figure();figprep();
    subplot('Position', PositionTop);
        plot(f,(N/2)*A,'o');
        hold on;grid on;
        plot(fe,real(X),'r*');
        %plot(fg,real(L),'k-')
        plot(fg,real(Lx),'b-')
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
        %plot(fg,imag(L),'k-')
        plot(fg,imag(Lx),'b-');
        plot(f,-(N/2)*B,'o');
        box on;
        set(gca,'XLim',xrange);
        ylabel('Im$[X(f)]$')
        xlabel('$f$')
        legend({'CFT[$x(t)$]','FFT[$x(t)$]', 'DFT[$x(t)$]'});

    figsave(fullfile(scriptdir(), fname2),'-p10');
        
end