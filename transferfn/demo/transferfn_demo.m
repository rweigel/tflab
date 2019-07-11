% Compare two methods for estimating transfer function and impulse
% response between output u and input y.

% Time Domain Method: Directly estimate impulse response by solving for
% H = {h0, h1, ...} in
% ym(n)   = h0 + u(n-1)*h(1) + u(n-2)*h(2) + ...
% ym(n+1) = h0 + u(n-0)*h(1) + u(n-1)*h(2) + ...
% Transfer function is computed from DFT of H.

% Frequency domain method: Compute DFT coefficients U(f) and Y(f)
% for output u(t) and input y(t).
% Transfer function is H(f) = Y(f)/U(f).
% Impulse response is computed from inverse DFT of H(f).

% Note that in the following, only 50 coefficients are calculated in
% time domain.  To make a more direct comparison, the frequency domain
% estimates should be averaged so that the resolution in the frequency
% domain for both methods is the same.

clear;
addpath('~/git/m-rsw/time');
addpath('~/git/m-rsw/stats');
base = './figures/transfer_function_demo_';

N = 2000;

% Compute IRF for dx/dt + x/tau = delta(0), and IC of x_0 = 0 dx_0/dt = 0
% using forward Euler.
tau = 10;
dt = 1;
gamma = (1-dt/tau);
for i = 1:10*tau
    h(i,1) = gamma^(i-1);
    t(i,1) = dt*(i-1);
end

h = [0;h];

% Create an input signal
u = randn(N,1);

% Create an output that has a known relationship to input u via known
% impulse response function h:
% y(n) = 0*u(n) + h(1)*u(n-1) + h(2)*u(n-2) + ... h(N)*u(n-N)
y = filter(h,1,u);

% Add some noise
%y = y + 0.01*max(abs(y))*randn(N,1);

% Remove non-steady-state part of input and output
y  = y(length(h)+1:end);
u  = u(length(h)+1:end);

Na = 0;         % The number of acausal coefficients
Nc = length(h); % The number of causal coefficients
Ts = 0;         % Shift input with respect to output

% T is the output, X is the input.
% Each row of X contains time delayed values of u
[T,X] = time_delay(y,u,(Nc-1),Ts,Na,'pad');

% Compute model by solving for H = {h0,h1,...} in overdetermined set
% of equations
% ym(n)   = h0 + u(n-1)*h(1) + u(n-2)*h(2) + ...
% ym(n+1) = h0 + u(n-0)*h(1) + u(n-1)*h(2) + ...
% ...
% If any row above contains a NaN, it is omitted from the computation.
% Solve for Ym = X*H -> H = Ym\U 
LIN = basic_linear(X,T);

% Impulse response function using basic_linear() (last element is h0)
hbl = LIN.Weights(1:end-1);
hbl = [0;hbl]; % Zero is because Na = 0 -> h(t<=0) = 0.

% Predicted output using model coefficients H = LIN.Weights
% First length(LIN.Weights) values of yp are set to NaN
yp = basic_linear(X,LIN.Weights,'predict');

% Estimate transfer function for frequency domain method.
yaib = fft(y);
uaib = fft(u);

Rf   = yaib./uaib;
Nf   = length(Rf);
If   = abs(Rf);
phif = (180/pi)*atan2(imag(Rf),real(Rf));
ff   = [0:Nf/2-1]/length(Rf);
hf   = ifft(Rf);

% Estimate transfer function for time domain method.
Rt   = fft(hbl);
Nt   = length(Rt);
It   = abs(Rt);
phit = (180/pi)*atan2(imag(Rt),real(Rt));
ft   = [0:Nt/2]/Nt;
ht   = hbl;

figure(1);clf;
    %subplot('position',[0.1 0.55 0.8 0.38]);
        grid on;hold on;    
        plot(y(1)+10,'k','LineWidth',3);
        plot(u(1),'g','LineWidth',3);
        plot(u-5,'g');
        plot(y+5,'k');
        %th = title('$\dot{y}(t)+y(t)/\tau = u(t)$');
        %set(th,'Interpreter','Latex')

        lh = legend(' Output u',' Input y');
        set(lh,'FontSize',14);
        axis([0 length(u) -15 15]);
        set(gca,'FontSize',14);
        set(gca,'YTickLabel','');
        set(gca,'XTickLabel','');
if (0)
    subplot('position',[0.1 0.1 0.8 0.44])
        grid on;hold on;
        % Plot nothing with large LineWidth to make line in legend
        % larger than that in plot.
        plot(NaN,'k','LineWidth',2);
        plot(NaN,'b','LineWidth',2);
        plot(NaN,'y','LineWidth',2);

        plot(y,'k','LineWidth',2);
        plot(yp,'b');
        plot(y-yp,'y','LineWidth',2);

        lh = legend(' Output',' Predicted Output',' Error');
        xlabel('Time','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'YTickLabel','');
end
        if (writeimgs)
            print('-depsc',[base,'IRF_timeseries.eps']);
            print('-dpng','-r150',[base,'IRF_timeseries.png']);
        end

figure(2);clf;
    grid on;hold on;
    plot(NaN,'k','LineWidth',3);
    plot(NaN,'g','LineWidth',2);

    plot(exp(-(t(1:end))/tau)/exp(-t(1)/tau),'g','LineWidth',3);
    plot(h(2:end),'k','LineWidth',2);

    title('Impulse Response Comparison');
    lh = legend(' Continuous exact solution',...
                ' Discrete approximate solution',...
                'Location','North');
    xlabel('Time Since Impulse','FontSize',14);
    set(lh,'FontSize',14);
    set(lh,'Box','off');
    set(gca,'FontSize',14);
    axis([0 tau*5 0 1.1]);
    if (writeimgs)
        print('-depsc',[base,'IRF.eps']);
        print('-dpng','-r150',[base,'IRF.png']);
    end

figure(3);clf;grid on;hold on;
    plot(ff(2:end),If(2:length(ff)),'b','LineWidth',2,'Marker','.','MarkerSize',5);
    plot(ft(2:end),It(2:length(ft)),'g','LineWidth',2,'Marker','.','MarkerSize',5);
    xlabel('f');
    legend('Frequency Domain Method','Time Domain Method');
    title('Transfer Function Magnitude Estimate Comparison');
    set(gca,'FontSize',14);

    if (writeimgs)
        print('-depsc',[base,'Magnitude.eps']);
        print('-dpng','-r150',[base,'Magnitude.png']);
    end

figure(4);clf;grid on;hold on;
    plot(ff(2:end),phif(2:length(ff)),'b','LineWidth',2,'Marker','.','MarkerSize',5);
    plot(ft(2:end),phit(2:length(ft)),'g','LineWidth',2,'Marker','.','MarkerSize',5)
    xlabel('f');
    ylabel('\phi [degrees]')
    legend('Frequency Domain Method','Time Domain Method');
    set(gca,'YLim',[-180 180]);
    title('Transfer Function Phase Estimate Comparison');
    set(gca,'FontSize',14);
    if (writeimgs)
        print('-depsc',[base,'Phase.eps']);
        print('-dpng','-r150',[base,'Phase.png']);    
    end

figure(5);clf;grid on;hold on;
    plot(hf(2:end),'b','LineWidth',3,'Marker','.','MarkerSize',5);
    plot(ht(2:end),'g','LineWidth',2,'Marker','.','MarkerSize',5);
    legend('Frequency Domain Method','Time Domain Method');
    xlabel('Time Since Impulse');
    title('Impulse Response Estimate Comparison');
    set(gca,'FontSize',14);
    set(gca,'XLim',[0 100])
    if (writeimgs)
        print('-depsc',[base,'IRF2.eps']);
        print('-dpng','-r150',[base,'IRF2.png']);    
    end

% Visual inspection to determine phase and magnitude
t    = [0:length(y)/10-1]';
f    = 1/100; 
uin  = sin(2*pi*t*f);    
yout = filter(ht,1,uin);
% From time domain method phase plot, phase ~ -32
% Shift on this plot is 360*(150-159)/100 = ~32.4
% From time domain method transfer function mag. plot, ratio is ~8.62
% Ratio on this plot is ~max(yout)/max(uin) ~ 8.7  

figure(6);clf;hold on;grid on;
    plot(t,uin,'r','LineWidth',3,'Marker','.','MarkerSize',20);
    plot(t,yout,'b','LineWidth',3,'Marker','.','MarkerSize',20);
    legend('Output','Input');
    set(gca,'XLim',[149 160])

