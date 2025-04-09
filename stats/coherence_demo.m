clear

% https://www.fieldtriptoolbox.org/example/coherence_snr/
Nr = 50;  % Number of repetitions

Ns = 59;   % Number of segments; must be odd (due to evalfreq)
n = 100;   % Length of segment
N = n*Ns;  % Length of full time series that is segmented
f = fftfreqp(N); % positive DFT frequencies for a signal of length N
t = (0:N-1)';

% Create input/output signals.

% Input and output frequency index
fidx = 32;
assert(fidx <= length(f),'fidx exceeds number of frequencies.')

Axo = 0; % Amplitude of input signal
Ayo = 0; % Amplitude of output signal
Anx = 1; % Amplitude of noise in input signal
Any = 0.5; % Amplitude of noise in output signal

for i = 1:Nr

    p = pi*(-1 + 2*rand(N,1)); 
    
    % Input
    xo = Axo*cos(2*pi*f(fidx)*t);
    % Output
    %yo = Ayo*cos(2*pi*f(i)*t + p + pi/2);
    %yo = Ayo*cos(2*pi*f(fidx)*t + p);
    yo = Ayo*cos(2*pi*f(fidx)*t);

    x = xo + Anx*randn(N,1);
    y = x + Any*randn(N,1);

    [Cxy{1}(:,i),fxy{1}] = mscohere(y,x,rectwin(N/Ns),0,N/Ns);
    fxy{1} = fxy{1}/(2*pi);

    [Cxy{2}(:,i),fxy{2}] = mscohere_(y,x,n,Ns);
    
    opts = tflab_options(0);
        opts.tflab.loglevel = 1;
        opts.fd.evalfreq.functionargs = {[Ns,(Ns-1)/2], 'linear'};
    [Cxy{3}(:,i),fxy{3}] = coherence(y,x,1,opts);

    % coherence() returns coherence, not magnitude squared coherence
    Cxy{3}(:,i) = Cxy{3}(:,i).^2; 
    
    %fprintf('%.2f %.2f %.2f\n', mean(Cxy{1}), mean(Cxy{2}), mean(Cxy{3}))
    
    for c = 1:3
        SNR{c}(:,i) = Cxy{c}(:,i)./(1-Cxy{c}(:,i));
    end

    % Alternative method of computing SNR
    SNRa(:,i) = signaltoerror(y, y-x, 1, opts);

end

Cxy{3} = Cxy{3}(2:end-1,:);
fxy{3} = fxy{3}(2:end-1);
SNR{3} = SNR{3}(2:end-1,:);
SNRa = SNRa(2:end-1,:);

mSNRa = mean(SNRa,2);

for c = 1:3
    mCxy{c} = mean(Cxy{c},2);
    sCxy{c} = std(Cxy{c},0,2)/sqrt(Nr);
    mSNR{c} = mean(SNR{c},2);
end

lso = {'MATLAB mscohere()','Manual mscohere()','TFLab coherence()'};
figure(1);figprep();clf
    plot(y);
    hold on;grid on;
    plot(x);
    legend({'$x$','$y$'})
    xlabel('$t$')
figure(2);clf;figprep()
    plot(fxy{1},mCxy{1},'b','linewidth',4);
    hold on;grid on;
    plot(fxy{2},mCxy{2},'g','LineWidth',2)
    plot(fxy{3},mCxy{3},'r','LineWidth',2)
    errorbar(fxy{3},mCxy{3},sCxy{3},sCxy{3},'Color','r')
    xlabel('$f$')
    templ = '$\\langle\\gamma^2_{xy}\\rangle_{N_r}$';
    ts = sprintf('$N_t = %d, N_r = %d, N_s = %d$', N, Nr, Ns);
    title(ts);
    ylabel(sprintf(templ, Nr));
    templ = '$\\overline{\\langle\\gamma^2_{xy}\\rangle}_{N_r}$';
    for i = 1:3
        ls{i} = sprintf('%s %s = %.3f',...
                        lso{i}, sprintf(templ, Nr), mean(Cxy{i}(:)));
    end
    legend(ls{:});
figure(3);clf;figprep()
    plot(fxy{1},mSNR{1},'b','linewidth',4);
    hold on;grid on;
    plot(fxy{2},mSNR{2},'g','LineWidth',2)
    plot(fxy{3},mSNR{3},'r','LineWidth',2)
    plot(fxy{3},mSNRa,'k','LineWidth',1)
    xlabel('$f$')
    ylabel('$\langle \mbox{SE} \rangle_{N_r}$')
    title('coherence-based SEs computed using $\langle \gamma^2_{xy}/(1-\gamma^2_{xy})\rangle_{N_r}$')
    legend(lso{:},'TFLab signaltoerror()');

function [Cxy, fxy] = mscohere_(y,x,n,Ns)
    % Manual version of mscohere with given arguments
    [dfty,~] = fftu(reshape(y,n,Ns));
    [dftx,fxy] = fftu(reshape(x,n,Ns));
    sxx = abs(sum(dfty.*conj(dfty),2));
    syy = abs(sum(dftx.*conj(dftx),2));
    sxy = abs(sum(dfty.*conj(dftx),2));
    Cxy = (sxy.^2)./(sxx.*syy);
end
