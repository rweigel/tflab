clear

addpath('../../time/')
addpath('../../stats/')
addpath('../misc/')

writeimgs = 0;

tau  = 10;  % Filter decay constant
Ntau = 10;  % Number of filter coefficients = Ntau*tau + 1
N    = 2e4; % Simulation length
Nc   = 51;  
df   = 50;  % Width of rectangualar window
Nss  = 4;   % Ntau*tau*Nss
nb   = 0.0; % Noise in B
ne   = 0.0; % Noise in E
ndb  = 0.0; % Noise in dB

paramstring = sprintf('_ne_%.1f',ne);

% IRF for dx/dt + x/tau = delta(0), and ICs
% x_0 = 0 dx_0/dt = 0 approximated using forward Euler.
dt = 1;
gamma = (1-dt/tau);
for i = 1:Ntau*tau
    h(i,1) = gamma^(i-1);
    t(i,1) = dt*(i-1);
end
hstr = sprintf('(1-1/%d)^{t}; t=1 ... %d; h_{xy}(0)=0;', tau, length(h));
% Add a zero because MATLAB filter function requires it.  Also,
% not having it causes non-physical phase drift with frequency.
h  = [0;h];
th = [0:length(h)-1]';

% Exact transfer function
Zxy = fft(h);
Nh  = length(Zxy);
Zxy = Zxy(1:floor(Nh/2)+1);
fh  = [0:floor(Nh/2)]'/Nh;

% Exact transfer Function Phase
Pxy = (180/pi)*atan2(imag(Zxy),real(Zxy));

% Add extra values to get nice length
% (because we cut off non-steady state).
N = N + Nss*length(h);

% Noise
NE  = [ne*randn(N,1),ne*randn(N,1)];
NB  = [nb*randn(N,1),nb*randn(N,1)];
NdB = [ndb*randn(N,1),ndb*randn(N,1)];

% Create signals
B(:,1) = randn(N,1);
B(:,2) = randn(N,1);

E(:,2) = filter(h,1,B(:,1)+NB(:,1)) + NE(:,2);
E(:,1) = filter(h,1,B(:,2)+NB(:,2)) + NE(:,1);

for i = 1:size(B,2)
    B(:,i)  = B(:,i) - mean(B(:,i));
end
for i = 1:size(E,2)
    E(:,i) = E(:,i) - mean(E(:,i));
end

dB = diff(B);
dB = [dB;dB(end,:)];

% Remove non-steady-state part of signals
B  = B(Nss*length(h)+1:end,:);
E  = E(Nss*length(h)+1:end,:);
dB = dB(Nss*length(h)+1:end,:);

NE  = NE(Nss*length(h)+1:end,:);
NB  = NB(Nss*length(h)+1:end,:);
NdB = NdB(Nss*length(h)+1:end,:);

N = size(B,1);

if (0) % Window in time domain.
    for i = 1:size(B,2)
	Bw(:,i) = B(:,i).*parzenwin(length(B));
    end
    for i = 1:size(E,2)
	Ew(:,i) = E(:,i).*parzenwin(length(E));
    end
else
    Bw = B;
    Ew = E;
end

if (0)
    for i = 1:size(B,2)
	B(:,i)  = B(:,i)  - mean(B(:,i));
	dB(:,i) = dB(:,i) - mean(dB(:,i));
    end
    for i = 1:size(E,2)
	E(:,i) = E(:,i) - mean(E(:,i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw periodograms

N = size(B,1);
f = [0:N/2]'/N;

X    = [B,dB,E,NB,NdB,NE];
ftX  = fft(X);
pX   = ftX;

ftB   = pX(:,1:2);
ftdB  = pX(:,3:4);
ftE   = pX(:,5:6);
ftNB  = pX(:,7:8);
ftNdB = pX(:,9:10);
ftNE  = pX(:,11:12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time domain
[ZBL,feBL,HBL,tBL,EpBL] = transferfnTD(B,E,Nc);

if (0)
    HBL2   = Z2H(feBL,ZBL);
    ZBLi   = Zinterp(feBL,ZBL,fh);
    EpBL   = Hpredict(tBL,HBL,B);
    EpBL2  = Zpredict(feBL,ZBL,B);
    ZxyBLi = ZBLi(:,2);
    EyBL   = EpBL(:,2);
    EyBL2  = EpBL2(:,2);

    EyBL(1:Nc-1) = NaN;
    peBL = pe_nonflag(E(:,2),EyBL);
end

EpBL2  = Zpredict(feBL,ZBL,B);
ZxyBL  = ZBL(:,2);
hBL    = HBL(:,2);

% Transfer Function Phase
PxyBL  = (180/pi)*atan2(imag(ZBL(:,2)),real(ZBL(:,2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency Domain Rectangular
[ZR,feR,HR,tR,EpR] = transferfnFD(Bw,Ew,1,'rectangular',df);

if (0)
    ZRi      = Zinterp(feR,ZR,fh);
    HR       = Z2H(feR,ZR);
    HR       = fftshift(HR,1);
    NR       = (size(HR,1)-1)/2;
    tR       = [-NR:NR]';
    EpR      = Hpredict(tR,HR,B);
    EpR2     = Zpredict(feR,ZR,B);

    ZxyR  = ZR(:,2);
    ZxyRi = ZRi(:,2);
    hR    = HR(:,2);
    EyR   = EpR(:,2);
    EyR2  = EpR2(:,2);
end
EpR2 = Zpredict(feR,ZR,B);
ZxyR = ZR(:,2);
hR   = HR(:,2);

% Transfer Function Phase
PxyR  = (180/pi)*atan2(imag(ZxyR),real(ZxyR));
%PxyRi = (180/pi)*atan2(imag(ZxyRi),real(ZxyRi));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency Domain Parzen
[ZP,feP] = transferfnFD(B,E,1,'parzen');
ZPi      = Zinterp(feP,ZP,fh);
HP       = Z2H(feP,ZP,fh);
HP       = fftshift(HP,1);
NP       = (size(HP,1)-1)/2;
tP       = [-NP:NP]';
EpP      = Hpredict(tP,HP,B);
EpP2     = Zpredict(feP,ZP,B);

ZxyP  = ZP(:,2);
ZxyPi = ZPi(:,2);
hP    = HP(:,2);
EyP   = EpP(:,2);
EyP2  = EpP2(:,2);

% Transfer Function Phase
PxyP  = (180/pi)*atan2(imag(ZxyP),real(ZxyP));
PxyPi = (180/pi)*atan2(imag(ZxyPi),real(ZxyPi));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')

figure(1);clf;hold on; grid on;
    plot(NaN*B(1:2,1),'r','LineWidth',3)
    plot(NaN*E(1:2,2),'Color',[1,0.5,0],'LineWidth',3)
    plot(NaN*E(1:2,2),'k','LineWidth',3)
    plot(NaN*E(1:2,2),'Color',[0.5,0.5,0.5],'LineWidth',3)

    ts = sprintf('E_y = filter(h_{xy},1,B_x+\\deltaBx) + \\deltaE_y; \\deltaE_y = \\eta(0,%.2f); \\deltaB_x = \\eta(0,%.2f)',ne,nb);
    title(ts);
    plot(B(:,1)+15,'r')
    plot(NB(:,1)+5,'Color',[1,0.5,0])
    plot(E(:,2)-5,'k')
    plot(NE(:,2)-15,'Color',[0.5,0.5,0.5])
    xlabel('t (sample number-1)')
    set(gca,'YLim',[-30 30])
    legend('B_x+15 (input)','\deltaB_x+5 (noise)','E_y-5 (output)','\deltaE_y-15 (noise)')
    plotcmds(['timeseries',paramstring],writeimgs)

figure(2);clf;
    loglog(NaN*ftB(1:2,1),'r','LineWidth',3)
    hold on;grid on;

    loglog(NaN*ftE(1:2,2),'Color',[1,0.5,0],'LineWidth',3)
    loglog(NaN*ftNB(1:2,1),'k','LineWidth',3)
    loglog(NaN*ftNE(1:2,2),'Color',[0.5,0.5,0.5],'LineWidth',3)
    %ts = sprintf('E_y = filter(h_{xy},1,B_x+\\deltaBx) + \\deltaE_y; \\deltaE_y = \\eta(0,%.2f); \\deltaB_x = \\eta(0,%.2f)',ne,nb);
    %title(ts);
    title('Raw Periodograms')
    loglog(f(2:end),abs(ftB(2:N/2+1,1)),'r')
    loglog(f(2:end),abs(ftNB(2:N/2+1,1)),'Color',[1,0.5,0])
    loglog(f(2:end),abs(ftE(2:N/2+1,2)),'k')
    loglog(f(2:end),abs(ftNE(2:N/2+1,2)),'Color',[0.5,0.5,0.5])
    xlabel('f')
    legend('B_x (input)','\deltaB_x (noise)','E_y (output)','\deltaE_y (noise)')
    plotcmds(['rawperiodograms',paramstring],writeimgs)

figure(3);clf;grid on;
    me = mean(E(:,2));
    mb = mean(B(:,1));
    xc = xcorr(E(:,2)-me,B(:,1)-mb,'unbiased');
    tl = [-N+1:N-1];
    xc = fftshift(xc);
    plot(tl,xc,'Color','r','Marker','.','MarkerEdgeColor','k');grid on;
    set(gca,'XLim',[-3*length(h) 3*length(h)]);
    title('Raw Cross Correlation')
    xlabel('lag')
    legend('E_y,B_x')
    plotcmds(['crosscorrelation',paramstring],writeimgs)

figure(4);clf;
    hold on;grid on;
    plot(th,h,'k','Marker','+','MarkerSize',5,'LineWidth',5)
    plot(tBL,hBL,'m','LineWidth',4)
    plot(tR,hR,'b','LineWidth',2)
    plot(tP,hP,'g','LineWidth',2)
    xlabel('t')
    title('Impulse Responses')
    legend( sprintf('h_{xy} = %s',hstr),...
            sprintf('h_{xy} time domain; n_T = %d',length(hBL)),...
            sprintf('h_{xy} freq. domain Rectangular; n_R = %d', df),...
            sprintf('h_{xy} freq. domain Parzen; n_P = %d',length(feP))...
           )
   plotcmds(['impulse_responses',paramstring],writeimgs)

figure(5);clf;
    % Create padded impulse responses.
    tp = [min([tP(1),tR(1),th(1)]):max([tP(end),tR(end),th(end)])]';

    hp = interp1(th,h,tp);
    hp(isnan(hp)) = 0;
    for i = 1:size(HBL,2)
        hBLp(:,i) = interp1(tBL,HBL(:,i),tp);
        hRp(:,i)  = interp1(tR,HR(:,i),tp);
        hPp(:,i)  = interp1(tP,HP(:,i),tp);
    end

    hold on;grid on;
    plot(tp,hBLp(:,2)-hp,'m','LineWidth',4)
    plot(tp,hRp(:,2)-hp,'b','LineWidth',2)
    plot(tp,hPp(:,2)-hp,'g','LineWidth',2)
    xlabel('t')
    title('Impulse Response Errors')
    legend(...
            sprintf('\\deltah_{xy} time domain; n_T = %d',length(hBL)),...
            sprintf('\\deltah_{xy} freq. domain Rectangular; n_f = %d', df),...
            sprintf('\\deltah_{xy} freq. domain Parzen; n_P = %d',length(feP))...
            )
   plotcmds(['impulse_response_errors',paramstring],writeimgs)

figure(6);clf;
    hold on;grid on;
    plot(fh,abs(Zxy),'k','Marker','+','MarkerSize',10,'LineWidth',5)
    plot(feBL,abs(ZxyBL),'m','Marker','.','MarkerSize',25,'LineWidth',3);
    plot(feR,abs(ZxyR),'b','Marker','.','MarkerSize',15,'LineWidth',2);
    plot(feP,abs(ZxyP),'g','Marker','.','MarkerSize',10,'LineWidth',1);
    xlabel('f')
    title('Transfer Functions')
    legend(...
            'Z_{xy}',...
            sprintf('Z_{xy} time domain; n_T = %d',length(hBL)),...
            sprintf('Z_{xy} freq. domain rectangular (n_R = %d)', df),...
            sprintf('Z_{xy} freq. domain parzen; n_P = %d',length(feP))...
            )
   plotcmds(['transfer_functions',paramstring],writeimgs)

if (0)   
figure(7);clf;
    hold on;grid on;
    plot(fh,abs(Zxy),'k','Marker','+','MarkerSize',10,'LineWidth',5)
    plot(fh,abs(ZxyBLi),'m','Marker','.','MarkerSize',25,'LineWidth',3);
    plot(fh,abs(ZxyRi),'b','Marker','.','MarkerSize',15,'LineWidth',2);
    plot(fh,abs(ZxyPi),'g','Marker','.','MarkerSize',10,'LineWidth',1);
    xlabel('f')
    title('Interpolated Transfer Functions')
    legend(...
            'Z_{xy}',...
            sprintf('Z_{xy} time domain; n_T = %d',length(hBL)),...
            sprintf('Z_{xy} freq. domain rectangular (n_R = %d)', df),...
            sprintf('Z_{xy} freq. domain parzen; n_P = %d',length(feP))...
            )
   plotcmds(['transfer_functions_interpolated',paramstring],writeimgs)

figure(8);clf;
    hold on;grid on;
    plot(fh,abs(ZxyBLi)-abs(Zxy),'m','Marker','.','MarkerSize',20,'LineWidth',2);
    plot(fh,abs(ZxyRi)-abs(Zxy),'b','Marker','.','MarkerSize',20,'LineWidth',2);
    plot(fh,abs(ZxyPi)-abs(Zxy),'g','Marker','.','MarkerSize',20,'LineWidth',2);
    xlabel('f')
    title('Transfer Function Error')
    legend(...
            sprintf('\\deltaZ_{xy} time domain; n_T = %d',length(hBL)),...
            sprintf('\\deltaZ_{xy} freq. domain Rectangular; n_R = %d', df),...
            sprintf('\\deltaZ_{xy} freq. domain Parzen; n_P = %d',length(feP)),...
            'Location','SouthEast')
   plotcmds(['transfer_function_errors',paramstring],writeimgs)
end
figure(9);clf;
    hold on;grid on;
    plot(fh,abs(Pxy),'k','Marker','+','MarkerSize',10,'LineWidth',5)
    plot(feBL,abs(PxyBL),'m','Marker','.','MarkerSize',25,'LineWidth',3);
    plot(feR,abs(PxyR),'b','Marker','.','MarkerSize',25,'LineWidth',2);
    plot(feP,abs(PxyP),'g','Marker','.','MarkerSize',25,'LineWidth',1);
    xlabel('f');
    ylabel('degrees')
    title('Transfer Function Phases')
    legend(...
        '\phi_{xy}',...
        sprintf('\\phi_{xy} time domain; n_T = %d',length(hBL)),...
        sprintf('\\phi_{xy} freq. domain Rectangular; n_R = %d', df),...
        sprintf('\\phi_{xy} freq. domain Parzen; n_P = %d',length(feP)),...
        'Location','SouthEast'...
    )
   plotcmds(['transfer_function_phases',paramstring],writeimgs)

if (0)
figure(10);clf;
    hold on;grid on;
    plot(fh,abs(PxyBLi)-abs(Pxy),'m','Marker','.','MarkerSize',20,'LineWidth',2);
    plot(fh,abs(PxyRi)-abs(Pxy),'b','Marker','.','MarkerSize',20,'LineWidth',2);
    plot(fh,abs(PxyPi)-abs(Pxy),'g','Marker','.','MarkerSize',20,'LineWidth',2);
    xlabel('f')
    ylabel('degrees')
    title('Transfer Function Phase Errors')
    legend(...
            sprintf('\\delta\\phi_{xy} time domain; n_T = %d',length(hBL)),...
            sprintf('\\delta\\phi_{xy} freq. domain Rectangular; n_R = %d', df),...
            sprintf('\\delta\\phi_{xy} freq. domain Parzen; n_P = %d',length(feP)),...
            'Location','NorthEast')
   plotcmds(['transfer_function_phase_errors',paramstring],writeimgs)
end

figure(11);clf;
    hold on;grid on;
    plot(E(:,2),'k','LineWidth',3)
    plot(EpBL(:,2),'m')
    plot(EpR(:,2),'b')
    plot(EyP,'g')
    xlabel('t')
    title('Predictions (using H)')
    legend('E_y',...
            'E_y time domain',...
            'E_y freq. domain Rectangular',...
            'E_y freq. domain Parzen'...
            )
   plotcmds(['predictions_H',paramstring],writeimgs)

figure(12);clf;
    hold on;grid on;
    plot(E(:,2),'k','LineWidth',3)
    plot(EpBL2(:,2),'m')
    plot(EpR2(:,2),'b')
    plot(EyP2,'g')
    xlabel('t')
    title('Predictions (using Z)')
    legend('E_y',...
            'E_y time domain',...
            'E_y freq. domain Rectangular',...
            'E_y freq. domain Parzen'...
            )
   plotcmds(['predictions_Z',paramstring],writeimgs)

figure(13);clf;
    hold on;grid on;
    plot(E(:,2)-EpBL(:,2)+10,'m')
    plot(E(:,2)-EpR(:,2),'b')
    plot(E(:,2)-EyP-10,'g')
    peBL = pe_nonflag(E(:,2),EpBL(:,2));
    peR  = pe_nonflag(E(:,2),EpR(:,2));
    peP  = pe_nonflag(E(:,2),EyP);

    xlabel('t')
    title('Prediction Errors (using H)')
    set(gca,'YLim',[-20 20])
    legend(...
        sprintf('\\DeltaE_y+10 time domain; PE = %.3f',peBL),...        
        sprintf('\\DeltaE_y freq. domain Rectangular; PE = %.3f',peR),...
        sprintf('\\DeltaE_y-10 freq. domain Parzen; PE = %.3f',peP)...
        );
    plotcmds(['prediction_H_errors',paramstring],writeimgs)

figure(14);clf;
    hold on;grid on;
    plot(E(:,2)-EpBL2(:,2)+10,'m')
    plot(E(:,2)-EpR2(:,2),'b')
    plot(E(:,2)-EyP2-10,'g')
    peBL = pe_nonflag(E(:,2),EpBL2(:,2));
    peR = pe_nonflag(E(:,2),EpR2(:,2));
    peP = pe_nonflag(E(:,2),EyP2);
    xlabel('t')
    title('Prediction Errors (using Z)')
    set(gca,'YLim',[-20 20])
    legend(...
        sprintf('\\DeltaE_y+10 time domain; PE = %.3f',peBL),...        
        sprintf('\\DeltaE_y freq. domain Rectangular; PE = %.3f',peR),...
        sprintf('\\DeltaE_y-10 freq. domain Parzen; PE = %.3f',peP)...
        );
    plotcmds(['prediction_Z_errors',paramstring],writeimgs)
