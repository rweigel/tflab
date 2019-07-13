clear;

close all;
set(0,'defaultFigureWindowStyle','docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regression test. Compare OLS_REGRESS() with ROBUSTFIT() when no noise.

N = [99,100];
for n = N
    B = randn(n,1);
    E = B;

    opts = transferfnFD_options(1);
    opts.fd.evalfreq.function = @evalfreq;
    opts.fd.evalfreq.functionstr  = '1 DFT point per window';
    opts.fd.evalfreq.functionargs = {[1,1],'linear'};
    
    %opts.fd.regression.function = @robust_robustfit;
    %opts.fd.regression.functionstr = 'Robust regression using robustfit() function';
    %opts.fd.regression.functionargs = {[],[],'off'};
    %S2 = transferfnFD(B,E,opts);

    opts.fd.regression.function = @robust_v1;
    opts.fd.regression.functionstr = 'Robust regression using robustfit() function';
    %opts.fd.regression.functionargs = {};
    S2 = transferfnFD(B,E,opts);

    % Uses default regression function OLS_REGRESS().
    S1 = transferfnFD(B,E,opts);
    
    assert(S1.Metrics.PE - S2.Metrics.PE < eps);
    assert(S1.Metrics.CC - S2.Metrics.CC < eps);
    assert(S1.Metrics.MSE - S2.Metrics.MSE < eps);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation test
%
% B = randn(), E = B. With evalfreqs = DFT frequencies, should produce
% perfect predictions b/c # of free parameters in fitted Z equals number of
% data points.

N = [99,100];
for n = N
    B = randn(n,1);
    E = B;

    opts = transferfnFD_options(1);
    opts.fd.evalfreq.functionstr  = '1 DFT point per window';
    opts.fd.evalfreq.functionargs = {[1,0],'linear'};
    S = transferfnFD(B,E,opts);

    % Somewhat arbitrary threshold as we don't have an exact expectation
    % value for what it should be.
    assert(1-S.Metrics.PE < 10*eps);
    assert(S.Metrics.MSE < 10*eps);
    assert(1-S.Metrics.CC < 10*eps);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation Test
% B = cos(w*t), E = A(w)*cos(w*t + phi(w)). No leakage

clear E B
N = 101;
f = fftfreqp(N);
t = (0:N-1)';
for i = 1:length(f)
    A(i)   = (i-1);
    Phi(i) = -2*pi*f(i);
    B(:,i) = cos(2*pi*f(i)*t);
    E(:,i) = A(i)*cos(2*pi*f(i)*t+Phi(i));
end

B = sum(B,2);
E = sum(E,2);

opts = transferfnFD_options(0);
S2 = transferfnFD(B,E,opts);

if 0
    figure(1);clf
        plot(S2.fe,abs(S2.Z),'.');
        grid on;
    figure(2);clf
        plot(S2.fe,(180/pi)*S2.Phi,'.');
        grid on;
end

assert(max(abs(abs(S2.Z)-A')) < 1e-12);
assert(max(S2.Phi-Phi') < 1e-12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% API Test - Multiple Ouputs

N = 1000;
B = randn(N,2);
E = B;

opts = transferfnFD_options(1);
opts.td.window.width = N;
opts.td.window.shift = N;

%%% 1 input, one or two outputs
S1 = transferfnFD(B(:,1),E(:,1),opts);
S2 = transferfnFD(B(:,1),[E(:,1),E(:,1)],opts);

assert(all(S1.Predicted == S2.Predicted(:,1)));
assert(all(S1.Predicted == S2.Predicted(:,2)));

%%% 2 inputs, one or two outputs
S1 = transferfnFD(B,E(:,1),opts);
S2 = transferfnFD(B,E(:,2),opts);
S3 = transferfnFD(B,E,opts);

assert(all(S1.Predicted == S3.Predicted(:,1)));
assert(all(S2.Predicted == S3.Predicted(:,2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% API Test - Segmenting
% E and B are split into segments and transfer functions are computed for
% each segment.

N = 1000;
B = randn(N,1);
E = B;

opts = transferfnFD_options(1);
opts.td.window.width = N;
opts.td.window.shift = N;

S1 = transferfnFD(B,E,opts);
S2 = transferfnFD([B;B],[E;E],opts);

% Results for two segments in S2 should be identical to single
% segment in S1.
assert(all(S1.Predicted == S2.Segment.Predicted(:,1,1)));
assert(all(S1.Predicted == S2.Segment.Predicted(:,1,2)));
assert(all(S1.Z == S2.Segment.Z(:,:,1)))
assert(all(S1.Z == S2.Segment.Z(:,:,2)))
assert(all(S1.Z(:) == S2.Z(:)));

%%%

N = 1000;
B = randn(N,2);
E = B;

opts = transferfnFD_options(1);
opts.td.window.width = N;
opts.td.window.shift = N;

S1 = transferfnFD(B(:,1),E(:,1),opts);
S2 = transferfnFD([B(:,1);B(:,1)],[E;E],opts);
assert(all(S1.Predicted == S2.Segment.Predicted(:,1,1)));
assert(all(S1.Predicted == S2.Segment.Predicted(:,1,2)));

S3 = transferfnFD(B,E,opts);
S4 = transferfnFD([B;B],[E;E],opts);
assert(all(S3.Predicted(:,1) == S4.Segment.Predicted(:,1,1)));
assert(all(S3.Predicted(:,2) == S4.Segment.Predicted(:,2,1)));
assert(all(S3.Predicted(:,1) == S4.Segment.Predicted(:,1,2)));
assert(all(S3.Predicted(:,2) == S4.Segment.Predicted(:,2,2)));
assert(all(S3.Z(:) == S4.Z(:)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% API Test - Intervals
% When there are gaps in time in the input/data, one can pass a cell array
% of intervals and then the transfer function is computed on each interval.
% The intervals may be segemented by specifying a window width and window
% shift that is less than the interval length.

N = 1000;
B = randn(N,1);
E = B;

opts = transferfnFD_options(1);
opts.td.window.width = N;
opts.td.window.shift = N;

S1 = transferfnFD([B;B],[E;E],opts);
S2 = transferfnFD({B;B},{E;E},opts);
assert(all(S1.Z(:) == S2.Z(:)));
assert(all(S1.Segment.Predicted(:) == S2.Segment.Predicted(:)))

S1 = transferfnFD(B,E,opts);
S2 = transferfnFD({B,[B;B]},{E,[E;E]},opts);

assert(all(S1.Predicted == S2.Segment.Predicted(:,:,1)))
assert(all(S1.Predicted == S2.Segment.Predicted(:,:,2)))
assert(all(S1.Predicted == S2.Segment.Predicted(:,:,3)))

% Average TF should be 1.0 for all fe, same as S1.Z.
S3 = transferfnFD({B,[B;B]},{0.5*E,[1.0*E;1.5*E]},opts);
assert(all(S1.Z(:)-S3.Z(:) < 2*eps))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% API Test - Stack Regression
% When intervals and/or segments are used, the default is to compute a
% transfer function that is the average of each segment. 

N = 1000;
B = randn(N,2);
E = B;

% 1 input/1 output. When using 1 segment, non-stack average should be same
% as stack average result
opts = transferfnFD_options(1);
opts.transferfnFD.loglevel = 1;
opts.td.window.width = N; % Not needed as this is default when width = NaN.
opts.td.window.shift = N; % Not needed as this is default when window = NaN.

S1 = transferfnFD(B(:,1),E(:,1),opts);
opts.fd.stack.average.function = ''; % Don't compute stack average.
S2 = transferfnFD(B(:,1),E(:,1),opts);
assert(all(S1.Z(:) == S2.Z(:)))

% 2 inputs/2 outputs. When using 1 segment, stack regression should be
% same as stack average result
opts = transferfnFD_options(1);
opts.td.window.width = N; % Not needed as this is default.
opts.td.window.shift = N; % Not needed as this is default.

S1 = transferfnFD(B,E,opts);
opts.fd.stack.average.function = ''; % Don't compute stack average.
S2 = transferfnFD(B,E,opts);
assert(all(S1.Z(:) == S2.Z(:)))

% Compare stack average Z to stack regression Z. Results not expected to be
% identical. For the stack average method, Z for each segment in a given
% frequency band is computed by regressing on the DFTs in that band segment
% Z values are averaged. For the stack regression method, the DFTs in a
% given frequency band are computed for each segment and then the segment
% frequency band DFTs are combined and a single regression is performed.
% DFTs.
N = 1000;
B = randn(N,1);
E = B;

opts = transferfnFD_options(1);
opts.transferfnFD.loglevel = 1;
opts.td.window.width = N; % Will result in two intervals.
opts.td.window.shift = N; % Will result in two intervals.

S3 = transferfnFD([B;B],[E;E],opts);
opts.fd.stack.average.function = '';
S4 = transferfnFD([B;B],[E;E],opts); % window.width and window.shift ignored.
assert(all(abs(S3.Z(:) - S4.Z(:)) < 10*eps))

% Expect identical result because intervals are identical to segments and
% both compute non-stack Z.
opts = transferfnFD_options(1);
opts.td.window.width = N;
opts.td.window.shift = N;

opts.fd.stack.average.function = '';
S3 = transferfnFD([B;B],[E;E],opts);
S4 = transferfnFD({B,B},{E,E},opts);
assert(all(S3.Z(:) == S4.Z(:)))

%assert(all(S1.Z == S2.Z));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('transferfnFD_test.m: All tests passed.\n');


break

if size(Z_FDR,2) == 4
    Zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};
    Hstrs = {'h_{xx}','h_{xy}','h_{yx}','h_{yy}'};
end
if size(Z_FDR,2) == 2
    Zstrs = {'Z_{x}','Z_{y}'};
    Hstrs = {'h_{x}','h_{y}'};
end
if size(Z_FDR,2) == 1
    Zstrs = {'Z'};
    Zstrs = {'h'};    
end

fn = 0;

fn = fn+1;figure(fn);clf;hold on;grid on;box on;
    plot(B(:,1)+15,'r')
    plot(NB(:,1)+5,'Color',[1,0.5,0])
    if dim > 1
        plot(B(:,2)+30,'r');
        plot(NB(:,2)+20,'Color',[1,0.5,0]);
    end
    plot(E(:,1)-5,'k')
    plot(NE(:,1)-15,'Color',[0.5,0.5,0.5])
    if dim == 4
        plot(E(:,2)-30,'k')
        plot(NE(:,2)-20,'Color',[0.5,0.5,0.5])
    end

    xlabel('t (sample number-1)')
    %set(gca,'YLim',[-30 30])
    if dim == 1
        ts = sprintf('E = filter(h,1,B+\\deltaB) + \\deltaE; \\deltaE = \\eta(0,%.2f); \\deltaB = \\eta(0,%.2f)',ne,nb);
        title(ts);
        ls = {'B_x+15 (input)','\deltaB_x+5 (noise)',...
              'E-5 (output)','\deltaE-15 (noise)'};
    end
    if dim == 2
        ts = sprintf('E_x = filter(h_{x},1,B_x+\\deltaBx) + filter(h_{y},1,B_y+\\deltaBy) + \\deltaE_x');
        title(ts);
        ls = {'B_x+15 (input)','\deltaB_x+5 (noise)',...
              'B_y+30 (input)','\deltaB+20 (noise)',...
              'E_x-5 (output)','\deltaE_x-15 (noise)',...
              'E_y-30 (output)','\deltaE_x-20 (noise)',...
              };
    end
    legend(ls,'Location','Best');
    plotcmds(['timeseries',paramstring],writeimgs)

fn = fn+1;figure(fn);clf
    if dim >=1
        loglog(f(2:end),abs(ftB(2:end,1)),'r')
        hold on;grid on;box on;
        loglog(f(2:end),abs(ftNB(2:end,1)),'Color',[0.5,0,0])
        loglog(f(2:end),abs(ftE(2:end,1)),'k')
        loglog(f(2:end),abs(ftNE(2:end,1)),'Color',[0.5,0.5,0.5])
        ls = {'B (input)','\deltaB (noise)','E (output)','\deltaE (noise)'};
    end
    if dim == 2
        loglog(f(2:end),abs(ftB(2:end,2)),'g')
        hold on;grid on;box on;
        loglog(f(2:end),abs(ftNB(2:end,2)),'Color',[0,0.5,0])
        ls = {...
            'B_x (input)','\deltaB_x (noise)',...
            'E_x (output)','\deltaE_x (noise)',...
            'B_y (input)','\deltaB_y (noise)',...
            };
    end
    legend(ls,'Location','Best');
    title('Raw Periodograms');
    xlabel('f')
    plotcmds(['rawperiodograms',paramstring],writeimgs)

fn = fn+1;figure(fn);clf;hold on;grid on;
    xc = xcorr(E(:,1),B(:,1),'unbiased');
    tl = [-N+1:N-1];
    xc = fftshift(xc);
    plot(tl,xc,'Color','r','Marker','.','MarkerEdgeColor','r');
    set(gca,'XLim',[-3*length(h) 3*length(h)]);
    title('Raw Cross Correlation')
    xlabel('lag')
    legend('E,B')
    plotcmds(['crosscorrelation',paramstring],writeimgs)

for j = 1:size(H_FDR,2)
    fn = fn+1;figure(fn);clf;hold on;grid on;box on;;
    plot(th,H(:,j),'k','Marker','+','MarkerSize',5,'LineWidth',5);
    plot(t_TD,H_TD(:,j),'m','LineWidth',4);
    for i = 1:size(H_FDR,3)
        plot(t_FDR,H_FDR(:,j,i),'LineWidth',2);
    end
    xlabel('t')
    %set(gca,'XLim',[-3*length(h) 3*length(h)]);
    set(gca,'XLim',[-20,40]);    
    title(sprintf('%s',Hstrs{j}));
    legend(...
            'Exact',...
            sprintf('TD; n_T = %d',length(H_TD)),...
            sprintf('FD; %s; OLS; min E',windowfn),...
            sprintf('FD; %s; OLS; min E using regress()',windowfn),...
            sprintf('FD; %s; Robust; min E using robustfit()',windowfn),...
            sprintf('FD; %s; OLS; min B',windowfn),...
            sprintf('FD; %s; OLS; min B using regress()',windowfn),...
            sprintf('FD; %s; Robust; min B using robustfit()',windowfn),...
            'Location','Best')
   plotcmds(['impulse_response_functions',paramstring],writeimgs)
end

for j = 1:size(Z_FDR,2)
    fn = fn+1;figure(fn);clf;
        loglog(fh,abs(Z(:,j)),'k','Marker','+','MarkerSize',10,'LineWidth',5)
        hold on;grid on;box on;
        loglog(fe_TD,abs(Z_TD(:,j)),'m','Marker','.','MarkerSize',25,'LineWidth',3);
        for i = 1:size(Z_FDR,3)
            loglog(fe_FDR,abs(Z_FDR(:,j,i)),'Marker','.','MarkerSize',15,'LineWidth',2);
         end
        xlabel('f')
        title(sprintf('%s',Zstrs{j}));
        legend(...
                'Exact',...
                sprintf('TD; n_T = %d',length(H_TD)),...
                sprintf('FD; %s; OLS; min E',windowfn),...
                sprintf('FD; %s; OLS; min E using regress()',windowfn),...
                sprintf('FD; %s; Robust; min E using robustfit()',windowfn),...
                sprintf('FD; %s; OLS; min B',windowfn),...
                sprintf('FD; %s; OLS; min B using regress()',windowfn),...
                sprintf('FD; %s; Robust; min B using robustfit()',windowfn),...
                'Location','Best')
       plotcmds(['transfer_functions',paramstring],writeimgs)
end
   
break
if (0)
figure(5);clf;
    % Create padded impulse responses.
    tp = [min([t_FDP(1),t_FDR(1),th(1)]):max([t_FDP(end),t_FDR(end),th(end)])]';
    hp = interp1(th,h,tp);
    hp(isnan(hp)) = 0;
    for i = 1:size(H_TD,2)
        H_TDp(:,i)  = interp1(t_TD,H_TD(:,i),tp);
        H_FDRp(:,i) = interp1(t_FDR,H_FDR(:,i),tp);
        H_FDPp(:,i) = interp1(t_FDP,H_FDP(:,i),tp);
    end

    hold on;grid on;
    plot(tp,H_TDp(:,2)-hp,'m','LineWidth',4)
    plot(tp,H_FDRp(:,2)-hp,'b','LineWidth',2)
    plot(tp,H_FDPp(:,2)-hp,'g','LineWidth',2)
    xlabel('t')
    title('Impulse Response Errors')
    legend(...
            sprintf('\\deltah_{xy} time domain; n_T = %d',length(H_TD)),...
            sprintf('\\deltah_{xy} freq. domain Rectangular; n_f = %d', df),...
            sprintf('\\deltah_{xy} freq. domain Parzen; n_P = %d',length(fe_FDP))...
            )
   plotcmds(['impulse_response_errors',paramstring],writeimgs)
end

figure(6);clf;
    loglog(fh,abs(Z(:,2)),'k','Marker','+','MarkerSize',10,'LineWidth',5)
    hold on;grid on;
    loglog(fe_TD,abs(Z_TD(:,2)),'m','Marker','.','MarkerSize',25,'LineWidth',3);
    loglog(fe_FDR,abs(Z_FDR(:,2)),'b','Marker','.','MarkerSize',15,'LineWidth',2);
    loglog(fe_FDP,abs(Z_FDP(:,2)),'g','Marker','.','MarkerSize',10,'LineWidth',1);
    xlabel('f')
    title('Transfer Functions')
    legend(...
            'Z_{xy}',...
            sprintf('Z_{xy} time domain; n_T = %d',length(H_TD)),...
            sprintf('Z_{xy} freq. domain rectangular (n_R = %d)', df),...
            sprintf('Z_{xy} freq. domain parzen; n_P = %d',length(fe_FDP))...
            )
   plotcmds(['transfer_functions',paramstring],writeimgs)

figure(9);clf;
    hold on;grid on;
    plot(fh,abs(P(:,2)),'k','Marker','+','MarkerSize',10,'LineWidth',5)
    plot(fe_TD,abs(P_TD(:,2)),'m','Marker','.','MarkerSize',25,'LineWidth',3);
    plot(fe_FDR,abs(P_FDR(:,2)),'b','Marker','.','MarkerSize',25,'LineWidth',2);
    plot(fe_FDP,abs(P_FDP(:,2)),'g','Marker','.','MarkerSize',25,'LineWidth',1);
    xlabel('f');
    ylabel('degrees')
    title('Transfer Function Phases')
    legend(...
        '\phi_{xy}',...
        sprintf('\\phi_{xy} time domain; n_T = %d',length(H_TD)),...
        sprintf('\\phi_{xy} freq. domain Rectangular; n_R = %d', df),...
        sprintf('\\phi_{xy} freq. domain Parzen; n_P = %d',length(fe_FDP)),...
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

if (0)
figure(11);clf;
    hold on;grid on;
    plot(E(:,2),'k','LineWidth',3)
    plot(E_TD(:,2),'m')
    plot(E_FDR_wH(:,2),'b')
    plot(E_FDP_wH(:,2),'g')
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
    plot(E_TD_wZ(:,2),'m')
    plot(E_FDR(:,2),'b')
    plot(E_FDP(:,2),'g')
    xlabel('t')
    title('Predictions (using Z)')
    legend('E_y',...
            'E_y time domain',...
            'E_y freq. domain Rectangular',...
            'E_y freq. domain Parzen'...
            )
   plotcmds(['predictions_Z',paramstring],writeimgs)
end

figure(13);clf;
    hold on;grid on;
    plot(E(:,2)-E_TD(:,2)+10,'m')
    plot(E(:,2)-E_FDR_wH(:,2),'b')
    plot(E(:,2)-E_FDP_wH(:,2)-10,'g')
    peTD = pe_nonflag(E(:,2),E_TD(:,2));
    peFDR  = pe_nonflag(E(:,2),E_FDR_wH(:,2));
    peFDP  = pe_nonflag(E(:,2),E_FDP_wH(:,2));

    xlabel('t')
    title('Prediction Errors (using H)')
    set(gca,'YLim',[-20 20])
    legend(...
        sprintf('\\DeltaE_y+10 time domain; PE = %.3f',peTD),...        
        sprintf('\\DeltaE_y freq. domain Rectangular; PE = %.3f',peFDR),...
        sprintf('\\DeltaE_y-10 freq. domain Parzen; PE = %.3f',peFDP)...
        );
    plotcmds(['prediction_H_errors',paramstring],writeimgs)

figure(14);clf;
    hold on;grid on;
    plot(E(:,2)-E_TD_wZ(:,2)+10,'m')
    plot(E(:,2)-E_FDR(:,2),'b')
    plot(E(:,2)-E_FDP(:,2)-10,'g')
    peTD = pe_nonflag(E(:,2),E_TD_wZ(:,2));
    peFDR = pe_nonflag(E(:,2),E_FDR(:,2));
    peFDP = pe_nonflag(E(:,2),E_FDP(:,2));
    xlabel('t')
    title('Prediction Errors (using Z)')
    set(gca,'YLim',[-20 20])
    legend(...
        sprintf('\\DeltaE_y+10 time domain; PE = %.3f',peTD),...        
        sprintf('\\DeltaE_y freq. domain Rectangular; PE = %.3f',peFDR),...
        sprintf('\\DeltaE_y-10 freq. domain Parzen; PE = %.3f',peFDP)...
        );
    plotcmds(['prediction_Z_errors',paramstring],writeimgs)
