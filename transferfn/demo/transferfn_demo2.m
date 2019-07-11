% Comparison of two approaches to determining transfer function for
% and infinite half-plane.

clear;
addpath('../../stats')
addpath('../../time')
writeimgs = 0;

base = 'transfer_function_demo2'; % Output files will be named ./figures/base_...

%for alpha = [-1/2,1/2]
for alpha = [-1/2]

    if (alpha == -0.5)  
        vs = 'Ex_By';
    end
    if (alpha == +0.5)
        vs = 'Ex_dBydt';
    end

    N  = 10^4;
    t  = [0:N-1]';
    dp = -pi*(alpha+1)/2;
    f  = [1:N/2]/N;

    Ey = 0;
    Bx = 0;
    for i = 1:length(f)
        phi(i) = 2*pi*(rand(1)-0.5); % Random phase in [-pi,pi].
        Ey     = Ey + (1/f(i))               *cos(2*pi*t*f(i)   +phi(i));
        Bx     = Bx - (1/f(i))*(f(i)^(alpha))*cos(2*pi*t*f(i)+dp+phi(i));    
    end

    % Subtract off mean (not needed if random phase is used and large N).
    Ey = Ey-mean(Ey);
    Bx = Bx-mean(Bx);

    N  = length(Ey);
    t  = [0:N-1]';
    f  = [0:N/2]/N;

    % Time domain estimate of transfer function
    Nl    = 800;
    Na    = 800;
    %[T,X] = time_delay(Ey,Bx,Nl,0,Nl);
    [T,X] = time_delay(Ey,Bx,Nl,0,Na);
    LIN   = basic_linear(X,T);
    tbl   = [-Na+1:1:Nl];
    hbl   = LIN.Weights(1:end-1);    
    Rbl   = fft(hbl);
    phibl = (180/pi)*atan2(imag(Rbl),real(Rbl));
    Nbl   = length(hbl);
    frbl  = [0:Nbl/2]/Nbl;

    % Needed to get phase that is comparable to frequency domain method's.
    tmp    = fft([hbl(end/2:end);hbl(1:end/2-1)]);
    phibl2 = (180/pi)*atan2(imag(tmp),real(tmp))

    % Frequency domain estimate of transfer function
    Eyfft = fft(Ey);
    Bxfft = fft(Bx);
    Rft   = Eyfft./Bxfft;
    phift = (180/pi)*atan2(imag(Rft),real(Rft));
    Nft   = length(Ey);
    frft  = [0:Nft/2]/Nft;

    hft   = ifft(Eyfft./Bxfft);
    hft   = fftshift(hft);
    tft   = [-N/2:1:N/2-1];
    %tmp = fft(hft);
    %phift = (180/pi)*atan2(imag(tmp),real(tmp))

    figure(1);clf;
        plot(t,t*NaN,'k-','LineWidth',2); % To make legend lines are larger.
        hold on;grid on;
        plot(t,t*NaN,'g-','LineWidth',2);

        plot(t,Ey,'k');
        plot(t,Bx,'g');
        if (alpha == -0.5)  
            lh = legend('$E_y$','$B_x$');
            vs = 'Ex_By';
        end
        if (alpha == +0.5)
            lh = legend('$E_y$','$B_x''$');
            vs = 'Ey_dBxdt';
        end
        set(lh,'Interpreter','Latex');
        xlabel('t');
        if (writeimgs)
            print('-dpng','-r150',sprintf('figures/%s_timeseries_%s.png',base,vs));
            print('-depsc',sprintf('figures/%s_timeseries_%s.eps',base,vs));
            fprintf('Wrote figures/%s_timeseries_%s.{png,eps}\n',base,vs)
        end
        
    figure(2);clf;
        loglog(frft(2:Nft/2),abs(Bxfft(2:Nft/2)),'k','LineWidth',2,'Marker','.','MarkerSize',10);
        hold on;grid on;
        loglog(frft(2:Nft/2),abs(Eyfft(2:Nft/2)),'g','LineWidth',2,'Marker','.','MarkerSize',10);
        %loglog(f(2:N/2),abs(R(2:N/2)),'m','LineWidth',2,'Marker','.','MarkerSize',10);
        xlabel('f');
        if (alpha == -0.5)
            lh = legend('$\|\widetilde{E}_y\|$','$\|\widetilde{B}_x\|$');
        end
        if (alpha == +0.5)
            lh = legend('$\|\widetilde{E}_y\|$','$\|\widetilde{B}_x''\|$');
        end
        set(lh,'Interpreter','Latex');
        if (writeimgs)
            print('-dpng','-r150',sprintf('figures/%s_dft_%s.png',base,vs));
            print('-depsc',sprintf('figures/%s_dft_%s.eps',base,vs));
            fprintf('Wrote figures/%s_dft_%s.{png,eps}\n',base,vs);
        end

    figure(3);clf;
        plot(tbl,hbl,'b','LineWidth',3,'Marker','.','MarkerSize',30);
        hold on;grid on;
        plot(tft,hft,'r','LineWidth',2,'Marker','.','MarkerSize',20);
        set(gca,'XLim',[tbl(1)-5 tbl(end)+5]);
        xlabel('t');
        if (alpha == -0.5)
            th = title('Response of $E_y$ to $B_x = \delta(t)$');
        end
        if (alpha == 0.5)
            th = title('Response of $E_y$ to $B_x'' = \delta(t)$');
        end
        set(th,'Interpreter','Latex');    
        legend('Time domain method','Freq. domain method');
        if (writeimgs)
            print('-dpng','-r150',sprintf('figures/%s_irf_%s.png',base,vs));
            print('-depsc',sprintf('figures/%s_irf_%s.eps',base,vs));
            fprintf('Wrote figures/%s_irf_%s.{png,eps}\n',base,vs)
        end

    figure(4);clf;
        loglog(frbl(2:Nbl/2),Rbl(2:Nbl/2),'b','Marker','.','MarkerSize',10);
        hold on;grid on;
        loglog(frft(2:Nft/2),Rft(2:Nft/2),'r','Marker','.','MarkerSize',10);
        hold on;grid on;
        xlabel('f');
        if (alpha == -0.5)
            th = title('$\|\widetilde{E_y}\|/\|\widetilde{B}_x\|$');
        end
        if (alpha == +0.5)
            th = title('$\|\widetilde{E}_y\|/\|\widetilde{B}''_x\|$');
        end
        set(th,'Interpreter','Latex');
        legend('Time domain method','Freq. domain method');
        if (writeimgs)
            print('-dpng','-r150',sprintf('figures/%s_transferfn_%s.png',base,vs));
            print('-depsc',sprintf('figures/%s_transferfn_%s.eps',base,vs));
            fprintf('Wrote figures/%s_transferfn_%s.{png,eps}\n',base,vs)
        end

    figure(5);clf;
        plot(frbl(2:Nbl/2),phibl(2:Nbl/2),'b.','Marker','.','MarkerSize',10);
        hold on;grid on;
        plot(frft(2:Nft/2),phift(2:Nft/2),'r.','Marker','.','MarkerSize',10);
        plot(frbl(2:Nbl/2),phibl2(2:Nbl/2),'k.','Marker','.','MarkerSize',10);
        xlabel('f');
        if (alpha == -0.5)
            th = title('Phase $\|\widetilde{E_y}\|/\|\widetilde{B}_x\|$');
        end
        if (alpha == +0.5)
            th = title('Phase $\|\widetilde{E}_y\|/\|\widetilde{B}''_x\|$');
        end
        set(th,'Interpreter','Latex');
        legend('Time domain method','Freq. domain method','Alt Time domain calc.');
        if (writeimgs)
            print('-dpng','-r150',sprintf('figures/%s_phase_%s.png',base,vs));
            print('-depsc',sprintf('figures/%s_phase_%s.eps',base,vs));
            fprintf('Wrote figures/%s_phase_%s.{png,eps}\n',base,vs);
        end

end % alpha    
    
