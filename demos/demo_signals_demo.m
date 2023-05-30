S = demo_signals('powerlaw')

dock on;figure(1);close all;

figure();
    tsplot(S);
    % Or, manually
    if 0
        subplot(2,1,1)
        plot(S.In)
        subplot(2,1,2)
        plot(S.Out)
    end

figure()
    % Populate S with data needed for psdplot().
    [PSDIn,~,f] = psd(S.In);
    S.Metrics.PSD.Raw.fe = f;
    S.Metrics.PSD.Raw.In = PSDIn;

    PSDOut = psd(S.Out);
    S.Metrics.PSD.Raw.Out = PSDOut;
    psdplot(S, struct('type','raw'));
    
    % Or, manually
    if 0
        subplot(2,1,1)
            loglog(f, PSDIn,'.');
            grid on;
            xlabel('f')
            ylabel('$|\widetilde{\mbox{In}}|^2$','Interpreter','latex')
        subplot(2,1,2)
            loglog(f, PSDOut,'.');
            grid on;
            xlabel('f')
            ylabel('$|\widetilde{\mbox{Out}}|^2$','Interpreter','latex')
    end

