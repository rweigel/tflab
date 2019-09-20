function transferfnZ_plot(S1,S2)

figprep();

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

Zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};

if nargin == 1
    subplot('Position', PositionTop);
        loglog(S1.fe,abs(S1.Z),'k','marker','.','markersize',3);
        grid on;box on;hold on;
        unitstr = sprintf('[%s/(%s)]',...
                            S1.Options.info.inunit,...
                            S1.Options.info.outunit);
        title(S1.Options.description,'FontWeight','Normal');
        if size(S1.Z,2) == 1
            legend(sprintf('$|Z|$ %s', unitstr),'Location','NorthEast');
        else
            % TODO.
        end
        set(gca,'XTickLabel',[]);
        %adjust_exponent('y');
    subplot('Position', PositionBottom);
        semilogx(S1.fe, (180/pi)*(atan2(imag(S1.Z),real(S1.Z))),...
                'k','marker','.','markersize',3);
        grid on;box on;hold on;
        %ylabel('[degrees]');
        set(gca,'YLim',[-185,185])
        set(gca,'YTick',[-180:60:180]);
        xlabel(sprintf('$f$ [1/%s]', S1.Options.info.timeunit));
        legend('$\phi$ [degrees]','Location','NorthEast');
        %adjust_exponent();
end

if nargin == 2
    %assert(all(size(S1.Z) == size(S2.Z)), 'Required: size(S1.Z) == size(S2.Z)');
    % Assumes units are the same. TODO: Allow different unit labels.
    subplot('Position', PositionTop);
        loglog(S1.fe, abs(S1.Z));
        grid on;box on;hold on;
        loglog(S2.fe, abs(S2.Z));
        adjust_exponent('y');
        ylabel(sprintf('[%s/(%s)]',...
                    S1.Options.info.inunit,...
                    S1.Options.info.outunit));
        legend(['$|Z|$ $\,$ ', S1.Options.description],...
               ['$|Z|$ $\,$ ', S2.Options.description],...
               'Location','NorthEast');
        set(gca,'XTickLabel',[]);
    subplot('Position', PositionBottom);
        semilogx(S1.fe, (180/pi)*unwrap(atan2(imag(S1.Z),real(S1.Z))));
        grid on;box on;hold on;
        semilogx(S2.fe, (180/pi)*unwrap(atan2(imag(S2.Z),real(S2.Z))));
        ylabel('[degrees]');
        set(gca,'YLim',[-180,180])
        set(gca,'YTick',[-180:60:180]);
        adjust_exponent();
        xlabel(sprintf('$f$ [1/%s]', S1.Options.info.timeunit));
        legend(['$\phi$ $\,$ ', S1.Options.description],...
               ['$\phi$ $\,$ ', S2.Options.description],...
               'Location','NorthEast');
end