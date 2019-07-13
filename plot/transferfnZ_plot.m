function transferfnZ_plot(S1,S2)

figure;
figprep();

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

if nargin == 1
    loglog(S1.fe,abs(S1.Z));
    grid on;box on;hold on;
    exponent_relabel('y');
    unitstr = sprintf('[%s/%s]',...
                S1.Options.info.inunit,...
                S1.Options.info.outunit);
    title(sprintf('Method: %s',...
                S1.Options.description),...
                'FontWeight','Normal');
    legend(sprintf('$Z$ %s',unitstr),...
                'Location','NorthEast');
end

if nargin == 2
    % Assumes units are the same. TODO: Allow different unit labels.
    subplot('Position',PositionTop);
        loglog(S1.fe,abs(S1.Z));
        grid on;box on;hold on;
        loglog(S2.fe,abs(S2.Z));
        adjust_exponent('y');
        ylabel(sprintf('[%s/%s]',...
                    S1.Options.info.inunit,...
                    S1.Options.info.outunit));
        legend(['$Z$ Method: ',...
                    S1.Options.description],...
                ['$Z$ Method: ',...
                    S2.Options.description],...
                'Location','NorthEast');
        set(gca,'XTickLabel',[]);
    subplot('Position',PositionBottom);
        semilogx(S1.fe,(180/pi)*unwrap(atan2(imag(S1.Z),real(S1.Z))));
        grid on;box on;hold on;
        semilogx(S2.fe,(180/pi)*unwrap(atan2(imag(S2.Z),real(S2.Z))));
        adjust_exponent('x');
        ylabel('[degrees]');
        xlabel(sprintf('$f$ [1/%s]', S1.Options.info.timeunit));
        legend(['$\phi$ Method: ',...
                    S1.Options.description],...
                ['$\phi$ Method: ',...
                    S2.Options.description],...
                'Location','NorthEast');
end