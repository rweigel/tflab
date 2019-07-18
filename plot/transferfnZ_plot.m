function transferfnZ_plot(S1,S2)

figure;
figprep();

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

Zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};

if nargin == 1
    loglog(S1.fe, abs(S1.Z));
    grid on;box on;hold on;
    exponent_relabel('y');
    unitstr = sprintf('[%s/%s]',...
                        S1.Options.info.inunit,...
                        S1.Options.info.outunit);
    title(sprintf(S1.Options.description),'FontWeight','Normal');
    if size(S1.Z,2) == 1
        legend(sprintf('$Z$ %s', unitstr),'Location','NorthEast');
    else
        % TODO.
    end
end

if nargin == 2
    %assert(all(size(S1.Z) == size(S2.Z)), 'Required: size(S1.Z) == size(S2.Z)');
    % Assumes units are the same. TODO: Allow different unit labels.
    subplot('Position', PositionTop);
        loglog(S1.fe, abs(S1.Z));
        grid on;box on;hold on;
        loglog(S2.fe, abs(S2.Z));
        adjust_exponent('y');
        ylabel(sprintf('[%s/%s]',...
                    S1.Options.info.inunit,...
                    S1.Options.info.outunit));
        legend(['$Z$ $\,$ ', S1.Options.description],...
               ['$Z$ $\,$ ', S2.Options.description],...
               'Location','NorthEast');
        set(gca,'XTickLabel',[]);
    subplot('Position', PositionBottom);
        semilogx(S1.fe, (180/pi)*unwrap(atan2(imag(S1.Z),real(S1.Z))));
        grid on;box on;hold on;
        semilogx(S2.fe, (180/pi)*unwrap(atan2(imag(S2.Z),real(S2.Z))));
        adjust_exponent();
        ylabel('[degrees]');
        xlabel(sprintf('$f$ [1/%s]', S1.Options.info.timeunit));
        legend(['$\phi$ $\,$ ', S1.Options.description],...
               ['$\phi$ $\,$ ', S2.Options.description],...
               'Location','NorthEast');
end