function transferfnZ_plot(S1,S2)

figprep();

PositionTop = [0.1300 0.5400 0.7750 0.4];
PositionBottom = [0.1300 0.1100 0.7750 0.4];

if size(S1.Z,2) > 1
    Zstrs = {'Z_{xx}','Z_{xy}','Z_{yx}','Z_{yy}'};
    Phistrs = {'\phi_{xx}','\phi_{xy}','\phi_{yx}','\phi_{yy}'};
else
    Zstrs = {'Z'};
    Phistrs = {'Z'};
end

if nargin == 1
    subplot('Position', PositionTop);
        loglog(S1.fe,abs(S1.Z),'marker','.','markersize',10,'linewidth',2);
        grid on;box on;hold on;
        unitstr = sprintf('[%s/(%s)]',...
                            S1.Options.info.inunit,...
                            S1.Options.info.outunit);
        title(S1.Options.description,'FontWeight','Normal');
        if size(S1.Z,2) == 1
            legend(sprintf('$|Z|$ %s', unitstr),'Location','NorthEast');
        else
            for j = 1:size(S1.Z,2)
                ls{j} = sprintf('$|%s|$ %s',Zstrs{j},unitstr);
            end
            legend(ls,'Location','NorthEast');
        end
        set(gca,'XTickLabel',[]);
        %adjust_exponent('y');
    subplot('Position', PositionBottom);
        semilogx(S1.fe, (180/pi)*(atan2(imag(S1.Z),real(S1.Z))),...
                'linewidth',2,'marker','.','markersize',10);
        grid on;box on;hold on;
        %ylabel('[degrees]');
        set(gca,'YLim',[-185,185])
        set(gca,'YTick',[-180:60:180]);
        xlabel(sprintf('$f$ [1/%s]', S1.Options.info.timeunit));
        if size(S1.Z,2) == 1
            legend(sprintf('$|\\phi|$ %s', unitstr),'Location','NorthEast');
        else
            for j = 1:size(S1.Z,2)
                ls{j} = sprintf('$%s$ [$^\\circ$]',Phistrs{j});
            end
            ls
            legend(ls,'Location','NorthEast');
        end
        %legend('$\phi$ [degrees]','Location','NorthEast');
        %adjust_exponent();
end

if nargin == 2
    %assert(all(size(S1.Z) == size(S2.Z)), 'Required: size(S1.Z) == size(S2.Z)');
    % Assumes units are the same. TODO: Allow different unit labels.
    for j = 1:size(S1.Z,2)
        if j > 1; figure(); figprep(); end
        subplot('Position', PositionTop);
            loglog(S1.fe, abs(S1.Z(:,j)));
            grid on;box on;hold on;
            loglog(S2.fe, abs(S2.Z(:,j)));
            adjust_exponent('y');
            ylabel(sprintf('[%s/(%s)]',...
                        S1.Options.info.inunit,...
                        S1.Options.info.outunit));
            legend({sprintf('$|%s|$ %s',Zstrs{j},S1.Options.description),...
                    sprintf('$|%s|$ %s',Zstrs{j},S2.Options.description)},...
                   'Location','NorthEast');
            set(gca,'XTickLabel',[]);
        subplot('Position', PositionBottom);
            semilogx(S1.fe, (180/pi)*unwrap(atan2(imag(S1.Z(:,j)),real(S1.Z(:,j)))));
            grid on;box on;hold on;
            semilogx(S2.fe, (180/pi)*unwrap(atan2(imag(S2.Z(:,j)),real(S2.Z(:,j)))));
            ylabel('[degrees]');
            set(gca,'YLim',[-180,180])
            set(gca,'YTick',[-180:60:180]);
            %adjust_exponent();
            xlabel(sprintf('$f$ [1/%s]', S1.Options.info.timeunit));
            legend({sprintf('$|%s|$ $\\,$ %s',Phistrs{j},S1.Options.description),...
                    sprintf('$|%s|$ $\\,$ %s',Phistrs{j},S2.Options.description)},...
                   'Location','NorthEast');
    end
end