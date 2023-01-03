%% 2-D plot

colors = colormap(jet(2*length(dir_list)));

S1 = Zxx + Zyy;
S2 = Zxy + Zyx;
D1 = Zxx - Zyy;
D2 = Zxy - Zyx;

% Eqn 5.10, p83
sigma = (abs(D1).^2 + abs(S2).^2)./(abs(D2).^2);

% Here we do not need to compute alpha b/c found to be zero.
% Eqn 5.16, p 85
kappa = abs(S1)./abs(D1);

% Equation 5.19, p86
mu = sqrt(abs(imag(S2.*conj(D1)) + imag(D2.*conj(S1))))./abs(D2);

figure(1);clf;
pcolor([1:74],log10(1./fi),log10(abs(Zxy)));
%title('$\log_{10}(|Z_{xy}|)$','Interpreter','Latex');
xlabel('Site #');
ylabel('$\log_{10}(T)$ [s]','Interpreter','Latex');
hcb = colorbar;
set(get(hcb,'Title') ,'String',...
    '$\log_{10}(|Z_{xy}|)$','Interpreter','Latex');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
print -dpng samtex_Zxx.png

figure(2);clf;
sigma(sigma > 1) = NaN;
pcolor([1:74],log10(1./fi),sigma);
xlabel('Site #');
ylabel('$\log_{10}(T)$ [s]','Interpreter','Latex');
title('$\Sigma > 0.05$ suggests not 1-D layered',...
      'Interpreter','Latex');
hcb = colorbar;
set(get(hcb,'Title') ,'String',...
    '$\Sigma$','Interpreter','Latex');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
print -dpng samtex_Sigma.png

figure(3);clf;
kappa(kappa > 1) = NaN;
pcolor([1:74],log10(1./fi),kappa);
xlabel('Site #');
ylabel('$\log_{10}(T)$ [s]','Interpreter','Latex');
title('Swift Skew $\kappa > 0.1$ suggests not two quarter spaces',...
      'Interpreter','Latex');
hcb = colorbar;
set(get(hcb,'Title') ,'String',...
    '$\kappa$','Interpreter','Latex');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
print -dpng samtex_kappa.png


figure(4);clf;
mu(mu > 1) = NaN;
pcolor([1:74],log10(1./fi),mu);
xlabel('Site #');
ylabel('$\log_{10}(T)$ [s]','Interpreter','Latex');
title('$\mu > 0.05$ suggests not regional 1-D w/ galvanic',...
      'Interpreter','Latex');
hcb = colorbar;
set(get(hcb,'Title') ,'String',...
    '$\mu$','Interpreter','Latex');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
print -dpng samtex_mu.png
