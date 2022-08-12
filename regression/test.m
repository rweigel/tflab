clear
%addpath('m-statrobust');

opts = struct();
    opts.weightfn = 'huber';
    opts.stepmax = 50;
    opts.zeps = sqrt(eps);
    opts.hardcut = Inf; % 2.8
    opts.snstop = 1000;
    opts.verbose = 0;
    
Nk = 10;
ey = 0.1;
ex = 0.1;
a = 2;
N = 10;
x = [0:N-1]'/N;
for k = 1:Nk
    y = a*(x+ex) + ey*randn(size(x));
    %y(end) = 10*y(end); % Add outlier

    [a1(1,k),W] = robust_v1(x,y,opts);
    W1(:,k) = W;

    [a2(1,k),stats2] = robustfit(x,y,'huber',[],'off');
    W2(:,k) = stats2.w;

    a3(1,k) = regress(y,x);

    X = [x,y];
    [coeff,score,roots] = pca(X);
    a4(1,k) = coeff(2)/coeff(1);
    % Weighted residuals for each fit type and each of Nk fits.

    r1(:,k) = W1(:,k).*(y - a1(1,k)*x);
    r2(:,k) = W2(:,k).*(y - a2(1,k)*x);
    r3(:,k) = y - a3(1,k)*x;
    r4(:,k) = y - a4(1,k)*x;
end

fprintf('Parameter estimates\n');
fprintf('robust_v1 robustfit regress orthogonal\n');
fprintf('%5.2f    %5.2f     %5.2f   %5.2f',mean([a1;a2;a3;a4],2))
fprintf('\n');
fprintf('Average of weighted residuals across all %d fits\n',Nk);
fprintf('%5.2f    %5.2f     %5.2f   %5.2f',mean([r1(:),r2(:),r3(:),r4(:)])');
fprintf('\n');


figprep();figure(1);clf;
    plot(x,y,'.','MarkerSize',20);
    grid on;hold on;
    plot(x,a*x,'k-','LineWidth',2);
    plot(x,mean(a1)*x,'b-','LineWidth',4);
    plot(x,mean(a2)*x,'g-','LineWidth',2);
    plot(x,mean(a3)*x,'r-','LineWidth',2);
    plot(x,mean(a4)*x,'m-','LineWidth',2);
    xlabel('$x$');
    ylabel('$y$');
    legend('Measurements','Actual','robust\_v1','robustfit','regress','orthogonal','Location','NorthWest')

figprep();figure(2);clf;
    plot(x,y,'.','MarkerSize',20);
    grid on;hold on;
    plot(x,a*x,'k-','LineWidth',2);
    plot(x,mean(a1)*x,'b-','LineWidth',4);
    plot(x,mean(a2)*x,'g-','LineWidth',2);
    plot(x,mean(a3)*x,'r-','LineWidth',2);
    set(gca,'XLim',[0,x(end-1)])
    xlabel('$x$');
    ylabel('$y$');
    legend('Measurements','Actual','robust\_v1','robustfit','regress','Location','NorthWest')
    
figprep();figure(3);clf;
    subplot(3,1,1)
        qqplot(a1(1,:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','$a$');
        set(get(gca,'title'),'String','robust\_v1');
    subplot(3,1,2)
        qqplot(a2(1,:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','$a$');
        set(get(gca,'title'),'String','robustfit');
    subplot(3,1,3)
        qqplot(a3(1,:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','$a$');
        set(get(gca,'title'),'String','regress');

figprep();figure(4);clf
    subplot(3,1,1)
        qqplot(r1(:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','Weighted residuals');
        set(get(gca,'title'),'String','robust\_v1');
    subplot(3,1,2)
        qqplot(r2(:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','Weighted residuals');
        set(get(gca,'title'),'String','robustfit');
    subplot(3,1,3)
        qqplot(r3(:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','Weighted residuals');
        set(get(gca,'title'),'String','regress');

if 0
N = 4;

Z = 2-1i;

B = fft(rand(N,1));

E = B*transpose(Z);%+randn(size(B))+sqrt(-1)*randn(size(B));

%E(10) = 100;

Zo = regress(E,B)

[Zr,stats] = robustfit(B,E,[],[],'off');
Zr

Uo = regress([real(E);imag(E)],[real(B),-imag(B);imag(B),real(B)]);
Uo'
R = [1 1;1 -1];
% U = RZ
Z = R\Uo;
Z = Z(1) + sqrt(-1)*Z(2)
end
