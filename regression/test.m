clear
addpath('m-statrobust');

opts = struct();
    opts.weightfn = 'huber';
    opts.stepmax = 50;
    opts.zeps = sqrt(eps);
    opts.hardcut = Inf; % 2.8
    opts.snstop = 1000;
    opts.verbose = 0;
    
Nk = 100;
ex = 2;
a = 1;
N = 10;
for k = 1:Nk
    x = [0:N-1]';
    y = a*x + ex*randn(size(x));
    y(end) = 10*y(end);

    [a1(1,k),stats1] = robust_v1(y,x,opts);
    W1(:,k) = stats1.w;
    [a2(1,k),stats2] = robustfit(y,x,'huber',[],'off');
    W2(:,k) = stats2.w;
    a3(1,k) = regress(x,y);

    % Weighted residuals for each fit type and each of Nk fits.
    r1(:,k) = W1(:,k).*(y - a1(1,k)*x);
    r2(:,k) = W2(:,k).*(y - a2(1,k)*x);
    r3(:,k) = y - a3(1,k)*x;
end

fprintf('Parameter estimates\n');
fprintf('%.2f\t%.2f\t%.2f',mean([a1;a2;a3],2))
fprintf('\n');
fprintf('Average of weighted residuals across all %d fits\n',Nk);
fprintf('%.2f\t%.2f\t%.2f',mean([r1(:),r2(:),r3(:)])');
fprintf('\n');

figure(1);clf;
    plot(x,y,'.','MarkerSize',10);
    grid on;
    xlabel('$x$');
    ylabel('$y$');

figure(2);clf;
    subplot(3,1,1)
        qqplot(a1(1,:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','$a$');
        set(get(gca,'title'),'String','robust_v1()');
    subplot(3,1,2)
        qqplot(a2(1,:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','$a$');
        set(get(gca,'title'),'String','robustfit()');
    subplot(3,1,3)
        qqplot(a3(1,:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','$a$');
        set(get(gca,'title'),'String','regress()');

figure(3);clf
    subplot(3,1,1)
        qqplot(r1(:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','Weighted residuals');
        set(get(gca,'title'),'String','robust_v1()');
    subplot(3,1,2)
        qqplot(r2(:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','Weighted residuals');
        set(get(gca,'title'),'String','robustfit()');
    subplot(3,1,3)
        qqplot(r3(:));
        grid on;box on;
        set(get(gca,'ylabel'),'String','Weighted residuals');
        set(get(gca,'title'),'String','regress()');

break

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
