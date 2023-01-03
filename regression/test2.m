clear

opts = struct();
    opts.weightfn = 'huber';
    opts.stepmax = 50;
    opts.zeps = sqrt(eps);
    opts.hardcut = Inf; % 2.8
    opts.snstop = 1000;
    opts.verbose = 0;

Nk = 1;
N = 1000;
mx = 1;
my = 1;
ex = 0.2*randn(N,1);
ey = 0.2*randn(N,1);
ez = 0.1*randn(N,1);
x = randn(N,1);
y = randn(N,1);
for k = 1:Nk
    z = mx*x + my*y;
    x = x + ex;
    y = y + ey;
    z = z + ez;
    %y(end) = 10*y(end); % Add outlier

    [a1(:,k),W] = robust_v1([x,y],z,opts);
    W1(:,k) = W;

    [a2(:,k),stats2] = robustfit([x,y],z,'huber',[],'off');
    W2(:,k) = stats2.w;

    a3(:,k) = regress(z,[x,y]);

    X = [x,y,z];
    [coeff,score,roots] = pca(X);
    a4(:,k) = -coeff(1:end-1,end)/coeff(end,end);
   
    % Weighted residuals for each fit type and each of Nk fits.

    r1(:,k) = W1(:,k).*(y - a1(1,k)*x);
    r2(:,k) = W2(:,k).*(y - a2(1,k)*x);
    r3(:,k) = y - a3(1,k)*x;
    r4(:,k) = y - a4(1,k)*x;
end
a1
a2
a3
a4

clf
normal = coeff(:,2);
meanX = mean([x,y]);
[xgrid,ygrid] = meshgrid(linspace(min(x),max(x),5), ...
                         linspace(min(y),max(y),5));
s4 = mean(a4,2);
zgrid = xgrid*s4(1) + ygrid*s4(2);
%zgrid = (xgrid*mx + ygrid*my);
h = mesh(xgrid,ygrid,zgrid,'EdgeColor',[0 0 0],'FaceAlpha',0);
hold on;
plot3(x,y,z,'*')
xlabel('x')
ylabel('y')
zlabel('z')