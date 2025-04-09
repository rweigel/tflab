clear;

robustfit_matlab_options = {'huber',[],'off'};

robustfit_tflab_options = struct();
    robustfit_tflab_options.weightfn = 'huber';
    robustfit_tflab_options.stepmax  = 50;
    robustfit_tflab_options.zeps     = sqrt(eps);
    robustfit_tflab_options.hardcut  = 2.8;
    robustfit_tflab_options.snstop   = 1000;
    robustfit_tflab_options.verbose  = 0;

Nr = 10;   % Number of regressions to perform.
Np = 100;  % Number of points per regression.

% Regression is performed on z = mx*x + my*y.
mx = 1; % Slope in input x 
my = 1; % Slope in input y

% Anonymous function zao adds outliers to z
% Set last value of z to 10
zao = @(z) [z(1:end-2); 10; -10];
% For no outliers, use
%zao = @(z) z;

% Standard deviation of noise on inputs and output
exo = 0.5;
eyo = 0.5;
ezo = 0.5;

% Input values
xo = randn(Np,1);
yo = randn(Np,1);

for k = 1:Nr
    
    % Noise terms
    ex = exo*randn(Np,1);
    ey = eyo*randn(Np,1);
    ez = ezo*randn(Np,1);

    % Add noise to inputs
    x = xo + ex;
    y = yo + ey;
    
    % Compute output and add noise
    z = mx*x + my*y + ez;

    % Add outlier
    z = zao(z);

    [m1(:,k), info1] = regress_ols(z, [x,y], 'backslash');
    [m2(:,k), info2] = regress_ols(z, [x,y], 'regress');
    [m3(:,k), info3] = regress_tls(z, [x,y]);

    [m4(:,k), info4] = regress_robustfit_tflab( z,[x,y],robustfit_tflab_options);
    [m5(:,k), info5] = regress_robustfit_matlab(z,[x,y],robustfit_matlab_options{:});
    
    r1(:,k) = info1.Residuals;
    r2(:,k) = info2.Residuals;
    r3(:,k) = info3.Residuals;
    r4(:,k) = info4.Residuals;
    r5(:,k) = info5.Residuals;
end
m1a = mean(m1,2);
m2a = mean(m2,2);
m3a = mean(m3,2);
m4a = mean(m4,2);
m5a = mean(m5,2);

% TODO: Use table() function
hline = [repmat('-',1,57),'\n'];
columns = 'backslash regress  tls  robust_tflab robust_matlab';
fprintf(hline);
fprintf('Parameter estimates\n');
fprintf(hline);
fprintf('fit #  param  %s\n',columns);
dimlabels = {'mx','my'};
for k = 1:Nr
    for dim = 1:2
        ms = [m1(dim,k); m2(dim,k); m3(dim,k); m4(dim,k); m5(dim,k)];
        fprintf('%4d    %s     %5.2f    %5.2f   %5.2f     %5.2f       %5.2f\n',...
                k,dimlabels{dim},ms);
    end
end

msa = [m1a(dim); m2a(dim); m3a(dim); m4a(dim); m5a(dim)];

fprintf('\n');
for dim = 1:2
    fprintf('        %s     %5.2f    %5.2f   %5.2f     %5.2f       %5.2f\n',...
            dimlabels{dim},msa);
end
fprintf('\n');
fprintf(hline);
fprintf('Weighted residuals (average of across %d fits)\n', Nr);
fprintf(hline);
fprintf('              %s\n',columns);
fprintf('               %6.3f    %6.3f   %6.3f     %6.3f       %6.3f\n',...
         mean([r1(:),r2(:),r3(:),r4(:),r5(:)])');
fprintf('\n\n');

% Create a plane based on slopes mx and my.
[xgrid,ygrid] = meshgrid(linspace(min(x),max(x),5), ...
                         linspace(min(y),max(y),5));
zgrid = xgrid*mx + ygrid*my;                     

zgridr = real(zgrid);
zr = real(z);

% Create plane based on average of TLS (orthogonal least squares) slopes.
%m2avg = mean(m2,2);
%zgrid = xgrid*m2avg(1) + ygrid*m2avg(2);

figure(1);clf;
    mesh(xgrid,ygrid,zgridr,'EdgeColor',[0 0 0],'FaceAlpha',0);
    hold on;
    plot3(x,y,zr,'*');
    xlabel('x')
    ylabel('y')
    if isreal(z)
        zlabel('z');
    else
        zlabel('Re(z)');
    end

if ~isreal(z)    
figure(2);clf;
    mesh(xgrid,ygrid,imag(zgrid),'EdgeColor',[0 0 0],'FaceAlpha',0);
    hold on;
    plot3(x,y,imag(z),'*');
    xlabel('x')
    ylabel('y')
    zlabel('Im(z)')
end    