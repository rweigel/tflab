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
zao = @(z) [z(1:end-1); 10];
% For no outliers, use
% zao = @(z) z;

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

    % Add outliers
    z = zao(z);

    [m1(:,k), r1(:,k)] = regress_ols(z, [x,y], 'backslash');
    [m2(:,k), r2(:,k)] = regress_ols(z, [x,y], 'regress');
    [m3(:,k), r3(:,k)] = regress_tls(z, [x,y]);

    [m4(:,k), r4(:,k), w4(:,k)] = regress_robustfit_tflab( z,[x,y],robustfit_tflab_options);
    [m5(:,k), r5(:,k), w5(:,k)] = regress_robustfit_matlab(z,[x,y],robustfit_matlab_options{:});

end

% TODO: Use table() function
hline = [repmat('-',1,57),'\n'];
columns = 'backslash regress  tls  robust_tflab robust_matlab';
fprintf(hline);
fprintf('Parameter estimates\n');
fprintf(hline);
fprintf('fit #  param  %s\n',columns);
for k = 1:Nr
    dimlabels = {'mx','my'};
    for dim = 1:2
        m_means = mean([m1(dim,k); m2(dim,k); m3(dim,k); m4(dim,k); m5(dim,k)], 2);
        fprintf('%4d    %s     %5.2f    %5.2f   %5.2f     %5.2f       %5.2f\n',...
                k,dimlabels{dim},m_means);
    end
end

fprintf('\n');
fprintf(hline);
fprintf('Weighted residuals (average of across %d fits)\n', Nr);
fprintf(hline);
fprintf('              %s\n',columns);
fprintf('               %5.2f    %5.2f   %5.2f     %5.2f       %5.2f\n',...
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