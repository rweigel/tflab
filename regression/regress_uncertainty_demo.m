clear
Np = 100;
p = 1-.05/2;
sigma = 0.1;

x = [1:Np]'/Np;
err = sigma*randn(Np,1);
%err = err-mean(err);
y = x + err;

[b,bint,r,rint] = regress(y,[x,ones(Np,1)]);

bo = sum((x-mean(x)).*(y-mean(y)))/sum((x-mean(x)).^2);
ao = mean(y) - bo*mean(x);
bo = [bo,ao];

ymo = bo(1)*x + bo(2);
ym  = b(1)*x(:,1) + b(2);

% https://en.wikipedia.org/wiki/Simple_linear_regression
S = sum((ymo-y).^2); 
s2 = S/(Np-2); % Unbiased estimate of sigma^2
sb = sqrt(s2/sum((x-mean(x)).^2)); % Estimate of std of b
sa = sb*sqrt(sum(x.^2)/Np);        % Estimate of std of a
t = tinv(p,Np-2);
bintt = [bo(1) - sb*t,bo(1) + sb*t;...
         ao(1) - sa*t,ao(1) + sa*t];

% Theoretical variances in estimates (Bulmer, Principles of Statistics, p214)
Vb = sigma^2/sum((x-mean(x)).^2);
Va = sigma^2*(1/Np + mean(x)^2/sum((x-mean(x)).^2));
d = norminv(p,0,1); % 1.96 for p = 1-0.05/2.
binto = [bo(1) - d*sqrt(Vb),bo(1) + d*sqrt(Vb);...
         ao(1) - d*sqrt(Va),ao(1) + d*sqrt(Va)];

fprintf('regress: slope = %.5f [%.5f,%.5f]; intercept = %.5f [%.5f,%.5f]\n',...
    b(1),bint(1,:),b(2),bint(2,:));
fprintf('manual:  slope = %.5f [%.5f,%.5f]; intercept = %.5f [%.5f,%.5f]\n',...
    bo(1),bintt(1,:),bo(2),bintt(2,:));
fprintf('exact conf. intervals    [%.5f,%.5f];                     [%.5f,%.5f]\n',...
    binto(1,:),binto(2,:));

% Sample output:
% regress: slope = 0.97611 [0.91069,1.04153]; intercept = 0.01730 [-0.02075,0.05535]
% manual:  slope = 0.97611 [0.91069,1.04153]; intercept = 0.01730 [-0.02075,0.05535]
% exact conf. intervals    [0.90821,1.04401];                     [-0.02220,0.05679]

clf;
grid on;hold on;
plot(x,y,'.','MarkerSize',20)
plot(x,ymo,'k','LineWidth',2)
