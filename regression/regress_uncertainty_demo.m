clear
Np = 1000;
x = [1:Np]'/Np;
sigma = 0.1;
err = sigma*randn(Np,1);
y = x + err-mean(err);

[b,bint,r,rint] = regress(y,[x,ones(Np,1)],0.6827);

clf;
grid on;hold on;
plot(x,y,'.','MarkerSize',20)

bo = sum((x-mean(x)).*(y-mean(y)))/sum((x-mean(x)).^2);
ao = mean(y) - bo*mean(x);

bo = [bo,ao];
[b';bo]

Vb = sigma^2/sum((x-mean(x)).^2);
Va = sigma^2*(1/Np + mean(x)^2)/sum((x-mean(x)).^2);
binto = [bo(1) - sqrt(Vb),bo(1) + sqrt(Vb);...
         ao(1) - sqrt(Va),ao(1) + sqrt(Va)];

bint
binto

[bint(1,2)-bint(1,1),binto(2,2)-binto(2,1)]

ym = b(1)*x(:,1) + b(2);
ymo = bo(1)*x + bo(2);

plot(x,ymo,'k','LineWidth',2)
