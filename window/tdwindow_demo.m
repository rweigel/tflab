% See also
% https://ccrma.stanford.edu/~jos/sasp/Kaiser_DPSS_Windows_Compared.html

N = 51;
X = ones(N,1);

Xp = tdwindow(X, @parzenwin);
Xk = tdwindow(X, @kaiser, pi);
Xd = tdwindow(X, @dpss, 1, 1);
Xd = Xd/max(Xd);

figure();
    plot(X,'k','linewidth',2);
    hold on;grid on
    plot(Xp,'linewidth',2)
    plot(Xd,'linewidth',2)
    plot(Xk,'linewidth',2)
    legend('Data',...
            sprintf('parzen(%d)',N),...
            sprintf('dpss(%d,1,1)',N),...
            sprintf('kaiser(%d,pi)',N));
