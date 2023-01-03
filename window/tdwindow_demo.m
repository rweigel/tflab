N = 51;
X = ones(N,1);
clf;
plot(X,'k','linewidth',2);
hold on;
d = tdwindow(X, @dpss, 1, 1);
d = d/max(d);
plot(tdwindow(X, @parzenwin),'linewidth',2)
plot(d,'linewidth',2)
plot(tdwindow(X, @kaiser, pi),'linewidth',2)
legend('Data',...
        sprintf('parzen(%d)',N),...
        sprintf('dpss(%d,1,1)',N),...
        sprintf('kaiser(%d,pi)',N));
