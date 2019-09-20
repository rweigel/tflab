clear;

N = 1024;
n = 1024/4;

% Raw
F = fftfreqp(N)';
f = fftfreqp(n)';
t = (0:n-1)';
for i = 1:length(F)
    p(i)  = 2*pi*2*(rand(1)-0.5);
    %p(i)  = 0;
    Fu(i) = F(i);
    a(i)  = F(i);
    y(:,i) = a(i)*cos(2*pi*F(i)*t + p(i));
end

yr = sum(y,2);
%yr = yr/sum(yr);

% Parzen
wp = parzenwin(n);
wp = wp-mean(wp);
yp = yr.*wp;
%yp = yp/sum(yp);

% DPSS
NW = 2; % 2, 5/2, 3, 7/2, and 4 are typical vaues
wd = dpss(n,NW);
wd = wd(:,1)-mean(wd(:,1));
yd = yr.*wd;
%yd = yd/sum(yd);

Y = [yr,yp,yd];
P = (2/size(Y,1))*sqrt(fft(Y).*conj(fft(Y)));
P = P(1:length(f),:);

figure(1);clf;
    plot(t,yr);

figure(2);clf;axis square
    plot(Fu,a,'.','MarkerSize',20);
    hold on;grid on;
    plot(f(2:end),P(2:end,1),'r','Marker','.');
    plot(f(2:end),P(2:end,2),'g','Marker','.')
    plot(f(2:end),P(2:end,3),'b','Marker','.')
    legend('Exact','Raw','Parzen','DPSS');
break

window = ones(floor(N/4),1);
noverlap = 0;

sf = N*N/2;

opts = struct();
opts.fd = struct();
opts.fd.evalfreq = struct();
opts.fd.evalfreq.method = 'linear';
opts.fd.evalfreq.options = 0;

opts.fd.window = struct();
opts.fd.window.function = @parzenwin;

[Perr,ferr] = pwelch(err,window,noverlap);
[Perr2,ferr2] = smoothSpectra(err,opts);
[Perr3,ferr3] = pmtm(err);
[Perr4,ferr4] = periodogram(err);
Perr5 = fft(err).*conj(fft(err))/sf;
Perr5 = Perr5(1:length(err)/2+1);
ferr5 = [0:length(err)/2]/length(err);

figure();
clf;
loglog(2*pi./ferr(2:end),Perr(2:end),'LineWidth',2);
hold on;
loglog(1./ferr2(2:end),Perr2(2:end)/sf,'LineWidth',4);
loglog(2*pi./ferr3(2:end),Perr3(2:end),'LineWidth',2);
loglog(2*pi./ferr4(2:end),Perr4(2:end),'LineWidth',2);
loglog(1./ferr5(2:end),Perr5(2:end),'y','LineWidth',2);
grid on;
[lh,lo] = legend('pwelch','smoothSpectra','pmtm','periodogram','fft','Location','Best');
set(lo,'LineWidth',2);
