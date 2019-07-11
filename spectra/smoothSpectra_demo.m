N = 1023;
window = ones(floor(N/4),1);
noverlap = 0;

wd = dpss(N,1,64);
wd = wd(:,1)/mean(wd(:,1));

wp = parzenwin(N);
wp = wp/mean(wp);

y = 2*sin(128*pi*[0:N-1]'/N);
y = y + 2*sin(64*pi*[0:N-1]'/N);
yp = y.*wp;
yd = y.*wd;
[Pp,fp] = periodogram([y,yp,yd]);

clf;
semilogy(fp(2:end)/(2*pi),P(2:end,1),'r')
hold on;
semilogy(f(2:end)/(2*pi),P(2:end,2),'g')
semilogy(f(2:end)/(2*pi),P(2:end,3),'b')

legend('Raw','Parzen','DPSS');
semilogy(fp(2:end)/(2*pi),Pp(2:end,1),'r.')
semilogy(f(2:end)/(2*pi),P(2:end,2),'g.')
semilogy(f(2:end)/(2*pi),P(2:end,3),'b.')

plot([128/(2*N),128/(2*N)],[eps,1],'k');

break

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
