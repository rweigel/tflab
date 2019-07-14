clear;
N = 101;
f = fftfreqp(N);
t = (0:10*N-1)';
for i = 1:length(f)
    A(i)   = (i-1);
    phio(i) = 2*pi*rand(1);
    Phi(i) = phio(i)-2*pi*f(i);
    B(:,i) = cos(2*pi*f(i)*t+phio(i));
    E(:,i) = A(i)*cos(2*pi*f(i)*t+Phi(i)) + 0.1*A(i)*randn(length(t),1);
end

B = sum(B,2);
E = sum(E,2);

opts = transferfnFD_options(0);
%opts.fd.stack.average.function = '';
opts.td.window.width = N;
opts.td.window.shift = N;
S2 = transferfnFD(B,E,opts)

if 1
    figure(1);clf
        plot(S2.fe,abs(S2.Z),'.');
        grid on;
        break
    figure(3);clf
        plot(S2.fe,sqrt(S2.PSD.Out),'.');
        hold on;
        
end
break
    figure(2);clf
        plot(S2.fe,(180/pi)*S2.Phi,'.');
        grid on;
