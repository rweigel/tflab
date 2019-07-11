clear
N = 101;
f = fftfreqp(N);
t = (0:N-1)';
for i = 1:length(f)
    A(i)   = (i-1);
    Phi(i) = -2*pi*f(i);
    B(:,i) = cos(2*pi*f(i)*t);
    E(:,i) = A(i)*cos(2*pi*f(i)*t+Phi(i));
end

B = sum(B,2);
E = sum(E,2);

opts = transferfnFD_options(0);
S2 = transferfnFD(B,E,opts)

figure(1)
    plot(S2.PSD.Out);
break
figure(2);clf
    plot(S2.fe,(180/pi)*S2.Phi,'.');
    grid on;
figure(1);clf
    plot(S2.fe,abs(S2.Z),'.');
    grid on;
    