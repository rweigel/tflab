N = 100;
B = ones(N,1);
Z = zeros(N,1);
Z(1) = 1;

% B has only DC component. Output should be constant.
assert(all(Zpredict(Z,B) == ones(N,1)));

% B has only DC component. Output should be constant and
% independent on non-zero Z frequencies.
Z = ones(N,1);
assert(all(Zpredict(Z,B) == ones(N,1)));

assert(all(Zpredict(0*Z,B) == zeros(N,1)));

T = 5;
t = [0:N-1]';
B = sin(2*pi*t/T);
Z = abs(fft(B))/(N/2); 
Ep = Zpredict(Z,B);

% Z is non-zero and real at same freq. as non-zero freq. in B.
% Output should be same as input.
assert(all(abs(Ep - B) < 100*eps));
 
Bs = sin(2*pi*t/T-pi/2); % Shifted B
Z = fft(B)/(N/2); % Z has complex at same freq. B is non-zero.
Ep = Zpredict(Z,B);

assert(all(abs(Zpredict(Z,B) - Bs) < 100*eps));

Z = repmat(Z,1,4);
Bs = repmat(Bs,1,4);
Ep4 = Zpredict(Z,B);
tmp = Ep4-repmat(Ep,1,4);
assert(all(tmp(:) == 0))

clf;hold on;grid on;
plot(Bs,'-','Marker','.','MarkerSize',10);
plot(Ep);
grid on;
legend('B','Ep');

