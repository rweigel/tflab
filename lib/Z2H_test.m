%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One input argument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[h,t] = Z2H([1,1,1]);
assert(all(h == [1,0,0]));
assert(all(t == [0,1,-1]));

[h,t] = Z2H([1,1,1]');
assert(all(h == [1,0,0]'));
assert(all(t == [0,1,-1]'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Two input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constant in frequency domain is impulse in time domain.
% Even # of points
N = 6;
[~,fu] = fftfreq(N);
Z = ones(1,N);
Zp = Z(1:length(fu)); % Make z correspond to positive frequencies (with f=-0.5 mapped to f=0.5).
assert(all(Z2H(Zp,fu,N) == [1,0,0,0,0,0]));
assert(all(Z2H(Zp,fu,N) == Z2H(Z)));

% Odd # of points
N = 7;
[~,fu] = fftfreq(N);
Z = ones(1,N);
Zp = Z(1:length(fu));
% With odd # of points, ifft(ones(1,N)) has machine precision error.
assert(all(abs(Z2H(Zp,fu,N) - [1,0,0,0,0,0,0]) < eps));
assert(all(Z2H(Zp,fu,N) == Z2H(Z)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
