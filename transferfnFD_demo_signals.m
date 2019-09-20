function S = transferfnFD_demo_signals(tn, opts)

addpath([fileparts(mfilename('fullpath')),'/lib']);

assert(nargin == 2,'Two inputs are required');

if tn == -2
    N = opts.N;
    addpath([fileparts(mfilename('fullpath')),'/fft']);
    addpath([fileparts(mfilename('fullpath')),'/lib']);

    [~,F] = fftfreq(opts.N);
    F = F';
    t = (0:opts.n-1)';
    for i = 2:length(F)
        A(i,1)   = F(i);
        Phi(i,1) = -2*pi*F(i);
        B(:,i) = cos(2*pi*F(i)*t);
        E(:,i) = A(i)*cos(2*pi*F(i)*t + Phi(i));
        %E(:,i) = A(i)*cos(2*pi*F(i)*t);
    end

    B = sum(B,2);
    E = sum(E,2);

    %B = B/max(B);
    %E = E/max(E);

    [~,f] = fftfreq(opts.n);
    f = f';
    for i = 2:length(f)
        A(i,1)   = f(i);
        Phi(i,1) = -2*pi*f(i);
        Z(i,1) = A(i)*(cos(Phi(i)) + sqrt(-1)*sin(Phi(i)));
    end
    
    S.In  = B;
    S.Out = E;
    S.Z  = Z;
    S.fe = f;
    S.H = z2h(zinterp(f,Z,opts.n));
    S.tH = (0:opts.n-1)';
    return
end

if tn == -1
    N = opts.N;
    addpath([fileparts(mfilename('fullpath')),'/fft']);
    addpath([fileparts(mfilename('fullpath')),'/lib']);
    [~,f] = fftfreq(N);
    H = opts.H;
    f = f';
    [z,w] = freqz(opts.H,1,f,1);
    [Z,fi] = zinterp(f,z,opts.N);    

    B = randn(opts.N,1);
    E = zpredict(Z,B);

    S.In  = B;
    S.Out = E;
    S.Z  = z;
    S.fe = f;
    S.H = opts.H;
    S.tH = (0:size(opts.H,1)-1)';
    return
end

if tn == 0 % Prescribed impulse response
    description = '';

    H = opts.H;
    N = opts.N;
    
    B = randn(opts.N+length(H),1);
    E = filter(H,1,B);
    
    % Remove non-steady-state
    B = B(length(H)+1:end); 
    E = E(length(H)+1:end);
end

if tn == 1 % Low pass filter impulse responses

    dim = opts;
    
    description = '';
    
    tau  = 10;  % Filter decay constant
    Ntau = 100; % Number of filter coefficients = Ntau*tau + 1
    N    = 2e4; % Simulation length
    nR   = 2;   % Width of rectangualar window is 2*nR+1
    Nss  = 4;   % Will remove Nss*Ntau*tau from start of all time series 
    nb   = 0.0; % Noise in B
    ne   = 0.0; % Noise in E
    ndb  = 0.0; % Noise in dB

    % IRF for dx/dt + x/\tau = delta(0), and ICs x(t=0) = 0; dx/dt|_{t=0} = 0
    % approximated using forward Euler.
    dt = 1;
    gamma = (1-dt/tau);
    for i = 1:Ntau*tau
        h(i,1) = gamma^(i-1);
        t(i,1) = dt*(i-1);
    end
    hstr = sprintf('(1-1/%d)^{t}; t=1 ... %d; h_{xy}(0)=0;', tau, length(h));
    % Add a zero because MATLAB filter function requires it.
    % (Not having it also causes phase drift with frequency.)
    h  = [0;h];            % Exact IRF

    % Add extra values so length is same after cutting off non-steady state
    % part of E and B.
    N = N + Nss*length(h);

    % Noise
    NE  = [ne*randn(N,1),ne*randn(N,1)];
    NB  = [nb*randn(N,1),nb*randn(N,1)];
    NdB = [ndb*randn(N,1),ndb*randn(N,1)];

    % Create signals
    if dim == 1
        H(:,1) = h;
        B(:,1) = randn(N,1);
        E(:,1) = NE(:,1) + filter(H(:,1),1,B(:,1) + NB(:,1));
    end
    if dim == 2
        H = zeros(length(h),dim);
        if tn == 1 % Doing all equal will lead to rank deficient warnings.
            H(:,1) = 0.1*h;
            H(:,2) = 0.2*h;
        end
        if tn == 2
            H(:,1) = 0*h;
            H(:,2) = h;
        end
        B(:,1) = randn(N,1);
        B(:,2) = randn(N,1);
        E(:,1) = NE(:,1) + ...
                    + filter(H(:,1),1,B(:,1) + NB(:,1)) ...
                    + filter(H(:,2),1,B(:,2) + NB(:,2));
    end
    if dim == 4
        H = zeros(length(h),dim);
        if tn == 1 
            % Doing all equal will lead to rank deficient warnings.
            H(:,1) = 0.1*h;
            H(:,2) = 0.2*h;
            H(:,1) = 0.3*h;
            H(:,2) = 0.4*h;
        end
        if tn == 2
            % Need to explain why Zxx and Zyy are not zero for this case.
            H(:,1) = 0*h; 
            H(:,2) = h;
            H(:,3) = 0*h;
            H(:,4) = h;
        end
        B(:,1) = randn(N,1);
        B(:,2) = randn(N,1);
        E(:,1) = NE(:,1) + filter(H(:,1),1,B(:,1) + NB(:,1)) ...
                         + filter(H(:,2),1,B(:,2) + NB(:,2));
        E(:,2) = NE(:,2) + filter(H(:,3),1,B(:,1) + NB(:,1)) ...
                         + filter(H(:,4),1,B(:,2) + NB(:,2));
    end

    % Remove non-steady-state part of signals
    B  = B(Nss*length(h)+1:end,:);
    E  = E(Nss*length(h)+1:end,:);

    NE  = NE(Nss*length(h)+1:end,:);
    NB  = NB(Nss*length(h)+1:end,:);

end

% Compute exact transfer function
[Z,f] = h2z(H);
Z(Z==0) = eps; % So points show up on loglog plot.

S.In  = B;
S.Out = E;
S.Z  = Z;
S.fe = f';
S.H = H;
S.tH = (0:size(H,1)-1)';
%S.H  = [H ; 0*H(2:end)];         % Convert H to standard form
%S.tH = [tH ; -fliplr(tH(2:end));]; % Create negative time values
S.Options.description = description;
