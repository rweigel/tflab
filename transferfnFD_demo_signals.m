function S = transferfnFD_demo_signals(tn, dim)

assert(nargin == 2,'Two inputs (test number and dimension) are required');

if tn == 0 % Delta function impulse responses
    description = '';

    H = zeros(100,1);
    H(1) = 1;
    B = rand(2*length(H)+1,1);
    E = filter(H,1,B);
    % Remove non-steady-state (technically, could just do 2:end b/c only H(1)
    % is non-zero)
    B = B(length(H):end); 
    E = E(length(H):end);
end

if tn == 1 % Low pass filter impulse responses
    
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
[Z,f] = H2Z(H);
Z(Z==0) = eps; % So points show up on loglog plot.

tH = (0:size(H,1)-1)'; % Exact IRF time lags

S.In  = B;
S.Out = E;
S.Z  = Z;
S.fe = f;
S.H  = [H ; NaN*H(2:end)];         % Convert H to standard form
S.tH = [tH ; -fliplr(tH(2:end));]; % Create negative time values
S.Options.description = description;
