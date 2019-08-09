addpath([fileparts(mfilename('fullpath')),'/../misc/']);

for N = [9,10,99,100,999,1000];
    B = ones(N,1);
    Z = zeros(N,1);
    Z(1) = 1; % Set DC value to 1.
    
    for k = 3:N
        % Find first period > 2 that is an integer multiple of N.
        if mod(N,k) == 0
            T = k;
            break;
        end
    end

    % B has only DC component. Output should be constant.
    logmsg('Test 1; N = %4d\n',N);
        Ep = zpredict(Z,B);
        mae = max(abs(Ep - Z(1)));
        assert(mae <= eps);
        logmsg('  max abs(error) = %.1e\n',mae);
        
    logmsg('Test 2 N = %4d\n',N);
        % B has only DC component. Output should be constant even though Z is 
        % non-zero at non-DC frequencies.
        Z = ones(N,1);
        Ep = zpredict(Z,B);
        mae = max(abs(Ep - Z(1)));
        assert(mae <= eps);
        logmsg('  max abs(error) = %.1e\n',mae);
        
    logmsg('Test 3; N = %4d\n',N);
        % If Z is all zeros, output should be all zeros.
        Ep = zpredict(zeros(size(B)),B);
        mae = max(abs(Ep - zeros(size(B))));
        logmsg('  max abs(error) = %.1e\n',mae);        
        assert(mae <= eps);
        
    logmsg('Test 4; N = %4d\n',N);
        t = [0:N-1]';
        B = sin(2*pi*t/T);
        Z = abs(fft(B))/(N/2);
        Ep = zpredict(Z,B);
        mae = max(abs(Ep - B));
        logmsg('  T = %d; max abs(error) = %.1e\n',T,mae);
        % TODO: Document why mae is expected to scale with N.
        assert(mae < 10*N*eps);

    logmsg('Test 5; N = %4d\n',N);
        Bs = sin(2*pi*t/T-pi/2); % Shifted B
        Z = fft(B)/(N/2); 
        Ep = zpredict(Z,B);
        mae = max(abs(Ep - Bs));
        logmsg('  T = %d; max abs(error) = %.1e\n',T,mae);
        % TODO: Document why mae is expected to scale with N.        
        assert(mae < 10*N*eps);
        
    logmsg('Test 6; N = %4d\n',N);    
        % Test API. Z and B have same number of columns.
        Z = repmat(Z,1,4);
        Bs = repmat(Bs,1,4);
        Ep4 = zpredict(Z,B);
        tmp = Ep4-repmat(Ep,1,4);
        mae = max(tmp(:));
        logmsg('  max abs(error) = %.1e\n',mae);        
        assert(mae <= eps);
end
