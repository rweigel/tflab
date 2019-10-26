clear;

% NB: Tests are intentially not deterministic. TODO: Make deterministic by
% setting random number seed.

addpath([fileparts(mfilename('fullpath')),'/misc']);

close all;
set(0,'defaultFigureWindowStyle','docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic calculation test 1.
% B = randn(), E = B. With evalfreqs = DFT frequencies, should produce
% perfect predictions b/c # of free parameters in fitted Z equals number of
% data points.
logmsg(['Basic calculation; Test 1.1 - '...
        'B = randn(), E = B. 1 DFT point per freq. band.\n']);

N = [99,100];
for n = N
    B = randn(n,1);
    E = B;

    opts = transferfnFD_options(0);
    S = transferfnFD(B,E,opts);

    % TODO: Justify 10*eps.
    assert(1-S.Metrics.PE < 10*eps);
    assert(S.Metrics.MSE < 10*eps);
    assert(1-S.Metrics.CC < 10*eps);
end
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic calculation test 2.
% B = cos(w*t), E = A(w)*cos(w*t + phi(w)). No leakage
logmsg(['Basic calculation; Test 1.2. - '...
        'B = cos(w*t), E ~ A(w)*cos(w*t + phi(w)). No leakage.\n']);

clear E B
N = 101;
f = fftfreqp(N);
t = (0:N-1)';
for i = 2:length(f)
    A(i,1)   = (i-1);
    Phi(i,1) = -2*pi*f(i);
    B(:,i) = cos(2*pi*f(i)*t);
    E(:,i) = A(i)*cos(2*pi*f(i)*t + Phi(i));
end

B = sum(B,2);
E = sum(E,2);

% Compute DC component
A(1) = mean(E)/mean(B);
if A(1) >= 0
    Phi(1) = 0;
else
    Phi(1) = pi;
end
A(1) = abs(A(1));

opts = transferfnFD_options(0); 
S2 = transferfnFD(B,E,opts);

% Complex version of A
Ac = A.*(cos(Phi) + sqrt(-1)*sin(Phi));

% TODO: Justify 1e-12.
assert(max( real(S2.Z) - real(A) ) < 1e-12 );
assert(max( imag(S2.Z) - imag(A) ) < 1e-12 );
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic calculation test 3.
% B = randn(), H = [1]. evalfreqs = DFT frequencies.

logmsg(['Basic calculation; Test 1.3. - '...
                'H = [1,0,...] with varying # of zeros. '...
                '1 DFT point per freq. band.\n']);

for i = 1:3     
    H = zeros(i+1,1);
    H(1) = 1;
    S0 = transferfnFD_demo_signals(0, struct('H',H,'N',100));

    opts = transferfnFD_options(0);
    S1 = transferfnFD(S0.In, S0.Out, opts);
    
    % Computed H should match used H and be zero for lags longer than
    % used H.
    L = length(S0.H);
    assert(max(abs(S0.H - S1.H(1:L))) <= eps);
    assert(max(abs(S1.H(L+1:end))) <= eps);
    
    % Analytically, real part of Z is 1, imaginary part is 0.    
    re = real(S1.Z)-1; 
    assert(max(abs(re)) <= 1000*eps);
    assert(max(abs(imag(S1.Z))) <= 1000*eps);
    fprintf('---\n');
end
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regression test 1. OLS_REGRESS() using real and complex arguments
% Expect results to be identical to within machine precision.
logmsg(['Basic calculation; Test 2.1. - '...
                'ols_regress() using real and complex arguments.\n']);

B = randn(n,1);
E = B;

opts = transferfnFD_options(0);
opts.fd.regression.functionargs = {struct('realvalued',0)};

% Use regress() with complex matrices (default).
S1 = transferfnFD(B,E,opts);

% Use regress() with only real matrices.
opts.fd.regression.functionargs = {struct('realvalued',1)};
S2 = transferfnFD(B,E,opts);

% TODO: Justify 4*eps.
assert(all(abs(S1.Z(:) - S2.Z(:)) <= 4*eps),'');
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regression test 2. - Compare OLS_REGRESS() with ROBUSTFIT() when no noise.

logmsg(['Regression comparison; Test 2.2 - Compare ols_regress() w/ '...
                'robustfit() and no noise.\n']);
     
N = [99,100];
for n = N
    B = randn(n,1);
    E = B;

    opts = transferfnFD_options(0);

    opts.fd.evalfreq.function = @evalfreq;
    % Can't use 1 DFT point per freq. band because regression needs more
    % points to not have rank deficiency.
    opts.fd.evalfreq.functionstr  = '3 DFT points per freq. band';
    opts.fd.evalfreq.functionargs = {[1,1],'linear'};

    % Uses default regression function OLS_REGRESS().
    S1 = transferfnFD(B,E,opts);

    opts.fd.regression.function = @robust_v1;
    opts.fd.regression.functionargs = {};
    opts.fd.regression.functionstr = ...
                            'Robust regression using robust_v1() function';
    S2 = transferfnFD(B,E,opts);

    assert(S1.Metrics.PE - S2.Metrics.PE <= 2*eps);
    assert(S1.Metrics.CC - S2.Metrics.CC <= 2*eps);
    assert(S1.Metrics.MSE - S2.Metrics.MSE <= 2*eps);

    opts.fd.regression.function = @robust_robustfit;
    opts.fd.regression.functionargs = {[],[],'off'};
    opts.fd.regression.functionstr = ...
                            'Robust regression using robustfit() function';
    S3 = transferfnFD(B,E,opts);

    assert(S1.Metrics.PE - S3.Metrics.PE <= 2*eps);
    assert(S1.Metrics.CC - S3.Metrics.CC <= 2*eps);
    assert(S1.Metrics.MSE - S3.Metrics.MSE <= 2*eps);
end
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% API Test - Multiple Inputs and Ouputs
logmsg('API I/O Test; Test 3.1. - One or Two Outputs, One Input.\n');

N = 1000;
B = randn(N,2);
E = B;

opts = transferfnFD_options(1);
opts.td.window.width = N;
opts.td.window.shift = N;

% 1 input, one or two outputs
S1 = transferfnFD(B(:,1),E(:,1),opts);
fprintf('---\n');
S2 = transferfnFD(B(:,1),[E(:,1),E(:,1)],opts);

assert(all(S1.Predicted == S2.Predicted(:,1)));
assert(all(S1.Predicted == S2.Predicted(:,2)));

%%%
fprintf('\n');
%%%
logmsg('API I/O Test; Test 3.2. - Two Outputs, Two Inputs\n');

% 2 inputs, one or two outputs
S1 = transferfnFD(B,E(:,1),opts);
fprintf('---\n');
S2 = transferfnFD(B,E(:,2),opts);
fprintf('---\n');
S3 = transferfnFD(B,E,opts);

assert(all(S1.Predicted == S3.Predicted(:,1)));
assert(all(S2.Predicted == S3.Predicted(:,2)));
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% API Test - Segmenting
% E and B are split into segments and transfer functions are computed for
% each segment.
logmsg('API Segmenting; Test 3.3.\n');

N = 1000;
B = randn(N,1);
E = B;

opts = transferfnFD_options(1);
opts.td.window.width = N;
opts.td.window.shift = N;

S1 = transferfnFD(B,E,opts);
fprintf('---\n');
S2 = transferfnFD([B;B],[E;E],opts);

% Results for two segments in S2 should be same a single segment in S1.
assert(all(S1.Predicted == S2.Segment.Predicted(:,1,1)));
assert(all(S1.Predicted == S2.Segment.Predicted(:,1,2)));
assert(all(S1.Z == S2.Segment.Z(:,:,1)))
assert(all(S1.Z == S2.Segment.Z(:,:,2)))
assert(all(S1.Z(:) == S2.Z(:)));

%%%
fprintf('\n');
%%%
logmsg('API Segmenting; Test 3.4.\n');

N = 1000;
B = randn(N,2);
E = B;

opts = transferfnFD_options(1);
opts.td.window.width = N;
opts.td.window.shift = N;

S1 = transferfnFD(B(:,1),E(:,1),opts);
fprintf('---\n');
S2 = transferfnFD([B(:,1);B(:,1)],[E;E],opts);

assert(all(S1.Predicted == S2.Segment.Predicted(:,1,1)));
assert(all(S1.Predicted == S2.Segment.Predicted(:,1,2)));

fprintf('\n');
%%%
logmsg('API Segmenting; Test 3.5.\n');

S3 = transferfnFD(B,E,opts);
S4 = transferfnFD([B;B],[E;E],opts);
assert(all(S3.Predicted(:,1) == S4.Segment.Predicted(:,1,1)));
assert(all(S3.Predicted(:,2) == S4.Segment.Predicted(:,2,1)));
assert(all(S3.Predicted(:,1) == S4.Segment.Predicted(:,1,2)));
assert(all(S3.Predicted(:,2) == S4.Segment.Predicted(:,2,2)));
assert(all(S3.Z(:) == S4.Z(:)));
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% API Intervals
% When there are gaps in time in the input/data, one can pass a cell array
% of intervals and then the transfer function is computed on each interval.
% The intervals may be segemented by specifying a window width and window
% shift that is less than the interval length.
logmsg('API Intervals; Test 3.6.\n');

N = 1000;
B = randn(N,1);
E = B;

opts = transferfnFD_options(1);
opts.td.window.width = N;
opts.td.window.shift = N;

S1 = transferfnFD([B;B],[E;E],opts);
fprintf('---\n');
S2 = transferfnFD({B;B},{E;E},opts);
assert(all(S1.Z(:) == S2.Z(:)));
assert(all(S1.Segment.Predicted(:) == S2.Segment.Predicted(:)))

%%%
fprintf('\n');
%%%
logmsg('API Intervals; Test 3.7.\n');

S1 = transferfnFD(B,E,opts);
fprintf('---\n');
S2 = transferfnFD({B,[B;B]},{E,[E;E]},opts);
fprintf('---\n');
S3 = transferfnFD({B,[B;B]},{0.5*E,[1.0*E;1.5*E]},opts);

assert(all(S1.Predicted == S2.Segment.Predicted(:,:,1)))
assert(all(S1.Predicted == S2.Segment.Predicted(:,:,2)))
assert(all(S1.Predicted == S2.Segment.Predicted(:,:,3)))
% DC value will be different for S3, so omit from test:
assert(all(S1.Z(2:end)-S3.Z(2:end) < 10*eps));

% Average TF should be 1.0 for all fe, same as S1.Z.
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% API Test - Stack Regression
% When intervals and/or segments are used, the default is to compute a
% transfer function that is the average of each segment. 
logmsg('API Stack Regression; Test 3.8.\n');

N = 1000;
B = randn(N,2);
E = B;

% 1 input/1 output. When using 1 segment, non-stack average should be same
% as stack average result
opts = transferfnFD_options(1);
% The following two lines are not needed as this is default for behavior when
% transferfnFD_options(1)
opts.td.window.width = N; 
opts.td.window.shift = N; 

S1 = transferfnFD(B(:,1),E(:,1),opts);
fprintf('---\n');
opts.fd.stack.average.function = ''; % Don't compute stack average.
S2 = transferfnFD(B(:,1),E(:,1),opts);
assert(all(S1.Z(:) == S2.Z(:)))

fprintf('\n');
%%%
logmsg('API Stack Regression; Test 3.9.\n');

% 2 inputs/2 outputs. When using 1 segment, stack regression should be
% same as stack average result
opts = transferfnFD_options(1);
% The following two lines are not needed as this is default for behavior when
% transferfnFD_options(1)
opts.td.window.width = N;
opts.td.window.shift = N;

S1 = transferfnFD(B,E,opts);
opts.fd.stack.average.function = ''; % Don't compute stack average.
fprintf('---\n');
S2 = transferfnFD(B,E,opts);
assert(all(S1.Z(:) == S2.Z(:)))

fprintf('\n');
%%%
logmsg('API Stack Regression; Test 3.10.\n');

% Compare stack average Z to stack regression Z. Results not expected to be
% identical. For the stack average method, Z for each segment in a given
% frequency band is computed by regressing on the DFTs in that band segment
% Z values are averaged. For the stack regression method, the DFTs in a
% given frequency band are computed for each segment and then the segment
% frequency band DFTs are combined and a single regression is performed.
% DFTs.
N = 1000;
B = randn(N,1);
E = B;

opts = transferfnFD_options(1);
opts.td.window.width = N; % Will result in two intervals.
opts.td.window.shift = N; % Will result in two intervals.

S3 = transferfnFD([B;B],[E;E],opts);
fprintf('---\n');
opts.fd.stack.average.function = '';
S4 = transferfnFD([B;B],[E;E],opts); % window.width and window.shift ignored.
assert(all(abs(S3.Z(:) - S4.Z(:)) < 10*eps))

fprintf('\n');
%%%
logmsg('API Stack Regression; Test 3.11.\n');

% Verify that get same answer when continuous and discontinuous segments
% are used. Expecet identical results.
opts = transferfnFD_options(1);
opts.td.window.width = N;
opts.td.window.shift = N;
opts.fd.stack.average.function = '';

S3 = transferfnFD([B;B],[E;E],opts);
fprintf('---\n');

% Each element of cell array is treated as having gap in time stamps.
S4 = transferfnFD({B,B},{E,E},opts);
assert(all(S3.Z(:) == S4.Z(:)))

%assert(all(S1.Z == S2.Z));
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logmsg('transferfnFD_test.m: All tests passed.\n');
