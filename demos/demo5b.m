%% Transfer function is time shift
% E(t+1) = B(t)

clear;
addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
tflab_setpaths();

N = 11;     % Number of input/output points
H = [1,0];   % Should get exact H back from Z

N = 1001;     % Number of input/output points
H = [0,1];   % Will not get exact H back from Z, but better w/ incr. N.

In = randn(N+length(H),1);
Out = filter(H,1,In);

% Remove non-steady-state
In = In(length(H):end);
Out = Out(length(H):end); 

opts = tflab_options(0);
S1 = tflab(In,Out,opts);

dock on;figure(1);close all;

figure()
    tsplot(S1,struct('type','original'));

figure()
    tsplot(S1,struct('type','error'));

figure()
    dftplot(S1,struct('type','original'));

figure()
    zplot(S1);

figure()
    hplot(S1,[-5,5]);
