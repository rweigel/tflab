

%S = demo_signals('powerlaw')
D = demo_signals('simple')

dock on;figure(1);close all;

opts = tflab_options(0);
opts.td.window.width = 100;
opts.td.window.shift = 100;
opts.td.window.function = @tdwindow;
opts.td.window.functionstr = 'Rectangle';
opts.td.window.functionargs = {@rectwin};

S = struct('In',D.In,'Out',D.Out,'Options',opts);

S = tflab_preprocess(S);
S.Z = D.Z;
S.fe = D.fe;

S = tflab_metrics(S);

figure();
    tsplot(S);

figure();
    tsplot(S, struct('type','error'));
    
figure();
    dftplot(S, struct('type','original'));

figure();
    dftplot(S, struct('type','original-averaged'));

figure();
    dftplot(S, struct('type','error'));

figure();
    dftplot(S, struct('type','error-averaged'));
    
figure();
    snplot(S);

