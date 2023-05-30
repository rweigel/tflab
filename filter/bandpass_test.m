T  = 8;
w  = 2*pi/T;
t  = 0:T-1;
x  = sin(w*t);
fb = 1/T;

%% 1. Keep a single frequency
[x_bp,aib_bp] = bandpass_(x,fb);
assert(max(abs(x_bp-x))/eps <= 2);

figure(1);clf
    title('Test 1');
    plot(x,'k','LineWidth',4);hold on;
    plot(x_bp,'y','LineWidth',2);
    plot((x_bp-x)/eps,'k');
    legend('Original','Bandpassed','Difference/eps');


%% 2. Keep two frequencies using range
T2  = 4;
w2  = 2*pi/T2;
t  = 0:T-1;
x  = x + sin(w2*t);
fb = [1/T-eps,1/T2+eps];

[x_bp,aib_bp] = bandpass_(x,fb);
assert(max(abs(x_bp-x))/eps < 2);

figure(2);clf
    title('Test 2');
    plot(x,'k','LineWidth',4);hold on;
    plot(x_bp,'y','LineWidth',2);
    plot((x_bp-x)/eps,'k');
    legend('Original','Bandpassed','Difference/eps');


% 3. Power law
S = demo_signals('powerlaw')

S = tflab_preprocess(S);

dock on;figure(1);close all;

figure();
    tsplot(S);

figure();
    % Populate S with data needed for psdplot().
    [PSDIn,~,f] = psd(S.In);
    S.Metrics.PSD.Raw.fe = f;
    S.Metrics.PSD.Raw.In = PSDIn;

    PSDOut = psd(S.Out);
    S.Metrics.PSD.Raw.Out = PSDOut;
    psdplot(S, struct('type','raw'));
    