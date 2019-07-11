T  = 8;
w  = 2*pi/T;
t  = 0:T-1;
x  = sin(w*t);
fb = 1/T;

[x_bp,aib_bp] = bandpass(x,fb);

figure(1);clf
    plot(x,'b','LineWidth',2);hold on;
    plot(x_bp,'g');
    plot((x_bp-x)/eps,'k');
    legend('Original','Bandpassed','Difference/eps');

[x_bp,aib_bp] = bandpass(x,[fb-eps,fb+eps]);

assert(max(abs(x_bp-x))/eps < 2);

T2  = 4;
w2  = 2*pi/T2;
t  = 0:T-1;
x  = x + sin(w2*t);
fb = [1/T-eps,1/T2+eps];

[x_bp,aib_bp] = bandpass(x,fb);
    
assert(max(abs(x_bp-x))/eps < 2);
    
