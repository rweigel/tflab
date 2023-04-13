figure(1)
subplot(2,1,1)
    y = [1,1e5];
    plot(y,y,'k*');
    grid on;
subplot(2,1,2)
    plot(y,y,'k*');
    adjust_exponent();
    grid on;

if 0    
figure(2)
subplot(2,1,1)
    y = [1,1e5];
    semilogy(y,y,'k*');
    grid on;
subplot(2,1,2)
    semilogy(y,y,'k*');
    adjust_exponent();
    grid on;

figure(3)
subplot(2,1,1)
    y = [1,1e5];
    semilogx(y,y,'k*');
    grid on;
subplot(2,1,2)
    semilogx(y,y,'k*');
    adjust_exponent();
    grid on;
    
figure(3)
subplot(2,1,1)
    y = [-0.000001, 0,0.000001];
    plot([0,1,2],y,'k*');
    grid on;
subplot(2,1,2)
    plot([0,1,2],y,'k*');
    adjust_exponent();
    grid on;

figure(4)
subplot(2,1,1)
    y = [-1,-1.000001,-1.000002];
    plot([0,1,2],y,'k*');
    grid on;
subplot(2,1,2)
    plot([0,1,2],y,'k*');
    adjust_exponent();
    grid on;

figure(5)
subplot(2,1,1)
    y = [-0.000001, 0,0.000001];
    plot([0,1,2],y,'k*');
    grid on;
    set(gca,'YLim',[1,1+40e-7])    
subplot(2,1,2)
    plot([0,1,2],y,'k*')
    grid on;
    adjust_exponent();
    set(gca,'YLim',[1,1+40e-7])
end