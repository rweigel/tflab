figure(1);clf;
y = [1+eps,1.00000,1-eps];
subplot(2,1,1)
    plot([0,1,2],y,'k.','MarkerSize',20);
    grid on;
subplot(2,1,2)
    plot([0,1,2],y,'k.','MarkerSize',20);
    adjust_ylim();
    grid on;

figure(2);clf;
y = [1+eps,1.00000,1-eps];
subplot(2,1,1)
    semilogy([0,1,2],y,'k.','MarkerSize',20);
    grid on;
subplot(2,1,2)
    semilogy([0,1,2],y,'k.','MarkerSize',20);
    adjust_ylim();
    grid on;
    
figure(3);clf;
y = [1,1e1,1e3];
subplot(2,1,1)
    plot([0,1,2],y,'k.','MarkerSize',20);
    grid on;
subplot(2,1,2)
    plot([0,1,2],y,'k.','MarkerSize',20);
    adjust_ylim();
    grid on;

figure(4);clf;
y = [-eps,0.00000,eps];
subplot(2,1,1)
    plot([0,1,2],y,'k.','MarkerSize',20);
    grid on;
subplot(2,1,2)
    plot([0,1,2],y,'k.','MarkerSize',20);
    adjust_ylim();
    grid on;

figure(5);clf;
y = [1e-16, 1e-10, 1];
subplot(2,1,1)
    semilogy([0,1,2],y,'k.','MarkerSize',20);
    grid on;
subplot(2,1,2)
    semilogy([0,1,2],y,'k.','MarkerSize',20);
    adjust_ylim();
    grid on;
