figure(1)

subplot(2,1,1)
    y = [1,1.000001,1.000002];
    plot([0,1,2],y,'*');
    grid on;
subplot(2,1,2)
    plot([0,1,2],y,'*')
    adjust_yticks()
    grid on;

    
figure(2)
subplot(2,1,1)
    y = [-1,-1.000001,-1.000002];
    plot([0,1,2],y,'*');
    grid on;
subplot(2,1,2)
    plot([0,1,2],y,'*')
    adjust_yticks()
    grid on;

figure(3)
subplot(2,1,1)
    y = [-0.000001, 0,0.000001];
    plot([0,1,2],y,'*');
    grid on;
subplot(2,1,2)
    plot([0,1,2],y,'*')
    adjust_yticks()
    grid on;

figure(4)
subplot(2,1,1)
    y = [-1,-1.000001,-1.000002];
    plot([0,1,2],y,'*');
    grid on;
subplot(2,1,2)
    plot([0,1,2],y,'*')
    adjust_yticks()
    grid on;

figure(5)
subplot(2,1,1)
    y = [-0.000001, 0,0.000001];
    plot([0,1,2],y,'*');
    grid on;
    set(gca,'YLim',[1,1+40e-7])    
subplot(2,1,2)
    plot([0,1,2],y,'*')
    grid on;
    adjust_yticks()
    set(gca,'YLim',[1,1+40e-7])
    
    