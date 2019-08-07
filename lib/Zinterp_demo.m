%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments passed to interp1
args = {'linear',0};
if ischar(args{2})
    ts = sprintf('; interp1 args: ''%s'', ''%s''',args{:});
else
    ts = sprintf('; interp1 args: ''%s'', %d',args{:});
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,fe] = fftfreq(7);
Z = (1+1j)*[1:length(fe)]';
Z(1) = 3; % Make f = 0 value real (as it must be).
[Zg,fg] = Zinterp(fe,Z,20,args{:});

figure();
    subplot(2,1,1)    
    plot(fe,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['Zinterp(f,Z,N) usage with N even and f w/o 0.5 value',ts]);
    
    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,fe] = fftfreq(7);
Z = (1+1j)*[1:length(fe)]';
Z(1) = 3; % Make f = 0 value real (as it must be).
[Zg,fg] = Zinterp(fe,Z,19,args{:});

figure();
    subplot(2,1,1)    
    plot(fe,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['Zinterp(f,Z,N) usage with N odd and f w/o 0.5 value',ts]);
    
    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,fe] = fftfreq(8);
Z = (1+1j)*[1:length(fe)]';
Z(1) = 3;   % Make f = 0 value real (as it must be).
Z(end) = 2; % Make fe(end) value real (as it must be when fe(end)= 0.5).
[Zg,fg] = Zinterp(fe,Z,20,args{:});

figure();
    subplot(2,1,1)    
    plot(fe,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['Zinterp(f,Z,N) usage with N even and f w/ 0.5 value',ts]);
    
    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,fe] = fftfreq(8);
Z = (1+1j)*[1:length(fe)]';
Z(1) = 3; % Make f = 0 value real (as it must be).
Z(end) = 2; % Make fe(end) value real (as it must be when fe(end)= 0.5).
[Zg,fg] = Zinterp(fe,Z,19,args{:});

figure();
    subplot(2,1,1)    
    plot(fe,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['Zinterp(f,Z,N) usage with N odd and f w/ 0.5 value',ts]);
    
    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');
    
break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f does not have a 0.5 value
% fi has a f = 0.5 value
[~,fe] = fftfreq(7);
[~,fg] = fftfreq(20);

Z = (1+1j)*[1:length(fe)]';
Z(1) = 3; % Make f = 0 value real (as it must be).

Zg = Zinterp(fe,Z,fg,args{:});

figure();
    subplot(2,1,1)    
    plot(fe,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['f does not have 0.5, fi does',ts]);
    
    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fe and fg have a f = 0.5 value
[~,fe] = fftfreq(6);
[~,fg] = fftfreq(20);

Z = (1+1j)*[1:length(fe)]';
Z(1) = 3; % Make fe = 0 value real.
Z(end) = real(Z(end)); % Make fe = 0.5 value real.

Zg = Zinterp(fe,Z,fg,args{:});

figure();
    subplot(2,1,1)
    plot(fe,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['f and fi have 0.5 values',ts]);

    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fe has f = 0 value, fg does not
[~,fe] = fftfreq(7);
[~,fg] = fftfreq(21);
fg = fg(2:end); % Remove fg = 0 value.

Z = (1+1j)*[1:length(fe)]';
Z(1) = 10; % Make fe = 0 value real.

Zg = Zinterp(fe,Z,fg,args{:});

figure();
    subplot(2,1,1)
    plot(fe,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['f has 0 value, fi does not',ts]);

    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f and fi have no 0 value
[~,fe] = fftfreq(7);
[~,fg] = fftfreq(21);
fe = fe(2:end); % Remove fe = 0
fg = fg(2:end); % Remove fg = 0

Z = (1+1j)*[1:length(fe)]';

Zg = Zinterp(fe,Z,fg,args{:});

figure();
    subplot(2,1,1)
    plot(fe,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['f and fi have no 0 value',ts]); 
    
    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f has no 0 value, fi does
[~,fe] = fftfreq(7);
[~,fg] = fftfreq(21);
fe = fe(2:end); % Remove fe = 0 value.

Z = (1+1j)*[1:length(fe)]';

Zg = Zinterp(fe,Z,fg,args{:});

figure();
    subplot(2,1,1)
    plot(fe,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['f has no 0 value, fi does',ts]); 
    
    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');