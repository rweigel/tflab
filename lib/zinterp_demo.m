%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arguments passed to interp1
%args = {'pchip','extrap'};
%args = {'linear','extrap'};
args = {'linear',0}; % The default
if ischar(args{2})
    ts = sprintf('; interp1 args: ''%s'', ''%s''',args{:});
else
    ts = sprintf('; interp1 args: ''%s'', %d',args{:});
end
opts = struct('interp1args',{args});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,f] = fftfreq(7);
Z = (1+1j)*[1:length(f)]';
Z(1) = 3; % Make f = 0 value real (as it must be).
[Zg,fi] = zinterp(f',Z,20,opts);

figure();
    subplot(2,1,1)    
    plot(f,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['zinterp(f,Z,N) usage with N even and f w/o +0.5 value',ts]);
    
    subplot(2,1,2)    
    plot(f,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,f] = fftfreq(7);
Z = (1+1j)*[1:length(f)]';
Z(1) = 3; % Make f = 0 value real (as it must be).
[Zg,fi] = zinterp(f',Z,19,opts);

figure();
    subplot(2,1,1)    
    plot(f,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['zinterp(f,Z,N) usage with N odd and f w/o +0.5 value',ts]);
    
    subplot(2,1,2)    
    plot(f,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,f] = fftfreq(8);
Z = (1+1j)*[1:length(f)]';
Z(1) = 3;   % Make f = 0 value real (as it must be).
Z(end) = 2; % Make fe(end) value real (as it must be when fe(end)= +0.5).
[Zg,fi] = zinterp(f',Z,20,opts);

figure();
    subplot(2,1,1)    
    plot(f,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['zinterp(f,Z,N) usage with N even and f w/ 0.5 value',ts]);
    
    subplot(2,1,2)    
    plot(f,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,f] = fftfreq(8);
Z = (1+1j)*[1:length(f)]';
Z(1) = 3;   % Make f = 0 value real (as it must be).
Z(end) = 2; % Make fe(end) value real (as it must be when fe(end)= 0.5).
[Zg,fi] = zinterp(f',Z,19,opts);

figure();
    subplot(2,1,1)    
    plot(f,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['zinterp(f,Z,N) usage with N odd and f w/ 0.5 value',ts]);
    
    subplot(2,1,2)    
    plot(f,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f does not have a 0.5 value
% fi has a f = 0.5 value
[~,f] = fftfreq(7);
[~,fi] = fftfreq(20);

Z = (1+1j)*[1:length(f)]';
Z(1) = 3; % Make f = 0 value real (as it must be).

Zg = zinterp(f',Z,fi',opts);

figure();
    subplot(2,1,1)    
    plot(f,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['zinterp(f,Z,N) usage; f does not have +0.5, fi does',ts]);
    
    subplot(2,1,2)    
    plot(f,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f and fi have a f = 0.5 value
[~,f] = fftfreq(6);
[~,fi] = fftfreq(20);

Z = (1+1j)*[1:length(f)]';
Z(1) = 3; % Make fe = 0 value real.
Z(end) = real(Z(end)); % Make fe = 0.5 value real.

Zg = zinterp(f',Z,fi',opts);

figure();
    subplot(2,1,1)
    plot(f,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['zinterp(f,Z,N) usage; f and fi have +0.5 values',ts]);

    subplot(2,1,2)    
    plot(f,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fe has f = 0 value, fi does not
[~,f] = fftfreq(7);
[~,fi] = fftfreq(21);
fi = fi(2:end); % Remove fi = 0 value.

Z = (1+1j)*[1:length(f)]';
Z(1) = 10; % Make f = 0 value real.

Zg = zinterp(f',Z,fi',opts);

figure();
    subplot(2,1,1)
    plot(f,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['zinterp(f,Z,N) usage; f has 0 value, fi does not',ts]);

    subplot(2,1,2)    
    plot(f,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f and fi have no 0 value
[~,f] = fftfreq(7);
[~,fi] = fftfreq(21);
f = f(2:end); % Remove fe = 0
fi = fi(2:end); % Remove fg = 0

Z = (1+1j)*[1:length(f)]';

Zg = zinterp(f',Z,fi',opts);

figure();
    subplot(2,1,1)
    plot(f,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['zinterp(f,Z,N) usage; f and fi have no 0 value',ts]); 
    
    subplot(2,1,2)    
    plot(f,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f has no 0 value, fi does
[~,f] = fftfreq(7);
[~,fi] = fftfreq(21);
f = f(2:end); % Remove fe = 0 value.

Z = (1+1j)*[1:length(f)]';

Zg = zinterp(f',Z,fi',opts);

figure();
    subplot(2,1,1)
    plot(f,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('f','fi','Location','Best');
    title(['zinterp(f,Z,N) usage; f has no 0 value, fi does',ts]); 
    
    subplot(2,1,2)    
    plot(f,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fi,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('f','fi','Location','Best');