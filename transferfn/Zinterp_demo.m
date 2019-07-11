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
% fe does not have a f = 0.5 value
% fg has a f = 0.5 value
[~,fe] = fftfreq(7);
[~,fg] = fftfreq(20);

Z = (1+1j)*[1:length(fe)]';
Z(1) = 3; % Make fe = 0 value real.

Zg = Zinterp(fe,Z,fg,args{:});

figure();
    subplot(2,1,1)    
    plot(fe,real(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,real(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Re(Zg)');
    legend('Given','Interpolated','Location','Best');
    title(['fe does not have 0.5, fg does',ts]);
    
    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('Given at fe','Interpolated to fg','Location','Best');

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
    legend('Given','Interpolated','Location','Best');
    title(['fe and fe have 0.5 values',ts]);

    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('Given at fe','Interpolated to fg','Location','Best');

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
    legend('Given','Interpolated','Location','Best');
    title(['fe has 0 value, fg does not',ts]);

    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('Given at fe','Interpolated to fg','Location','Best');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fe and fg have no f = 0 value
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
    legend('Given','Interpolated','Location','Best');
    title(['fe and fg have no 0 value',ts]); 
    
    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('Given at fe','Interpolated to fg','Location','Best');
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fg has no f = 0 value, fe does
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
    legend('Given','Interpolated','Location','Best');
    title(['fe has no 0 value, fg does',ts]); 
    
    subplot(2,1,2)    
    plot(fe,imag(Z),'k.','MarkerSize',30);
    hold on;grid on;box on;
    plot(fg,imag(Zg),'g.','MarkerSize',20);
    xlabel('f');
    ylabel('Im(Zg)');
    legend('Given at fe','Interpolated to fg','Location','Best');