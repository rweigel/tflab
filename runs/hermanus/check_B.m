
id = 'KAP103';
% Read input/output data
[B,E,Header] = data_KAP03(id);
start = [strrep(Header.STARTTIME,' ','T'),'.000'];
stop  =  [strrep(Header.ENDTIME,' ','T'),'.000'];
timedelta = str2num(Header.DELTA_T);
dno = datenum([start(1:10),' ',start(12:end-4)],31);
tm = dno + timedelta*[0:size(B,1)-1]'/86400;
Header.STARTTIME
Header.ENDTIME
if 1
    Hermanus = load('check_B_Hermanus');
    t = datetime(Hermanus.data.DateTimeVector(:,1:6));
    X = Hermanus.data.x_parameter2;
    Y = Hermanus.data.x_parameter3;
    Z = Hermanus.data.x_parameter4;
end

addpath([fileparts(mfilename('fullpath')),'/../plot']);

tm = datenum(tm);

figure(1);clf;
    figprep();
    plot(tm,B(:,1),'r');
    hold on;
    plot(datenum(t),X-mean(X),'b');
    grid on;
    legend(sprintf('$B_x$ %s',id),'$B_x$ Hermanus');
    ylabel('[nT]');
    datetick;
    xl = cellstr(get(gca,'XTickLabel'));
    xl{1} = [xl{1},'/',start(1:4)];
    set(gca,'XTickLabel',xl);

    set(gcf,'color','w');
    set(gcf,'defaultFigureColor', [1,1,1]); % Background color to white.
    export_fig('compare_Bx.png','-r200');

figure(2);clf;
    figprep();
    plot(tm,B(:,2),'r');
    hold on;
    plot(datenum(t),Y-mean(Y),'b');
    grid on;
    legend(sprintf('$B_y$ %s',id),'$B_y$ Hermanus');
    ylabel('[nT]');
    datetick;
    xl = cellstr(get(gca,'XTickLabel'));
    xl{1} = [xl{1},'/',start(1:4)];
    set(gca,'XTickLabel',xl);

    set(gcf,'color','w');
    set(gcf,'defaultFigureColor', [1,1,1]); % Background color to white.
    export_fig('compare_By.png','-r200');

figure(3);clf;
    figprep();
    plot(tm,-B(:,1)/20,'r');    
    hold on;
    plot(datenum(tm),E(:,1),'b');

    set(gcf,'color','w');
    set(gcf,'defaultFigureColor', [1,1,1]); % Background color to white.
    export_fig('compare_By.png','-r200');
    