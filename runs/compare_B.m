
id = 'KAP103';
% Read input/output data
[B,E,H] = KAP03_data(id);
start = [strrep(H.STARTTIME,' ','T'),'.000'];
stop  =  [strrep(H.ENDTIME,' ','T'),'.000'];
timedelta = str2num(H.DELTA_T);
dno = datenum([start(1:10),' ',start(12:end-4)],31);
tm = dno + timedelta*[0:size(B,1)-1]'/86400;

if 0
    [D,M]=hapi('http://localhost:9998/INTERMAGNET/hapi','minute/definitive/her','','2003-11-08','2003-12-06');

    fid = fopen('data/Hermanus/id-minute_definitive_her_parameters-undefined_time.min-2003-11-08Z_time.max-2003-12-06.csv','r');
    D = textscan(fid,'%s%d%f%f%f%f','CollectOutput',1,'Delimiter',',');
    fclose(fid);
    th = datenum(double(D.DateTimeVector(:,1:6)));    
end

addpath([fileparts(mfilename('fullpath')),'/../plot']);

figure(1);clf;
figprep();
plot(tm,B(:,1),'r');
hold on;
plot(th,D.H-mean(D.H(1:1440)),'b');
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
plot(th,D.D-mean(D.D(1:1440)),'b');
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
