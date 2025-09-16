addpath([fileparts(mfilename('fullpath')),'/m_map']);

M = m_shaperead('countries/ne_110m_admin_0_countries')

clf
%m_proj('miller','long',[12 32],'lat',[-35 -25]);
m_proj('miller','long',[14 30],'lat',[-35 -28]);
m_etopo2('shadedrelief','gradient',3);
hold on
%m_coast('patch',[1 .85 .7]);
%m_elev('contourf',[500:500:6000]);
m_grid('box','fancy','tickdir','in');
for k=1:length(M.ncst)
     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color','k','linewidth',1);
     %m_patch(M.ncst{k}(:,1),M.ncst{k}(:,2),'y')
end

hold on

cape_fold = load('Nguuri_2001_Cape_Fold.txt');

m_line(20 + cape_fold(:,1),-32 + cape_fold(:,2),...
        'Color','k','LineWidth',2,'LineStyle','--')
m_text(22,-33.3,'Cape Fold Belt','FontSize',16,'FontWeight','Bold')
% kap103.edi has (no indication of geomagnetic or geographic)
%  REFLAT=  -32:07:48 (-32.1300)
%  REFLONG= 020:27:35 ( 46.4597)
% However, this must be geographic b/c if the above is geomagnetic, then
% geographic is 41.440S	49.597W, which is in the Atlantic Ocean.
%
%https://www.mtnet.info/data/kap03/kap03.html has (and "The data are in
%GEOMAGNETIC co-ordinates")
% lat = -32.14, long = 20.47
% Which differs from kap103.edi by 0.01 degrees in both lat and long.
kap103 = [-32 - 8/60  - 21/3600, 20 + 28/60 + 4/3600];
middelpos = [-31.906, 20.234];

m_plot(kap103(2), kap103(1), 'k.','MarkerSize',16)
m_text(kap103(2), kap103(1)-0.4, 'KAP 103','FontSize',16)

m_plot(middelpos(2), middelpos(1), 'k.','MarkerSize',16)
m_text(middelpos(2), middelpos(1)+0.4,'Middelpos','FontSize',16)

m_text(22, -30, 'South Africa', 'FontWeight', 'Bold','FontSize',18)
%m_text(16, -27, 'Namibia', 'FontWeight', 'Bold','FontSize',18)

% White background color
set(gcf,'color','w');
set(gcf,'defaultFigureColor',[1,1,1]); 
exportgraphics(gcf, 'map.pdf')
%print -dpdf map.pdf
%addpath([fileparts(mfilename('fullpath')),'/../../plot/export_fig/'])
%export_fig('map.pdf');