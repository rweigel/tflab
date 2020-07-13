addpath([fileparts(mfilename('fullpath')),'/m_map']);
addpath([fileparts(mfilename('fullpath')),'/../../plot/export_fig/'])

M = m_shaperead('countries/ne_110m_admin_0_countries')

clf
m_proj('miller','long',[12 32],'lat',[-35 -25]);
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
kap103 = [-32 - 8/60  - 21/3600, 20 + 28/60 + 4/3600];
middelpos = [-31.5436,20.1405];

m_plot(kap103(2), kap103(1), 'k.','MarkerSize',16)
m_text(kap103(2)+0.2, kap103(1)+0.2, 'KAP 103','FontSize',16)

m_plot(middelpos(2), middelpos(1), 'k.','MarkerSize',16)
m_text(middelpos(2)+0.2, middelpos(1)+0.2,'Middelpos','FontSize',16)

m_text(22, -28, 'South Africa', 'FontWeight', 'Bold','FontSize',18)
m_text(16, -27, 'Namibia', 'FontWeight', 'Bold','FontSize',18)

% White background color
set(gcf,'color','w');
set(gcf,'defaultFigureColor',[1,1,1]); 
export_fig('map.pdf');