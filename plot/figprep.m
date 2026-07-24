function figprep(fontsize)

if nargin == 0
    fontsize = 16;
end

set(gcf,'color','w');
set(gcf,'defaultFigureColor',[1,1,1]); 

set(gcf,'DefaultAxesXTickLabelRotationMode','manual')
set(gcf,'DefaultLegendAutoUpdate', 'off');
set(gcf,'DefaultTextInterpreter','latex');
set(gcf,'DefaultLegendInterpreter','latex');
set(gcf,'DefaultAxesTickLabelInterpreter','latex');
set(gcf,'DefaultAxesTitleFontWeight','normal');
set(gcf,'DefaultTextFontSize',18);
set(gcf,'DefaultAxesFontSize',18);
