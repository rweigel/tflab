function figprep(writefile,w,h)

set(gcf,'DefaultTextInterpreter','latex')
set(gcf,'DefaultLegendInterpreter','latex')
set(gcf,'DefaultAxesTickLabelInterpreter','latex');
set(gcf,'DefaultTextFontSize',14);
set(gcf,'DefaultAxesFontSize',14);

if nargin == 0
  writefile = 0;
end

if writefile
    % Open each figure in new window instead of docking.
    set(0,'defaultFigureWindowStyle','normal');
    set(0,'defaultFigurePosition', [0 0 w h]);
    set(0,'defaultFigureColor', [1,1,1]); % Background color to white.
else
    % Dock figure windows
    set(0,'defaultFigureWindowStyle','docked');
end
