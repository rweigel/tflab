function errorbars_(h,x,y,dyl,dyu,lw)
%ERRORBARS_  Errorbars using same colors as markers

if isempty(dyl),return,end
colors = get(h,'Color');
if ~iscell(colors)
    colors = {colors};
end
for comp = 1:size(y,2)
    args = {'y','Color',colors{comp},'LineWidth',lw};
    errorbar_(x,y(:,comp),dyl(:,comp),dyu(:,comp),args{:});
end
