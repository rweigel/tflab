function errorbars_(h,x,y,dyl,dyu,lw)
%ERRORBARS_  Errorbars using same colors as markers

if isempty(dyl)
    return
end

if isscalar(lw)
    lw = repmat(lw, size(y,2));
end

ymin = NaN;
scaley = get(gca,'YScale');
if strcmp(scaley,'log')
    Ip = find(y-dyl > 0);        % Positive lower error bar value
    ymin = min(y(Ip) - dyl(Ip));
end

for j = 1:size(y,2)
    color = get(h(j),'Color');
    args = struct('Color',color,'LineWidth',lw(j));
    errorbar_(x,y(:,j),dyl(:,j),dyu(:,j),'y',ymin,args);
end
