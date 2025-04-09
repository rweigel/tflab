function errorbars_(h,x,y,dyl,dyu,lw)
%ERRORBARS_  Errorbars using same colors as markers

if isempty(dyl)
    return
end

if iscell(x)
    % TODO: Check that x,y,dyl, and dyu have same dimensionality and type
end

if iscell(x)
    xc = x;
    yc = y;
    dylc = dyl;
    dyuc = dyu;
else
    for j = 1:size(x,2)
        xc{j} = x(:,j);
        yc{j} = y(:,j);
        dylc{j} = dyl(:,j);
        dyuc{j} = dyu(:,j);
    end
end

lw = repmat(lw, length(yc));

ymin = NaN;
scaley = get(gca,'YScale');
if strcmp(scaley,'log')
    y = cat(1,yc{:});
    dyl = cat(1,dylc{:});
    Ip = find(y-dyl > 0);        % Positive lower error bar value
    ymin = min(y(Ip) - dyl(Ip));
end

for j = 1:length(yc)
    color = get(h(j),'Color')
    args = struct('Color',color,'LineWidth',lw(j));
    errorbar_(xc{j},yc{j},dylc{j},dyuc{j},'y',ymin,args);
end
