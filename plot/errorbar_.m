function errorbar_(x,y,dyl,dyu,dir,varargin)
%ERRORBAR_  Similar to ERRORBAR except handles loglog and semilogx

if nargin < 4
    dyu = dyl;
end
if nargin < 5
    % Only dir = 'y' is handled.
    dir = 'y';
end
if nargin < 6
    linestyle = {'k'};
else
    linestyle = varargin;
end

scalex = get(gca,'XScale');
scaley = get(gca,'YScale');

if strcmp(scalex,'linear') && strcmp(scaley,'linear')
    for i = 1:length(x)
        plot([x(i),x(i)],[y(i)+dyu(i),y(i)-dyl(i)],linestyle{:});
    end
end
if strcmp(scalex,'log') && strcmp(scaley,'log')
    Ip = find(y-dyl > 0);        % Positive lower error bar value
    ymin = min(y(Ip) - dyl(Ip)); % Minimum positive lower error bar value
    for i = 1:length(x)
        if isnan(y(i)) || isnan(dyl(i)) || isnan(dyu(i))
            continue;
        end
        if y(i) - dyl(i) > 0
            add_head = 0;
            ylow = y(i) - dyl(i);
        else
            add_head = 1;
            ylow = ymin;
        end
        loglog([x(i),x(i)],[y(i)+dyu(i),ylow],linestyle{:});
        if add_head
            loglog(x(i),ylow,'ko','MarkerSize',10);
        end
    end
end
if strcmp(scalex,'log') && strcmp(scaley,'linear')
    for i = 1:length(x)
        semilogx([x(i),x(i)],[y(i)+dyu(i),y(i)-dyl(i)],linestyle{:});
    end
end
if strcmp(scalex,'linear') && strcmp(scaley,'log')
    for i = 1:length(x)
        semilogy([x(i),x(i)],[y(i)+dyu(i),y(i)-dyl(i)],linestyle{:});
    end
end