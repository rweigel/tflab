function errorbar_(x,y,dyl,dyu,dir,ymin,args)
%ERRORBAR_  Similar to ERRORBAR except handles loglog and semilogx

if nargin < 4
    dyu = dyl;
end
if nargin < 5
    % Only dir = 'y' is handled.
    dir = 'y';
end
if nargin < 6
    ymin = NaN;
end
if nargin < 7
    args = struct('Color', 'k');
end

fields = fieldnames(args);
linestyle = {};
for f = 1:length(fields)
    linestyle{end+1} = fields{f};
    linestyle{end+1} = args.(fields{f});
end

scalex = get(gca,'XScale');
scaley = get(gca,'YScale');

if strcmp(scaley,'log')
    Ip = [];
    if isnan(ymin)
        Ip = find(y-dyl > 0);        % Positive lower error bar value
        ymin = min(y(Ip) - dyl(Ip)); % Minimum positive lower error bar value
    end
    for i = 1:length(x)
        if isnan(y(i)) || isnan(dyl(i)) || isnan(dyu(i))
            continue;
        end
        if y(i) - dyl(i) > 0
            add_head = 0;
            ylow = y(i) - dyl(i);
        else
            add_head = 1;
            % Place head ymin/5 below lowest y.
            ylow = ymin/5;
        end
        color = args.Color;
        plot([x(i),x(i)],[y(i)+dyu(i),ylow],linestyle{:});
        if add_head
            plot(x(i),ylow,'o','MarkerSize',10, 'Color', color);
        end
    end
else
    for i = 1:length(x)
        plot([x(i),x(i)],[y(i)+dyu(i),y(i)-dyl(i)],linestyle{:});
    end
end
