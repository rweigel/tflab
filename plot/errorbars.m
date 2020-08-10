function errorbars(x,y,yl,yu,dir,linestyle)
% Similar to ERRORBAR except handles loglog and semilogx

    if nargin < 4
        yu = yl;
    end
    if nargin < 5
        dir = 'y';
    end    
    if nargin < 6
        linestyle = 'k';
    end

    scalex = get(gca,'XScale');
    scaley = get(gca,'YScale');

    scale = scaley;
    if strcmp(dir,'x')
        scale = scalex;
    end
    
    if strmatch(scalex,'linear','exact') & strmatch(scaley,'linear','exact')
        for i = 1:length(x)
            plot([x(i),x(i)],[y(i)+yu(i),y(i)-yl(i)],linestyle);
        end
    end
    if strmatch(scalex,'log','exact') & strmatch(scaley,'log','exact')
        for i = 1:length(x)
            add_head = 0;
            if y(i) - yl(i) <= 0
                add_head = 1;
                ylims = get(gca,'YLim');
                ylow = ylims(1);
            else
                ylow = y(i)-yl(i);
            end
            loglog([x(i),x(i)],[y(i)+yu(i),ylow],linestyle);
            if add_head
                loglog(x(i),ylow,'kx');
            end
        end
    end
    if strmatch(scalex,'log','exact') & strmatch(scaley,'linear','exact')
        for i = 1:length(x)
            semilogx([x(i),x(i)],[y(i)+yu(i),y(i)-yl(i)],linestyle);
        end
    end
    if strmatch(scalex,'linear','exact') & strmatch(scaley,'log','exact')
        for i = 1:length(x)
            semilogy([x(i),x(i)],[y(i)+yu(i),y(i)-yl(i)],linestyle);
        end
    end    

end