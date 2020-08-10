function period_lines(m)

    set(gcf, 'DefaultLegendAutoUpdate', 'off')
    yl = get(gca,'YLim');
    plot([60,60],yl,'--','Color',[0.5,0.5,0.5]);
    text(60,yl(1),'1m','VerticalAlignment','bottom');
    hr = 3600;
    if m >= 60*10
        yl = get(gca,'YLim');
        loglog([hr,hr],yl,'--','Color',[0.5,0.5,0.5]);
        text(hr,yl(1),'1h','VerticalAlignment','bottom');
    end
    if m >= hr*2
        yl = get(gca,'YLim');
        loglog([6*hr,6*hr],yl,'--','Color',[0.5,0.5,0.5]);
        text(6*hr,yl(1),'6h','VerticalAlignment','bottom');
    end            
    if m >= hr*7
        yl = get(gca,'YLim');
        loglog([12*hr,12*hr],yl,'--','Color',[0.5,0.5,0.5]);
        text(12*hr,yl(1),'12h','VerticalAlignment','bottom');
    end
    if m >= hr*13 
        yl = get(gca,'YLim');
        loglog([24*hr,24*hr],yl,'--','Color',[0.5,0.5,0.5]);
        text(24*hr,yl(1),'1d','VerticalAlignment','bottom');
    end    
end