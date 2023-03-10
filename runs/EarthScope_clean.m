function [B,E,t,infile,outfile] = EarthScope_clean(id,plot_,print_)

if nargin < 2
    plot_ = 1;
end
if nargin < 3
    print_ = 1;
end

addpath(fullfile(fileparts(mfilename('fullpath'))),'..');
tflab_setpaths();

infile  = fullfile(scriptdir(),'data',id,[id,'_raw.mat']); 
outfile = strrep(infile,'_raw','_clean');

if exist(outfile,'file')
    logmsg('Loading: %s\n',relpath(outfile));
    load(outfile);
    logmsg('Loaded:  %s\n',relpath(outfile));
    return;
end

logmsg('Cleaning data for %s\n',id);

[B,E,t,infile] = EarthScope_read(id);

Ir = (1123950:size(B,1))'; % Points to remove

% TODO: Creating extra matrices not needed if plot_ and print_ = 0.
E1 = E;
E2(Ir,:) = NaN;
E3 = despike(E2, 1e7, [5, 5]);
E4 = naninterp1(E3);

B1 = B;
B2(Ir,:) = NaN;
B2 = removemean(B2);
B3 = despike(B2, 1e4, [5, 5]);
B4 = naninterp1(B3);

E = E4;
B = B4;

logmsg('Saving: %s\n',relpath(outfile));
save(outfile,'B','E','t','infile','outfile');
logmsg('Saved:  %s\n',relpath(outfile));

if plot_
    figure(1);figprep();clf;

    subplot(4,1,1)
        plot(E1)
        title('Raw $\mathbf{E}$');
        grid on
        axis tight;
        xlim = get(gca,'XLim');
        legend('$E_x$','$E_y$');
        zoominfo('E raw');

    subplot(4,1,2)
        plot(E2)
        title('Chunk removed');
        grid on
        axis tight;
        set(gca,'XLim',xlim);        
        legend('$E_x$','$E_y$');
        zoominfo('E after chunk removal');
        
    subplot(4,1,3)
        logmsg('Despiking E\n')
        plot(E3);
        title('Despiked');
        grid on
        axis tight;
        set(gca,'XLim',xlim);
        legend('$E_x$','$E_y$');
        zoominfo('E after despiking');
        
    subplot(4,1,4)
        logmsg('Interpolating over spike windows in E\n')
        plot(E4)
        title('Interpolated');
        grid on
        axis tight;
        set(gca,'XLim',xlim);
        legend('$E_x$','$E_y$');
        zoominfo('E after interpolation');        
        xlabel('Measurement number');    
end

if print_
    fname = sprintf('%s_E.png',outfile(1:end-4));
    logmsg(sprintf('Writing: %s\n',relpath(fname)));
	print(fname,'-dpng','-r300');
    logmsg(sprintf('Wrote:   %s\n',relpath(fname)));
end


if plot_
    figure(2);figprep();clf;
    
    subplot(4,1,1)
        plot(B1(:,1:2))
        title('Raw $\mathbf{B}$');
        grid on
        axis tight;
        legend('$B_x$','$B_y$');
        zoominfo('B');
        xlim = get(gca,'XLim');

    subplot(4,1,2)
        plot(B2(:,1:2))
        title('Chunk and then mean removed');
        grid on
        axis tight;
        set(gca,'XLim',xlim);
        legend('$B_x$','$B_y$');
        zoominfo('E');
        
    subplot(4,1,3)
        logmsg('Despiking B\n'); 
        plot(B3(:,1:2));
        title('Despiked');
        grid on
        axis tight;
        set(gca,'XLim',xlim);
        legend('$B_x$','$B_y$');

    subplot(4,1,4)
        logmsg('Interpolating over spike windows in B\n')
        plot(B4(:,1:2))
        title('Interpolated');
        grid on
        axis tight;
        set(gca,'XLim',xlim);
        legend('$B_x$','$B_y$');
        xlabel('Measurement number')    
end

if print_
    fname = sprintf('%s_B.png',outfile(1:end-4));
    logmsg(sprintf('Writing: %s\n',relpath(fname)));
	print(fname,'-dpng','-r300');
    logmsg(sprintf('Wrote:   %s\n',relpath(fname)));
end

function zoominfo(varstr)
    zoom off;
    hB = zoom(gca);
    set(hB,'ActionPreCallback',  @msgpre);
    set(hB,'ActionPostCallback', @msgpost);
    function msgpost(obj,evd)
        xlims = round(evd.Axes.XLim);
        minx = max(0,xlims(1));
        fprintf('Showing %s over XLim = [%d,%d]\n', varstr, minx, xlims(2));
    end
    function msgpre(obj,evd), end
end

end

