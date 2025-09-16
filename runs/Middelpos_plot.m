function Middelpos_plot(rundir, filestr, print_figs)

logmsg('\n');
for tfn = 1:5
    fname{tfn} = fullfile(rundir, sprintf('%s-tf%d.mat',filestr,tfn));
    TFs{tfn} = loadtf(fname{tfn});
end

period_range = [1, 86400 + 3600];

copts = struct();
    copts.print = print_figs; % Set to 1 to print pdf of each figure created.
    copts.printOptions.printDir = fullfile(rundir,'figures');
    copts.printOptions.printFormats = {'pdf','png'};
    copts.title = '';

if print_figs
    dock off;close all;
    if exist(copts.printOptions.printDir, 'dir')
        rmdir(copts.printOptions.printDir,'s');
    end
    mkdir(copts.printOptions.printDir);
else
    dock on;figure(1);close all;
end

if 1
%% Time series plots
tsopts = copts;
if (1)
    figure();
        tsopts.type = 'original';
        tsplot(TFs{1},tsopts);
end

if (1)
    tsopts = copts;
    tsopts.type = 'error';
    for i = 1:length(TFs)
        figure();
            tsplot(TFs{i},tsopts);
    end
end

%% DFT plots
figure();
    dftopts = copts;
    dftopts.type = 'original-averaged';
    dftplot(TFs{1},dftopts);

figure();
    dftopts = copts;
    dftopts.type = 'original-averaged-phases';
    dftplot(TFs{1},dftopts);

%figure();
    %dftopts.type = 'error-averaged-magphase';
    %dftplot(TFs,dftopts);

%% SN plots
figure();
    snopts = copts;
    snopts.type = 1;
    snopts.period_range = period_range;
    snplot(TFs([1,3]),snopts);

figure();
    snopts = copts;
    snopts.type = 1;
    snopts.period_range = period_range;
    snplot(TFs,snopts);
    
%% Z plots
figure();
    zopts = copts;
    zopts.type = 1;
    zopts.period_range = period_range;
    zplot(TFs([1,3]),zopts);

figure();
    zopts = copts;
    zopts.type = 1;
    zopts.period_range = period_range;
    zplot(TFs,zopts);


%% qq plots    
figure();
    qqopts = copts;
    %qqopts.printOptions.printDir = fullfile(rundir,'figures','qqplot');
    fidx = 20; % frequency number
    comp = 2;  % component (x = 1, y = 2)
    qqplot_(TFs{1},qqopts,comp,fidx);

figure();
    qqopts = copts;
    qqopts.type = 'combined';
    %qqplot_(TFs,qqopts);
end

if 1
    %outdir = fullfile(scriptdir(),'data','KAP03','KAP103','tfs-20031108-20031118');
    %outdir = fullfile(scriptdir(),'data','KAP03','KAP103','tfs-20031121-20031204');
    outdir = fullfile(scriptdir(),'data','KAP03','KAP103','tfs-20031108-20031204');
    f = fullfile(outdir, 'KAP103-tf3.mat');
    TFs{6} = loadtf(f);

    f = fullfile(outdir, 'KAP103-tf1.mat');
    TFs{7} = loadtf(f);

    TFs{6} = tflab_preprocess(TFs{6});
    TFs{6} = tflab_metrics(TFs{6});
    TFs{7} = tflab_preprocess(TFs{7});
    TFs{7} = tflab_metrics(TFs{7});
    %TFs{6}.Z = TFs{6}.Z/10;
    %TFs{7}.Z = TFs{7}.Z/10;
    figure();
        zopts = copts;
        zopts.type = 1;
        zopts.period_range = [10, 30000];
        zplot({TFs{1}, TFs{6}, TFs{7}},zopts);

    %keyboard
    if 0
        figure();
            dftopts = copts;
            dftopts.type = 'original-averaged';
            dftopts.period_range = period_range;
            dftplot({TFs{1}, TFs{6}},dftopts);
            %dftplot(TFs{6},dftopts);
            %keyboard
    
        figure();
            dftopts = copts;
            dftopts.type = 'original-raw';
            % Check amplitudes
            if 0
                for tf = [1, 6]
                    t = 1:size(TFs{tf}.In, 1);
                    T = 60/TFs{tf}.Metadata.timedelta;
                    TFs{tf}.In(:,1) = tf*100*sin(2*pi*t/T);
                    TFs{tf}.In(:,2) = tf*1000*sin(2*pi*t/T);
                    if isfield(TFs{tf}, 'DFT')
                        TFs{tf} = rmfield(TFs{tf},'DFT');
                    end
                end
            end
            %dftopts.period_range = period_range;
            %dftplot({TFs{1}, TFs{6}},dftopts);
            dftplot(TFs{6},dftopts);
            %dftplot(TFs{1},dftopts);
    
    By_K = TFs{6}.In(:,2);
    By_M = TFs{1}.In(:,2);
    
    I = 1:5:size(By_M,1);
    By_M = By_M(I);
    By_M = By_M(1:length(By_K));
    By_M = By_M(end/2:end);
    By_K = By_K(end/2:end);
    figure();clf;hold on;
        plot(By_M);
        plot(By_K);
        legend('Middelpos By', 'KAP103 By')
    
    [dftBy_K, fK] = fftu(By_K);
    [dftBy_M, fM] = fftu(By_M);
    dftBy_K = abs(dftBy_K)/(length(By_K)/2);
    dftBy_M = abs(dftBy_M)/(length(By_M)/2);
    
    %dftBy_K = abs(fft(By_K))/length(By_K);
    %dftBy_M = abs(fft(By_M))/length(By_M);
    
    figure();clf;
        loglog(fM, dftBy_M);
        hold on;
        loglog(fK, dftBy_K);
        legend('Middelpos', 'KAP103')
    
        TFs{6}.In = TFs{1}.In;
        TFs{6}.Out = TFs{1}.Out;
        TFs{6}.Z = -TFs{6}.Z;
        
        TFs{7}.In = TFs{1}.In;
        TFs{7}.Out = TFs{1}.Out;
        TFs{7}.Z = -TFs{7}.Z;    
        
    
            
        figure();
            snopts = copts;
            snplot(TFs,snopts);
    end
end
    
if print_figs == 1
    figHTML(copts.printOptions.printDir)
end
