
clear
if 0
    stationids = {'KAP103', 'KAP106', 'KAP109', 'KAP112'};
    
    figure();
    for id = 1:length(stationids)
        stationid = stationids{id};
        EDI = read_edi([scriptdir(),'/data/KAP03/edi/',lower(stationid),'.edi']);
        
        Z = EDI.Z;
        fe = EDI.fe;
        Z(Z > 1e31) = NaN; % TODO: Get from EMPTY= line in file; kap003 file has EMPTY= 0.1000000E+33
        loglog(1./fe, abs(Z(:,2)), 'o', 'LineWidth',2)
        Midd 
        ZVAR = EDI.ZVAR;
        ZVAR(ZVAR > 1e31) = NaN;
    end
    legend(stationids)
end

%range = '20031108-20031204';
range = '20031108-20031118';
%range = '20031121-20031204';
outdir = fullfile(scriptdir(), 'data', 'KAP03_Middelpos');
fname_tmp = fullfile(outdir, sprintf('KAP03_Middelpos_compare_%s.mat',range));

start   = '20120712';
stop    = '20121107';

if ~exist(fname_tmp, 'file')

    filestr = 'Middelpos';
    dirstr  = sprintf('tfs-%s-%s',start,stop);
    rundir  = fullfile(scriptdir(),'data',filestr,dirstr);
    for tfn = 1:1
        fname = fullfile(rundir, sprintf('%s-tf%d.mat',filestr,tfn));
        TFs{tfn} = loadtf(fname);
    end
    
    
    rundir = fullfile(scriptdir(),'data','KAP03','KAP103',sprintf('tfs-%s',range));
    f = fullfile(rundir, 'KAP103-tf3.mat');
    TFs{6} = loadtf(f);
    
    f = fullfile(rundir, 'KAP103-tf1.mat');
    TFs{7} = loadtf(f);
    
    TFs{6} = tflab_preprocess(TFs{6});
    TFs{6} = tflab_metrics(TFs{6});
    TFs{7} = tflab_preprocess(TFs{7});
    TFs{7} = tflab_metrics(TFs{7});
    save(fname_tmp,'TFs')
else
    fprintf('Reading %s\n',fname_tmp)
    load(fname_tmp)
end

period_range = [10, 30000];

copts = struct();
    copts.print = 0; % Set to 1 to print pdf of each figure created.
    copts.printOptions.printDir = fullfile(outdir,'figures');
    copts.printOptions.printFormats = {'pdf','png'};
    copts.title = '';

%TFs{6}.Z = TFs{6}.Z/10;
%TFs{7}.Z = TFs{7}.Z/10;

if 1
figure();
    zopts = copts;
    zopts.type = 1;
    zopts.period_range = period_range;
    zplot({TFs{1}, TFs{6}, TFs{7}},zopts);
end

figure();
    TFs{1}.Options.description = sprintf('Middelpos %s-%s',start,stop);
    TFs{6}.Options.description = sprintf('KAP103 %s',range);
    dftopts = copts;
    dftopts.type = 'original-averaged';
    dftopts.period_range = period_range;
    dftplot({TFs{1}, TFs{6}},dftopts);
    title(range)
    %dftplot(TFs{6},dftopts);
    %keyboard

figure();
    TFs{1}.Options.description = sprintf('Middelpos %s-%s',start,stop);
    TFs{6}.Options.description = sprintf('KAP103 %s',range);
    dftopts = copts;
    dftopts.type = 'original-raw';
    dftopts.period_range = period_range;
    dftplot({TFs{1}, TFs{6}},dftopts);
    title(range)
    %dftplot(TFs{6},dftopts);
    %keyboard

if 0

     figure();
        dftopts = copts;
        dftopts.type = 'original-raw';
        % Check amplitudes by replacing data with known dft amplitudes
        for tf = [1, 6]
            t = 1:size(TFs{tf}.In, 1);
            T = 60/TFs{tf}.Metadata.timedelta;
            % Multiply by tf to prevent overlap
            TFs{tf}.In(:,1) = tf*100*sin(2*pi*t/T);
            TFs{tf}.In(:,2) = tf*1000*sin(2*pi*t/T);
            if isfield(TFs{tf}, 'DFT')
                TFs{tf} = rmfield(TFs{tf},'DFT');
            end
        end
        %dftopts.period_range = period_range;
        % Peak in Bx for TFs{1} should be 100
        % Peak in By for TFs{1} should be 1000
        % Peak in Bx for TFs{6} should be 600
        % Peak in By for TFs{6} should be 6000
        dftplot({TFs{1}, TFs{6}},dftopts);
        %dftplot(TFs{6},dftopts);
        %dftplot(TFs{1},dftopts);
end

if 0
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
end

if 0
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
if 0
clc
clear all
fc=0.1;% cut off frequency
w=2*pi*fc;% convert to radians per second
fn=25; %nyquivst frequency = sample frequency/2;
order = 6; %6th order filter, high pass
[b14, a14]=butter(order,(w/fn),'high');
fvtool(b14,a14);
end