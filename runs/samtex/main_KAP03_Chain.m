addpath('zread');
addpath('../lib');
addpath('../fft');

set(0,'defaultAxesFontSize',20)

% Get absolute path of this demo file
script_dir = fileparts(mfilename('fullpath'));

%% Read files from mtnet.info.
% Read a file in https://www.mtnet.info/data/samtex/samtex.zip
% (Backup of samtex.zip stored at http://mag.gmu.edu/git-data/zread in
% case the mtnet.info copy is no longer available.)
% Download samtex.zip if not found. Extract to ./mtnet.info/samtex
ddir = fullfile(script_dir,'data','mtnet.info','samtex');
if ~exist(ddir,'dir')
    url = 'https://www.mtnet.info/data/samtex/samtex.zip';
    fprintf('Downloading and extracting %s\n',url);
    mkdir(ddir);
    unzip(url,ddir);
end

kap03_dir = fullfile(ddir,'Final_with_DeBeers','edi_files','KAP03');
edi_file = fullfile(kap03_dir,'kap003.edi');

S = read_edi(edi_file);
S

%% Read and interpolate KAP03 EDI files
% Need to download https://github.com/rweigel/transferfn

%% Read all kap* files
edi_files = fullfile(kap03_dir,'kap*.edi');
dir_list = dir(edi_files);

for i = 1:length(dir_list)
    edi_file = fullfile(dir_list(i).folder, dir_list(i).name);
    fprintf('Reading %s\n',edi_file);
    F{i} = read_edi(edi_file);
    Ns = 86400*10/5; % Number of samples in 10 days with dt = 5 seconds.
    Npd = 7;         % Number of frequency points per decade.
    fi = evalfreq(Ns,Npd);
    %[Zi,fi] = zinterp(F{i}.fe,F{i}.Z,Ns); % Use logrithmically spaced fi.
    F{i}.Z(F{i}.Z > 1e30) = nan;
    [Zi,fi] = zinterp(F{i}.fe,F{i}.Z,fi');  % Use uniformly spaced fi.
    F{i}.Zi = Zi;
    F{i}.fi = fi;
    if i == 1
        clf;
        grid on;
        hold on;
        xlabel('f');
        title('$\log_{10}(|Z_{xy}|)$','Interpreter','Latex')
    end
    %loglog(fi,abs(Zi(:,2)),'Color',colors(i,:));
    loglog(fi,abs(Zi(:,2)));
    hold on;
    Zxx(:,i) = Zi(:,1);
    Zxy(:,i) = Zi(:,2);
    Zyx(:,i) = Zi(:,3);
    Zyy(:,i) = Zi(:,4);
end

