% Report software bugs/issues/feature requests at
% https://github.com/hapi-server/client-matlab/issues
% Report data server issues to rweigel@gmu.edu

%% Set-up
% Download hapi.m if not found in path.
u = 'https://raw.githubusercontent.com/hapi-server/client-matlab/master';
if exist('hapi','file') ~= 2
    u = sprintf('%s/hapi.m',u);
    urlwrite(u,'hapi.m');
end
% Download hapi_plot.m if not found in path.
if exist('hapi','file') ~= 2
    u = sprintf('%s/hapi_plot.m',u);
    urlwrite(u,'hapi_plot.m');
end

%% Get data and metadata
server     = 'http://rweigel.dynu.net/servers/INTERMAGNET/hapi';
dataset    = 'her/definitive/minute'; % 
% Use parameters='' to request all data. Multiple parameters
% can be requested using a comma-separated list, e.g., parameters='DOY,Component 1'
parameters = ''; % 
start      = '2003-11-08T16:30:00Z';
stop       = '2003-12-05T08:54:30Z';
opts       = struct('logging',1);

[data,meta] = hapi(server,dataset,parameters,start,stop,opts);

%% Display data and metadata
data
meta
fprintf('meta.parameters = ');
meta.parameters{:}

save check_B_Hermanus data meta
%% Plot data
%hapiplot(data,meta)