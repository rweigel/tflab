%% Configure

close all
set(0,'defaultFigureWindowStyle','docked');
%set(0,'defaultFigureWindowStyle','normal');

lemi = 0; % Set to 1 to use lemimt code to compute Z.
if lemi == 1
    lemimt_dir = '../../lemimt/';
end

%% Run
S = run('KAP103',lemi);
%S = {S{1},S{3}};

S{1}.Metrics.Segment = struct();

S{1}.ZCL = struct();

S{1}.ZCL.Bootstrap = struct();

for i = 1:size(S{1}.Segment.Z,2)
    tmp = squeeze(S{1}.Segment.Z(:,i,:));
    S{1}.ZCL.Bootstrap.x_95(:,:,i) = bootcl(tmp, @mean, 1000, 0.025);
    S{1}.ZCL.Bootstrap.x_1sigma(:,:,i) = bootcl(tmp, @mean, 1000, 0.341);
    [abs(S{1}.ZCL.Bootstrap.x_95(:,1,i)),...
     abs(S{1}.ZCL.Bootstrap.x_1sigma(:,1,i)),...
     abs(S{1}.Z(:,i)),...
     abs(mean(tmp,2)),...
     abs(S{1}.ZCL.Bootstrap.x_1sigma(:,2,i)),...
     abs(S{1}.ZCL.Bootstrap.x_95(:,2,i))]
end
%S{1}.ZVAR = std(abs(S{1}.Segment.Z),0,3)...
%                /sqrt(size(S{1}.Segment.Z,3));
S{1}.Metrics.Segment.fe = S{1}.Segment.fe;
S{1}.Metrics.Segment.SN = mean(S{1}.Segment.Metrics.SN,3);

S{1}.Metrics.Segment.SN_1sigma_normal = ...
    std(S{1}.Segment.Metrics.SN,0,3)/sqrt(size(S{1}.Segment.Metrics.SN,3));

Sr = {S{1},S{3}};
plots(Sr,0) % Or plots('KAP03',0) 

return;

S = run('Middelpos',lemi);
plots(S) % Or plots('Middelpos')

compare % Compares Middelpos and KAP103
