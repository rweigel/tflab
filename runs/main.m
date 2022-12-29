%% Configure

close all
set(0,'defaultFigureWindowStyle','docked');
%set(0,'defaultFigureWindowStyle','normal');

lemimt = 1; % Set to 1 to use lemimt code to compute Z, 0 otherwise.
if lemimt == 1
    lemimt = '../../lemimt/';
end

%% Run

S = run('KAP103',lemimt);
Sr = {S{1},S{3}};
plots(Sr,0) % Or plots('KAP103',0)

S = run('Middelpos',lemimt);
plots(S,0) % Or plots('Middelpos',0)

compare % Compares Middelpos and KAP103
