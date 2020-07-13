%% Configure

close all
set(0,'defaultFigureWindowStyle','docked');
%set(0,'defaultFigureWindowStyle','normal');

lemi = 0; % Set to 1 to use lemimt code to compute Z.
if lemi == 1
    lemimt_dir = '../../lemimt/';
end

%% Run
S = run('Middelpos',lemi);
plots(S) % Or plots('Middelpos')

S = run('KAP103',lemi);
plots(S) % Or plots('KAP103') 

compare % Compares Middelpos and KAP103
