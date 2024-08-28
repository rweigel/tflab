function [t_simple, t_latex] = evalfreq_log(N,varargin)
%EVALFREQ_LOG  Display table of evaluation frequencies and periods
%
%   [fe,Ic,Ne] = EVALFREQ_LOG(N) takes the same inputs as evalfreq(),
%   returns the same outputs, and prints a table.
%
%   See also EVALFREQ, EVALFREQ_DEMO, EVALFREQ_PLOT.

if nargin > 1
    [fe,Ic,Ne] = evalfreq(N,varargin{:});
    if length(varargin) == 1
        if length(varargin{1}) == 1
            fprintf('evalfreq_log.m: evalfreq(%d,%d) returns:\n',N,varargin{1});
        else
            fprintf('evalfreq_log.m: evalfreq(%d,[%d,%d]) returns:\n',N,varargin{:});            
        end
    else
        if length(varargin{1}) == 1
            fprintf('evalfreq_log.m: evalfreq(%d,%d,%s) returns:\n',N,varargin{:});
        else
            fprintf('evalfreq_log.m: evalfreq(%d,[%d,%d],''%s'') returns:\n',N,varargin{:});            
        end        
    end
else
    [fe,Ic,Ne] = evalfreq(N);
    fprintf('evalfreq_log.m: evalfreq(%d) returns:\n',N);
end

[~,f] = fftfreq(N);

fprintf('evalfreq_log.m: %d DFT frequencies\n',N);
t_simple = sprintf('   fe        fl        fu        Te        Tl        Th           il         ih      Ne\n');
    
t_latex = '';
t_latex = append(t_latex,'\\begin{tabular}{c c c c c c c c c}',newline());
t_latex = append(t_latex,'\\hline',newline());
t_latex = append(t_latex,'$f_e$ [Hz] & $f_l$ [Hz] & $f_u$ [Hz] & $T_e$ [s] & $T_l$ [s] & $T_u$ [s] & $i_l$ & $i_u$ & $N_e$ \\\\',newline());
t_latex = append(t_latex,'\\hline',newline());

io = 1;
if fe(1) == 0
    io = 2;
end
for i = io:length(fe)
    tmp = sprintf('$%7.2e$ & $%7.2e$ & $%7.2e$ & $%7.2e$ & $%7.2e$ & $%7.2e$ & $%10d$ & $%10d$ & $%7d$',...
                    fe(i),f(Ic(i)-Ne(i)),f(Ic(i)+Ne(i)),...
                    1/fe(i),1/f(Ic(i)+Ne(i)),1/f(Ic(i)-Ne(i)),...
                    Ic(i)-Ne(i),Ic(i)+Ne(i),... 
                    2*Ne(i)+1);
    tmp = append(tmp,'\\\\',newline());

    % Regex is brittle.
    tmpl = regexprep(tmp,'e(\+|-)([0-9].?)','\\\\cdot 10^{$1$2}');
    tmpl = regexprep(tmpl,'\{-0','{-');
    tmpl = regexprep(tmpl,'\{\+0','{');
    t_latex = append(t_latex,tmpl);

    tmp = replace(tmp,'\\','');
    tmp = replace(tmp,'$','');
    tmp = replace(tmp,' &','');
    t_simple = append(t_simple,tmp);
   
end
t_latex = append(t_latex,'\\hline \\\\',newline());
t_latex = append(t_latex,'\\end{tabular}',newline());
fprintf(t_simple);
%fprintf(t_latex);
