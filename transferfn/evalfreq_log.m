function evalfreq_log(N,varargin)
%EVALFREQ_LOG
%
%   [fe,Ic,Ne] = evalfreq_log(...) takes the same inputs as EVALFREQ(),
%   returns the same outputs, and prints the frequencies.
%
%   See also EVALFREQ, EVALFREQ_PLOT.

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

fprintf('evalfreq_log.m: __________________________________________\n');
for i = 1:length(fe)
    fprintf('evalfreq_log.m: fe = %.2g, flow = %.2g, fhigh = %.2g, Ne = %3d\n',...
            fe(i),f(Ic(i)-Ne(i)),f(Ic(i)+Ne(i)),Ne(i));
end
fprintf('evalfreq_log.m: __________________________________________\n');
for i = 1:length(fe) 
    fprintf('evalfreq_log.m: Te = %.2f, Tlow = %.2f, Thigh = %.2f, Ne = %3d\n',...
            1/fe(i),1/f(Ic(i)+Ne(i)),1/f(Ic(i)-Ne(i)),Ne(i));
end
fprintf('evalfreq_log.m: __________________________________________\n');