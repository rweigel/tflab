function [X,a,b] = prewhiten(X,varargin)

if strcmp(varargin{1},'diff')
    a = [1];
    b = [1,-1];
    X = filter(b,a,X);
end

if strcmp(varargin{1},'yulewalker')
    % Not implemented
    [S,f] = spectrogram(X);
    [b,a] = yulewalk(N,f.'/f(end),1./mean(abs(S.')));
    X = filter(b,a,X);
end
