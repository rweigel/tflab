function [fe,Ic,Ne] = evalfreq(N,Npd,scale)
%EVALFREQ Logarithmically spaced center frequencies and windows
%
%  [fe,Ic,Ne] = evalfreq(N) returns approximately 7 evaluation frequencies
%  per decade in frequency range [2/N,1/4]. Ic are the center indices and
%  Ne is the number of points to left and right of center frequency to use
%  as a inclusion window. For fe = 0, Ne = 0. 
%
%  The left boundary of a frequency window i (for i > 2) starts at fe(i-1)
%  and the number of DFT points between fe(i) and fe(i-1) is used to
%  compute the window width. See EVALFREQ_DEMO for a visualization of the
%  frequencies and windows.
% 
%  [fe,Ic,Ne] = evalfreq(N,Npd) returns approximately Npd frequencies per
%  decade.
%
%  [fe,Ic,Ne] = evalfreq(N,dN,'linear') returns evaluation frequencies
%  separated by dN points and windows of size 2*dN+1. For fe = 0, Ne = 0.
%
%  [fe,Ic,Ne] = evalfreq(N,[dN,w],'linear') returns evaluation frequencies
%  separated by dN points and windows of size 2*w+1.
%
%  See also EVALFREQ_DEMO.

if nargin < 2
    Npd = 7;
end
if nargin < 3
    scale = 'logarithmic';
end

assert(any(strcmp(scale,{'linear','logarithmic'})),...
        'scale must be ''linear'' or ''logarithmic''');

% TODO: Allow fl and fu to be given as opposed to 2/N and 1/4.
% TODO: Allow Npd to be array [Npd,No], where No is overlap

[~,f] = fftfreq(N); % Unique DFT frequencies

if strcmp(scale,'linear')
    if length(Npd) > 1
        w = Npd(2);
        Npd = Npd(1);
        assert(w >= 0,'w >= 0 is required.');
    else
        w = Npd;
    end
    assert(Npd > 0,'dN must be greater than zero.');
    
    % Linearly spaced center frequencies
    Ic = [1,2+w:Npd:length(f)-w]; % Indices of center points
    % e.g., if Npd = 3, first window center at 5 and
    % window will extend from 2 through 8.
    % second window center at 12 and window from 9 through 15.
    Ne = w*ones(size(f));
    Ne(1) = 0; % f = 0 
    fe = f(Ic);
end

if strcmp(scale,'logarithmic')
    assert(Npd > 0,'Npd must be greater than zero.');
    
    Nd = log10(1/4)-log10(2/N); % Number of decades

    % Nominal evaluation frequencies. 
    fen = logspace(log10(2/N),log10(1/4),round(Nd*Npd));

    if isempty(fen) % Will happen when N < 10.
        fen = 0.25;
    end

    % fea will be an array of evaluation frequencies that are also actual
    % frequencies. (Needed so that when windowing is used, the same
    % number to left and right of eval freq are used).
    fea = [f(2),fen,f(end)];
    for i = 1:length(fea)
        [~,Ic(i)] = min(abs(f-fea(i)));
    end
    fea = f(Ic);
    % Actual evaluation frequencies, fea, are now equal to an actual
    % frequency.

    % Compute number of points to left and right to use for window.
    % Window will be symmetric in linear space (if asymmetric window is
    % used, one would need to change location of evaluation frequencies
    % to be at center of window).
    for i = 2:length(fea)-1
        Il = find(f >= fea(i-1),1); % Find frequency nearest above or equal to previous eval freq.
        Iu = Ic(i)+(Ic(i)-Il); % Upper frequency is determined by how many frequencies below were used.
        fe(i-1) = fea(i);
        N(i-1) = (Iu(1)-Il(end))/2;
    end
    Ic = Ic(2:end-1);

    fe = [0,fe];
    Ic = [1,Ic];
    Ne  = [0,N];
    [fe,Iu] = unique(fe);
    Ic = Ic(Iu);
    Ne = Ne(Iu);
end