function [fe,Ic,Ne] = evalfreq(N,Npd,scale)
%EVALFREQ Linear or logarithmicaly spaced center frequencies and bands
%
%  [fe,Ic,Ne] = EVALFREQ(N) returns approximately 7 evaluation frequencies
%  per decade in frequency range [2/N,1/4]. Ic are the center indices and
%  Ne is the number of points to left and right of center frequency to use
%  as a inclusion band. For fe = 0, Ne = 0. 
%
%  Use evalfreq_plot(N) to visualize the evaluation frequencies and their
%  associated frequency bands.
%
%  The left boundary of a frequency window i (for i > 2) starts at fe(i-1)
%  and the number of DFT points between fe(i) and fe(i-1) is used to
%  compute the window width. See EVALFREQ_DEMO for a visualization of the
%  frequencies and windows.
% 
%  [fe,Ic,Ne] = EVALFREQ(N, Npd) returns approximately Npd frequencies per
%  decade.
%
%  [fe,Ic,Ne] = EVALFREQ(N,w,'linear') returns evaluation frequencies
%  separated by w points and windows of size 2*w+1. For fe = 0, Ne = 0.
%
%  [fe,Ic,Ne] = EVALFREQ(N,[dN,w],'linear') returns evaluation frequencies
%  separated by dN points and windows of size 2*w+1.
%
%  See also FFTFREQ, FFTFREQP, EVALFREQ_DEMO.

if nargin < 2
    Npd = 7; % Number per decade
end
if nargin < 3
    scale = 'logarithmic';
end

assert(any(strcmp(scale,{'linear','logarithmic'})),...
        'scale must be ''linear'' or ''logarithmic''');

% TODO: Allow fl and fu to be given as opposed to 2/N and 1/4.
% TODO: Allow Npd to be array [Npd,No], where No is overlap

[~,f] = fftfreq(N); % Unique DFT frequencies (includes 0.5 if N even).

if strcmp(scale,'linear')
    if length(Npd) > 1
        dN = Npd(1);
        w = Npd(2);
        assert(w >= 0,'w >= 0 is required.');
    else
        w = Npd;
        dN = Npd;
    end
    assert(dN > 0,'dN must be greater than zero.');
    
    % Linearly spaced center frequencies
    Ic = [1,2+w:dN:length(f)-w]; % Indices of window centers
    % e.g., if w = 1, first window center at 3 and
    % window will extend from 2 through 4; second window center at 3+dN.
    if mod(N,2) == 0
        Ic = [Ic(1:end-1),length(f)];
        % f = 0.5 is always real and so imaginary part of Z cannot be
        % regressed with other frequencies. We could allow f = 0.5
        % to be regressed with real parts for f < 0.5, but this would
        % require a significant amount of extra code in the regression
        % functions. As a result, we always treat f = 0.5 in the same
        % ways as f = 0 - Z for these frequencies is always based on one
        % DFT point, so the result of regression (for 1-D case) is
        %   Z(f=0) = E(f=0)/B(f=0)
        %   Z(f=0.5) = E(f=0.5)/B(f=0.5)
        % In case where B has more than one component,
        %   Z(f=0,:) = NaN
        %   Z(f=0.5,:) = NaN
    end

    Ne = w*ones(size(Ic));
    Ne(1) = 0; % f = 0
    if mod(N,2) == 0
        Ne(end) = 0; 
    end
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
    if mod(N,2) == 0
        % Remove f = 0.5, which is always real. See note above.
        f = f(1:end-1); 
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