function X = despike(X,thresh,window,fid)
%DESPIKE  Threshold depspike
%
%  Xd = DESPIKE(X,thresh,window) despikes each column of X by finding
%  indices i where abs(X(i+1) - X(i)) >= thresh and replacing
%  X(i - window(1):i + window(2)) with NaN.
%
%  Example X = [1, 1, 2, 2]'
%  despike(X,1,[0,0]) => [1, NaN, 2, 2]'
%  despike(X,1,[1,0]) => [NaN, NaN, 2, 2]'
%  despike(X,1,[1,1]) => [NaN, NaN, NaN, 2]'
%
%  See also naninterp1.

% TODO: Allow interpolation arguments to be passed as is done in
%       naninterp1?
%       If single row, despike it.

Nr = size(X,1); % Number of rows.

for c = 1:size(X,2)                         % Loop over columns
    I = find(abs(diff(X(:,c))) >= thresh);  % Find spikes
    for i = 1:length(I)                     % Loop over spike indices
        a = window(1);
        b = window(2);
        if (I(i) - a < 1)
            % Window extends to indices <= 0
            a = 0;
        end
        if (I(i) + b >= size(X,1))
            % Window extends to indices > length of column.
            b = 0;
        end
        % Replace values in window around spike with NaN.
        X(I(i)-a:I(i)+b,c) = NaN;
    end
    Nb = sum(isnan(X(:,c)));
    msg = sprintf('Removed %d possible spikes in column %d. (%.2f)\n',Nb,c,100*Nb/Nr);
    logmsg('%s',sprintf('Removed %d possible spikes in column %d. (%.2f)\n',Nb,c,100*Nb/Nr));
    if nargin > 3 && fid > 0
        fprintf(fid,msg);
    end
end