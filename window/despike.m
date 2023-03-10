function X = despike(X,thresh,window)
%DESPIKE - Basic depspiking algorithm
%
%  Xd = DESPIKE(X,thresh,window) despikes each column of X by finding
%  indices i where the change is >= thresh and replacing data from
%  i - window(1) to i + window(2) with NaN.
%
%  See also naninterp1.

% TODO?: Allow interpolation arguments to be passed as is done in
%        naninterp1.

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
    logmsg(['despike.m: Removed %d possible spikes in column %d.',...
             ' (%.2f%%)\n'],Nb,c,100*Nb/Nr);
end