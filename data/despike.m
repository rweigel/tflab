function E = despikeE(E,threshold,a,b)

if nargin < 2,threshold = 1;end
if nargin < 3,a = 2;end
if nargin < 4,b = 4;end

t = [1:size(E,1)]';

for c = 1:size(E,2) % Columns
    I = find(abs(diff(E(:,c))) >= threshold); 
    for i = 1:length(I)
        ax = a;
        bx = b;
        if (I(i) - a < 1),ax = 0;end
        if (I(i) + b >= size(E,1)),bx = 0;end
        E(I(i)-ax:I(i)+bx,c) = NaN;
    end
    Ig = ~isnan(E(:,c));
    x = t(Ig);
    y = E(Ig,c);
    if length(x) < 2
        warning('despikeE: Despiking only left 1 element\n');
        return
    end
    E(:,c) = interp1(x,y,t);
    f = (size(E,1)-length(find(Ig == 1)))/size(E,1);
    fprintf('despikeE.m: Removed %d possible spikes in E(:,%d). %.2f%% of data modified\n',length(I),c,100*f);
end