function exponent_relabel(dir)
%EXPONENT_RELABEL - Relabel axes number with exponents
%
%  exponent_relabel() Relabels x and y axis numbers on axis returned by
%  gca()
%
%  exponent_relabel(dir) % Relabels dir axis numbers only, where dir = x,
%  y, or z.
%  
%  Relabels 10^{-1} to 0.1, 10^{0} to 1, 10^{1} to 10, and 10^{2} to 100.
%

ax = gca();
if nargin == 0
    exponent_relabel('x');
    exponent_relabel('y');
    exponent_relabel('z');
    return;
end

assert(any(strcmp(dir,{'x','y','z'})),'dir must be x, y, or z');

drawnow;

l = get(ax,[dir,'TickLabel']);

for i = 1:length(l)
    if strcmp(l{i},'10^{-1}')
        l{i} = '0.1';
    end
    if strcmp(l{i},'$10^{-1}$')
        l{i} = '$0.1$';
    end
    if strcmp(l{i},'10^{0}')
        l{i} = '1';
    end
    if strcmp(l{i},'$10^{0}$')
        l{i} = '$1$';
    end
    if strcmp(l{i},'10^{1}')
        l{i} = '10';
    end
    if strcmp(l{i},'$10^{1}$')
        l{i} = '$10$';
    end
    if strcmp(l{i},'10^{2}')
        l{i} = '100';
    end
    if strcmp(l{i},'$10^{2}$')
        l{i} = '$100$';
    end
end

set(ax,[dir,'TickLabel'],l);