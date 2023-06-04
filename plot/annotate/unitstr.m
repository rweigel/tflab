function s = unitstr_(unit, padstart)
    prepad = ' ';
    if nargin < 2 || padstart == 0
        prepad = '';
    end
    s = '';
    if ~isempty(unit)
        s = sprintf('%s[%s]', prepad, unit);
    end
end
