fe = fftfreqp(5);
fg = fftfreqp(5);

Z = repmat(1+1j,length(fe),1);
Z(1) = 1;

% If fe == fg, no interpolation done.
Zg = Zinterp(fe,Z,fg);
assert(all(Zg(:) == Z(:)));

% No zero frequency given for fe. Zg(1) should be zero.
Zg = Zinterp(fe(2:end),Z(2:end,:),fg);
assert(Zg(1) == 0 + 0j);
assert(all(Zg(2:end) == Z(2:end)));

% No zero frequency given for fe and fg.
Zg = Zinterp(fe(2:end),Z(2:end,:),fg(2:end));
assert(all(Zg(2:end) == Z(2:end)));

% No zero frequency given for fg.
Zg = Zinterp(fe,Z,fg(2:end));
assert(all(Zg == Z(2:end)));
