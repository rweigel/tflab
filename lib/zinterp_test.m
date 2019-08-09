% See zinterp_demo for extensive visual verification tests.
% The tests here are primarily API/documentation tests.

addpath([fileparts(mfilename('fullpath')),'/../misc/']);

fe = fftfreqp(5)';
fg = fftfreqp(5)';

Z = repmat(1+1j,length(fe),1);
Z(1) = 1;

% zinterp(f,Z,N) usage
logmsg(dbstack, 'Test 1\n');
% Check that positive Zi values match given Z values when no interp
% performed.
[Zi,fi] = zinterp(fe,Z,5);
[~,I] = intersect(fi,fg);
assert(all(Zi(I) == Z));
% Check that f > 0 values of Z are conj of f < 0 values of Z
assert(all(Zi(2:3) == conj(Zi(4:5))));

logmsg(dbstack, 'Test 2\n');
% If fe == fg, no interpolation done.
Zg = zinterp(fe,Z,fg);
assert(all(Zg(:) == Z(:)));

logmsg(dbstack, 'Test 3\n');
% No zero frequency given for fe. Zg(1) should be zero.
Zg = zinterp(fe(2:end),Z(2:end,:),fg);
assert(Zg(1) == 0 + 0j);
assert(all(Zg(2:end) == Z(2:end)));

logmsg(dbstack, 'Test 4\n');
% No zero frequency given for fe and fg.
Zg = zinterp(fe(2:end),Z(2:end,:),fg(2:end));
assert(all(Zg(2:end) == Z(2:end)));

logmsg(dbstack, 'Test 5\n');
% No zero frequency given for fg.
Zg = zinterp(fe,Z,fg(2:end));
assert(all(Zg == Z(2:end)));
