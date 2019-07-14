fe = [1/16,1/8,1/4];
Z  = [1,2,4]';
Z  = Z+sqrt(-1)*Z;
fg = linspace(0,0.5,10);
Z1 = Zinterp(fe,Z,fg,1);
[fe',Z]
[fg',Z1]

Iz = find(fg >= fe(end) | fg <= fe(1));
if ~all(Z1(Iz) == 0)
    error('All Z values at frequencies outside of fe range should be zero.')
end

fe = [0,1/16,1/8,1/4];
Z  = [0,1,2,4]';
Z  = Z+sqrt(-1)*Z;
fg = linspace(0,0.5,10);
Z1 = Zinterp(fe,Z,fg,1);
[fe',Z]
[fg',Z1]

Iz = find(fg >= fe(end) | fg <= fe(1));
if ~all(Z1(Iz) == 0)
    error('All Z values at frequencies outside of fe range should be zero.')
end

fe = [0,1/8,1/4,1/2];
Z  = [10,1,2,4]';
Z  = Z+sqrt(-1)*Z;
fg = linspace(0,0.5,10);
Z1 = Zinterp(fe,Z,fg,1);
[fe',Z]
[fg',Z1]

fe = transpose([0,1/8,1/4,1/2]);
Z  = [10,1,2,4]';
Z  = Z+sqrt(-1)*Z;
fg = linspace(0,0.5,10);
Z2 = Zinterp(fe,Z,fg,1)
[fg',Z2]

if Z1 ~= Z2
    error('Z1 and Z2 should be equal.')
end

fe = transpose([0,1/8,1/4,1/2]);
Z  = [10,1,2,4]';
Z  = Z+sqrt(-1)*Z;
fg = transpose(linspace(0,0.5,10));
Z3 = Zinterp(fe,Z,fg,1)
[fg,Z3]

if Z2 ~= Z3
    error('Z2 and Z3 should be equal.')
end

fe = transpose([0,1/8,1/4,1/2]);
Z  = [10,1,2,4]';
Z  = Z+sqrt(-1)*Z;
fg = linspace(0,0.5,10);
Z4 = Zinterp(fe,Z,fg,1)
[fg',Z4]

if Z3 ~= Z4
    error('Z3 and Z4 should be equal.')
end

