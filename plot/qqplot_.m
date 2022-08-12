function qqplot_(S,idx)

if nargin == 1
    idx = 3;
end

Rediduals = S.Regression.Residuals{idx};

PositionTop = [0.1300 0.5400 0.7750 0.38];
PositionBottom = [0.1300 0.1100 0.7750 0.38];
figprep();

f = S.fe(idx);
subplot('Position',PositionTop);
    qqplot(real(Rediduals));
    grid on;box on;
    title(sprintf('$f=%g; T=%.1f$',f,1/f));
    set(gca,'XLabel',[]);
    ylabel('Re[Error$(f)]$')
    adjust_exponent();
subplot('Position',PositionBottom);
    qqplot(imag(Rediduals));
    grid on;box on;
    title('')
    ylabel('Im[Error$(f)]$')
    adjust_exponent()

end