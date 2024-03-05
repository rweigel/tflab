function line = markeropts(N, s)
    if N == 1 % Single point
        ms = 30;
        line = {'.','markersize', ms};
    else
        ms = max(22-2*s,1);
        line = {'.','markersize',ms};
    end
end
