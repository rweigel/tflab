clc;

sigma = 10;
offset = 1;

[Zc,Zx,Rc,Rx,Rcm,Rxm] = compute(offset,sigma);
sqrtN = sqrt(2*length(Rc));
fprintf('max(Zx-Zc)/(sqrt(N)*eps)  = %.1f\n', (max(Zx - Zc))/(sqrtN*eps));
fprintf('max(Rx-Rc)/(sqrt(N)*eps)  = %.1f\n', (max(Rx - Rc))/(sqrtN*eps));
fprintf('max(Rxm-Rcm)/(sqrt(N)*eps) = %.1f\n', (max(Rxm - Rcm))/(sqrtN*eps));

function [Zc,Zx,Rc,Rx,Rcm,Rxm] = compute(offset,sigma)
    
    % Estimate
    %   real(E) = real(B*Z) + real(Noise)
    %   imag(E) = imag(B*Z) + imag(Noise)
    %
    %   Er = Zr*Br - Zi*Bi + dZr + Noiser
    %   Ei = Zr*Bi + Zi*Br + dZi + Noisec
    % using regress() where r and i subscripts are the real and imaginary
    % parts. If offset = 0, dZr = dZi = 0.
    %   
    % The result regress() using this form is written in complex form using
    %   Zx = Zr + 1j*Zi and the residuals are 
    %   Rx = [Rxr + 1j*Rxi] 
    % Rxx are the resisudals computed manually instead of using the output
    % of regress();
    %
    % Estimate 
    %   E = B*Z + Noise
    % using [Zc,~,Rc] = regress(E,B) where E, B, Z, and Noise are complex.
    %
    % Zx and Zc should match to order of (machine precision)/sqrt(N).
    %

    Zc = [2*(1+1j); 3*(1+1j)];
    Nd = length(Zc);
    N = 2*length(Zc) + offset; % Number of regression values.
    N = 2*N;                   % Make over-determined

    Noisec = sigma*randn(N,1) + 1j*sigma*randn(N,1);

    ftBo = randn(N,Nd) + 1j*randn(N,Nd);

    dZc = 4*(1+1j);
    
    if offset
        ftB = [ftBo, ones(size(ftBo,1),1)];
        Zc = [Zc; dZc];
        ftE = ftB*Zc;
        ftE = ftE + Noisec;
        [Zc,~,Rc] = regress(ftE,ftB);
        Rcm = ftE - ftB*Zc; % Manual calculation of residuals
    else
        ftE = ftBo*Zc;
        ftE = ftE + Noisec;
        [Zc,~,Rc] = regress(ftE,ftBo);
        Rcm = ftE - ftBo*Zc ; % Manual calculation of residuals
    end

    % Do equivalent of above using only real values.

    E = [real(ftE); ...
         imag(ftE)];
    
    if offset
        B = [real(ftBo),  ones(size(ftBo,1),1), -imag(ftBo), zeros(size(ftBo,1),1);...
             imag(ftBo), zeros(size(ftBo,1),1),  real(ftBo),  ones(size(ftBo,1),1)];
    else
        B = [real(ftBo), -imag(ftBo);...
             imag(ftBo),  real(ftBo)];
    end
    
    [Zx,~,Rx] = regress(E,B);
    Rxm = E - B*Zx; % Manual calculation of residuals
    Rxm = Rxm(1:end/2) + sqrt(-1)*Rxm(end/2+1:end);

    Zx = Zx(1:end/2) + sqrt(-1)*Zx(end/2+1:end);
    Rx = Rx(1:end/2) + sqrt(-1)*Rx(end/2+1:end);
end
