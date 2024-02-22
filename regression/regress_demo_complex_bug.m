% Demonstrating an error in MATLAB's regress function. When inputs are
% complex, the confidence limits on the complex part of the computed
% slopes are the same as the computed slope. This error is indicated
% by a start in the table displayed when this script is executed.

N = 10;
sigma = 0.1;

if 0
    % The following code was posted at
    % https://www.mathworks.com/matlabcentral/answers/2080041-possible-bug-in-regress-when-y-and-x-are-complex
    % to demonstrate an error in regress() when both y and x are complex.
    % Model equation: y = x*b + error, where b = (1 + 1j).
    % The 95% confidence interval the imaginary part of b has zero width.
    rng(2);
    x = randn(N,1) + 1j*randn(N,1);
    y = x*(1+1j) + sigma*(randn(N,1) + 1j*randn(N,1));
    [b,bint] = regress(y,x)
    
    % Here we emulate how the regression is done internally. We get the
    % same result for b and the correct 95% confidence intervals on
    % both the real and imaginary parts of b. (I have verified via
    % simulation that the confidence intervals computed  for both components 
    % are consistent with a 95% confidence intervals.)
    y = [real(y);imag(y)];
    x = [real(x),-imag(x);imag(x),real(x)];
    [b,bint] = regress(y,x);
    b = b(1) + 1j*b(2)
    bint = [bint(1,1)+1j*bint(2,1),bint(1,2)+1j*bint(2,2)]
end   

if 0
    [x,y] = xy_(1,N,sigma);
    [Z,dZ,Info] = regress_ols(y,x,'regress');
    fprintf('\n1-D complex (2-D) using regress\n')
    print_(Z,Info,'regress')
    
    fprintf('1-D complex (2-D) using regress-real\n')
    [Z,dZ,Info] = regress_ols(y,x,'regress-real');
    print_(Z,Info,'regress-real')

end

if 0
    [x,y] = xy_(2,N,sigma);
    [Z,dZ,Info] = regress_ols(y,x,'regress');
    fprintf('\n2-D complex (4-D) using regress\n')
    print_(Z,Info,'regress')
    
    [Z,dZ,Info] = regress_ols(y,x,'regress-real');
    fprintf('2-D complex (4-D) using regress-real\n')
    print_(Z,Info,'regress-real')
end

if 0
    fprintf('2-D complex (4-D) using regress-real\n')
    Ne = 5000;
    for i = 1:Ne
        [x,y] = xy_(2,N,sigma);
        [Z(i,:),dZ(i,:),Info] = regress_ols(y,x,'regress-real');
        ZCL95u(i,:) = Info.ZCL95u;
        ZCL95l(i,:) = Info.ZCL95l;
        %print_(Z,Info,'regress-real')
    end
    % See Devore, Probability and Statistics for Engineering and the Sciences, 8th edition Figure 7.3.
    Nxr = sum(real(ZCL95u(:,1)) < 1 | real(ZCL95l(:,1)) > 1);
    Nxi = sum(imag(ZCL95u(:,1)) < 1 | imag(ZCL95l(:,1)) > 1);
    fprintf('\nTest of confidence intervals for 2-D complex. Excpected answer for Ne -> Inf is 0.05.\n')
    fprintf('Fraction of real(Z1) estimates for which its error bars do not overlap with real(Z1) actual: %0.2f\n',Nxr/Ne);
    fprintf('Fraction of imag(Z1) estimates for which its error bars do not overlap with imag(Z1) actual: %0.2f\n',Nxi/Ne);

    Nxr = sum(real(ZCL95u(:,2)) < 1 | real(ZCL95l(:,2)) > 1);
    Nxi = sum(imag(ZCL95u(:,2)) < 1 | imag(ZCL95l(:,2)) > 1);
    fprintf('Fraction of real(Z2) estimates for which its error bars do not overlap with real(Z2) actual: %0.2f\n',Nxr/Ne);
    fprintf('Fraction of imag(Z2) estimates for which its error bars do not overlap with imag(Z2) actual: %0.2f\n',Nxi/Ne);
end

function [x,y] = xy_(dim,N,sigma)

    er = sigma*randn(N,1);
    ei = sigma*randn(N,1);

    if dim == 1
        %x = randn(N,1);
        x = randn(N,1) + 1j*randn(N,1);
        y = x + er + 1j*ei;
        return;
    end
    
    x1 = randn(N,1);
    x1 = x1 + 1j*x1;
    x2 = randn(N,1);
    x2 = x2 + 1j*x2;
    
    x = [x1,x2];
    y = x1 + x2 + er + 1j*ei;

end

function print_(Z,Info,fun)
    star = '';
    if strcmp(fun,'regress')
        star = '*';
    end

    if length(Z) == 1
        data = [Info.ZCL95l,Z,Info.ZCL95u];
        fprintf('real([Z-, Z, Z+]) = [%5.2f, %5.2f, %5.2f]\n',real(data));
        fprintf('imag([Z-, Z, Z+]) = [%5.2f, %5.2f, %5.2f]%s\n',imag(data),star);
        return;
    end

    for i = 1:size(Z,2)
        data = [Info.ZCL95l(i),Z(i),Info.ZCL95u(i)];
        fprintf('real([Z%d-, Z%d, Z%d+]) = [%5.2f, %5.2f, %5.2f]\n',i,i,i,real(data));
        fprintf('imag([Z%d-, Z%d, Z%d+]) = [%5.2f, %5.2f, %5.2f]%s\n',i,i,i,imag(data),star);
    end
end
