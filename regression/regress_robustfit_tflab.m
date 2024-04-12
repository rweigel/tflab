function [Z,dZ,Info] = regress_robustfit_tflab(ftE,ftB,varargin)
%REGRESS_ROBUSTFIT_TFLAB Robust regression with hard cut-off
%
%   REGRESS_ROBUSTFIT_TFLAB(ftE,ftB) uses robust regression with a Huber
%   weight function, a maximum number of steps of 50, zeps = sqrt(eps)
%   a hard cut-off of 2.8, and a SN stop value of 1000.
%
%   This algorithm downweights both the real and imaginary parts of ftB
%   and ftE when there is a large complex residual. For example, for
%   size(ftB,2) = 1, the system of 2*N equations being solved for in ftE =
%   Z*ftB, where ftE and ftB are complex, is
%
%     ftEr(1) = Zr*ftBr(1) - Zi*ftBi(1)
%     ...
%     ftEr(N) = Zr*ftBr(N) - Zi*ftBi(N)
%     ftEi(1) = Zr*ftBi(1) + Zi*ftBr(1)
%     ...
%     ftEi(N) = Zr*ftBi(N) + Zi*ftBr(N)
%
%     where N = size(ftE,1), r = real, i = imaginary.
%
%   In robust regression, one typically downweights each row individually
%   and there are 2N weights. However, N weights are computed based on
%   the residual of abs(ftE - ftE_regression).
%
%   REGRESS_ROBUSTFIT_TFLAB(ftE,ftB,opts) uses options in structure opts,
%   which has defaults of
%
%      opts = struct();
%      opts.weightfn = 'huber';  % or 'bisquare'
%      opts.stepmax  = 50;
%      opts.zeps     = sqrt(eps);
%      opts.hardcut  = 2.8;
%      opts.snstop   = 1000;
%      opts.verbose  = 0;        % 1 to print information during execution
%
%   References
%
%   * Andrews 1972, "Robust Estimates of Location"
%   * Hill and Holland 1977, "Two Robust Alternatives to Least-Squares
%     Regression.
%   * Huber 1981, "Robust statistics"
%   * Street et al., 1988, "A Note on Computing Robust Regression Estimates
%     via Iteratively Reweighted Least Squares"
%   * DuMouchel and O'Brien 1989, "Integrating a Robust Option Into a
%     Multiple Regression Environment.
%   * Fox and Weisberg 2013, "Robust Regression", Accessed 02/01/2019
%     http://users.stat.umn.edu/~sandy/courses/8053/handouts/robust.pdf
%
%   See also
%
%   * https://wis.kuleuven.be/stat/robust/Programs/LIBRA/contents-20160628.pdf
%   * https://www.mathworks.com/help/stats/robustfit.html
%   * https://www.gnu.org/software/gsl/doc/html/lls.html#robust-linear-regression
%     (GSL code and notation very similar to MATLAB's)

optsd = struct();
optsd.weightfn = 'huber';
optsd.stepmax = 50;
optsd.zeps = sqrt(eps);
optsd.hardcut = 2.8;
optsd.snstop = 1000;
optsd.verbose = 0;

if nargin < 3
    opts = optsd;
else
    if ~isempty(varargin{1})
        opts = varargin{1};
        for f = fieldnames(optsd)
            if ~isfield(opts,f)
                opts.(f) = optsd.(f);
            end
        end
    else
        opts = optsd;
    end
end

stats = struct();

p = size(ftB,2); % Number of input variables.

% Leverage
hf = 1;                      % leverage factor
H  = ftB*inv(ftB'*ftB)*ftB'; % Hat matrix
h  = diag(H);                % leverage value (0 <= h <= 1)
stats.Leverage = abs(h);

hf = sqrt(1-abs(h));
% hf is the leverage factor from Huber 1981. Use abs() to account for
% complex h (roundoff makes it small but non-zero). Equation 9.10 of
% Huber 1981 is
%    s*hf_i*psi(r_i/(s*hf_i))
% and using
%    w(x) = psi(x)/x (Equation 8.32)
% gives
%    w_i = psi(r_i/(s*hf_i))

Z(1,:) = regress_ols(ftE,ftB); % Intial estimate.

step = 2;
laststep = 0;
so = std(abs(ftE)); % Standard deviation of |output|.

while 1

    % Residuals from previous Z
    R = ftE - ftB*transpose(Z(step-1,:));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Robust estimate of standard deviation using median absolute residual
    % (Equation 7.2 of Huber 1981) who notes the motivation is Andrews
    % 1972. The factor of 0.6745 makes it unbiased when when residuals are
    % gaussian. Excludes the smallest p values as recommended by Hill and
    % Holland 1977, which was cited by DuMouchel and O'Brien 1989. However,
    % instead of using a factor of 2.1 as in Hill and Holland 1977, we use
    % a factor of 1.483 = 1/0.6745.
    %
    % TODO: Compare with madc() in
    % https://wis.kuleuven.be/stat/robust/Programs/LIBRA/contents-20160628.pdf
    % which includes an additional finite sample correction and uses
    % the same formula as https://ijpam.eu/contents/2014-91-3/7/7.pdf:
    % mad = median(abs(x-median(x)))/0.6745.
    Rsrt = sort(abs(R)); % Sorted absolute value of residuals
    s = median(Rsrt(p:end))/0.6745;
    %s = median(abs(R-median(R)))/0.6745; % More modern method?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if so/s > opts.snstop && laststep == 0
        if opts.verbose
            fprintf('Stopping at step %d because std(abs(ftE))/mad(abs(residuals)) > 1000\n',step);
        end
        if step == 2
            % No weights have been computed.
            W = ones(size(R));
            stats.W = W;
        end
        break;
    end

    if opts.verbose
        fprintf('Step %d: std(abs(ftE))/mad(abs(residuals)) = %.2g\n',step,so/s);
    end

    if strcmp(opts.weightfn,'bisquare')
        % W = (1-(R/const)^2)^2 for R <= const
        % W = 0 otherwise.
        const = 4.685;    % 95% efficiency when the errors are gaussian
        Rs = R/(s*const); % Normalize residuals by MAD then scale by const
        Rs = Rs./hf;      % As per Equation 9.10 of Huber 1981; see note above.
        W = zeros(size(R));
        W = (abs(Rs) < 1).*(1 - Rs.^2).^2; % Bi-square weights
        if laststep
            % Hard cut-off for adjfactor*|residual|/MAD > opts.hardcut
            Iz = abs(Rs*const) > opts.hardcut;
            W(Iz) = 0;
        end
    elseif strcmp(opts.weightfn,'huber')
        % W = 1 for |R| <= const
        % W = const/|R| for |R| > const
        const = 1.345;    % 95% efficiency when the errors are normal
        Rs = R/(s*const); % Normalize residuals by MAD then scale by const
        Rs = Rs./hf;      % As per Equation 9.10 of Huber 1981; see earlier note about Eqn. 9.10
        W = ones(size(Rs));
        W(abs(Rs) > 1) = 1./abs(Rs(abs(Rs) > 1)); % Huber weights
        if laststep
            % Hard cut-off for adjfactor*|residual|/MAD > opts.hardcut
            Iz = abs(Rs*const) > opts.hardcut;
            W(Iz) = 0;
        end
    else
        error('Weight function %s not known.',opts.weightfn);
    end

    if all(W == 0)
        warning('All weights are zero');
    end

    Z(step,:) = regress_ols(ftE.*sqrt(W),ftB.*repmat(sqrt(W),1,p));

    if laststep
        if opts.verbose
            fprintf('Step %d: Hard cut-off last-step zeroed %d values (%.2f%%)\n',step,sum(Iz),100*sum(Iz)/length(R));
        end
        break;
    end

    if step > 2
        ratio = abs(Z(step,:)-Z(step-1,:))./abs(Z(step-1,:));
        if opts.verbose
            fprintf('Step %d: ratio = abs(Z(step,:)-Z(step-1,:))./abs(Z(step-1,:)) = %.2g\n',step,ratio);
        end
        stop = any(ratio < opts.zeps);
        if step == opts.stepmax || stop == 1
            if stop && opts.verbose
                fprintf('Step %d: Last step because any(ratio) < opts.zeps\n',step);
            end
            if step == opts.stepmax && opts.verbose
                fprintf('Step %d: Last step because step == opts.stepmax = %d\n',step,opts.stepmax);
            end
            if step == opts.stepmax
                warning('Iteration limit of %d reached',opts.stepmax);
            end
            laststep = 1;
        end
    end
    stats.W(step,:) = W;

    step = step + 1;

end

Info = struct();
Info.Residuals_ = R;
Info.Weights = stats.W(end,:);

dZ = [];
Z = Z(end,:);

% Un-weighted residuals
R = ftE - ftB*transpose(Z);

