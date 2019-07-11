function [Z,f,H,t,Ep] = transferfnTD(B,E,Nc,Na)
%TRANSFERFNTD - Compute transfer function using time domain method.
%
%   [Z,f,H,t,Ep] = TRANSFERFNTD(B,E,Nc,Na) - B and E are time series with
%   time varying with row number.  The * given below is used to indicate
%   convolution. H is the impulse response at time lags t and Z is its DFT
%   at frequency values f.
% 
%   Solves
%       E(:,1) = Hx*B(:,1)
%   or
%       E(:,1) = Hxx*B(:,1) + Hxy*B(:,2)
%   or
%       E(:,1) = Hxx*B(:,1) + Hxy*B(:,2)
%       E(:,2) = Hyx*B(:,2) + Hyy*B(:,2)
%   giving
%       Hxx = H(:,1)
%       Hxy = H(:,2)
%       Hyx = H(:,3)
%       Hyy = H(:,4)
%
%   The rows of H correspond to t = [-Na+1:Nc]';

if ~exist('Na')
    Na = 0; % The number of acausal coefficients
end
Ts = 0; % Shift input with respect to output (non-zero not implemented)

if size(B,2) == 1 && size(E,2) == 1
    % Ex = Z*Bx
    [T,X] = time_delay(E(:,1),B(:,1),Nc,Ts,Na,'pad');
    LIN   = basic_linear(X,T);
    Ep(:,1) = basic_linear(X,LIN.Weights,'predict');
    H = LIN.Weights(1:Nc+Na);
elseif (size(B,2) == 2 && size(E,2) == 1) || (size(B,2) == 2 && size(E,2) == 2)
    % Ex = Zxx*Bx + Zxy*By
    [T,Xxx] = time_delay(E(:,1),B(:,1),Nc,Ts,Na,'pad');
    [T,Xxy] = time_delay(E(:,1),B(:,2),Nc,Ts,Na,'pad');
    LINr1   = basic_linear([Xxx,Xxy],T); % r1 = Row 1 of transfer function
    Ep(:,1) = basic_linear([Xxx,Xxy],LINr1.Weights,'predict');
    H(:,1) = LINr1.Weights(1:Nc+Na);
    H(:,2) = LINr1.Weights(Nc+Na+1:end-1);    
    if size(B,2) == 2 && size(E,2) == 2
        % Ex = Zxx*Bx + Zxy*By
        % Ey = Zyx*Bx + Zyy*By
        [T,Xyx] = time_delay(E(:,2),B(:,1),Nc,Ts,Na,'pad');
        [T,Xyy] = time_delay(E(:,2),B(:,2),Nc,Ts,Na,'pad');
        LINr2   = basic_linear([Xyx,Xyy],T); % r2 = Row 2 of transfer function
        Ep(:,2) = basic_linear([Xyx,Xyy],LINr2.Weights,'predict');
        H(:,3) = LINr2.Weights(1:Nc+Na);
        H(:,4) = LINr2.Weights(Nc+Na+1:end-1);
    end
else
    error('Number of columns for input and output not supported.');
end

t = [-Na+1:Nc]';

if (Na == 0)
    % Zero because Na = 0 -> h(t<=0) = 0.
    H = [zeros(1,size(H,2));H];
    t = [0;t];
end

% Transfer Function
Z = fft(H);
N = size(Z,1);
Z = Z(1:floor(N/2)+1,:);
f = [0:floor(N/2)]'/N;