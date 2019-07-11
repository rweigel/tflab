function S = transferfnFD2(So,opts)

weightfn = 'bisquare';
for i = 2:size(So.In_FT,1) % Loop over evaluation frequencies
    for j = 1:2 % Loop over outputs
        cols = [1,2];
        if j == 2
            cols = [3,4];
        end
        x = So.In_FT{i,j,1};
        y = So.Out_FT{i,j,1};
        if 1
            Zk(i,cols,1) = regress(y,x);
            %Zk(i,cols,1) = mvregress(x,y);
            for k = 2:size(So.In_FT,3) % Loop over sections
                xk = So.In_FT{i,j,k};
                yk = So.Out_FT{i,j,k};
                x = cat(1,x,xk);
                y = cat(1,y,yk);
                Zk(i,cols,k) = regress(yk,xk);
                % Sk should match So.Mean.Z as this is doing the same
                % calculation as is done in transferfnFD.
                %Zk(i,cols,1) = mvregress(x,y);            
            end
        end
        
        % x and y now have all input and output FTs
        Zols1(i,cols) = x\y;
        Zrob1(i,cols)  = robustfit(x,y,[],[],'off');

        if 1
            Zols2(i,cols) = (ctranspose(x)*x)^(-1)*ctranspose(x)*y;
            Zols3(i,cols) = regress(y,x);


            Zrob1(i,cols)  = robustfit(x,y,[],[],'off');
            %Zmvr(i,cols)  = mvregress(x,y);
            %addpath('m/m-statrobust'); % So one can put breakpoints in MATLAB's code.

            Eols1{i,1} = y - x*Zols1(i,cols).';
            Erob1{i,1} = y - x*Zrob1(i,cols).';
            SN{i,1}    = mean(Eols1{i,1}.*conj(Eols1{i,1}))/mean(y.*conj(y));

            Zrob2(i,cols) = regress(y,x); % Start with OLS
            Zo = Zrob2(i,cols);
            stepmax = 50;
            step = 1;

            % MATLAB's robust fit seems to do this without accounting for
            % fact that residuals can be complex.
            [Q,R] = qr(x,0);
            E = x/R;
            % MATLAB's version:
            % h = min(.9999, sum(E.*E,2));
            % adjfactor = 1 ./ sqrt(1-h);
            % Correct version:
            hr = min(.9999, sum(real(E).*real(E),2));
            hi = min(.9999, sum(imag(E).*imag(E),2));
            adjfactor = 1./ sqrt(1-hr) + sqrt(-1)./ sqrt(1-hi);
            adjfactor = 1; % No leverage correction.
            laststep = 0;
            while 1
                % stop = 1 if either
                % |Zx-Zx_last|/|Zx_last| < eps
                % or
                % |Zy-Zy_last|/|Zy_last| < eps
                stop = any( abs(Zrob2(i,cols)-Zo) < eps*abs(Zo) );
                if step > 1 && (step == stepmax || stop == 1)
                    fprintf('Stopping at step %d for Te = %8d\n',step,round(1/So.fe(i)));
                    laststep = 1;
                end
                step = step + 1;
                Zo = Zrob2(i,cols);            
                E = y - x*transpose(Zrob2(i,cols));    % Residuals (Errors)
                Es = sort(abs(E));  % Sorted residuals
                p = 2; % Number of input variables.
                s = median(Es(p:end))/0.6745; % Median abs redisual estimation of standard deviation
                if strmatch(weightfn,'bisquare');
                    const = 4.685;
                    E = adjfactor.*E/(s*const);
                    w = (abs(E)<1) .* (1 - E.^2).^2; % Bi-square weights
                    if laststep
                        w(abs(E/(s*const)) > 2.8) = 0;
                    end
                elseif strmatch(weightfn,'huber');
                    const = 1.345;
                    E = adjfactor.*E/(s*const);
                    w = (abs(E)>1)./abs(E); % Huber weights                            
                    if laststep
                        w(abs(E/(s*const)) > 2.8) = 0;
                    end
                end
                Zrob2(i,cols) = regress(y.*sqrt(w),x.*sqrt([w,w]));
                if laststep
                    break;
                end
            end
        end
    end
end

S = struct();
fprintf('transferfnFD2: OLS\n');
S.OLS = createStruct(Zols1);

fprintf('transferfnFD2: Robust 1\n');
S.Robust1 = createStruct(Zrob1);

fprintf('transferfnFD2: Robust 2\n');
S.Robust2 = createStruct(Zrob2);

function Sx = createStruct(Z);

    Sx = struct();
    Sx.Z = Z;

    Nint = size(So.In,3);
    %Nint = size(Z,3);
    for k = 1:Nint

        Sx.fe = So.fe;

        N = size(So.In,1);
        f = [0:N/2]'/N; % Assumes N is even.
        Sx.H = Z2H(So.fe,Z,f);

        Sx.Predicted(:,:,k) = Zpredict(So.fe,Z,So.In(:,:,k));

        Sx.PE(1,:,k)  = pe_nonflag(So.Out(:,:,k),Sx.Predicted(:,:,k));
        Sx.MSE(1,:,k) = mse_nonflag(So.Out(:,:,k),Sx.Predicted(:,:,k));
        Sx.CC(1,:,k)  = cc_nonflag(So.Out(:,:,k),Sx.Predicted(:,:,k));

        Sx.In_PSD(:,:,k)     = smoothSpectra(So.In(:,:,k),opts);
        Sx.Out_PSD(:,:,k)    = smoothSpectra(So.Out(:,:,k),opts);
        Sx.Error_PSD(:,:,k)  = smoothSpectra(So.Out(:,:,k) - Sx.Predicted(:,:,k),opts);
        Sx.SN(:,:,k)         = Sx.Out_PSD(:,:,k)./Sx.Error_PSD(:,:,k);
        Sx.Coherence(:,:,k)  = smoothCoherence(So.Out(:,:,k),Sx.Predicted(:,:,k),opts);

        Sx.Phi(:,:,k) = atan2(imag(Sx.Z),real(Sx.Z));

        fprintf('transferfnFD2.m: Interval %d/%d: PE/CC/MSE of In_x = %.2f/%.2f/%.3f\n',k,Nint,Sx.PE(1,1,k),Sx.CC(1,1,k),Sx.MSE(1,1,k)); 
        fprintf('transferfnFD2.m: Interval %d/%d: PE/CC/MSE of In_y = %.2f/%.2f/%.3f\n',k,Nint,Sx.PE(1,2,k),Sx.CC(1,2,k),Sx.MSE(1,2,k)); 
    end
end

end