function Savg = transferfnAverage(S,opts,Ik)

Savg        = struct();
Savg.Mean   = struct();
%Savg.Huber  = struct();
%Savg.Median = struct();

Method = struct();
Method.Mean.Function   = @cmean;
%Method.Median.Function = @cmedian;
%Method.Huber.Function  = @cmlochuber;
%Method.MeanPEWeighted.Function  = @cmean;    
%Method.MeanPEWeighted.WeightStr = 'PE';        

% TODO: Account for this when computing error spectra.
a = opts.td.Ntrim;
%b = opts.td.window.width-opts.td.Ntrim+1;
b = size(S.In,1)-opts.td.Ntrim+1;

if nargin < 3
    Ik = [1:size(S.In_PSD,3)];
else
    if isfield(S,'ao')
        S.ao = S.ao(:,:,Ik);
        S.bo = S.bo(:,:,Ik);
    else
        S.Z  = S.Z(:,:,Ik);
    end        
end


function z = cmedian(z,dim)
    z = median(real(z),dim) + sqrt(-1)*median(imag(z),dim);
    % There seems to be a bug in MATLAB so that the median of a complex
    % number is not always correct. If
    % u = [-1+2*sqrt(-1),0+2*sqrt(-1),1+2*sqrt(-1),2+2*sqrt(-1)];
    %
    % Returns true:  real(mean(u,2)) == mean(real(u),2)
    % Returns false: real(median(u,2)) == median(real(u),2))
    %
    % So the median is calculated by splitting z into real and
    % imaginary components. For consistency, this is done for mean() even
    % though it is not needed.
end

function z = cmean(z,dim)
    % Not actually needed, but used for consistent notation.
    z = mean(real(z),dim) + sqrt(-1)*mean(imag(z),dim);
end

function z = cmlochuber(z,dim)
    if dim == 1
        z = ( mlochuber(real(z)) + sqrt(-1)*mlochuber(imag(z)) );
    end
    if dim == 2
        % .' is non-conjugate transpose. mlochuber assumes observations are
        % rows and variables are columns.
        z = ( mlochuber(real(z.')) + sqrt(-1)*mlochuber(imag(z.')) ).';
    end
    if dim == 3
        for i = 1:size(z,2)
            za(:,i) = cmlochuber(squeeze(z(:,i,:)),2);
        end
    end
end        

function w = weights(S,i,method)
    if strmatch('PE',method,'exact')
        if i < 3
            tmp = squeeze(S.PE(1,1,:))';
        else
            tmp = squeeze(S.PE(1,2,:))';
        end
        w = repmat(tmp,size(S.Z,1),1);
        w = w/mean(tmp);
    else
        w = ones([size(S.Z,1),size(S.Z,3)]);
    end
end

keys = fieldnames(Method);
for f = 1:length(keys)
    
    Savg.(keys{f}).In_PSD  = Method.(keys{f}).Function(S.In_PSD,3); 
    Savg.(keys{f}).Out_PSD = Method.(keys{f}).Function(S.Out_PSD,3); 

    if isfield(S,'ao')
        
        Savg.(keys{f}).ao(1,:) = Method.(keys{f}).Function(squeeze(S.ao),2);
        Savg.(keys{f}).bo(1,:) = Method.(keys{f}).Function(squeeze(S.bo),2);
        for i = 1:2
            Savg.(keys{f}).ao_CI95(:,i) = boot(squeeze(S.ao(1,i,:)),@(x) Method.(keys{f}).Function(x,1),1000,50);
            Savg.(keys{f}).bo_CI95(:,i) = boot(squeeze(S.bo(1,i,:)),@(x) Method.(keys{f}).Function(x,1),1000,50);        
        end

    else            
        tmp = getfield(Method,keys{f});
        afunc = tmp.Function;
        wstr = '';
        if isfield(tmp,'WeightStr')
            wstr = tmp.WeightStr;
        end
        
        T = struct();
        T.fe = S.fe;

        for i = 1:size(S.Z,2)
            z = squeeze(S.Z(:,i,:));
            if length(wstr) > 0
                W = weights(S,i,wstr);
                z = z.*W;
            end

            % Z
            T.Z(:,i)        = afunc(z,2);
            T.Z_StdErr(:,i) = std(z,0,2)/sqrt(size(z,2));
            T.Z_CI95(:,i,:) = boot(transpose(z),@(x) cmean(x,1),1000,159);
            
            % |Z|
            T.Zabs(:,i)  = afunc(abs(z),2);
            T.Zabs2(:,i) = abs(T.Z(:,i));

            T.Zabs_StdErr(:,i)  = std(abs(z),0,2)/sqrt(size(z,2));
            T.Zabs_CI95(:,i,:)  = boot(abs(z)',@(x) mean(x,1),1000,159);
            T.Zabs2_StdErr(:,i) = T.Zabs_StdErr(:,i);

            % Phi
            phi = atan2(imag(z),real(z));
            phi = unwrap(phi,[],2);

            T.Phi(:,i)   = afunc(phi,2);
            T.Phi2(:,i)  = atan2(imag(T.Z(:,i)),real(T.Z(:,i)));
            T.Phi_StdErr(:,i) = std(phi,0,2)/sqrt(size(z,2));            
            T.Phi_CI95(:,i,:)   = boot(phi',@(x) mean(x,1),1000,159);
            T.Phi2_StdErr(:,i) = T.Phi_StdErr(:,i);

            % H
            h = squeeze(S.H(:,i,:));
            if strcmp(keys{f},'Huber')
                % mlochuber too slow. TODO: Consider only a restricted range of h.
                T.H(:,i) = repmat(NaN,size(h,1),1);
                T.H_StdErr(:,i) = repmat(NaN,size(h,1),1);
            else
                T.H(:,i) = afunc(h,2);
                T.H_StdErr(:,i) = std(h,0,2)/sqrt(size(h,2));
                % Bootstrap is too slow for this.
            end
        end
        Savg = setfield(Savg,keys{f},T);
    end
end


for k = 1:size(S.In,3)

    keys = fieldnames(Savg);
    for f = 1:length(keys)
        T = getfield(Savg,keys{f});
        if isfield(S,'ao')
            Savg.(keys{f}).Predicted(:,1,k) = Savg.(keys{f}).ao(1)*S.In(:,1,k) + Savg.(keys{f}).bo(1)*S.In(:,2,k);
            Savg.(keys{f}).Predicted(:,2,k) = Savg.(keys{f}).ao(2)*S.In(:,1,k) + Savg.(keys{f}).bo(2)*S.In(:,2,k);
        else
            %keys{f}
            %[Savg.(keys{f}).Predicted(:,:,k),Savg.(keys{f}).Zi] = ...
            %        Zpredict(Savg.(keys{f}).fe,Savg.(keys{f}).Z,S.In(:,:,k));
            Savg.(keys{f}).Predicted(:,:,k) = Zpredict(Savg.(keys{f}).fe,Savg.(keys{f}).Z,S.In(:,:,k));
        end
        Savg.(keys{f}).PE(1,:,k)  = pe(S.Out(a:b,:,k),Savg.(keys{f}).Predicted(a:b,:,k));
        Savg.(keys{f}).MSE(1,:,k) = mse(S.Out(a:b,:,k),Savg.(keys{f}).Predicted(a:b,:,k));
        Savg.(keys{f}).CC(1,:,k)  = cc(S.Out(a:b,:,k),Savg.(keys{f}).Predicted(a:b,:,k));

        Savg.(keys{f}).Error_PSD(:,:,k) = smoothSpectra(S.Out(:,:,k)-Savg.(keys{f}).Predicted(:,:,k),opts);    
        Savg.(keys{f}).Coherence(:,:,k) = smoothCoherence(S.Out(:,:,k),Savg.(keys{f}).Predicted(:,:,k),opts);
        Savg.(keys{f}).SN(:,:,k)        = S.Out_PSD(:,:,k)./Savg.(keys{f}).Error_PSD(:,:,k);    
    end
    
    cstr = ['x','y'];
    for c = 1:size(Savg.(keys{1}).PE,2)
        fprintf('transferfnAverage.m: Interval %02d: %s PE in-sample: %6.3f; using: | ', k,cstr(c),S.PE(1,c,k));
        for f = 1:length(keys)
            fprintf('%s: %6.3f | ',keys{f},Savg.(keys{f}).PE(1,c,k));
        end
        fprintf('\n');
    end
    
end

keys = fieldnames(Savg);
for f = 1:length(keys)
    for i = 1:size(Savg.(keys{1}).PE,2)
        Savg.(keys{f}).PE_CI95(:,i) = boot( squeeze( Savg.(keys{f}).PE(1,i,:) ),@(x) mean(x,1),1000,50);
        Savg.(keys{f}).CC_CI95(:,i) = boot( squeeze( Savg.(keys{f}).CC(1,i,:) ),@(x) mean(x,1),1000,50);
        Savg.(keys{f}).MSE_CI95(:,i) = boot( squeeze( Savg.(keys{f}).MSE(1,i,:) ),@(x) mean(x,1),1000,50);
    end
end

end


