function S = tflab_uncertainty(S)

if iscell(S)
    for j = 1:length(S)
        S{j} = tflab_uncertainty(S{j});
    end
    return
end

if ~isfield(S,'Segment') 
    return;
end
if ~isfield(S.Segment,'Z')
    return;
end

% S.Segment.Z has size = [length(fe), Nc, Ns]
% where Nc is number of components of Z and Ns is number of time series
% segments used to estimate Z

% Loop over components of Z
for j = 1:size(S.Segment.Z,2) 

    % Extract all segment estimates of Z for component i
    tmp = squeeze(S.Segment.Z(:,j,:));

    re = real(tmp);
    im = imag(tmp);

    eb = (std(re,0,2)+1i*std(im,0,2))/sqrt(size(tmp,2)); 

    ZCL.Z.Normal.x_1sigma(:,j,:) = [S.Z(:,j)-eb/2,S.Z(:,j)+eb/2];
    ZCL.Z.Bootstrap.x_95(:,j,:) = bootcl(tmp, @mean, 1000, 0.025);
    ZCL.Z.Bootstrap.x_1sigma(:,j,:) = bootcl(tmp, @mean, 1000, 0.341);

    ZCL.Magnitude.Normal.x_1sigma(:,j,:) = ...
        propagateError(ZCL.Z.Normal.x_1sigma(:,j,:), S.Z(:,j),'magnitude');

    ZCL.Magnitude.Bootstrap.x_95(:,j,:) = ...
        propagateError(ZCL.Z.Bootstrap.x_95(:,j,:), S.Z(:,j),'magnitude');
        %bootcl(tmp, @(x) mean(abs(x)), 1000, 0.025);

    ZCL.Magnitude.Bootstrap.x_1sigma(:,j,:) = ...
        propagateError(ZCL.Z.Bootstrap.x_1sigma(:,j,:), S.Z(:,j),'magnitude');
        %bootcl(tmp, @(x) mean(abs(x)), 1000, 0.341);

end

for j = 1:size(S.Segment.Metrics.SN,2) 
    tmp = squeeze(S.Segment.Metrics.SN.Smoothed(:,j,:));
    SNCL.Normal.x_1sigma_normal(:,j,:) = std(tmp,0,2)/sqrt(size(tmp,2));
    SNCL.Bootstrap.x_1sigma(:,j,:) = bootcl(tmp, @mean, 1000, 0.341);
    SNCL.Bootstrap.x_95(:,j,:) = bootcl(tmp, @mean, 1000, 0.025);
end

S.ZCL = ZCL;
S.Metrics.SNCL = SNCL;

function E = propagateError(U,Z,type)
    U = squeeze(U);
    U = U(:,2) - U(:,1);
    reZ = abs(real(Z));
    reU = real(U);
    imZ = abs(imag(Z));
    imU = imag(U);
    if strcmp(type,'magnitude')
        dE = (reZ.*reU + imZ.*imU)./abs(Z);
        %dEx = abs(U); % Naive estimate
        E = [abs(Z)-dE/2,abs(Z)+dE/2];
    end
end
end