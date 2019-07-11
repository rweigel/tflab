function transferfnSummary(S,Savg,desc)

fprintf('___________________________________________________________________________\n');
fprintf('%s\n',desc);
fprintf('___________________________________________________________________________\n');

keys = fieldnames(Savg);

if isfield(S,'ao')

    % Get padding for fprintf.
    for f = 1:length(keys)
        L(f) = length(keys{f});
    end
    mL = max(L);
    for f = 1:length(keys)
        s{1,f} = repmat(' ',1, mL - length(keys{f}));
    end
    %
    
    fprintf('Ave ao: Using:\n');
    keys = fieldnames(Savg);
    for f = 1:length(keys)
        fprintf('%s %s %6.3f [%.3f,%.3f]\n',keys{f},s{1,f},Savg.(keys{f}).ao(1,2),Savg.(keys{f}).ao_CI95(:,2));
    end

    fprintf('Ave bo: Using:\n');
    for f = 1:length(keys)
        fprintf('%s %s %6.3f [%.3f,%.3f]\n',keys{f},s{1,f},Savg.(keys{f}).bo(1,2),Savg.(keys{f}).bo_CI95(:,2));
    end

end

metric = {'PE','CC','MSE'};

Savg.InSample = S;
keys = fieldnames(Savg);

% Get padding
for m = 1:length(metric)
    for j = 1:length(keys)
        L(m,j) = length(metric{m}) + length(keys{j});
    end
end
l = max(L(:));
for m = 1:length(metric)
    for f = 1:length(keys)
        s{m,f} = repmat(' ',1,l-L(m,f));
    end
end
%

comps = ['x','y'];
for m = 1:length(metric)
    for j = 1:size(S.(metric{m}),2)
        keys = fieldnames(Savg);
        %fprintf('Ave %s %s in-sample: %6.3f [%.3f,%.3f]\n',...
        %    comps(j),metric{m},mean(S.(metric{m})(1,j,:)),boot(S.(metric{m})(1,j,:),@mean,1000,50));
        for f = 1:length(keys)
            cistr = [metric{m},'_CI95'];
            if isfield(Savg.(keys{f}),cistr)
                fprintf('Ave %s %s %s: %s %7.3f [%.3f,%.3f]\n',...
                    comps(j),metric{m},keys{f},s{m,f},...
                    mean(Savg.(keys{f}).(metric{m})(1,j,:)),...
                    Savg.(keys{f}).(cistr)(:,j));
            else
                fprintf('Ave %s %s %s: %s %7.3f\n',...
                    comps(j),metric{m},keys{f},s{m,f},...
                    mean(Savg.(keys{f}).(metric{m})(1,j,:)));
            end
        end
    end
end

fprintf('___________________________________________________________________________\n')