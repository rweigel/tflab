function tmp = combineBootstrap(Bootstrap)
% COMBINEBOOTSTRAP  E.g., Bootstrap(j).ZCL95l => Bootstrap.ZCL95l(j,:)

for j = 1:length(Bootstrap)
  names = fieldnames(Bootstrap(j));
  for f = 1:length(names)
      if ~isempty(Bootstrap(j).(names{f}))
        tmp.(names{f})(j,:) = Bootstrap(j).(names{f});
        Nc = size(Bootstrap(j).(names{f}),2);
        isReal.(names{f}) = isreal(Bootstrap(j).(names{f}));
      end
  end
end
for j = 1:length(Bootstrap)
    for f = 1:length(names)
        if isempty(Bootstrap(j).(names{f}))
            if isReal.(names{f})
                one = 1;
            else
                one = (1+1j);
            end
            tmp.(names{f})(j,:) = one*nan*ones(1,Nc);
        end
    end
end
