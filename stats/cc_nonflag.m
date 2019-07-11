function c = cc_nonflag(x,y)
%CC_NONFLAG

for k = 1:size(x,2)
    c(k) = corr(x(:,k),y(:,k),'rows','complete');
end

end

