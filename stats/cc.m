function c = cc(x,y)
%CC

for k = 1:size(x,2)
    c(k) = corr(x(:,k),y(:,k));
end

end
