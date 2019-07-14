function Ep = Hpredict(t,H,B)
%HPREDICT
%

a = find(t == 0,1)-1;

B = [B;zeros(a,size(B,2))];

Ep(:,1) = filter(H(:,2),1,B(:,2));
Ep(:,2) = filter(H(:,3),1,B(:,1));
%Ep(1:size(Hi,1)-1,:) = NaN;

a = find(t == 0,1);

Ep = Ep(a:end,:);
%Ep = [Ep; NaN*ones(a-1,2)];


