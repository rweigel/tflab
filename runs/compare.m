matfile1 = [scriptpath,sprintf('/data/Middelpos/Middelpos')];
matfile2 = [scriptpath,sprintf('/data/KAP03/KAP103')];

F1 = load(matfile1);
F2 = load(matfile2);

[B1,E1] = Middelpos_data();
[B2,E2] = KAP03_data('KAP103');

mean(B1)

ti = [1:size(B2,1)]';
for i = 1:size(B2,2)
    tg = find(~isnan(B2(:,i)));
    B2(:,i) = interp1(tg,B2(tg,i),ti);
end
mean(B2)

F1.S{1}.Options.description = 'Middelpos';
F2.S{1}.Options.description = 'KAP103';

transferfnZ_plot({F1.S{1},F2.S{1}},...
    struct('period_range',period_range,'title',ptitle,...
           'filename',filename,'savefmt',savefmt));

