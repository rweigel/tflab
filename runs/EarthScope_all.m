
print_figs = 1;

ids = {'VAQ58','ORF03','ORG03'};
for i = 1:length(ids)
    if 0
      EarthScope_main(ids{i})
    end
    if 1
      EarthScope_plot(ids{i},print_figs)
    end
end

