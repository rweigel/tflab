function figsave_(popts,ext)
    if popts.print == 0
        return;
    end
    if nargin == 1
        printName = popts.printOptions.printName;
    else
        % Remove TeX
        ext = regexprep(ext,'\{|\}','');
        ext = regexprep(ext,'\$','');
        printName = sprintf('%s-%s',popts.printOptions.printName,ext);
    end

    figsave(fullfile(popts.printOptions.printDir,printName),...
            popts.printOptions.export_fig,...
            popts.printOptions.printFormats);

end
