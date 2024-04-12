function figureHTML(imgDir)

%title = 'Figures';
fid = fopen(fullfile(imgDir,'figures.html'),'w');

header = fileread(fullfile(scriptdir(),'figHTML.html'));
header = strrep(header,'__TITLE__','figures');
fprintf(fid,header);

dlist = dir(imgDir);
[~,idx] = sort([dlist.datenum]);
dlist = dlist(idx);
k = 1;
for i = 1:length(dlist)
    if dlist(i).isdir == 1 || ~endsWith(dlist(i).name,'.png')
        continue;
    end
    fprintf(fid,'<figure>\n');
    fprintf(fid,'  <a href="%s" target="_blank"><img src="%s"></a>\n',...
            dlist(i).name, dlist(i).name);
    fprintf(fid,'  <figcaption>Figure %d</figcaption>\n',k);
    fprintf(fid,'</figure>\n');
    k = k + 1;
end
fprintf(fid,"</body>");
fclose(fid);
logmsg("Wrote %s\n",fullfile(imgDir,'figures.html'));
end