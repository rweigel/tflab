function [B,E,t,files] = LEMI_read(inpath,ext,start,stop,D)
% LEMI_READ  Read LEMI text files
%
%   [B,E,t,files] = LEMI_READ(inpath,ext,start,stop,D)
%
%   Assumes inpath contains LEMI text files with names that start with
%   YYYYMMDD and end with .ext where .ext is typically .tnn with nn indicating the 
%   location of the MT device e.g. 20210920000000.t07
%
%   start and stop are strings that specify the day range of files to read
%   (inclusive) and must have the format YYYYMMDD.
%
%   The E-field components are derived from the E-field voltages using
%
%   Ex = (E1-E3)/Dx
%   Ey = (E2-E4)/Dy
%
%   where Dx and Dy are from the input D = [Dx, Dy]. Dx and Dy are the distances
%   between the ends of each E-field sensors. The values can be derived from
%   the lengths L1, L2, L3, and L4 of the arms of each of the E-field sensors.
%   This data is available in the corresponding *.INF files. For example, if
%
%   L1 = 25, L2 = 50, L3 = 50, and L4 = 50, then
%
%     Dx = L1 + L3 = 75
%     Dy = L2 + L4 = 100
%
%   Output:
%   B = B-field data [Bx,By,Bz] [nT]
%   E = E-field data [E1,E2,E3,E4] [microvolt/m]
%   t = time [MATLAB datenumber]
%
%   Creates .mat version of data for each text files in the same directory.
    dlist = dir(inpath);
    Data = [];
    data = [];
    B = [];
    E = [];
    t = [];
    for i = 1:length(dlist)
        fnametxt = fullfile(dlist(i).folder, dlist(i).name);
        fnamemat = [fnametxt,'.mat'];
        if ~endsWith(dlist(i).name, ext)
            continue
        end
        try
            if str2num(dlist(i).name(1:8)) < str2num(start) || str2num(dlist(i).name(1:8)) > str2num(stop)
                fprintf('Skipping %s b/c out of time range.\n',fnametxt);
                continue
            end
        catch
            fprintf('Failed read; Error reading file %s\n',dlist(i).name);
            continue
        end
        if exist(fnamemat,'file')
            fprintf('Reading: %s\n',fnamemat);
            load(fnamemat);
        else
            try
                fprintf('Reading: %s\n',fnametxt);
                data = load(fnametxt);
                fprintf('Writing: %s\n',fnamemat);
                save(fnamemat,'data');
            catch ex
                fprintf('Failed read; skipping. Exception:\n');
                ex
            end

        end
        files{i} = fnametxt;
        if ~isempty(data)
            Data = [Data; data];
        end
    end
    if isempty(Data)
        return
    end

    B = Data(:,7:9);
    V = Data(:,12:15);
    Ex = (V(:,1)-V(:,3))/D(1);
    Ey = (V(:,2)-V(:,4))/D(2);
    E = [Ex,Ey];
    t = datenum(Data(:,1:6));
end