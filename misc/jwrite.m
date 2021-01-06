%function jwrite(S)
% http://www.mtnet.info/docs/jformat.txt
File = load('data/Middelpos/Middelpos.mat');
S = File.S;
S{1}.Options.info.longitude = 20.1405;
S{1}.Options.info.latitude = -31.5436;
S{1}.Options.info.elevation = 1142;
S{1}.Options.info.azimuth = 0;

fe = S{1}.fe;

Z = S{1}.Z;

rho = abs(Z).^2./(5*fe/S{1}.Options.info.timedelta);
phi = (180/pi)*S{1}.Phi;
fname = [S{1}.Options.info.stationid, '.dat'];

comment = '# rho in ohms/m, Z in mV/km/nT, phi in degrees';
fid = fopen(fname,'w');
fprintf(fid,'%s\n',comment);
fprintf(fid,'>STATION   = %s\n',S{1}.Options.info.stationid);
fprintf(fid,'>MISDAT    = 0.100000E+33\n');
fprintf(fid,'>AZIMUTH   = %.2f\n',S{1}.Options.info.azimuth);
fprintf(fid,'>LATITUDE  = %.2f\n',S{1}.Options.info.latitude);
fprintf(fid,'>LONGITUDE = %.4f\n',S{1}.Options.info.longitude);
fprintf(fid,'>ELEVATION = %.1f\n',S{1}.Options.info.elevation);
fprintf(fid,'%s\n',S{1}.Options.info.stationid);
comp = {'RXX','RXY','RYX','RYY'};
for c = 1:length(comp)
    fprintf(fid,'%s\n%d\n',comp{c}, length(fe));
    for i = 2:size(rho,1)
        fprintf(fid,'%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n',...
                -fe(i),...
                rho(i,c), phi(i,c),...
                0.100000E+33,0.100000E+33,...
                0.100000E+33,0.100000E+33, 1, 1);
    end
end
comp = {'ZXX','ZXY','ZYX','ZYY'};
for c = 1:length(comp)
    fprintf(fid,'%s\n%d\n',comp{c}, length(fe));
    for i = 2:size(rho,1)
        fprintf(fid,'%10.3e %10.3e %10.3e %10.3e %10.3e\n',...
                -fe(i),...
                real(Z(i,c)), imag(Z(i,c)),...
                0.100000E+33, 1);
    end
end
fclose(fid);
type(fname)

