function [COORDS, LATT] = Read_VASP_Structure()
[fid,message] = fopen('CONTCAR');
tmp = fgetl(fid); 
tmp = fgetl(fid); 
scale_factor = str2num(tmp);
LATT = zeros(3);
for i = 1 : 3
tmp = fgetl(fid);
LATT(i,:) = str2num(tmp);
end
LATT = LATT*scale_factor; 
tmp = fgetl(fid); 
ntype = str2num(tmp);
if isempty(ntype)   
tmp = fgetl(fid); 
ntype = str2num(tmp);
end
natom = sum(ntype);
tmp_mode = fgetl(fid);
if (tmp_mode(1) == 's') || (tmp_mode(1) == 'S')      # modified 210701
#if (tmp_mode(1) == 's') | (tmp_mode(1) == 'S')      # modified 210701
tmp_mode = fgetl(fid);
sss = fscanf(fid,'%g %g %g %s %s %s',[6,natom]);
else
sss = fscanf(fid,'%g %g %g',[3,natom]);
end
ss=sss';
COORDS = ss(:,1:3);
fclose(fid);
if (tmp_mode(1) == 'c') || (tmp_mode(1) == 'C') || (tmp_mode(1) == 'k') || (tmp_mode(1) == 'K')   # modified 210701 
#if (tmp_mode(1) == 'c') | (tmp_mode(1) == 'C') | (tmp_mode(1) == 'k') | (tmp_mode(1) == 'K')    # modified 210701

COORDS = COORDS*scale_factor;
COORDS = COORDS/LATT;
end
COORDS = COORDS - floor(COORDS);
