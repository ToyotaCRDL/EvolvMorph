function [sgrp_name] = anasym_stokes(lat, coord, numIons, atomType, tolerance)
global ORG_STRUC
getPy    =[ORG_STRUC.USPEXPath '/FunctionFolder/getInput.py'];
spgBINDIR=[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup'];
i = 1;
while i <= length(numIons)
if numIons(i) == 0
numIons(i) = [];
atomType(i) = [];
i = i - 1;
end
i = i + 1;
end
if (size(lat,1) == 3) & (size(lat,2) == 3)
lat = latConverter(lat);
end
lat(4:6) = lat(4:6)*180/pi; 
fp=fopen('sym.in','w');
fprintf(fp, 'testing\n');
fprintf(fp, '%6.4f accuracy\n', tolerance);
fprintf(fp, '2  form of lattice parameters: to be entered as lengths and angles\n');
fprintf(fp, '%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f a,b,c,alpha,beta,gamma\n', lat(1:6));
fprintf(fp, '1  form of primitive lattice vectors\n');
fprintf(fp, '1 0 0 primitive lattice vectors\n');
fprintf(fp, '0 1 0\n');
fprintf(fp, '0 0 1\n');
fprintf(fp, '%4d number of atoms in the primitive unit cell\n', sum(numIons));
for ii = 1 : length(numIons)
for jj = 1 : numIons(ii)
fprintf(fp, '%2d',ii);
end
end
fprintf(fp, ' type of each atom\n');
for coordLoop = 1 : size(coord,1)
fprintf(fp, '%12.6f %12.6f %12.6f\n', coord(coordLoop,:));
end
fclose(fp);
[a, b] = unix([spgBINDIR '/findsym_new < sym.in > sym.out']);
sgrp_name=python_uspex(getPy, '-f sym.out -b Space -c 3');
sgrp_name = regexprep(sgrp_name,'\n',' ');
[no_bombed, nothing] = unix('grep bombed sym.out');
if (no_bombed == 0) || isempty(sgrp_name)
a = 1;
end
fp1=fopen('symmetrized.cif','w');
if a ~= 0
status = ['Error while determining the crystal symmetry'];
sgrp_name = 0;
else
[a,b]=unix('echo end_of_file >> sym.out'); % to easily identify the end of file
toWrite = 0;
readNormal = 1;
currentType = '0';
typeA = 0;
handle = fopen('sym.out');
try
while 1
tmp = fgetl(handle);
if ~isempty(findstr(tmp, 'end_of_file'))
break; 
end
if readNormal == 0
coords = sscanf(tmp,'%*s %s %g %g %g %g');
if isempty(coords)
break; 
end
if isempty(findstr(coords(1), currentType))
typeA = typeA + 1;
currentType = coords(1);
end
s = megaDoof(atomType(typeA));
fprintf(fp1, '%4s %4s %8.4f %8.4f %8.4f %8.4f\n', s, s, coords(2:5));
end
if toWrite
fprintf(fp1, '%-40s\n', tmp);
end
if ~isempty(findstr(tmp, 'data_')) || ~isempty(findstr(tmp, 'loop_'))   # modified 210709
#if ~isempty(findstr(tmp, 'data_')) | ~isempty(findstr(tmp, 'loop_'))   # modified 210709
toWrite = 1 - toWrite;
end
if ~isempty(findstr(tmp, '_atom_site_occupancy'))
readNormal = 0;
toWrite = 0;
end
end
catch
status = ['Error while determining the crystal symmetry'];
sgrp_name = 0;
end
status = fclose(handle);
end
if sgrp_name == 0 
fprintf(fp1, '_symmetry_space_group_name_H-M \"P 1\" \n');
fprintf(fp1, '_symmetry_Int_Tables_number 1\n');
fprintf(fp1, '\n');
fprintf(fp1, '_cell_length_a %6.3f\n', lat(1));
fprintf(fp1, '_cell_length_b %6.3f\n', lat(2));
fprintf(fp1, '_cell_length_c %6.3f\n', lat(3));
fprintf(fp1, '_cell_angle_alpha %6.3f\n', lat(4));
fprintf(fp1, '_cell_angle_beta  %6.3f\n', lat(5));
fprintf(fp1, '_cell_angle_gamma %6.3f\n', lat(6));
fprintf(fp1, '\n');
fprintf(fp1, 'loop_\n');
fprintf(fp1, '_atom_site_label\n');
fprintf(fp1, '_atom_site_type_symbol\n');
fprintf(fp1, '_atom_site_fract_x\n');
fprintf(fp1, '_atom_site_fract_y\n');
fprintf(fp1, '_atom_site_fract_z\n');
fprintf(fp1, '__atom_site_occupancy\n');
iC = 1;
for typeA = 1 : length(numIons)
for i = 1 : numIons(typeA)
s = megaDoof(atomType(typeA));
tmp = [s ' ' s ' ' num2str(coord(iC,1),6) ' ' num2str(coord(iC,2),6) ' ' num2str(coord(iC,3),6) ' 1.0000'];
iC = iC + 1;
end
end
sgrp_name = '1';
end
fclose(fp1);
sgrp_name = str2num(sgrp_name);
