function exp_cos_similarity_xrd_volchanged = exp_cos_similarity_xrd_with_target_volchanged(Ind_No, bodyCount)

### newly added on 210416
### modified 210802

global ORG_STRUC
atomType = ORG_STRUC.atomType;
tempFolder = [ORG_STRUC.homePath '/CalcFoldTemp-allraw'];
codeFolder = [ORG_STRUC.USPEXPath '/FunctionFolder/Tool'];
specificFolder = [ORG_STRUC.homePath, '/', ORG_STRUC.specificFolder];

global POP_STRUC
lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
coor = POP_STRUC.POPULATION(Ind_No).COORDINATES;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
symg = POP_STRUC.POPULATION(Ind_No).symg;
order = POP_STRUC.POPULATION(Ind_No).order;
count = POP_STRUC.POPULATION(Ind_No).Number;

Write_POSCAR_arb(atomType, count, symg, numIons, lattice, coor, tempFolder);
cd(tempFolder);

### print latticevol 210806 ###
x = fopen('vol_default','w');
#disp(['latvol-------',num2str(ORG_STRUC.latVolume)]);
fprintf(x, num2str(ORG_STRUC.latVolume));
fclose(x);

### input copy 220913 ###
[nothing,nothing] = unix(['cp ' ORG_STRUC.homePath '/additional_INPUT.csv .']);

###############################

[nothing,nothing] = unix(['cp CONTCAR CONTCAR-' num2str(bodyCount)]);

[nothing,nothing] = unix(['python ' codeFolder '/poscar-to-cif-mpd2021.py CONTCAR-' num2str(bodyCount) ' CONTCAR-' num2str(bodyCount)]);

[nothing,nothing] = unix(['python ' codeFolder '/exp_cos_similarity_xrd_volchanged_uspex_220919.py '  specificFolder '/xrd_target.csv CONTCAR-' num2str(bodyCount) '.cif']);
a = fopen('exp_cos_similarity_xrd_volchanged.csv', 'r');
exp_cos_similarity_xrd_volchanged = str2num(fgets(a));
[nothing,nothing] = unix(['cp exp_cos_similarity_xrd_volchanged.csv exp_cos_similarity_xrd_volchanged.csv-' num2str(bodyCount)]);
[nothing,nothing] = unix(['cp CONTCAR-volchanged.cif CONTCAR-volchanged-' num2str(bodyCount) '.cif']);
[nothing,nothing] = unix(['rm CONTCAR CONTCAR-volchanged.cif exp_cos_similarity_xrd_volchanged.csv']);
cd ..;
