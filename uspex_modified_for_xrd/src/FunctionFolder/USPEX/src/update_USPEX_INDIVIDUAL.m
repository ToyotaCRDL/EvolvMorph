function update_USPEX_INDIVIDUAL(POPULATION, resFolder, gen, atomType)
global USPEX_STRUC
global FORCE_STRUC
ID = POPULATION.Number;
USPEX_STRUC.POPULATION(ID).gen         =  gen; 
USPEX_STRUC.POPULATION(ID).ToCount     =  1; 
USPEX_STRUC.POPULATION(ID).Fitness     = 1000; 
Common_var = {'COORDINATES','LATTICE','numIons','FINGERPRINT','order',...
'Parents','S_order','symg','K_POINTS','Enthalpies','struc_entr','howCome'};
for i = 1:length(Common_var)
eval(['USPEX_STRUC.POPULATION(ID).' Common_var{i} ' = POPULATION.' Common_var{i} ';']);
end

### modified on 210415 ###
#Private_var = {'cell', 'numMols', 'MOLECULES', 'numBlocks', 'Surface_numIons',...
#'gap','hardness','powerfactor', 'mag_moment','magmom_ions', 'magmom_ini', 'ldaU','dielectric_tensor','elasticMatrix','elasticProperties'};
Private_var = {'cell', 'numMols', 'MOLECULES', 'numBlocks', 'Surface_numIons',...
'gap','hardness','powerfactor', 'mag_moment','magmom_ions', 'magmom_ini', 'ldaU','dielectric_tensor','elasticMatrix','elasticProperties',...
'cos_similarity_xrd_volchanged','exp_cos_similarity_xrd_volchanged'};

############# modified 210718 210802 210820 210830 211227 #############

for i = 1:length(Private_var)
if isfield(POPULATION, Private_var{i})
eval(['USPEX_STRUC.POPULATION(ID).' Private_var{i} ' = POPULATION.' Private_var{i} ';']);
end
end
if isfield(POPULATION, 'Vol')
volume = POPULATION.Vol;
else
volume = det(POPULATION.LATTICE);
end
USPEX_STRUC.POPULATION(ID).Vol = volume;
USPEX_STRUC.POPULATION(ID).density = calcDensity( POPULATION.numIons, atomType, volume);
safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
if isfield(POPULATION, 'RelaxStep')
FORCE_STRUC.POPULATION(ID).atomType     = atomType;
FORCE_STRUC.POPULATION(ID).numIons      = POPULATION.numIons;
FORCE_STRUC.POPULATION(ID).RelaxStep    = POPULATION.RelaxStep;
safesave([resFolder '/FORCE.mat'], FORCE_STRUC);
end
