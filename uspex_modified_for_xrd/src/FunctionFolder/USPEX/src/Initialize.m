function Initialize()
global ORG_STRUC
global POP_STRUC
global USPEX_STRUC
global POOL_STRUC
[nothing,homePath] = unix('pwd');
homePath(end) = [];
USPEXPath= homePath;
[nothing, uspexmode]=unix('echo -e $UsPeXmOdE');
if findstr(uspexmode,'exe')
[nothing,USPEXPath]= unix('echo -e $USPEXPATH');
USPEXPath(end)=[];
end
if ~exist('Seeds','dir')
[nothing, nothing] = unix('mkdir Seeds');
end
if ~exist('AntiSeeds','dir')
[nothing, nothing] = unix('mkdir AntiSeeds');
end
if ORG_STRUC.fixRndSeed > 0
rng( ORG_STRUC.fixRndSeed+1, 'twister' );
end
if (ORG_STRUC.dimension==3)  &&  (ORG_STRUC.varcomp==1)
N_T = size(ORG_STRUC.numIons,1);
splitting = zeros(1,N_T);
findSplit_VC(N_T, 0, ORG_STRUC.minAt, ORG_STRUC.maxAt, splitting);
end
calcType_str = [num2str(ORG_STRUC.dimension) '' num2str(ORG_STRUC.molecule) '' num2str(ORG_STRUC.varcomp)];
if str2num(calcType_str) < 0
calcType = -1*str2num(calcType_str);
calcType_str = ['M' num2str(calcType)];
end
path(path,[USPEXPath '/FunctionFolder/USPEX/' calcType_str]);
eval(['initialize_POP_STRUC_' calcType_str '()']);
POP_STRUC.resFolder = ORG_STRUC.resFolder;
USPEX_STRUC = struct('POPULATION',{}, 'GENERATION',{}, 'SYSTEM', {});
USPEX_STRUC(1).GENERATION = struct('quasiEntropy', {}, 'convex_hull', {}, 'composEntropy',{}, 'BestID',{}, 'ID',{});
USPEX_STRUC.GENERATION(1) = QuickStart(USPEX_STRUC.GENERATION);

### modified on 210415 ###
#USPEX_STRUC(1).POPULATION = struct('COORDINATES', {}, 'LATTICE', {}, 'numIons',{}, 'MOLECULES',{},...
#'symg',{}, 'howCome',{},'order',{}, 'Fitness', {}, 'FINGERPRINT', {}, 'K_POINTS', {}, ...
#'ToCount',{},'S_order',{}, 'gen',{}, 'struc_entr',{}, 'Enthalpies', {}, 'Parents',{},'Vol',{},...
#'cell',{},'Surface_numIons',{},'numBlocks',{},'gap',{},'hardness',{},'mag_moment',{},'magmom_ions',{}, 'magmom_ini',{}, 'ldaU',{}, 'dielectric_tensor',{});

USPEX_STRUC(1).POPULATION = struct('COORDINATES', {}, 'LATTICE', {}, 'numIons',{}, 'MOLECULES',{},...
'symg',{}, 'howCome',{},'order',{}, 'Fitness', {}, 'FINGERPRINT', {}, 'K_POINTS', {}, ...
'ToCount',{},'S_order',{}, 'gen',{}, 'struc_entr',{}, 'Enthalpies', {}, 'Parents',{},'Vol',{},...
'cell',{},'Surface_numIons',{},'numBlocks',{},'gap',{},'hardness',{},'mag_moment',{},'magmom_ions',{}, 'magmom_ini',{}, 'ldaU',{}, 'dielectric_tensor',{},...
'cos_similarity_xrd_volchanged', {}, 'exp_cos_similarity_xrd_volchanged', {});
#######################   modified 210713 210718 210802 210820  210830 211227 #######


USPEX_STRUC(1).SYSTEM    = struct( 'atomType', {}, 'atomType_symbol',{}, 'Fp_weight', {});
USPEX_STRUC.POPULATION(1) = QuickStart(USPEX_STRUC.POPULATION);
atomType = ORG_STRUC.atomType;
USPEX_STRUC.SYSTEM(1).atomType = atomType;
for i = 1:length(atomType)
USPEX_STRUC.SYSTEM(1).atomType_symbol{i} = megaDoof(atomType(i));
end
USPEX_STRUC.SYSTEM(1).Fp_weight = ORG_STRUC.weight;
safesave('Current_POP.mat', POP_STRUC);
safesave([POP_STRUC.resFolder '/USPEX.mat'], USPEX_STRUC);
if exist('POOL_STRUC','var')
safesave([POP_STRUC.resFolder '/POOL.mat'], POOL_STRUC);
end
[nothing, nothing] = unix(['echo "Gen   ID    Origin   Composition    Enthalpy   Volume  Density   Fitness   KPOINTS  SYMM  Q_entr A_order S_order Magmoment-Type" >>' ORG_STRUC.resFolder '/Individuals']);
if (ORG_STRUC.dimension==-4)  &&  (ORG_STRUC.molecule==0)  &&  (ORG_STRUC.varcomp==0)
[a,b]= unix(['echo "                                   (kcal/mol)  (A^3)   (g/cm^3)" >>' ORG_STRUC.resFolder '/Individuals']);
[a,b]= unix(['echo "  ID   Origin     Composition  Enthalpy      Volume(A^3)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/OUTPUT.txt']);
[a,b]= unix(['echo "                              (kcal/mol)" >>' ORG_STRUC.resFolder '/OUTPUT.txt']);
else
[a,b]= unix(['echo "                                      (eV)     (A^3)   (g/cm^3)" >>' ORG_STRUC.resFolder '/Individuals']);
[a,b]= unix(['echo "  ID   Origin     Composition  Enthalpy(eV)  Volume(A^3)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/OUTPUT.txt']);
end
[a,b]=unix(['echo " ID    Origin    Enthalpy   Parent-E   Parent-ID" >>' ORG_STRUC.resFolder '/origin']);
WriteGenerationStart();
