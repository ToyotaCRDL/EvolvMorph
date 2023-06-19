function initialize_POP_STRUC_300()
global ORG_STRUC
global POP_STRUC
POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, 'resFolder', {},'generation',{}, ...
'DoneOrder',{}, 'bodyCount', {}, 'ranking',{},'bad_rank',{}, 'fitness', {}, 'convex_hull',{});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, 'LATTICE', {}, 'INIT_LAT', {}, 'numIons', {}, 'INIT_numIons', {},...
'struc_entr',{}, 'order',{},'dielectric_tensor',{}, 'gap',{}, 'hardness',{}, 'powerfactor',{}, 'FINGERPRINT',{}, 'symg',{}, 'K_POINTS',{},...
'mag_moment',{}, 'magmom_ions',{}, 'magmom_ini',{}, 'ldaU', {}, ...
'S_order',{}, 'Step', {}, 'Enthalpies', {}, 'Error',{},'Done',{},'ToDo',{},'Parents',{},'howCome',{},'JobID',{},'Folder',{}, 'Number',{},...   # added on 210405
'cos_similarity_xrd_volchanged', {}, 'exp_cos_similarity_xrd_volchanged', {})    # added on 210713   210718 210802

POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);
POP_STRUC(1).SOFTMUTATED = struct('FINGERPRINT',{}, 'mutatedAt', {}, 'fallBack', {});
POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'fingerprint',{},'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{},'numIons',{});
POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;
createCompostion_SingleBlock();
nSplit = length( ORG_STRUC.firstGeneSplit );
if nSplit == 0
error('ERROR:  The minAt and maxAt is too small, please check your input file ...');
end
ORG_STRUC.wrong_spacegroups = 0;
for i = 1: ORG_STRUC.initialPopSize
numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand(1)*nSplit) );
numIons = ORG_STRUC.numIons*numBlocks;
numMols = ORG_STRUC.numMols*numBlocks;
if ~isempty(ORG_STRUC.numMols)  
[candidate, lat] = Random_Init_310(i, numMols);
else
[candidate, lat] = Random_Init_300(i, numIons);
end
POP_STRUC.POPULATION(i).LATTICE = lat;
POP_STRUC.POPULATION(i).COORDINATES = candidate;
POP_STRUC.POPULATION(i).howCome = '  Random  ';
POP_STRUC.POPULATION(i).numIons = numIons;
end
if ORG_STRUC.wrong_spacegroups > 0
disp([' ']);
disp(['ATTENTION! In ' num2str(ORG_STRUC.wrong_spacegroups) ' / ' ...
num2str(ORG_STRUC.initialPopSize) ' cases actually generated symmetry was different.']);
disp([' ']);
end
pick_Seeds();
if ORG_STRUC.doFing
pickAntiSeeds();
end
Start_POP_300();
function createCompostion_SingleBlock()
global ORG_STRUC
if isempty(ORG_STRUC.minAt) | isempty(ORG_STRUC.maxAt)
ORG_STRUC.minAt = sum(ORG_STRUC.numIons);
ORG_STRUC.maxAt = sum(ORG_STRUC.numIons);
end
N_T = size(ORG_STRUC.numIons,1);
splitting = zeros(1,N_T);
findSplit_VC(N_T, 0, ORG_STRUC.minAt, ORG_STRUC.maxAt, splitting);
IPS = ORG_STRUC.initialPopSize;
fp = fopen('Seeds/compositions', 'w');
if exist('Seeds/Anti-compositions')
[nothing, nothing] = unix('mv Seeds/Anti-compositions Seeds/Anti-compositions-back');
end
for i=1:size(ORG_STRUC.firstGeneSplit,1)
for j=1:size(ORG_STRUC.firstGeneSplit,2)
fprintf(fp, '%4d', ORG_STRUC.firstGeneSplit(i,j));
end
fprintf(fp, '\n');
end
fclose(fp);
