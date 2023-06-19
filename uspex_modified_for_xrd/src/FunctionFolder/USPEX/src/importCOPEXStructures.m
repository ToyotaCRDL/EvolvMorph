function [addon_copex,OFF_STRUC]  = importCOPEXStructures(OFF_STRUC, generation)
global ORG_STRUC
addon_copex = 0;
if      ORG_STRUC.pluginType == 0
return;
elseif  mod( generation-1, ORG_STRUC.pluginType )>0
return;
end
calcType = str2num( [num2str(ORG_STRUC.dimension),num2str(ORG_STRUC.molecule),num2str(ORG_STRUC.varcomp)] );
resFolder = ORG_STRUC.resFolder;
if exist( [resFolder, '/COPEX.mat'], 'file' )
load([resFolder, '/COPEX.mat']);
else
disp('Warning : no COEPX.mat found!')
return
end
load([resFolder, '/COPEX.mat']);
addon_copex= length(COPEX_POP);
for addon = 1 : addon_copex
dummy = OFF_STRUC.POPULATION(end); 
OFF_STRUC.POPULATION(end+1) = dummy;
POP = createPOPULATION(COPEX_POP(addon), calcType);
if ~isempty( POP )
OFF_STRUC.POPULATION(end)   = POP;
OFF_STRUC.POPULATION(end).Parents = [];
info_parents = struct('parent',{}, 'enthalpy', {});
info_parents(1).parent = POP.Parents.parent;  
info_parents.enthalpy  = POP.Enthalpies(end);
OFF_STRUC.POPULATION(end).Parents = info_parents;
OFF_STRUC.POPULATION(end).howCome = ' COPEX ';
if ORG_STRUC.reoptOld
OFF_STRUC.POPULATION(end).Step = max(length(ORG_STRUC.abinitioCode)-1,1);
else
OFF_STRUC.POPULATION(end).Step = length(ORG_STRUC.abinitioCode)+1;
end
end
end
disp([num2str(addon_copex) ' structures were imported from COPEX']);
function POPULATION = createPOPULATION(COPEX_POP, calcType)
if calcType ==300 | calcType==301

### modified on 210415 ###
#POP = struct('COORDINATES',{},'INIT_COORD',{},'LATTICE',{},'INIT_LAT',{},'numIons',{},'INIT_numIons',{}, ...
#'struc_entr',{},'order',{},'S_order',{},'dielectric_tensor',{}, 'gap',{}, 'hardness',{}, 'mag_moment',{}, 'symg',{}, ...
#'FINGERPRINT',{}, 'K_POINTS',{}, 'Step', {}, 'Enthalpies', {}, 'Error',{},'Done',{},...
#'ToDo',{},'Parents',{},'howCome',{},'JobID',{},'Folder',{}, 'numBlocks', {},'Number',{});
POP = struct('COORDINATES',{},'INIT_COORD',{},'LATTICE',{},'INIT_LAT',{},'numIons',{},'INIT_numIons',{}, ...
'struc_entr',{},'order',{},'S_order',{},'dielectric_tensor',{}, 'gap',{}, 'hardness',{}, 'mag_moment',{}, 'symg',{}, ...
'FINGERPRINT',{}, 'K_POINTS',{}, 'Step', {}, 'Enthalpies', {}, 'Error',{},'Done',{},...
'ToDo',{},'Parents',{},'howCome',{},'JobID',{},'Folder',{}, 'numBlocks', {},'Number',{},...
'cos_similarity_xrd_volchanged',{},'exp_cos_similarity_xrd_volchanged',{});
###################### modified 210718 210802 210820 210830 #####

POPULATION(1) = QuickStart( POP );
Names = fieldnames(POP);
for i = 1:size(Names, 1);
if isfield(COPEX_POP, Names{i})
value = getfield(COPEX_POP, Names{i});
POPULATION(1) = setfield(POPULATION(1), Names{i}, value);
end
end
else
POP=[];
end
