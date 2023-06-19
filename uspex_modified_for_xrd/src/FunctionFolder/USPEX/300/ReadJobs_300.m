function ReadJobs_300()
global ORG_STRUC
global POP_STRUC
for indic = 1:ORG_STRUC.numParallelCalcs
whichInd = find([POP_STRUC.POPULATION(:).Folder]==indic);
if ~isempty (whichInd)
Step = POP_STRUC.POPULATION(whichInd).Step;
disp(['Structure' num2str(whichInd) ' step' num2str(Step) ' at CalcFold' num2str(indic) ]);
if POP_STRUC.POPULATION(whichInd).JobID
if (ORG_STRUC.platform > 0) || (ORG_STRUC.numParallelCalcs > 1)     # modified 210630
#if (ORG_STRUC.platform > 0) | (ORG_STRUC.numParallelCalcs > 1)     # modified 210630 
disp(['JobID=' num2str(POP_STRUC.POPULATION(whichInd).JobID) ]);
end
doneOr = checkStatusC(whichInd);
if doneOr
if (POP_STRUC.POPULATION(whichInd).JobID == 0.01) || (POP_STRUC.POPULATION(whichInd).JobID == 0.02)    # modified 210630
#if (POP_STRUC.POPULATION(whichInd).JobID == 0.01) | (POP_STRUC.POPULATION(whichInd).JobID == 0.02)       # modified 210630
POP_STRUC.POPULATION(whichInd).Step = length([ORG_STRUC.abinitioCode]) + 1;
else
Error = Reading(ORG_STRUC.abinitioCode(Step),whichInd, indic);
end
POP_STRUC.POPULATION(whichInd).JobID = 0;
if POP_STRUC.POPULATION(whichInd).Error > ORG_STRUC.maxErrors
POP_STRUC.POPULATION(whichInd).Done = 1;
POP_STRUC.POPULATION(whichInd).ToDo = 0;
POP_STRUC.POPULATION(whichInd).Folder=0;
elseif POP_STRUC.POPULATION(whichInd).Step > length ([ORG_STRUC.abinitioCode])
POP_STRUC.bodyCount = POP_STRUC.bodyCount + 1;
POP_STRUC.POPULATION(whichInd).Number = POP_STRUC.bodyCount;
LATTICE = POP_STRUC.POPULATION(whichInd).LATTICE;
numIons = POP_STRUC.POPULATION(whichInd).numIons;
COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES;
atomType = ORG_STRUC.atomType;
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, atomType);
[order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);
POP_STRUC.POPULATION(whichInd).order =  order;
POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
POP_STRUC.POPULATION(whichInd).struc_entr = structureQuasiEntropy(whichInd, atom_fing);
POP_STRUC.POPULATION(whichInd).S_order    = StructureOrder(FINGERPRINT, V, numIons, ORG_STRUC.deltaFing, ORG_STRUC.weight);

### added on 210416 ###
### added on 210718 ###
if ORG_STRUC.optType == 104 || ORG_STRUC.optType == -104
  POP_STRUC.POPULATION(whichInd).cos_similarity_xrd_volchanged = cos_similarity_xrd_with_target_volchanged(whichInd, POP_STRUC.bodyCount);
end

### added on 210802 ###
if ORG_STRUC.optType == 105 || ORG_STRUC.optType == -105
  POP_STRUC.POPULATION(whichInd).exp_cos_similarity_xrd_volchanged = exp_cos_similarity_xrd_with_target_volchanged(whichInd, POP_STRUC.bodyCount);
end

disp('Relaxation is done.')
disp(' ')
POP_STRUC.DoneOrder(whichInd) = POP_STRUC.bodyCount;
WriteIndividualOutput_300(whichInd);
POP_STRUC.POPULATION(whichInd).Folder=0;
POP_STRUC.POPULATION(whichInd).Done = 1;
POP_STRUC.POPULATION(whichInd).ToDo = 0;
end
safesave ('Current_POP.mat', POP_STRUC)
end
end
end
end
