function Error = Reading(code, Ind_No, indic)
global POP_STRUC
global ORG_STRUC
Error = 0;
cd (['CalcFold' num2str(indic)])
Gen = POP_STRUC.generation;
Step = POP_STRUC.POPULATION(Ind_No).Step;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
try
sumIons = sum(numIons);
catch
sumIons = 0;
end
maxErrors = ORG_STRUC.maxErrors;
Const_Lat = ORG_STRUC.constLattice; 
molecule  = ORG_STRUC.molecule; 
dimension = ORG_STRUC.dimension;
minDistMatrice = ORG_STRUC.minDistMatrice;
collectForces  = ORG_STRUC.collectForces;
ID = ['Gen' num2str(Gen) '-Ind' num2str(Ind_No) '-Step' num2str(Step)]; %For ERROR output
TotalStep = length([ORG_STRUC.abinitioCode]);
GoodBad = Read_AbinitCode(code, 0, ID);  
if GoodBad==1
if (dimension ~=-4)
[COORDINATES, LATTICE] = Read_Structure(code, Const_Lat);
else
PROTEINS_STRUC = Read_Tinker_Structure(POP_STRUC.backbone_atoms);
end
if (molecule ~=1) & (dimension ~=-4) 
if ~distanceCheck(COORDINATES, LATTICE, numIons, minDistMatrice)
POP_STRUC.POPULATION(Ind_No).Error = maxErrors + 1;
POP_STRUC.POPULATION(Ind_No).Enthalpies(end) = 100000;
disp('The structure after relaxation cannot satisfy distance constraint, please see the Warning ...');
newPOP.numIons     = numIons;
newPOP.LATTICE     = LATTICE;
newPOP.COORDINATES = COORDINATES;
brokenPOSCARStr = createPOSCARStr(POP_STRUC.POPULATION(Ind_No), newPOP, ORG_STRUC.atomType);
USPEXmessage(0, brokenPOSCARStr, 1);
Error = 1;
end
end
if POP_STRUC.POPULATION(Ind_No).Error <= maxErrors;
if (dimension ~=-4)
POP_STRUC.POPULATION(Ind_No).COORDINATES = COORDINATES;
POP_STRUC.POPULATION(Ind_No).LATTICE = LATTICE;
else
POP_STRUC.POPULATION(Ind_No).ANGLES      = PROTEINS_STRUC.angles;
POP_STRUC.POPULATION(Ind_No).RESIDUES    = PROTEINS_STRUC.residues;
POP_STRUC.POPULATION(Ind_No).SEC_STRUCT  = PROTEINS_STRUC.secondary_structure;
POP_STRUC.POPULATION(Ind_No).COORDINATES = PROTEINS_STRUC.backbone_crd_norm;
POP_STRUC.POPULATION(Ind_No).LATTICE     = PROTEINS_STRUC.lattice;
POP_STRUC.POPULATION(Ind_No).numIons     = PROTEINS_STRUC.numIons;
ORG_STRUC.weight   = 1;
ORG_STRUC.numIons  = PROTEINS_STRUC.numIons;
ORG_STRUC.atomType = PROTEINS_STRUC.atomType;
end
POP_STRUC.POPULATION(Ind_No).Enthalpies(Step) = Read_AbinitCode(code, 1, ID);
if molecule == 1
readMOL(Ind_No, code); 
end
if collectForces == 1
if ~isfield(POP_STRUC.POPULATION(Ind_No), 'RelaxStep')
POP_STRUC.POPULATION(Ind_No).RelaxStep(1).Energy = [];
POP_STRUC.POPULATION(Ind_No).RelaxStep(1).LATTICE= [];
POP_STRUC.POPULATION(Ind_No).RelaxStep(1).COORD = [];
POP_STRUC.POPULATION(Ind_No).RelaxStep(1).FORCE = [];
POP_STRUC.POPULATION(Ind_No).RelaxStep(1).STRESS= [];
end
allLattice=  Read_AbinitCode(code, 9, ID);
allCoords =  Read_AbinitCode(code, 8, ID);
allForces =  Read_AbinitCode(code, 7, ID);
allStess  =  Read_AbinitCode(code,-2, ID);
allEnergy =  Read_AbinitCode(code,-1, ID);
relaxStep = min( size(allCoords, 1), size(allForces, 1) )/sum(numIons);
for k = 1:relaxStep
if isempty(allLattice)
tmpLattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
elseif 3 == size(allLattice, 1)
tmpLattice = allLattice;
else
tmpLattice = allLattice( (k-1)*3+1:k*3, : );
end
POP_STRUC.POPULATION(Ind_No).RelaxStep(Step).Energy(k)                   = allEnergy(k);
POP_STRUC.POPULATION(Ind_No).RelaxStep(Step).LATTICE(1:3, 1:3, k)        = tmpLattice;
POP_STRUC.POPULATION(Ind_No).RelaxStep(Step).COORD  (1:sumIons, 1:3, k)  = allCoords( (k-1)*sumIons+1:k*sumIons, : );
POP_STRUC.POPULATION(Ind_No).RelaxStep(Step).FORCE  (1:sumIons, 1:3, k)  = allForces( (k-1)*sumIons+1:k*sumIons, : );
POP_STRUC.POPULATION(Ind_No).RelaxStep(Step).STRESS (1:3, 1:3, k)        = allStess(  (k-1)*3+1:k*3, : );
end
end
if ((ORG_STRUC.optType == 6) | (ORG_STRUC.optType == 8)) && (Step == TotalStep)   # modified 210630
#if ((ORG_STRUC.optType == 6) | (ORG_STRUC.optType == 8)) & (Step == TotalStep)   # modified 210630
POP_STRUC.POPULATION(Ind_No).dielectric_tensor = Read_AbinitCode(code, 3, ID);
end
if ((ORG_STRUC.optType == 7) & (Step == TotalStep)) || ((ORG_STRUC.optType == 8) & (Step == TotalStep-1))   # modified 210630
#if ((ORG_STRUC.optType == 7) & (Step == TotalStep)) | ((ORG_STRUC.optType == 8) & (Step == TotalStep-1))   # modified 210630
POP_STRUC.POPULATION(Ind_No).gap = Read_AbinitCode(code, 4, ID);
end


if ORG_STRUC.optType == 9
POP_STRUC.POPULATION(Ind_No).mag_moment = Read_AbinitCode(code, 5, ID);
end
if (ORG_STRUC.optType == 8) && (Step == TotalStep-1)   # modified 210630
#if (ORG_STRUC.optType == 8) & (Step == TotalStep-1)   # modified 210630
if (POP_STRUC.POPULATION(Ind_No).gap < 0.1) 
POP_STRUC.POPULATION(Ind_No).Step = Step + 2; 
end
end
if ((ORG_STRUC.optType == 14) && (Step == TotalStep))   # modified 210630
#if ((ORG_STRUC.optType == 14) & (Step == TotalStep))   # modified 210630
POP_STRUC.POPULATION(Ind_No).powerfactor = Read_AbinitCode(code, 2, ID);
end
if ((ORG_STRUC.optType > 1100) && (ORG_STRUC.optType<1112)) && (Step == TotalStep)   # modified 210630
#if ((ORG_STRUC.optType > 1100) & (ORG_STRUC.optType<1112)) & (Step == TotalStep)   # modified 210630
POP_STRUC.POPULATION(Ind_No).elasticMatrix=Read_AbinitCode(code, 6, ID);
elasticProperties=calcElasticProperties(POP_STRUC.POPULATION(Ind_No).elasticMatrix, ...
POP_STRUC.POPULATION(Ind_No).numIons, ORG_STRUC.atomType, det(POP_STRUC.POPULATION(Ind_No).LATTICE));
POP_STRUC.POPULATION(Ind_No).elasticProperties = elasticProperties;
end
if (ORG_STRUC.checkConnectivity == 1) && (Step == TotalStep)   # modified 210630
#if (ORG_STRUC.checkConnectivity == 1) & (Step == TotalStep)   # modified 210630
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
POP_STRUC.POPULATION(Ind_No).hardness = calcHardness(COORDINATES, LATTICE, numIons);
end
POP_STRUC.POPULATION(Ind_No).Step = Step + 1;
end
Clean_AbinitCode(code);
else
Error = 1;
if (code~=1) && (code~=8) && (code~=9) && (code~=14)   # modified 210630
#if (code~=1) & (code~=8) & (code~=9) & (code~=14)      # modified 210630
POP_STRUC.POPULATION(Ind_No).Error = maxErrors + 1;
else
POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 1; 
end
[a,b]=unix(['echo PROBLEM_reading Structure' num2str(Ind_No)]);
[nothing, nothing] = unix('pwd');
end
cd ..
function brokenPOSCARStr = createPOSCARStr(oldPOP, POP, atomType )
atomTypStr = [];
for i = 1:length(atomType)
atomTypStr = [atomTypStr, ' ', megaDoof( atomType(i) ) ' '];
end
brokenPOSCARStr= ['The structure after relaxation cannot satisfy distance constraint\n', '1.00\n'];
for i = 1:3
brokenPOSCARStr= [ brokenPOSCARStr, num2str(POP.LATTICE(i,:)) '\n'];
end
brokenPOSCARStr= [ brokenPOSCARStr, atomTypStr, '\n'];
brokenPOSCARStr= [ brokenPOSCARStr, num2str(POP.numIons),'\n', 'Direct\n'];
for i = 1:size( POP.COORDINATES, 1 )
brokenPOSCARStr= [ brokenPOSCARStr, num2str(POP.COORDINATES(i,:)) '\n' ];
end
brokenPOSCARStr= [brokenPOSCARStr, '\n\nStruture before relaxation, for checking\n', '1.00\n'];
for i = 1:3
brokenPOSCARStr= [ brokenPOSCARStr, num2str(oldPOP.LATTICE(i,:)) '\n'];
end
brokenPOSCARStr= [ brokenPOSCARStr, atomTypStr, '\n'];
brokenPOSCARStr= [ brokenPOSCARStr, num2str(oldPOP.numIons),'\n', 'Direct\n'];
for i = 1:size( oldPOP.COORDINATES, 1 )
brokenPOSCARStr= [ brokenPOSCARStr, num2str(oldPOP.COORDINATES(i,:)) '\n' ];
end
