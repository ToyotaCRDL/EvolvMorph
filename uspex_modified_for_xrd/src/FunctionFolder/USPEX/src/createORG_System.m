function createORG_System(inputFile)
global ORG_STRUC
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];
collectForces = python_uspex(getPy, ['-f ' inputFile ' -b collectForces -c 1']);
if ~isempty(collectForces)
ORG_STRUC.collectForces = str2num(collectForces);
end
fixRndSeed = python_uspex(getPy, ['-f ' inputFile ' -b fixRandSeed -c 1']);
if ~isempty(fixRndSeed)
ORG_STRUC.fixRndSeed = str2num(fixRndSeed);
end
optType = python_uspex(getPy, ['-f ' inputFile ' -b optType -c 1']);
if ~isempty(optType)
if optType(1) == '-'  
ORG_STRUC.opt_sign = -1;
optType(1) = [];
end
optType(end) = [];
if ~isempty(str2num(optType)) 
ORG_STRUC.optType=abs(str2num(optType));
else
if strcmpi(optType, 'enthalpy')
ORG_STRUC.optType = 1;
elseif strcmpi(optType, 'volume')
ORG_STRUC.optType = 2;
elseif strcmpi(optType, 'hardness')
ORG_STRUC.optType = 3;
elseif strcmpi(optType, 'struc_order')
ORG_STRUC.optType = 4;
elseif strcmpi(optType, 'aver_dist')
ORG_STRUC.optType = 5;
elseif strcmpi(optType, 'diel_sus')
ORG_STRUC.optType = 6;
elseif strcmpi(optType, 'gap')
ORG_STRUC.optType = 7;
elseif strcmpi(optType, 'diel_gap')
ORG_STRUC.optType = 8;
elseif strcmpi(optType, 'mag_moment')
ORG_STRUC.optType = 9;
elseif strcmpi(optType, 'quasientropy')
ORG_STRUC.optType = 10;
elseif strcmpi(optType, 'powerfactor')
ORG_STRUC.optType = 14;


### added on 210718 ###
elseif strcmpi(optType, 'cos_similarity_xrd_volchanged')
ORG_STRUC.optType = 104;

### added on 210802 ###
elseif strcmpi(optType, 'exp_cos_similarity_xrd_volchanged')
ORG_STRUC.optType = 105;




#######################
end
end
#optType_OK=[1:10,14,1101:1111]; 
optType_OK=[1:10,14,104,105,1101:1111];    # modified on 210713 210718 210802 210820 210830 211227
if isempty( find( ORG_STRUC.optType==optType_OK) )
fprintf(fp, 'you did not chose a valid optType. Program STOPS...\n');
disp(['you did not chose a valid optType. Program STOPS...']);
quit;
end
end
calculationType = python_uspex(getPy, ['-f ' inputFile ' -b calculationType -c 1']);
if ~isempty(calculationType)
if calculationType(1) == '-'  
calculationType(1) = [];
ORG_STRUC.dimension = -1*str2num(calculationType(1));
else
ORG_STRUC.dimension = str2num(calculationType(1));
end
ORG_STRUC.molecule  = str2num(calculationType(2));
ORG_STRUC.varcomp   = str2num(calculationType(3));
end
pluginType = python_uspex(getPy, ['-f ' inputFile ' -b pluginType -c 1']);
if ~isempty(pluginType)
ORG_STRUC.pluginType = str2num(pluginType);
if 1==ORG_STRUC.pluginType
ORG_STRUC.startNextGen = 0;
end
else
ORG_STRUC.pluginType = 0;
end
numIons = python_uspex(getPy, ['-f ' inputFile ' -b numSpeci -e EndNumSpeci'], 1);
ORG_STRUC.numIons = numIons;
if ORG_STRUC.molecule == 1
ORG_STRUC.numMols = ORG_STRUC.numIons;
end
atomType = python_uspex(getPy, ['-f ' inputFile ' -b atomType -e EndAtomType']);
N_type   = size(ORG_STRUC.numIons,2);
atomType1 = GetElement(N_type, atomType);
for i=1:length(atomType1)
if atomType1(i)==0
atomType1(i)=[];
end
end
ORG_STRUC.atomType = atomType1;
valences = python_uspex(getPy, ['-f ' inputFile ' -b valences -e EndValences']);
if isempty(valences) 
ORG_STRUC.valences = zeros(1,length(ORG_STRUC.atomType));
for i = 1 : length(ORG_STRUC.atomType)
ORG_STRUC.valences(i) = str2num(valence(ORG_STRUC.atomType(i)));
end
else
ORG_STRUC.valences = str2num(valences);
end
topologyRandom = python_uspex(getPy, ['-f ' inputFile ' -b topologyRandom -c 1']);
if isempty(topologyRandom)
ORG_STRUC.topologyRandom = 0;
else
ORG_STRUC.topologyRandom = str2num(topologyRandom);
end
coordinatingNumbers = python_uspex(getPy, ['-f ' inputFile ' -b coordinatingNumbers -e EndCoordinatingNumbers']);
if isempty(coordinatingNumbers) 
ORG_STRUC.coordinatingNumbers = 'NONE';
else
ORG_STRUC.coordinatingNumbers = str2num(coordinatingNumbers);
end
NvalElectrons = python_uspex(getPy, ['-f ' inputFile ' -b valenceElectr -e EndValenceElectr']);
if isempty(NvalElectrons) 
ORG_STRUC.NvalElectrons = zeros(1,length(ORG_STRUC.atomType));
for i = 1 : length(ORG_STRUC.atomType)
ORG_STRUC.NvalElectrons(i) = str2num(valenceElectronsNumber(ORG_STRUC.atomType(i)));
end
else
ORG_STRUC.NvalElectrons = str2num(NvalElectrons);
end
minAt = python_uspex(getPy, ['-f ' inputFile ' -b minAt -c 1']);
maxAt = python_uspex(getPy, ['-f ' inputFile ' -b maxAt -c 1']);
if     ORG_STRUC.varcomp==0 & ORG_STRUC.dimension==3
if isempty(minAt) && isempty(maxAt)
elseif   ~isempty(minAt) && ~isempty(maxAt)
ORG_STRUC.minAt = str2num(minAt);
ORG_STRUC.maxAt = str2num(maxAt);
elseif   isempty(minAt) || isempty(maxAt)
status = 'Please specify BOTH maximum and minimum amount of atoms in the structure for variable composition calculation'
quit;
end
elseif ORG_STRUC.varcomp==1 & ORG_STRUC.dimension==3
if isempty(minAt)
status = 'Please specify the minimum amount of atoms in the structure for variable composition calculation'
quit;
end
ORG_STRUC.minAt = str2num(minAt);
if isempty(maxAt)
status = 'Please specify the maximum amount of atoms in the structure for variable composition calculation'
quit;
end
ORG_STRUC.maxAt = str2num(maxAt);
end
ExternalPressure = python_uspex(getPy, ['-f ' inputFile ' -b ExternalPressure -c 1']);
if ~isempty(ExternalPressure)
ORG_STRUC.ExternalPressure = str2num(ExternalPressure);
else
ORG_STRUC.ExternalPressure = 0;
end
minVectorLength = python_uspex(getPy, ['-f ' inputFile ' -b minVectorLength -c 1']);
if isempty(minVectorLength)
minVectorLength = 0;
Vector = zeros(1,length(ORG_STRUC.atomType));
for i = 1 : length(ORG_STRUC.atomType)
Vector(i)= 1.8*str2num(covalentRadius(ceil(ORG_STRUC.atomType(i))));
end
if ORG_STRUC.varcomp == 0
minVectorLength = max(Vector);
else
minVectorLength = min(Vector);
end
ORG_STRUC.minVectorLength = minVectorLength; 
else
ORG_STRUC.minVectorLength = str2num(minVectorLength);
end
hardCore = python_uspex(getPy, ['-f ' inputFile ' -b IonDistances -e EndDistances']);
if isempty(hardCore)
tempFolder = [ORG_STRUC.homePath '/CalcFoldTemp'];
codeFolder = [ORG_STRUC.USPEXPath '/FunctionFolder/Tool'];
V_atom     = zeros(1,length(ORG_STRUC.atomType));
hardCore   = zeros(length(ORG_STRUC.atomType), length(ORG_STRUC.atomType));
for i = 1 : length(ORG_STRUC.atomType)
V_atom(i) = calcDefaultVolume(1, ORG_STRUC.atomType(i), ORG_STRUC.ExternalPressure, 0, tempFolder, codeFolder);
end
for i = 1 : length(ORG_STRUC.atomType)
for j = i:length(ORG_STRUC.atomType)
if ORG_STRUC.molecule == 0
hardCore(i,j) = min( 0.2*(V_atom(i)^(1/3) + V_atom(j)^(1/3)), 1.2 );
else
hardCore(i,j) = 0.45*(V_atom(i)^(1/3) + V_atom(j)^(1/3));
end
end
end
else
hardCore = str2num(hardCore);
end
if size(hardCore,1) == size(hardCore,2)   
for i = 1 : length(ORG_STRUC.atomType)
for j = i : length(ORG_STRUC.atomType)
ORG_STRUC.minDistMatrice(i,j) = hardCore(i,j);
ORG_STRUC.minDistMatrice(j,i) = hardCore(i,j);
end
end
end
ORG_STRUC.goodBonds = zeros(length(ORG_STRUC.atomType));
goodBonds = python_uspex(getPy, ['-f ' inputFile ' -b goodBonds -e EndGoodBonds']);
if isempty(goodBonds)
gB = calcDefaultGoodBonds(ORG_STRUC.atomType); 
else
gB = str2num(goodBonds);
end
if size(gB,1) == 1        
for i = 1 : length(ORG_STRUC.atomType)
for j = 1 : length(ORG_STRUC.atomType)
ORG_STRUC.goodBonds(i,j) = gB;
end
end
else
for i = 1 : length(ORG_STRUC.atomType)
for j = i : length(ORG_STRUC.atomType)
ORG_STRUC.goodBonds(i,j) = gB(i,j);
ORG_STRUC.goodBonds(j,i) = gB(i,j);
end
end
end
magRatioVal = python_uspex(getPy, ['-f ' inputFile ' -b magRatio -e EndMagRatio']);
if isempty(magRatioVal)
ORG_STRUC.magRatio(1:7)=[0.1, 0.9/4, 0.9/4, 0.9/4, 0.9/4, 0, 0];
else
ORG_STRUC.magRatio = str2num(magRatioVal);
ORG_STRUC.magRatio = abs(ORG_STRUC.magRatio(1:7))/sum(abs(ORG_STRUC.magRatio));
end
ldaU = python_uspex(getPy, ['-f ' inputFile ' -b ldaU -e EndLdaU']);
if isempty(ldaU)
ORG_STRUC.ldaU(1:length(ORG_STRUC.atomType))=0;
else
ORG_STRUC.ldaU = str2num(ldaU);
if size(ORG_STRUC.ldaU, 1) == 1;
ORG_STRUC.ldaU(2,:) = 0;
end
end
mlPeerFilter = python_uspex(getPy, ['-f ' inputFile ' -b mlPeerFilter -c 1']);
if isempty(mlPeerFilter)
ORG_STRUC.mlPeerFilter = 0;
else
ORG_STRUC.mlPeerFilter = str2num(mlPeerFilter);
end
mlMinPopulation = python_uspex(getPy, ['-f ' inputFile ' -b mlMinPopulation -c 1']);
if isempty(mlMinPopulation)
ORG_STRUC.mlMinPopulation = 128;
else
ORG_STRUC.mlMinPopulation = str2num(mlMinPopulation);
end
pickUpYN = python_uspex(getPy, ['-f ' inputFile ' -b pickUpYN -c 1']);
if ~isempty(pickUpYN)
ORG_STRUC.pickUpYN = str2num(pickUpYN);
end
pickUpGen = python_uspex(getPy, ['-f ' inputFile ' -b pickUpGen -c 1']);
if ~isempty(pickUpGen)
ORG_STRUC.pickUpGen = str2num(pickUpGen);
end
if ORG_STRUC.pickUpYN == 0
ORG_STRUC.pickUpGen = 1; 
end
pickUpFolder = python_uspex(getPy, ['-f ' inputFile ' -b pickUpFolder -c 1']);
if ~isempty(pickUpFolder)
ORG_STRUC.pickUpFolder = str2num(pickUpFolder);
end
