function Start()
global ORG_STRUC
global POP_STRUC
global POOL_STRUC
global ANTISEEDS
global USPEX_STRUC
global FORCE_STRUC
rand('state',sum(100*clock));
c = clock;
rand_init = round(c(5)+c(6));
for i = 1:rand_init
dummy = rand;
dummy = randperm(5);
end
disp('  ');
while 1
stillRunningInfo();
[nothing,homePath] = unix('pwd');
homePath(end) = [];
USPEXPath= homePath;
#[nothing, uspexmode]=unix('echo -e $UsPeXmOdE');
[nothing, uspexmode]=unix('echo $UsPeXmOdE');   # modiied on 210423
if findstr(uspexmode,'exe')
#[nothing,USPEXPath]= unix('echo -e $USPEXPATH');
[nothing,USPEXPath]= unix('echo $USPEXPATH');   # modified on 210423
USPEXPath(end)=[];
end
if exist('Current_POP.mat')
load Current_POP.mat
load Current_ORG.mat
checkCompatibility();
#if isempty(findstr(ORG_STRUC.homePath, homePath)) || length(ORG_STRUC.homePath)~=length(homePath)
if isempty(findstr(ORG_STRUC.homePath, homePath)) || length(ORG_STRUC.homePath)~=length(homePath)   # modified on 210423
mesg=[ORG_STRUC.homePath ' --> ' homePath ];
ORG_STRUC.homePath=homePath;
USPEXmessage(1011, mesg, 0);
end
ORG_STRUC.USPEXPath=USPEXPath;
cd (POP_STRUC.resFolder)
load USPEX.mat
if exist('FORCE.mat')
load FORCE.mat
end
if exist('POOL.mat')
load POOL.mat
end
cd ..
if exist('ANTISEEDS.mat')
load ANTISEEDS.mat
else
pickAntiSeeds();
end
if ORG_STRUC.fixRndSeed > 0 
rng( ORG_STRUC.fixRndSeed+POP_STRUC.generation, 'twister' );
end
else
inputFile = 'INPUT.txt';
createORGStruc(inputFile);
[a,b]=unix (['cp ' inputFile ' ' ORG_STRUC.resFolder '/Parameters.txt']);
path(path,[USPEXPath '/FunctionFolder/USPEX/' calcTypeStr]);
eval(['PreCheck_' calcTypeStr '()']);
ORG_STRUC.pickedUP = 0;
CreateCalcFolder(1);
if ORG_STRUC.pickUpYN
PickUp();
else
Initialize();
end
CreateCalcFolder();
end
safesave ('Current_POP.mat', POP_STRUC)
safesave ('Current_ORG.mat', ORG_STRUC)
if ((ORG_STRUC.dimension==3) && (ORG_STRUC.molecule==1)) || ...
((ORG_STRUC.dimension==1) && (ORG_STRUC.molecule==1) && (ORG_STRUC.varcomp==0))
[a,b]=unix (['cp MOL_* ' POP_STRUC.resFolder '/']);
elseif (ORG_STRUC.dimension==2) && (ORG_STRUC.molecule==0)
[a,b]=unix (['cp POSCAR_SUBSTRATE ' POP_STRUC.resFolder '/']);
end
path(path,[USPEXPath '/FunctionFolder/USPEX/' calcTypeStr]);
eval(['EA_' calcTypeStr '()']);
if ~isempty(ORG_STRUC.stopFitness)
disp('Run statistics analysis:');
statistics(ORG_STRUC.stopFitness);
end
if ORG_STRUC.repeatForStatistics > 1
RepeatRun();
else
if exist('still_running')
[a,b]=unix('rm still_running');
[a,b]=unix('rm POSCAR');
[a,b]=unix('rm POSCAR_order');
end
disp('USPEX IS DONE!');
[a,b]=unix('echo 1 > USPEX_IS_DONE');
disp(' ');
quit
end
end
function calcType_str = calcTypeStr()
global ORG_STRUC
calcType_str = [num2str(ORG_STRUC.dimension) '' num2str(ORG_STRUC.molecule) '' num2str(ORG_STRUC.varcomp)];
if str2num(calcType_str) < 0
calcType = -1*str2num(calcType_str);
calcType_str = ['M' num2str(calcType)];
end
