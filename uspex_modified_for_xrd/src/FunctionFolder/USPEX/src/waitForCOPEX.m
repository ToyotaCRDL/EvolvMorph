function waitForCOPEX()
global ORG_STRUC
if (ORG_STRUC.platform > 0) || (ORG_STRUC.numParallelCalcs > 1)   # modified 210630
#if (ORG_STRUC.platform > 0) | (ORG_STRUC.numParallelCalcs > 1)   # modified 210630
if ORG_STRUC.startNextGen==0
[nothing, nothing] = unix(['rm ', ORG_STRUC.homePath, '/still_running']);
fclose ('all');
quit();
end
else
while ~ORG_STRUC.startNextGen
load('Current_ORG.mat');
if ORG_STRUC.startNextGen==0
disp('Wait for COPEX')
pause(30);
end
end
end
