function Backupresult(fileName, bodyCount)

##### newly added file on 210414 #####
global ORG_STRUC
tempFolder = [ORG_STRUC.homePath 'CalcFoldTemp-allraw'];
global POP_STRUC

if ~exit(tempFolder,'dir')
[nothing] = unix(['mkdir -p ' tempFolder]);
end

[nothing,nothing] = unix(['cp ' fileName ' ' tempFolder '/' fileName '-' num2str(bodyCount+1)]);
[nothing,nothing] = unix('rm ' fileName]);


