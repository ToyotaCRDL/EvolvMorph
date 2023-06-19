function Backupvasprun(bodyCount)

##### newly added file on 210414 #####
global ORG_STRUC
tempFolder = [ORG_STRUC.homePath 'CalcFoldTemp-allraw'];
global POP_STRUC

if ~exit(tempFolder,'dir')
[nothing] = unix(['mkdir -p ' tempFolder]);
end

[nothing,nothing] = unix(['cp vasprun.xml ' tempFolder '/vasprun.xml-' num2str(bodyCount+1)]);

