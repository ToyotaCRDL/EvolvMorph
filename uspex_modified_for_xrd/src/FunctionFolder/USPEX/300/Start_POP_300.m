function Start_POP_300()
global ORG_STRUC
global POP_STRUC
N_Step = length([ORG_STRUC.abinitioCode]);
POP_STRUC.DoneOrder = zeros(1, length(POP_STRUC.POPULATION));
for i = 1:length(POP_STRUC.POPULATION)
POP_STRUC.POPULATION(i).INIT_COORD = POP_STRUC.POPULATION(i).COORDINATES;
POP_STRUC.POPULATION(i).INIT_LAT   = POP_STRUC.POPULATION(i).LATTICE;
POP_STRUC.POPULATION(i).INIT_numIons=POP_STRUC.POPULATION(i).numIons;
if isempty(POP_STRUC.POPULATION(i).Step)
POP_STRUC.POPULATION(i).Step = 1;
POP_STRUC.POPULATION(i).Enthalpies = 100000*ones(1, N_Step);
POP_STRUC.POPULATION(i).K_POINTS = ones(N_Step, 3);
if ORG_STRUC.opt_sign == -1
POP_STRUC.POPULATION(i).dielectric_tensor =100000*ones(1,6);
POP_STRUC.POPULATION(i).gap = 100000;
POP_STRUC.POPULATION(i).mag_moment = 0;
POP_STRUC.POPULATION(i).powerfactor = 100000;
POP_STRUC.POPULATION(i).cos_similarity_xrd_volchanged = 100000;   # added on 210718
POP_STRUC.POPULATION(i).exp_cos_similarity_xrd_volchanged = 100000;   # added on 210802
else
POP_STRUC.POPULATION(i).dielectric_tensor =zeros(1,6);
POP_STRUC.POPULATION(i).gap = 0;
POP_STRUC.POPULATION(i).hardness = 0;
POP_STRUC.POPULATION(i).mag_moment = 100000;
POP_STRUC.POPULATION(i).powerfactor = 0;
POP_STRUC.POPULATION(i).cos_similarity_xrd_volchanged = 0;   # added on 210718
POP_STRUC.POPULATION(i).exp_cos_similarity_xrd_volchanged = 0;   # added on 210802
end
end
POP_STRUC.POPULATION(i).Error = 0;
POP_STRUC.POPULATION(i).Folder = 0;
POP_STRUC.POPULATION(i).ToDo = 1;
POP_STRUC.POPULATION(i).Done = 0;
POP_STRUC.POPULATION(i).Number = 0;
end
ReRank();
