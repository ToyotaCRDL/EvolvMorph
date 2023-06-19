function fitness = CalcFitness_300()
global POP_STRUC
global ORG_STRUC
fitness = zeros(1,length(POP_STRUC.POPULATION));
if ORG_STRUC.optType == 1 
for fit_loop = 1:length(POP_STRUC.POPULATION)
factor = POP_STRUC.POPULATION(fit_loop).numIons/ORG_STRUC.numIons;
fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end)/factor ;
end
elseif ORG_STRUC.optType == 2 
for fit_loop = 1:length(POP_STRUC.POPULATION)
fitness(fit_loop) = det(POP_STRUC.POPULATION(fit_loop).LATTICE)/sum(POP_STRUC.POPULATION(fit_loop).numIons);
end
elseif ORG_STRUC.optType == 3 
for i = 1 : length(POP_STRUC.POPULATION)
fitness(i) = -1*POP_STRUC.POPULATION(i).hardness;
end
elseif ORG_STRUC.optType == 4 
for fit_loop = 1:length(POP_STRUC.POPULATION)
if ORG_STRUC.opt_sign == -1   # added on 210414
  fitness(fit_loop) = -1*POP_STRUC.POPULATION(fit_loop).S_order;   # added on 210414
elseif ORG_STRUC.opt_sign == 1   # added on 210414
  fitness(fit_loop) = 1*POP_STRUC.POPULATION(fit_loop).S_order;   # added on 210414
#fitness(fit_loop) = -1*POP_STRUC.POPULATION(fit_loop).S_order;
end   # added on 210414
end
elseif ORG_STRUC.optType == 5 
for i = 1 : length(POP_STRUC.POPULATION)-1
for j = i+1 : length(POP_STRUC.POPULATION)
dist_ij = cosineDistance(POP_STRUC.POPULATION(i).FINGERPRINT, POP_STRUC.POPULATION(j).FINGERPRINT, ORG_STRUC.weight);
fitness(i) = fitness(i) + dist_ij^2;
fitness(j) = fitness(j) + dist_ij^2;
end
end
fitness = -sqrt(fitness);
elseif ORG_STRUC.optType == 6 
for i = 1 : length(POP_STRUC.POPULATION)
fitness(i) = -1*sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3;
end
elseif ORG_STRUC.optType == 7 
for i = 1 : length(POP_STRUC.POPULATION)
fitness(i) = -1*POP_STRUC.POPULATION(i).gap;
end
elseif ORG_STRUC.optType == 8 
Egc = 4; 
for i = 1 : length(POP_STRUC.POPULATION)
if POP_STRUC.POPULATION(i).gap >= Egc
fitness(i) = -1*(sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(POP_STRUC.POPULATION(i).gap/Egc)^2; 
else
fitness(i) = -1*(sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(POP_STRUC.POPULATION(i).gap/Egc)^6; 
end
end
elseif ORG_STRUC.optType == 9
for fit_loop = 1 : length(POP_STRUC.POPULATION)
fitness(fit_loop) = -1*POP_STRUC.POPULATION(fit_loop).mag_moment;
end
elseif ORG_STRUC.optType == 10 
for fit_loop = 1 : length(POP_STRUC.POPULATION)
factor = POP_STRUC.POPULATION(fit_loop).numIons/ORG_STRUC.numIons;
fitness(fit_loop) = -1*POP_STRUC.POPULATION(fit_loop).struc_entr/factor;
end
elseif ORG_STRUC.optType == 14 
for fit_loop = 1 : length(POP_STRUC.POPULATION)
fitness(fit_loop) = -1*POP_STRUC.POPULATION(fit_loop).powerfactor;
end

### newly added on 210414 for similarity function ###
### newly added on 210718 for cos similarity with vol-changed ###
elseif ORG_STRUC.optType == 104
for fit_loop = 1:length(POP_STRUC.POPULATION)
if ORG_STRUC.opt_sign == -1
  fitness(fit_loop) = -1*POP_STRUC.POPULATION(fit_loop).cos_similarity_xrd_volchanged;
elseif ORG_STRUC.opt_sign == 1
  fitness(fit_loop) = 1*POP_STRUC.POPULATION(fit_loop).cos_similarity_xrd_volchanged;
end
end

### newly added on 210802 for exp cos similarity with vol-changed ###
elseif ORG_STRUC.optType == 105
for fit_loop = 1:length(POP_STRUC.POPULATION)
if ORG_STRUC.opt_sign == -1
  fitness(fit_loop) = -1*POP_STRUC.POPULATION(fit_loop).exp_cos_similarity_xrd_volchanged;
elseif ORG_STRUC.opt_sign == 1
  fitness(fit_loop) = 1*POP_STRUC.POPULATION(fit_loop).exp_cos_similarity_xrd_volchanged;
end
end

######################################################

elseif (ORG_STRUC.optType > 1100) & (ORG_STRUC.optType < 1112)  
whichPara= mod(ORG_STRUC.optType,110);
for i = 1 : length(POP_STRUC.POPULATION)
if isempty(POP_STRUC.POPULATION(i).elasticProperties) | (POP_STRUC.POPULATION(i).elasticProperties(end)==0)
fitness(i)=NaN;
else
%disp(['Individual ', num2str(i) ]);
fitness(i) = -1*POP_STRUC.POPULATION(i).elasticProperties(whichPara);
end
end

end
#fitness = ORG_STRUC.opt_sign*fitness;      # modified 210414 (erased)

if ORG_STRUC.abinitioCode(1) > 0
for i = 1 : length(fitness)
if POP_STRUC.POPULATION(i).Enthalpies(end) > 99999    
fitness(i) = 100000;
end
end
end
