function LatticeMutation_300(Ind_No)
global POP_STRUC
global ORG_STRUC
global OFF_STRUC
goodMutant = 0;
goodMutLattice = 0;
count = 1;
while goodMutant + goodMutLattice ~= 2
count = count + 1;
if count > 50
#%disp('failed to do LatMutation in 50 attempts, switch to Random');
disp('failed to do LatMutation in 50 attempts, switch to Random');     # modified on 210414
USPEXmessage(512,'',0);
Random_300(Ind_No);
break;
end
toMutate = find (ORG_STRUC.tournament>RandInt(1,1,[0,max(ORG_STRUC.tournament)-1]));
ind = POP_STRUC.ranking(toMutate(end));
COORDINATES=POP_STRUC.POPULATION(ind).COORDINATES;
numIons=POP_STRUC.POPULATION(ind).numIons;
LATTICE=POP_STRUC.POPULATION(ind).LATTICE;
order =POP_STRUC.POPULATION(ind).order;
[MUT_LAT,strainMatrix, new_Coord] = lattice_atom_Mutation(COORDINATES, numIons, LATTICE, order);
coord = new_Coord*MUT_LAT;
[coord, MUT_LAT] = optLattice(coord, MUT_LAT);
new_Coord = coord/MUT_LAT;
goodMutant = distanceCheck(new_Coord,MUT_LAT,numIons, ORG_STRUC.minDistMatrice);
goodMutant = ( goodMutant & (CompositionCheck(numIons/ORG_STRUC.numIons)) );
goodMutLattice = latticeCheck(MUT_LAT);
if goodMutant + goodMutLattice == 2
OFF_STRUC.POPULATION(Ind_No).COORDINATES = new_Coord;
OFF_STRUC.POPULATION(Ind_No).LATTICE = MUT_LAT;
OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
OFF_STRUC.POPULATION(Ind_No).howCome = 'LatMutate';
info_parents = struct('parent', {},'strainMatrix', {}, 'enthalpy', {});
info_parents(1).parent= num2str(POP_STRUC.POPULATION(ind).Number);
info_parents.enthalpy = POP_STRUC.POPULATION(ind).Enthalpies(end)/sum(numIons);
info_parents.strainMatrix=strainMatrix;
info_parents.strainMatrix=0;
OFF_STRUC.POPULATION(Ind_No).Parents = info_parents;
disp(['Structure ' num2str(Ind_No) ' generated by latmutation']);
end
end