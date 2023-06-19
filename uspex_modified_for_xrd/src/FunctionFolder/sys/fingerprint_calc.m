function [order, fing, atom_fing]=fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons)
global ORG_STRUC
Rmax  = ORG_STRUC.RmaxFing;
sigm  = ORG_STRUC.sigmaFing;
delta = ORG_STRUC.deltaFing;
species = size(numIons,2);
Ncell = sum(numIons);
Nfull = size(dist_matrix, 1);  
normaliser = 1;
if ORG_STRUC.dimension == 0
V = 1;
normaliser = 0;
end
if sum(Ncell)==0
order=[];
fing =[];
atom_fing=[];
else
numBins = round(Rmax/delta);
fing = zeros(species*species, numBins);
if ORG_STRUC.dimension ~= 0
fing(:,1) = -1;
end
order = zeros(1, Ncell);
atom_fing = zeros(Ncell, species, numBins);
sigm = sigm/sqrt(2*log(2));
interval=zeros(2,1);
for bins = 2:numBins
for i = 1:Ncell
for j = 1:Nfull
if (dist_matrix(j,i)>0) && (abs(dist_matrix(j,i)-delta*(bins-0.5))<4*sigm)   # modified 210630
#if (dist_matrix(j,i)>0) & (abs(dist_matrix(j,i)-delta*(bins-0.5))<4*sigm)   # modified 210630
R0 = dist_matrix(j,i);
interval(2)=5*sign(delta*bins-R0);
if abs((delta*bins-R0)/sqrt(2)) <= 5*sigm
interval(2) = (delta*bins-R0)/(sqrt(2)*sigm);
end
interval(1)=5*sign(delta*(bins-1)-R0);
if abs((delta*(bins-1)-R0)/sqrt(2)) <= 5*sigm
interval(1) = (delta*(bins-1)-R0)/(sqrt(2)*sigm);
end
my_erf=approx_erf(interval, ORG_STRUC.erf_table);
delt=0.5*(my_erf(2)-my_erf(1));
atom_fing(i, typ_j(j), bins) = atom_fing(i, typ_j(j), bins) + delt/(Ni(j)*R0^2); 
fing((typ_i(i)-1)*species+typ_j(j), bins) = fing((typ_i(i)-1)*species+typ_j(j), bins) + delt/R0^2; 
end
end
atom_fing(i, :, bins) = V*atom_fing(i, :, bins)/(4*pi*delta) - normaliser;
for j = 1:species
weight = numIons(j)/sum(numIons); 
order(i) = order(i) + weight * delta*atom_fing(i, j, bins)*atom_fing(i, j, bins)/power(V/Ncell,1/3);  
end
end
for i = 1:species
for j = 1:species
fing((i-1)*species+j, bins) = V*fing((i-1)*species+j, bins)/(4*pi*numIons(i)*numIons(j)*delta) - normaliser;
end;
end;
end
order = sqrt(order);
end
