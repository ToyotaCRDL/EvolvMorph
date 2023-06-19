function fitness = AntiSeedsCorrection(fitness)
global ORG_STRUC
global POP_STRUC
global ANTISEEDS
fit_len = round(length(fitness)*ORG_STRUC.bestFrac);
antiSeedsMax = ORG_STRUC.antiSeedsMax*((sum(fitness(1:fit_len))/fit_len)-min(fitness)); 
if (POP_STRUC.generation > ORG_STRUC.antiSeedsActivation + 1) && ~isempty(ANTISEEDS)    
if antiSeedsMax < 0.5*ANTISEEDS(end).Max
antiSeedsMax = 0.5*ANTISEEDS(end).Max;
end
end
count = 0;
antiSeedsSigma = 0;
for fit1 = 1 : fit_len - 1
f1 = POP_STRUC.POPULATION(fit1).FINGERPRINT;
for fit2 = fit1+1 : fit_len
f2 = POP_STRUC.POPULATION(fit2).FINGERPRINT;
if isempty(f1) || isempty(f2)   # modified 210630
#if isempty(f1) | isempty(f2)   # modified 210630
cos_dist = 1;
else
if ORG_STRUC.varcomp == 1
cos_dist = cosineDistance(f1, f2, 1);
else
cos_dist = cosineDistance(f1, f2, ORG_STRUC.weight);
end
end
count = count + 1;
antiSeedsSigma = antiSeedsSigma + cos_dist;
end
end
if (fit_len < 2) || (antiSeedsSigma < 0.000000000000001)
antiSeedsSigma = 0.001;
count = 1;
end
antiSeedsSigma = ORG_STRUC.antiSeedsSigma*antiSeedsSigma/count;
if POP_STRUC.generation == 1
for antiS_loop = 1 : length(ANTISEEDS)
ANTISEEDS(antiS_loop).Sigma = antiSeedsSigma;
ANTISEEDS(antiS_loop).Max   = antiSeedsMax;
end
end
if ~isempty(ANTISEEDS)
if ~isempty(ANTISEEDS(1).FINGERPRINT) 
for antiS_loop = 1 : length(ANTISEEDS)
f1 = ANTISEEDS(antiS_loop).FINGERPRINT;
for fit_loop = 1 : length(POP_STRUC.POPULATION)
f2 = POP_STRUC.POPULATION(fit_loop).FINGERPRINT;
if isempty(f1) | isempty(f2)
cos_dist = 1;
else
if ORG_STRUC.varcomp == 1
cos_dist = cosineDistance(f1, f2, 1);
else
cos_dist = cosineDistance(f1, f2, ORG_STRUC.weight);
end
end
fitness(fit_loop) = fitness(fit_loop) + ANTISEEDS(antiS_loop).Max*exp(-cos_dist^2/(2*ANTISEEDS(antiS_loop).Sigma^2));
end
end
end
end
if ORG_STRUC.antiSeedsActivation > 0
if POP_STRUC.generation >= abs(ORG_STRUC.antiSeedsActivation)
for i = 1 : length(POP_STRUC.POPULATION)
f1 = POP_STRUC.POPULATION(i).FINGERPRINT;
ANTISEEDS(end+1).FINGERPRINT = f1;
ANTISEEDS(end).Sigma = antiSeedsSigma;
ANTISEEDS(end).Max   = antiSeedsMax;
end
end
elseif ORG_STRUC.antiSeedsActivation < 0
if POP_STRUC.generation >= abs(ORG_STRUC.antiSeedsActivation)
[tmp,ind] = min(fitness);
f1 = POP_STRUC.POPULATION(ind).FINGERPRINT;
ANTISEEDS(end+1).FINGERPRINT = f1;
ANTISEEDS(end).Sigma = antiSeedsSigma;
ANTISEEDS(end).Max   = antiSeedsMax;
%status = ['AntiSeed added at generation' num2str(POP_STRUC.generation)]
end
end
if ~isempty(ANTISEEDS)
if isempty(ANTISEEDS(1).FINGERPRINT) 
ANTISEEDS(1) = [];
end
end
safesave ('ANTISEEDS.mat', ANTISEEDS)
