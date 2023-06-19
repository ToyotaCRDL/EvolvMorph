function [same] = sameComposition(c1,c2)
global ORG_STRUC
L = length(c1);
tolerance = 0.0001;
ratio = 0; 
for i = 1 : L
if (abs(c1(i)) < tolerance) && (abs(c2(i)) < tolerance)   # modified 210630
#if (abs(c1(i)) < tolerance) & (abs(c2(i)) < tolerance)   # modified 210630
continue;
elseif (abs(c1(i)) < tolerance) | (abs(c2(i)) < tolerance)   # modified 210630
#elseif (abs(c1(i)) < tolerance) | (abs(c2(i)) < tolerance)   # modified 210630
ratio = 0;
break;
elseif ratio == 0
ratio = c1(i)/c2(i);
elseif abs(c1(i)/c2(i) - ratio) > tolerance
ratio = 0;
break;
end
end
if ratio == 0
same = 0;
else
same = 1;
end
if (ORG_STRUC.dimension==0) && (ORG_STRUC.varcomp==1) && (same == 1)
for i = 1 : L
if (c1(i)~=c2(i))
same = 0;
end
end
end
