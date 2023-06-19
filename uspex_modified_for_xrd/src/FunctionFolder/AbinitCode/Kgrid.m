function [Kpoints, Error] = Kgrid(LATTICE, Kresol, dimension)
angLattice = latConverter(LATTICE); 
vol        = abs(det(LATTICE));
dist    = zeros(1,3);
dist(3) = vol/(angLattice(1)*angLattice(2)*sin(angLattice(6)));
dist(2) = vol/(angLattice(1)*angLattice(3)*sin(angLattice(5)));
dist(1) = vol/(angLattice(2)*angLattice(3)*sin(angLattice(4)));
Kpoints = ceil(1./(dist*Kresol));
Error   = 0;
if abs(dimension) == 2   
Kpoints(1, 3) = 1;
end
if (max(Kpoints) > 20) && (vol > 50)    # modified 210630 
#if (max(Kpoints) > 20) & (vol > 50)    # modified 210630
Error   = 1;
Kpoints = [1 1 1];
end
