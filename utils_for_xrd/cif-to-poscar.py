#!/usr/bin/env python
from pymatgen.io.vasp import Poscar
import sys
from pymatgen.io import cif

# example : python cif-to-poscar.py MgO.cif POSCAR-MgO

cifname = sys.argv[1]
struc = cif.CifParser(cifname).get_structures(primitive=False)

pos = Poscar(struc[0])
pos.write_file(sys.argv[2])
