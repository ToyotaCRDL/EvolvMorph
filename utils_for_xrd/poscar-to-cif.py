#!/usr/bin/env python
from pymatgen.io.vasp import Poscar
import sys
from pymatgen.io import cif

#### example : python poscar-to-cif.py POSCAR-MgO MgO-v2

posname = sys.argv[1]
struc = Poscar.from_file(posname,read_velocities=None).structure

w = cif.CifWriter(struc)
w.write_file(sys.argv[2]+'.cif')
