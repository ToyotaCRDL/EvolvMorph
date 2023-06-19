import os
import sys
import math
import numpy as np
import copy
import time
import itertools
import collections
from pymatgen.core import Structure, Lattice

homepath = os.path.join(os.path.dirname(__file__), '..')
# # sys.path.append(homepath)
# sys.path.append(os.path.join(homepath, "soap/_soap"))
# # print(homepath)
sys.path.append(os.path.join(homepath, "soap"))
# # from soap._soap import rfunction
# # print(sys.path)
from _operator_local import _get_ref_gvecs, _get_ref_chigs, _index_in_rcut

class CrystalStructure(Structure):
    """
    CrystalStructure object
    """

    def rtransform(self, sigma=0.5, multi_species=False):
        """
        Reciprocal-lattice transformation

        Paramters
        ---------
        sigma : float
            Width of gaussian allocated at atoms [A]
        multi_species : bool
            Switch to distinguish single element or multi element systems (default : False : single element)

        Returns
        -------
        _self : CrystalStructure object
            CrystalStructure object in reciprocal space
        """
        _self = copy.deepcopy(self)
        rjs = _self.cart_coords
        # rjs[:,:2] -= rjs[0,:2]
        rjs[:,:3] -= rjs[0,:3] # for 3d

        _self.rjs = rjs
        _self.multi_species=multi_species

        cs = _CrystalStructure(_self.lattice, _self.species, _self.frac_coords)
        _self.chis = cs.rtransform(sigma, multi_species)
        _self.rcut = _self.chis[0].rcut

        return _self


    def is_reciprocal(self):
        """
        True : object is in reciprocal space : true, vs in real space : false
        """
        if hasattr(self, 'chis'):
            varnames = np.array( list(self.chis[0].__dict__.keys()))
            if np.isin("si",varnames) == False: ## real space in case that there is no si
                return False
            else:
                return True
        else:
            varnames = np.array( list(self.__dict__.keys()))
            if np.isin("si",varnames) == False: ## real space in case that there is no si
                return False
            else:
                return True

    def nums_atoms_per_elements(self):
        """
        Returns
        -------
        natoms: (num elements) list
            Numbers of each elements
        elems: (num elements) list
            Symbol of each elements

        """
        if self.is_reciprocal()==False:
            print("Error: this method is only allowed for the real space structure.")
            exit()

        elems = []
        for site in self.sites:
            elems.append( site.specie.symbol )
        elems = list(dict.fromkeys(elems)) ## eliminate duplication
        natoms = []
        for elem in elems:
            count = 0
            for site in self.sites:
                if elem == site.specie.symbol:
                    count += 1
            natoms.append(count)

        return natoms, elems


    def write_POSCAR(self, poscar_name):
        """
        Parameters
        ----------
        poscar_name: string
            file name printed out

        """
        ## return here in case that there is no si for instance variable
        if self.is_reciprocal() == False:
            self.to(filename=poscar_name)
            return
        ## if there is si, becaue it is in reciprocal space, for convenience, add abs_chigs to the vx column of initial velocity in POSCAR format.
        ## ref: https://www.vasp.at/wiki/index.php/CONTCAR
        self.chis[0].to(filename=poscar_name)
        with open(poscar_name, 'a') as f:
            print('' , file = f)
            for abs_chig in self.chis[0].si:
                print(str(abs_chig) + '  0.0  0.0' , file = f)


    def write_CIF(self, cif_name, scaling=False, scale_to_chi=None, scaling_method="volume"):
        """
        Parameters
        ----------
        cif_name : str
            cif file name printed out
        scaling : bool
           True : scaling is applied
        scale_to_chi : CrystalStructure object
            target structure when scaling is applied
        scaling_method :
            "area": scaling based on the area of 2D lattice
            "volume": scaling based on the volume of lattice
        """
        print("\n                       Output CIF file")
        if self.is_reciprocal() == True:
            if hasattr(self, 'chis'):
                lattice = self.chis[0].lattice.reciprocal_lattice
                _chi_real = CrystalStructure(lattice, self.species, self.rjs, coords_are_cartesian=True)
            else:
                lattice = self.lattice.reciprocal_lattice
                _chi_real = CrystalStructure(lattice, self.species, self.frac_coords)
        else:
            _chi_real = copy.deepcopy(self)

        if scaling:
            _chi_real = _chi_real.scale_chi(scale_to_chi, scaling_method)

        # Output file
        _chi_real.to(filename=cif_name)


    def scale_density(self, const_natoms_density = 0.1):
        """
        normalize the lattice constant that density of number of atoms becomes const_natoms_density

        Parameters
        ----------
        const_natoms_density: flaot

        Returns
        -------
        _self: object
            normalized CrystalStructure object

        """
        _self = copy.deepcopy(self)
        newvolume = (_self.natoms_density() / const_natoms_density) * _self.lattice.volume
        # print( _self.natoms_density(), const_natoms_density,  _self.lattice.volume)
        _self.scale_lattice(newvolume)
        return _self


    def natoms_density(self):
        """
        Returns
        -------
        natoms_density: float
            density of number of atoms [1/A3^]

        """
        return self.num_sites / self.volume  # number / unitcell volume A^3


    def scale_chi(self, chi, scaling_method="volume"):
        """
        To avoid too much deformation 

        Parameters
        -----------
        chi : CrystalStructure object
            structure with the desired volume
        scaling_method : str
            "area": scaling based on the area of 2D lattice
            "volume": scaling based on the volume
        """
        _self = copy.deepcopy(self)
        lat = _self.lattice.matrix
        abc = tuple(np.sqrt(np.sum(lat ** 2, axis=1)).tolist())
        versors = lat / abc
        if scaling_method=="area":
            new_area = np.linalg.norm(np.cross(chi.lattice.matrix[0], chi.lattice.matrix[1]))
            _lat = copy.deepcopy(lat)
            ab = abc[:2]
            versors = versors[:2,:2]
            geo_factor = np.linalg.norm(np.cross(versors[0], versors[1]))
            ratios = np.array(ab) / ab[1]
            new_b = (new_area / (geo_factor * np.prod(ratios))) ** (1 / 2.0)
            _lat[:2,:2] = versors * (new_b * ratios)
            # print("\n area =", new_area, ",", np.linalg.norm(np.cross(_lat[0],_lat[1])), "\n")
        else:
            new_volume = chi.volume
            geo_factor = abs(np.dot(np.cross(versors[0], versors[1]), versors[2]))
            ratios = np.array(abc) / abc[2]
            new_c = (new_volume / (geo_factor * np.prod(ratios))) ** (1 / 3.0)
            _lat = versors * (new_c * ratios)
            # print("\n volume =", new_volume, ",", abs(np.dot(np.cross(_lat[0], _lat[1]), _lat[2])), "\n")

        _self.lattice = Lattice(_lat)
        return _self


class _CrystalStructure(Structure):
    """
    _CrystalStructure object
    """    
    def rtransform(self, sigma=0.5, multi_species=False):
        """
        Reciprocal-lattice transformation

        Paramters
        ---------
        sigma : float
            Width of gaussian allocated at atoms [A]
        multi_species : bool
            Switch to distinguish single element or multi element systems (default : False : single element)

        Returns
        -------
        r_struct : list of _CrystalStructure object
            list of CrystalStructure object in reciprocal space as element
        """
        
        cc = 6.0
        const_natoms_density = 0.1

        self.sigma = sigma
        self.rcut = (2.0*cc)**0.5 / sigma

        ## First: single species treatment
        if multi_species==False:
            print("Single species treatment")

            _struct = copy.deepcopy(self)
            lattice = _struct.lattice
            r_lattice = lattice.reciprocal_lattice

            num_m=20
            mrange = np.arange(-(num_m-1), num_m)  # if num_m=5, mrange = [-4 -3 -2 -1  0  1  2  3  4]
            ref_mvecs = []
            for m3 in mrange:
                for m2 in mrange:
                    for m1 in mrange:
                        ref_mvecs.append(np.array([m1, m2, m3]))
            ref_gvecs = _get_ref_gvecs(ref_mvecs, r_lattice.matrix)
            A_sign    = []
            for j, rj in enumerate(_struct.cart_coords):
                A_sign.append(1.0)
            ref_chigs = _get_ref_chigs(ref_gvecs, _struct.cart_coords, A_sign, self.sigma, _struct.volume)
            
            ##---- search gvecs < rcut
            mvecs = []
            gvecs = []
            chigs     = []
            abs_chigs = []
            iis = _index_in_rcut(ref_gvecs, self.rcut)
            for ii in iis:
                mvecs.append(ref_mvecs[ii])
                gvecs.append(ref_gvecs[ii])
                chigs.append(ref_chigs[ii])
                abs_chigs.append(abs(ref_chigs[ii]))
                                        
            if max( [max(xx) for xx in mvecs] ) > num_m:
                print("   Error: search gvecs within rcut failed in transform.")
                print("    Found max m is", max( [max(xx) for xx in mvecs] ) )
                print("    num_m defined as max m is ", num_m)
                exit()
            ##--- G indicates G=0 point, X's are other G points
            r_elems = []
            for gvec in gvecs:
                if np.linalg.norm(gvec) == 0.0:
                    r_elems.append("G")
                else:
                    r_elems.append("X")

            ##--- reciprocal structure object
            r_struct = _CrystalStructure(r_lattice, r_elems, mvecs)

            r_struct.r_elems = r_elems
            r_struct.rcut    = self.rcut
            r_struct.gvecs   = gvecs
            r_struct.si      = abs_chigs
            r_struct.chigs   = chigs
            r_struct.sigma   = self.sigma
            r_struct.ref_mvecs = ref_mvecs
            r_struct.ref_chigs = ref_chigs

            ## ---check
            # for abs_chigs, r_elem, mvec in zip(abs_chigs, r_elems, mvecs):
            #     print(abs_chigs, r_elem, mvec)

            return [r_struct]

        ## Other cases: take "multiple species" treatment
        else:
            print("Multiple species treatment")

            # get species on each site
            arr = [self.species[i].name for i in range(len(self.species))]
            # print(arr)

            # get combinations of species
            count_dict = collections.Counter(arr)
            if len(count_dict.keys())<2:
                print("\n Error: multiple species are not found, species =", arr)
                print("Please set species = \"single\"")
                exit()
            c_list = list(itertools.combinations(count_dict.keys(), 2))

            # for loop to make r_struct for each pair of species
            ## 
            r_struct_list=[]
            for pair in c_list:
                self.sigma = sigma

                # make _struct of species in pair
                _struct = copy.deepcopy(self)
                rm_sp_names = [s for s in arr if s not in pair]
                _struct.remove_species(rm_sp_names)

                lattice = _struct.lattice
                r_lattice = lattice.reciprocal_lattice

                ##---- search gvecs < rcut
                num_m=20
                mrange = np.arange(-(num_m-1), num_m)  # if num_m=5, mrange = [-4 -3 -2 -1  0  1  2  3  4]
                ref_mvecs = []
                for m3 in mrange:
                    for m2 in mrange:
                        for m1 in mrange:
                            ref_mvecs.append(np.array([m1, m2, m3]))
                ref_gvecs = _get_ref_gvecs(ref_mvecs, r_lattice.matrix)
                A_sign    = []
                for j, rj in enumerate(_struct.cart_coords):
                    if _struct.species[j].name==pair[0]:
                        A_sign.append(-1.0)
                        continue
                    A_sign.append(1.0)

                ref_chigs = _get_ref_chigs(ref_gvecs, _struct.cart_coords, A_sign, self.sigma, _struct.volume)
            
                ##---- search gvecs < rcut
                mvecs = []
                gvecs = []
                chigs     = []
                abs_chigs = []
                iis = _index_in_rcut(ref_gvecs, self.rcut)
                for ii in iis:
                    mvecs.append(ref_mvecs[ii])
                    gvecs.append(ref_gvecs[ii])
                    chigs.append(ref_chigs[ii])
                    abs_chigs.append(abs(ref_chigs[ii]))
                    
                    
                if max( [max(xx) for xx in mvecs] ) > num_m:
                    print("   Error: search gvecs within rcut failed in transform.")
                    print("    Found max m is", max( [max(xx) for xx in mvecs] ) )
                    print("    num_m defined as max m is ", num_m)
                    exit()

                ##--- G indicates G=0 point, X's are other G points
                r_elems = []
                for gvec in gvecs:
                    if np.linalg.norm(gvec) == 0.0:
                        r_elems.append("G")
                    else:
                        r_elems.append("X")

                ##--- reciprocal structure object
                r_struct = _CrystalStructure(r_lattice, r_elems, mvecs)

                r_struct.r_elems = r_elems
                r_struct.rcut    = self.rcut
                r_struct.gvecs   = gvecs
                r_struct.si      = abs_chigs
                r_struct.chigs   = chigs
                r_struct.sigma   = self.sigma
                r_struct.ref_mvecs = ref_mvecs
                r_struct.ref_chigs = ref_chigs
                r_struct.sp_names= arr
                r_struct.pair    = pair

                r_struct_list.append(r_struct)

            return r_struct_list
       

    def is_reciprocal(self):
        """
        True : object is in reciprocal space : true, vs in real space : false
        """
        varnames = np.array( list(self.__dict__.keys()))
        if np.isin("si",varnames) == False: ## real space in case that there is no "si"
            return False
        else:
            return True
        
    def have_power_spectrum(self):
        """
        True : object includes Power spectrum
        Cache in order to speed up calculation with get_1d_power_spectrum of soap.
        """
        varnames = np.array( list(self.__dict__.keys()))
        if np.isin("power_spectrum",varnames) == False: ## in case that no power_spectrum
            return False
        else:
            return True
        
    def have_dpower_spectrum_dt(self):
        varnames = np.array( list(self.__dict__.keys()))
        if np.isin("dpower_spectrum_dt",varnames) == False: 
            return False
        else:
            return True

    def have_dpower_spectrum_ds(self):
        varnames = np.array( list(self.__dict__.keys()))
        if np.isin("dpower_spectrum_ds",varnames) == False: 
            return False
        else:
            return True


    def nums_atoms_per_elements(self):
        """
        Returns
        -------
        natoms: (num elements) list
            Numbers of each elements
        elems: (num elements) list
            Symbol of each elements

        """
        if self.is_reciprocal()==False:
            print("Error: this method is only allowed for the real space structure.")
            exit()

        elems = []
        for site in self.sites:
            elems.append( site.specie.symbol )
        elems = list(dict.fromkeys(elems)) 
        natoms = []
        for elem in elems:
            count = 0
            for site in self.sites:
                if elem == site.specie.symbol:
                    count += 1
            natoms.append(count)

        return natoms, elems


    def get_si(self, element=""):
        """
        Parameters
        ----------
        element: string
           get si of selected element. empty case : get all
        Returns
        -------
        si: (num of spots) ndarray
            spot strength norm

        """
        if self.is_reciprocal() == False:
            si = np.array( [ 1. for _ in range(self.num_sites)] )
            return si
        else: ## in reciprocal space
            if element != "":
                indices = [i for i, x in enumerate(self.r_elems) if x == element]
                si = [self.si[i] for i, x in enumerate(self.r_elems) if x == element]
                return np.array(si), indices
            else:
                return np.array(self.si)

    def set_si(self, si):
        if self.is_reciprocal() == True:
            self.si = si
        else:
            print(" Error: set_si is not allowed if chi is real space.")
        return

    def get_chigs(self, element=""):
        """
        Parameters
        ----------
        element: string
            get chigs of selected element. empty case : get all

        Returns
        -------
        chigs: (num of spots) ndarray
            complex number spot strength 

        """
        if self.is_reciprocal() == False:
            chigs = np.array( [ 1. for _ in range(self.num_sites)] )
            return chigs
        else: ## in reciprocal space
            if element != "":
                indices = [i for i, x in enumerate(self.r_elems) if x == element]
                chigs = [self.chigs[i] for i, x in enumerate(self.r_elems) if x == element]
                return np.array(chigs), indices
            else:
                return np.array(self.chigs)

    def set_chigs(self, chigs):
        if self.is_reciprocal() == True:
            self.chigs = chigs
            self.set_si( [ abs(x) for x in chigs] )
        else:
            print(" Error: set_chigs is not allowed if chi is real space.")
        return



##-----------------------------------------------------------------------------##
if __name__ == "__main__":
    struct_file1 = './examination/CONTCAR_Ti_hcp'
    # struct_file1 = './examination/CONTCAR_Mo_bcc'
    # struct_file1 = './examination/CONTCAR_Cu_fcc'

    cstruct=CrystalStructure.from_file(struct_file1).scale_density()

    #cstruct2 = cstruct.scale_density()

    # Reciprocal transformation parameters
    sigma = 0.5
    rstruct1 = cstruct.rtransform(sigma=sigma)

    poscar1 = "./examination/POSCAR_r1-dev"
    rstruct1.write_POSCAR(poscar1)

    print(rstruct1.is_reciprocal())
    print(rstruct1.chis[0].get_si())
    # print(rstruct1)

    print("Normal ends")
