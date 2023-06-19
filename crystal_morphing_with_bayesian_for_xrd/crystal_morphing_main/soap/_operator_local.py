import re
import copy
import numpy as np
import math
from pymatgen.core import Lattice, Structure, Molecule
from scipy.special import sph_harm as _sph_harm

from functools import wraps
import time
def stop_watch(func) :
    @wraps(func)
    def wrapper(*args, **kargs) :
        start = time.time()
        result = func(*args,**kargs)
        elapsed_time =  time.time() - start
        print(f"Time : {func.__name__} took {round(float(elapsed_time),3)} seconds.")
        return result
    return wrapper

##------------------ local functions
# @stop_watch
def _get_from_functions(self,functions, ris):
    values = [[0 for j in range(self.nlmax+1)] for i in range(self.nbasis)]
    for _n in range(self.nbasis):
        nn = _n + 1
        for ll in range(self.nlmax+1):
            key = str([nn,ll])
            values[_n][ll] = list(functions[key](ris))
    return values

#------------------- main funcitons
# @stop_watch
def _nlmchi_arr(self, chi, xyzis, xyz_center=0.0):
    '''
    Returns
    -------
    nlmchi_val_arr: complex, nparray(nbasis,nlmax+1,2*nlmax+1)
    '''
    element_list = list(dict.fromkeys(chi.r_elems))
    # _, center_indices  = chi.get_si( element = self.center_element_name)
    count_element_list = [elem for elem in element_list if elem != self.center_element_name]
    
    if chi.is_reciprocal() == True:
        ris = np.linalg.norm(xyzis,axis=1)
        # coef = (2.0/math.pi)**0.5 / (self.sigma_soap)**3.0
        from _soap import rfunction
        
        for count_element in count_element_list:
            sis, indices = chi.get_si( element = count_element)
            xyzis_part = np.array([xyzis[i] for i in indices])
            ris_part = np.array([ris[i] for i in indices])
            if self.use_fast_rfunctions:
                r2values = self._get_from_functions(self.r2values_functions, ris_part)
                # exit()
            else:
                r2values = rfunction.r2_arr(ris_part, self.nlmax, self.rcut,
                        self.sigma_soap, self.drad,
                        self.rbasis, self.rgrid)
            nlmchi_val_arr = rfunction.nlmchi_arr(xyzis_part, ris_part, 
                                                      sis, r2values,
                                                      self.sigma_soap)

    else: 
        print("Error; real space SOAP will be implemented in future")
        exit()
    return nlmchi_val_arr

# @stop_watch
def _dnlmchi_dt_arr(self, chi, xyzis, xyz_center=0.0):
    '''
    Returns
    ------
    dnlmchi_dt_arr: complex, nparray(nbasis,nlmax+1,3,3,2*nlmax+1)

    '''
    bb   = chi.lattice.matrix
    mfrac = chi.frac_coords
    element_list = list(dict.fromkeys(chi.r_elems))
    # _, center_indices  = chi.get_si( element = self.center_element_name)
    # count_element_list = [elem for elem in element_list if elem != self.center_element_name]
    ris = np.linalg.norm(xyzis,axis=1)

    from _soap import rfunction
    sis  = chi.get_si()
    gaunt_arr_flatten  = self.gaunt_arr.reshape(-1, order='F')
    if self.use_fast_rfunctions:
        r2values = self._get_from_functions(self.r2values_functions, ris)
        r3values = self._get_from_functions(self.r3values_functions, ris)
    else:
        r2values = rfunction.r2_arr(ris, self.nlmax, self.rcut,
                        self.sigma_soap, self.drad,
                        self.rbasis, self.rgrid)
        r3values = rfunction.r3_arr(ris, self.nlmax, self.rcut,
                        self.sigma_soap, self.drad,
                        self.rbasis, self.rgrid)
    dnlmchi_dt_mat_arr = rfunction.dnlmchi_dt_arr(bb, xyzis, ris, sis,
                                                r2values, r3values,
                                                mfrac, gaunt_arr_flatten,
                                                self.sigma_soap)
    return dnlmchi_dt_mat_arr

# @stop_watch
def _dnlmchi_ds_arr(self, chi, xyzis, xyz_center=0.0):
    '''
    Returns
    -------
    dnlmchi_ds_arr: complex, nparray(nbasis,nlmax+1,num_i,2*nlmax+1)
    
    '''
    element_list = list(dict.fromkeys(chi.r_elems))
    # si, center_indices  = chi.get_si( element = self.center_element_name)
    # count_element_list = [elem for elem in element_list if elem != self.center_element_name]
    ris = np.linalg.norm(xyzis,axis=1)
    from _soap import rfunction
    if self.use_fast_rfunctions:
        r2values = self._get_from_functions(self.r2values_functions, ris)
    else:
        r2values = rfunction.r2_arr(ris, self.nlmax, self.rcut,
                        self.sigma_soap, self.drad,
                        self.rbasis, self.rgrid)
    dnlmchi_ds_arr = rfunction.dnlmchi_ds_arr(xyzis, ris,
                                      r2values,
                                      self.sigma_soap)
    return dnlmchi_ds_arr

### for speed up
def save_gaunt(self, bool_arr):
    from sympy.physics.wigner import gaunt
    if bool_arr==True:
        self.gaunt_arr = np.zeros((self.nlmax+1, 2*self.nlmax+1, 2, 3))
        for ll in range(self.nlmax+1):
            llds = _lld_indices_trialngle_relation(ll)
            for mm in range(-ll, ll+1, 1):
                for i, lld in enumerate(llds):
                    self.gaunt_arr[ll, mm+self.nlmax, i, 0]  = float( gaunt(ll,lld,1,-mm, mm+1,-1, prec=64) )
                    self.gaunt_arr[ll, mm+self.nlmax, i, 1]  = float( gaunt(ll,lld,1,-mm, mm-1,+1, prec=64) )
                    self.gaunt_arr[ll, mm+self.nlmax, i, 2]  = float( gaunt(ll,lld,1,-mm, mm,   0, prec=64) )
    else:
        self.gaunt = {}
        for ll in range(self.nlmax+1):
            llds = _lld_indices_trialngle_relation(ll)
            for mm in range(-ll, ll+1, 1):
                for lld in llds:
                    self.gaunt[str([ll,lld,1,-mm, mm+1,-1])] = float( gaunt(ll,lld,1,-mm, mm+1,-1, prec=64) )
                    self.gaunt[str([ll,lld,1,-mm, mm-1,+1])] = float( gaunt(ll,lld,1,-mm, mm-1,+1, prec=64) )
                    self.gaunt[str([ll,lld,1,-mm, mm,   0])] = float( gaunt(ll,lld,1,-mm, mm,   0, prec=64) )
    return

def _jgaunt(self, ll, mm, theta, phi):
    from sympy.physics.wigner import gaunt
    jml = np.zeros(3, dtype=np.complex)
    ## x
    llds = _lld_indices_trialngle_relation(ll)
    for lld in llds:
        # term1 = self.sph_harm(lld, mm+1, theta, phi).conjugate() * gaunt(ll,lld,1,-mm, mm+1,-1, prec=64)
        # term2 = self.sph_harm(lld, mm-1, theta, phi).conjugate() * gaunt(ll,lld,1,-mm, mm-1,+1, prec=64)
        # term3 = self.sph_harm(lld, mm,   theta, phi).conjugate() * gaunt(ll,lld,1,-mm, mm,   0, prec=64)
        ### fast version
        term1 = self.sph_harm(lld, mm+1, theta, phi).conjugate() * self.gaunt[str([ll,lld,1,-mm, mm+1,-1])]
        term2 = self.sph_harm(lld, mm-1, theta, phi).conjugate() * self.gaunt[str([ll,lld,1,-mm, mm-1,+1])]
        term3 = self.sph_harm(lld, mm,   theta, phi).conjugate() * self.gaunt[str([ll,lld,1,-mm, mm,   0])]

        ## x
        jml[0] += term1 - term2
        ## y
        jml[1] += 1j*(term1 + term2)
        ## z
        jml[2] += term3*(2.0)**0.5
    return jml

def _jgaunt_real(self, ll,mm, theta, phi):
    from sympy.physics.wigner import gaunt
    jml = np.zeros(3)
    ## x
    llds = _lld_indices_trialngle_relation(ll)
    for lld in llds:
        ## x
        jml[0] += self.sph_harm_real(lld, -mm-1, theta, phi) * gaunt(ll,lld,1,mm,-mm-1,1, prec=64)
        ## y
        jml[1] += self.sph_harm_real(lld, -mm+1, theta, phi) * gaunt(ll,lld,1,mm,-mm+1,-1, prec=64)
        ## z
        jml[2] += self.sph_harm_real(lld, -mm, theta, phi)   * gaunt(ll,lld,1,mm,-mm,0, prec=64)

    return jml


def _lld_indices_trialngle_relation(ll):
    llds = []

    if ll == 0:
        llds = [1]
    else:
        llds = [ll-1, ll+1]
    return llds

# spherical harmonics
def sph_harm(self, l, m, theta, phi):
    if l < 0:
        return 0.0
    if abs(m) > l:
        return 0.0
    return _sph_harm(m,l,phi, theta )

def sph_harm_real(self, l, m, theta, phi):
    if l < 0:
        return 0.0
    if abs(m) > l:
        return 0.0
    if m > 0:
        return np.sqrt(2)*(-1)**m*_sph_harm(m,l, phi, theta).real
    if m == 0:
        return _sph_harm(m,l,phi,theta).real
    else:
        return np.sqrt(2)*(-1)**m*_sph_harm(-m,l,phi,theta).imag

##--------------------- functions other than soap ----------------------------------##
def _get_ref_gvecs(ref_mvecs, bb):
    from _soap import rfunction
    return rfunction.get_ref_gvecs(ref_mvecs, bb)
    ## python ver
    ref_gvecs = []
    for mm in ref_mvecs:
        gvec = mm@bb
        ref_gvecs.append(gvec)
    return ref_gvecs
    
def _get_ref_chigs(ref_gvecs, rjs, Ajs, sigma, vu):
    from _soap import rfunction
    return rfunction.get_ref_chigs(rjs, ref_gvecs, Ajs, sigma, vu)
    ## python ver
    ref_chigs = []
    for gvec in ref_gvecs:
        chig = rfunction.chig(rjs, gvec, Ajs, sigma, vu)
        ref_chigs.append(chig)
    return ref_chigs


def _index_in_rcut(ref_gvecs, rcut):
    from _soap import rfunction
    iis = rfunction.index_in_rcut(ref_gvecs, rcut) # if r > rcut : -1 is needed. 
    iis = iis[iis!=-1]
    # print("iis=",len(iis))
    return iis
    ## python ver
    iis = []
    for ii, gvec in enumerate(ref_gvecs):
        if math.sqrt(gvec[0]*gvec[0]+gvec[1]*gvec[1]+gvec[2]*gvec[2]) < rcut:
            iis.append(ii)
    return iis
##---------------------------------------------------------------------------------##

