
import sys
import os
import copy
import numpy as np
import math
import cmath
import itertools
from pymatgen.core import Structure
from pymatgen.transformations import standard_transformations as stf
from scipy import interpolate
##
homepath = os.path.join(os.path.dirname(__file__), '..')
sys.path.append(os.path.join(homepath, "crystal_structure"))
from crystal_structure import CrystalStructure
##
# from tqdm import tqdm

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

class Soap():
    """
    Integration of SOAP, distance, power spectrum, and so on.
    """
    # from _operator_local import sph_harm_real, sph_harm
    # from _operator_local import _nlmchi, _dnlmchi_dt, _dnlmchi_ds
    from _operator_local import _nlmchi_arr, _dnlmchi_dt_arr, _dnlmchi_ds_arr
    from _operator_local import _get_from_functions
    from _operator_local import _jgaunt_real, _jgaunt, save_gaunt

    def __init__(self, nbasis, nlmax, rcut, sigma_soap, nrad, ltype, center_element_name, 
                 use_fast_rfunctions=True, cache_power_spectrum=True):
        """
        Parameters
        ----------
        nbasis : int
            Number of radial basis
        nlmax : int
            Max number of spherical harmonics indiex l
        rcut : float
            Cutoff radius of counting region
        sigma_soap : float
            Gaussian width of the atom model
        nrad : int
            Number of integral points of r2 and r3 functions
        ltype : int
            Kernel type
            1=Original SOAP, 2 = Alchemical SOAP**, 3 = unifrom element***
            *  Calculate integrals with the same elements
            ** For different type of atoms : Guassian with electron affinity
               Ref: S. De, A. P. Bartó, G. Csáyi, and M. Ceriotti, Phys. Chem.Chem. Phys. 18, 13754 (2016).
            *** to calculate a different element by considering it as another
        center_element_name : str
            "element symbol" or "ALL"
            In case for element symbol, SOAP is calculated with the center of the specific atom.
            In case for ALL, SOAP is calculated with the center of the all atoms.
        use_fast_rfunctions: bool
            If True, use R2 and R3 interpolated functions.
        cache_power_spectrum: bool
            If True, the power spectrum of the chi is cached, for speed up.

        """
        self.nbasis = nbasis
        self.nlmax  = nlmax
        self.sigma_soap  = sigma_soap
        self.nrad   = nrad
        self.rcut   = rcut
        self.ltype  = ltype
        self.center_element_name = center_element_name
        self.use_fast_rfunctions = use_fast_rfunctions
        self.cache_power_spectrum = cache_power_spectrum
        
        self.eps   = 1e-6

        ## Calculation of the gaunt coeffcients is bottleneck.
        ## Save them as instance variable to cut the computational cost.
        bool_arr=True
        self.save_gaunt(bool_arr)
        ## Radial basis sets of SOAP is bottleneck.
        ## Save them as instance variable to cut the computational cost.
        from _soap import rfunction
        self.drad, self.rbasis, self.rgrid =  rfunction.radial_basis_set(self.nbasis, self.nrad, 
                                                                         self.rcut, self.sigma_soap)

        ### fast R2, R3 prepartion
        ## Number of data samplings
        if self.use_fast_rfunctions:
            num_i = 100
            self.r2values_functions, self.r3values_functions = self._set_fast_rfunctions(num_i)        

    def _set_fast_rfunctions(self, num_i):
        '''
        Returns
        -------
        r2values_functions, r3values_functions: dictionary of functions
            r2values[key] where string of n and l; for example of key, "[1, 0]".
            Input of the function is list of r coodinates; 
                For example, r2values_functions['[1, 0]']([0.0, 0.25, 0.3, 2.0]). 
                This returns a list of the r2values.
        '''
        ## rcut*1.1 is used because this extention is considered to be good for the interpolation accuracy.
        dr = (self.rcut*1.1)/num_i
        sample_ris = [ float(elem)*dr for elem in range(num_i+1)]
        ## interpolate
        from _soap import rfunction
        r2values_functions = {}
        r3values_functions = {}
        for _n in range(self.nbasis):
            nn = _n + 1
            for ll in range(self.nlmax+1):
                key = str([nn,ll])
                r2values=[]
                r3values=[]
                for ri in sample_ris:
                    r2value = rfunction.r2(nn,ll, self.nlmax, self.rcut, self.sigma_soap, ri, 
                        self.drad, self.rbasis, self.rgrid)
                    r3value = rfunction.r3(nn,ll, self.nlmax, self.rcut, self.sigma_soap, ri, 
                        self.drad, self.rbasis, self.rgrid)
                    r2values.append(r2value)
                    r3values.append(r3value)
                # Quadratic spline interpolation
                ## ref: https://qiita.com/maskot1977/items/913ef108ff1e2ba5b63f
                ip_method = lambda x, y: interpolate.interp1d(x, y, kind="quadratic") # "quadratic" "cubic"
                r2values_functions[key] = ip_method(sample_ris, r2values)
                r3values_functions[key] = ip_method(sample_ris, r3values)
        return r2values_functions, r3values_functions

    # @stop_watch
    def power_spectrum(self,chi):
        """
        Parameters
        ----------
        chi: CrystalStructure obj

        Returns
        -------
        pnnl: (nbasis, nbasis, nlmax+1) numpy object
            Power spectrum = P_{n, n', l} (chi)

        """
        xyzis = chi.cart_coords
        nlmchi = {}
        _nlmchi = self._nlmchi_arr(chi, xyzis)

        for _nn in range(self.nbasis):
            nn = _nn + 1
            for ll in range(self.nlmax+1):
                key = str([nn,ll])
                nlmchi[key] = np.array(_nlmchi[_nn,ll,:])

        pnnl = np.zeros((self.nbasis, self.nbasis, self.nlmax+1), dtype="complex")
        for _nn in range(self.nbasis):
            nn = _nn + 1
            for _nnd in range(self.nbasis):
                nnd = _nnd + 1
                for ll in range(self.nlmax+1):
                    key1 = str([nn,ll])
                    key2 = str([nnd,ll])
                    res = np.dot(nlmchi[key1].conjugate(), nlmchi[key2])
                    pnnl[_nn, _nnd, ll] = 2.0*math.pi*(2./(2*ll+1))**0.5*res

        return pnnl

    # @stop_watch
    def get_1d_power_spectrum(self, chi):
        '''
        If chi has power_spectrum, the saved data is used.
        Otherwise, it is calculated and cached on the chi.
        
        NOTE: This procedure is developed to speed up. But this coding is a bit dangerous to induce a future but.
              Remind that the chi update should be performed only by CrystalStructure or _CrystalStructure (from scratch).
              If you update chi manually, this causes a bug.
        '''
        if self.cache_power_spectrum == True :
            if chi.have_power_spectrum() == False:
                ps1d = self.power_spectrum(chi).flatten() 
                chi.power_spectrum = ps1d
            else:
                ps1d = chi.power_spectrum
        elif self.cache_power_spectrum == False:
            ps1d = self.power_spectrum(chi).flatten() 
        return ps1d
    
    # @stop_watch
    def overlapInteg(self, chi1, chi2):
        """
        Parameters
        ----------
        chi1, chi2: CrystalStructure obj

        Returns
        -------
        valueL: float
            Overlap integral L =  int |<chi1|R|chi2>|^2 dR
        """
        ps1d_1 = self.get_1d_power_spectrum(chi1)
        ps1d_2 = self.get_1d_power_spectrum(chi2)
        valueL = np.dot(ps1d_1.conjugate(), ps1d_2).real
        return valueL

    def similarity(self, chi1, chi2):
        '''
        Returns
        -------
        float
            Similarity tilde L

        '''
        return (self.overlapInteg(chi1, chi2)/np.sqrt(self.overlapInteg(chi1, chi1)*self.overlapInteg(chi2, chi2)))

    def distance(self, chi1, chi2):
        """
        Returns
        -------
        float
            SOAP distance
        """
        similarity = self.similarity(chi1, chi2)
        if similarity > 1.0: ## Sometimes it's just a little bit of a numerical error
            similarity = 1.0
        return np.sqrt(2.0 - 2.0*similarity)

    ##--------------------------------------------------------------------------##
    # @stop_watch
    def dpower_spectrum_dt(self, chi):
        """
        Derivative of P with reciprocal basis

        Returns
        -------
        dpp_dt: (3; xi, 3; eta, nbasis; n, nbasis; n', nlmax+1; l) numpy object
            dP_{n,n',l}(delta chi(t_{xi, eta}))/dt_{xi, eta}
        """
        import time
        start = time.time()

        xyzis = chi.cart_coords

        if chi.is_reciprocal() == False:
            print("Error: derivative routines are not allowed for real space chi.")
            exit()
        nlmchi = {}
        dnlmchi_dt = {}
        _nlmchi         = self._nlmchi_arr(chi, xyzis)
        _dnlmchi_dt_arr     = self._dnlmchi_dt_arr(chi, xyzis)

        for _nn in range(self.nbasis):
            nn = _nn + 1
            for ll in range(self.nlmax+1):
                key_nl          = str([nn,ll])
                nlmchi[key_nl]  = _nlmchi[_nn,ll,:]
                dnlmchi_dt[key_nl] = _dnlmchi_dt_arr[_nn,ll]

        # elapsed_time = time.time() - start
        # print ("\n                dt-1, elapsed_time:{0}".format(elapsed_time) + "[sec]")

        start = time.time()
        dpp_dt = np.zeros((3, 3, self.nbasis, self.nbasis, self.nlmax+1), dtype="complex")
        for _nn in range(self.nbasis):
            nn = _nn + 1
            for _nnd in range(self.nbasis):
                nnd = _nnd + 1
                for ll in range(self.nlmax+1):
                    key_nl = str([nn,ll])
                    key_ndl = str([nnd,ll])
                    coef = 2.0*math.pi*(2./(2*ll+1))**0.5
                    for xi in range(3):
                        for eta in range(3):
                            _ndlmchi     = nlmchi[ key_ndl ]
                            _dnlmchi_dt  = dnlmchi_dt[ key_nl ][xi, eta]
                            res1         = np.dot(_ndlmchi.conjugate(), _dnlmchi_dt)
                            _dndlmchi_dt = dnlmchi_dt[ key_ndl ][xi, eta]
                            _nlmchi      = nlmchi[ key_nl ]
                            res2         = np.dot(_dndlmchi_dt.conjugate(), _nlmchi)
                            dpp_dt[xi, eta, _nn, _nnd, ll] =  coef*(res1+res2)
        # elapsed_time = time.time() - start
        # print ("\n                dt-2, elapsed_time:{0}".format(elapsed_time) + "[sec]")
        return dpp_dt

    ##--------------------------------------------------------------------------##
    # @stop_watch
    def dpower_spectrum_ds(self,chi):
        """
        Derivative of P with si

        Returns
        -------
        dpp_ds: (ii, nbasis; n, nbasis; n', nlmax+1; l) numpy object
            dP_{n,n',l}(delta chi(s_{ii}))/ds_{ii}
        """
        import time
        start = time.time()

        xyzis = chi.cart_coords

        if chi.is_reciprocal() == False:
            print("Error: derivative routines are not allowed for real space chi.")
            exit()
        nlmchi = {}
        dnlmchi_ds = {}
        num_ii = 0
        _nlmchi = self._nlmchi_arr(chi, xyzis)
        _dnlmchi_ds_arr = self._dnlmchi_ds_arr(chi, xyzis)
        num_ii = len(_dnlmchi_ds_arr[0,0,:,0])
        for _nn in range(self.nbasis):
            nn = _nn + 1
            for ll in range(self.nlmax+1):
                key_nl     = str([nn,ll])
                nlmchi[key_nl]     = _nlmchi[_nn,ll,:]
                dnlmchi_ds[key_nl] = _dnlmchi_ds_arr[_nn,ll]

        # elapsed_time = time.time() - start
        # print ("\n                ds-1, elapsed_time:{0}".format(elapsed_time) + "[sec]")
        
        start = time.time()
        from _soap import rfunction
        dpp_ds = np.zeros((num_ii, self.nbasis, self.nbasis, self.nlmax+1), dtype="complex")
        for _nn in range(self.nbasis):
            nn = _nn + 1
            for _nnd in range(self.nbasis):
                nnd = _nnd + 1
                for ll in range(self.nlmax+1):
                    key_nl = str([nn,ll])
                    key_ndl = str([nnd,ll])
                    _ndlmchi     = nlmchi[ key_ndl ]
                    _dnlmchi_ds  = dnlmchi_ds[ key_nl ]
                    _dndlmchi_ds = dnlmchi_ds[ key_ndl ]
                    _nlmchi      = nlmchi[ key_nl ]
                    res12 = rfunction.dpp_ds_res(_nlmchi, _ndlmchi, _dnlmchi_ds, _dndlmchi_ds)
                    dpp_ds[:, _nn, _nnd, ll] =  2.0*math.pi*(2./(2*ll+1))**0.5*(res12[:])

        # elapsed_time = time.time() - start
        # print ("\n                ds-2, elapsed_time:{0}".format(elapsed_time) + "[sec]")

        return dpp_ds


    ##--------------------------------------------------------------------------##
    # @stop_watch
    def doverlapInteg(self, chi1, chi2, variation, ps1d_1=None, dpp=None):
        """
        Returns
        -------
        if variation == "dt":
            dL_dt: (3; xi, 3; eta) numpy object
                dL(chi1, chi2)/d chi2(t_{xi, eta})
        if variation == "ds"
            dL_ds: (num_ii; num of points) numpy object
                dL(chi1, chi2)/d chi2(s_{ii})
        """
        if np.all(ps1d_1==None):
            ps1d_1 = self.get_1d_power_spectrum(chi1)

        if variation == "dt":
            if np.all(dpp==None):
                dpp_dt = self.dpower_spectrum_dt(chi2)
            else:
                dpp_dt = dpp

            dL_dt  = np.zeros((3, 3))
            for xi in range(3):
                for eta in range(3):
                    ps1d_2 = dpp_dt[xi, eta].flatten()
                    dL_dt[xi, eta] = np.dot(ps1d_1.conjugate(), ps1d_2).real

            return dL_dt

        elif variation == "ds":
            if np.all(dpp==None):
                dpp_ds = self.dpower_spectrum_ds(chi2)
            else:
                dpp_ds = dpp

            num_ii = len(dpp_ds)
            dL_ds  = np.zeros(num_ii)
            for ii in range(num_ii):
                ps1d_2 = dpp_ds[ii].flatten()
                dL_ds[ii] = np.dot(ps1d_1.conjugate(), ps1d_2).real

            return dL_ds

        else:
            print("Error: variation should be \"dt\" or \"ds\". Currently using ", variation)
            exit()

    # @stop_watch
    def dsimilarity(self, chi1, chi2, variation):
        """
        Returns
        -------
        if variation == "dt":
            dtildeL_dt: (3; xi, 3; eta) numpy object
                d tilde L(chi1, chi2)/d chi2(t_{xi, eta})
        if variation == "ds"
            dtildeL_ds: (num_ii; num of poins) numpy object
                d tilde L(chi1, chi2)/d chi2(s_{ii})
        """
        l11 = self.overlapInteg(chi1, chi1)
        l22 = self.overlapInteg(chi2, chi2)
        l12 = self.overlapInteg(chi1, chi2)
        # print(l11, l22, l12)
        # exit()
        denom = 1.0/(l11*l22 )**0.5
        coef  = -1.0 *l12/ l22
        dL12 = self.doverlapInteg(chi1, chi2, variation)
        dL22 = self.doverlapInteg(chi2, chi2, variation)
        if variation == "dt":
            dtildeL_dt = np.zeros((3, 3))
            for xi in range(3):
                for eta in range(3):
                    dtildeL_dt[xi,eta] = denom*(dL12[xi, eta] + coef*dL22[xi,eta] )
            return dtildeL_dt
        elif variation == "ds":
            num_ii = len(dL22)
            dtildeL_ds = np.zeros(num_ii)
            for ii in range(num_ii):
                # dtildeL_ds[ii] = 2./l11**2.0*dL12[ii]*(l12-l11)
                dtildeL_ds[ii] = denom*(dL12[ii] + coef*dL22[ii] )
            #     print(ii, dtildeL_ds[ii], dL12[ii], coef*dL22[ii], coef)
            # exit()
            return dtildeL_ds
        else:
            print("Error: variation should be \"dt\" or \"ds\". Currently using ", variation)
            exit()

    def ddistance_soap(self, chi1, chi2, variation, cost_type="linear", rjs=None, d_des=0.0):
        """
        differentiate soap distance 

        Paramters
        ----------
        chi1 : CrystalStructure object
            target structure
        chi2 : CrystalStructure object
            updated structure
        variation : str
            "dt": differentiation of the components of a deformation matrix
            "dr": Differentiation of atomic coordinates
        cost_type : str
            type of cost function
            "linear": distance d、dD=|dd^2|/|dd^2|^{1/2}
            others: distance d^2、dd^2
        rjs : ndarray
            coordinates of atoms
        d_des : float
            target distance

        Returns
        -------
        ddis :  if variation == "dt":
                    (3; xi, 3; eta) numpy object
                if variation == "dr":
                    (j; num of atoms, 3; eta) numpy object
        dsoap : float
        """
        c = d_des
        dsoap = self.distance(chi1,chi2)
        # print("  c = {},  d_soap = {:.5f}".format(c,dsoap))

        if variation == "dr":
            ddis_ds  = -2.0*self.dsimilarity(chi1, chi2, "ds")
            ds_dr = self.ds_dr(chi2, rjs)
            ddis = np.zeros((len(rjs), 3))
            for j in range(len(rjs)):
                for eta in range(3):
                    for ii in range(len(ddis_ds)):
                        ddis[j, eta] += ds_dr[j, ii, eta]*ddis_ds[ii]

        else:
            ddis  = -2.0*self.dsimilarity(chi1, chi2, variation)

        if dsoap==0.0:
            ddis = np.zeros_like(ddis)
        else:
            ddis = (1.0 - c/dsoap) * ddis

        if cost_type=="linear":
            denom = np.where(ddis == 0.0, 1.0, ddis)
            ddis  = ddis/abs(denom)**0.5

        # print(ddis)
        return ddis, dsoap


    # @stop_watch
    def ddistance_ps(self, chi1, chi2, variation, cost_type="linear", l11=None, rjs=None, d_des=0.0):
        """
        Derivative of distance of Power spectra
          (P1-P2)^2/P1^2 = (L11+L22-2L12)/L11

        Paramters
        ----------

        chi1 : CrystalStructure object
            target structure
        chi2 : CrystalStructure object
            updated structure
        variation : str
            "dt": differentiation of the components of a deformation matrix
            "dr": Differentiation of atomic coordinates
        cost_type : str
            type of cost function
            "linear": distance d、dD=|dd^2|/|dd^2|^{1/2}
            others: distance d^2、dd^2
        l11 : float
            integral between target structures
        rjs : ndarray
            coordinates of atoms
        d_des : float
            target distance

        Returns
        -------
        _dcost: if variation == "dt":
                    (3; xi, 3; eta) numpy object
                if variation == "dr":
                    (j; num of atoms, 3; eta) numpy object
        dps: float
        """
        c = d_des
        ps1d_1 = self.get_1d_power_spectrum(chi1)
        ps1d_2 = self.get_1d_power_spectrum(chi2)
        if variation == "dt":
            dpp = self.dpower_spectrum_dt(chi2)
        if variation == "dr":
            dpp = self.dpower_spectrum_ds(chi2)

        if l11==None:
            l11 = self.overlapInteg(chi1, chi1)
        # l12 = self.overlapInteg(chi1, chi2)
        # l22 = self.overlapInteg(chi2, chi2)
        l12 = np.dot(ps1d_1.conjugate(), ps1d_2).real
        l22 = np.dot(ps1d_2.conjugate(), ps1d_2).real
        if l11+l22-2.0*l12 < 0.0:
            dps = 0.0
        else:
            # dps = math.sqrt((l11+l22-2.0*l12)/l11)
            dps = math.sqrt((l11+l22-2.0*l12)/math.sqrt(l11*l22))
        # print("  c = {},  d_ps = {:.5f}".format(c,dps))

        # coef = -(l11-l12)/(2.0*l11*l12) * (1.0 + l11/l12)
        if variation == "dr":
            dL12 = self.doverlapInteg(chi1, chi2, "ds", ps1d_1, dpp)
            dL22 = self.doverlapInteg(chi2, chi2, "ds", ps1d_2, dpp)
        else:
            dL12 = self.doverlapInteg(chi1, chi2, variation, ps1d_1, dpp)
            dL22 = self.doverlapInteg(chi2, chi2, variation, ps1d_2, dpp)

        if variation == "dt":
            _dcost = np.zeros((3, 3))
            for xi in range(3):
                for eta in range(3):
                    # _dcost[xi,eta] = 2.0*( dL22[xi, eta] - dL12[xi, eta] ) /l11
                    _dcost[xi,eta] = 2.0*( dL22[xi, eta] - dL12[xi, eta] )/math.sqrt(l11*l22) \
                                    -(l11+l22-2.0*l12)/l22*dL22[xi, eta]/math.sqrt(l11*l22)
        elif variation == "ds":
            num_ii = len(dL22)
            _dcost = np.zeros(num_ii)
            for ii in range(num_ii):
                _dcost[ii] = 2.0*(dL22[ii] - dL12[ii]) /l11
        elif variation == "dr":
            sw_multi=False
            bool_list=[True]*len(rjs)
            try:
                pair=chi2.pair
            except:
                pass
            else:
                sw_multi=True
                sp_names=chi2.sp_names
                bool_list=[True if s in pair else False for s in sp_names]
            num_ii = len(dL22)
            ds_dr = self.ds_dr(chi2, rjs, multi=sw_multi)
            _dcost = np.zeros((len(rjs), 3))
            for j in range(len(rjs)):
                if bool_list[j]==False:
                    continue
                dL_ds = (2.0*(dL22 - dL12) - (l11+l22-2.0*l12)/l22*dL22) / math.sqrt(l11*l22)
                for eta in range(3):
                    _dcost[j, eta] = np.dot(ds_dr[j,:,eta], dL_ds)

        else:
            print("Error: variation should be \"dt\" or \"ds\". Currently using ", variation)
            exit()

        if dps==0.0:
            _dcost = np.zeros_like(_dcost)
        else:
            _dcost = (1.0 - c/dps) * _dcost

        if cost_type=="linear":
            denom = np.where(_dcost == 0.0, 1.0, _dcost)
            _dcost  = _dcost/abs(denom)**0.5

        # print(_dcost)
        return _dcost, dps

    def ds_dr(self, chi2, rjs, multi=False):
        """
        differentiation of atomic coordinates of point intensity

        Paramters
        ----------
        chi2 : _CrystalStructure object
            updated structure
        rjs : ndarray
            coordinate of atoms
        multi : bool
            multi-element case : True

        Returns
        --------
        ds_dr : (j; num of atoms, num_ii; num of G points, 3; eta) numpy object
        """
        from _soap import rfunction
        
        exp=math.e
        sigma=chi2.sigma
        vu=chi2.lattice.reciprocal_lattice.volume
        tol_F=1.e-20

        ds_dr = np.zeros((len(rjs), len(chi2.gvecs), 3))
        gvecs = chi2.gvecs
        ## Multi species treatment
        if multi==True:
            pair=chi2.pair
            sp_names=chi2.sp_names
            bool_list=[True if s in pair else False for s in sp_names]
            comp_sp_names=list(itertools.compress(sp_names, bool_list))
            _rjs=rjs[bool_list]

            Ajs=[]
            for i, rj in enumerate(_rjs):
                if comp_sp_names[i]==pair[0]:
                    Ajs.append(-1.0)
                    continue
                Ajs.append(1.0)
            
            ds_dr[bool_list,:,:] = rfunction.ds_dr(_rjs, gvecs, Ajs, sigma, vu, tol_F)

        ## Single species treatment
        else:
            # for j, rj in enumerate(rjs):
            #     for i, gvec in enumerate(chi2.gvecs):
            #         F=0
            #         for rj_for_F in rjs:
            #             F += exp**( -1.0j*np.dot(gvec, rj_for_F) )
            #         if abs(F) < tol_F:
            #             # print("-----Chack-----", abs(F))
            #             dchig = -gvec*math.sin(np.dot(gvec, rj))
            #         else:
            #             dchig = 1.0j*gvec/2./abs(F)*(exp**(1.0j*np.dot(gvec, rj))*F
            #                                     - exp**(-1.0j*np.dot(gvec, rj))*F.conjugate())
            #             dchig = dchig.real
            #         dchig = dchig*exp**(-1.0*np.dot(gvec, gvec)*sigma*sigma/2.0) / vu
            #         ds_dr[j, i] = dchig
            ###
            Ajs = [1.0 for elem in rjs]
            ds_dr = rfunction.ds_dr(rjs, gvecs, Ajs, sigma, vu, tol_F)
        return ds_dr



if __name__ == "__main__":

    struct_path1 = './examination/MgO_1.cif'
    struct_path2 = './examination/MgO_2.cif'

    # SOAP parameters
    nbasis = 10
    nlmax  = 4
    sigma_soap  = 0.4
    nrad   = 60
    ltype  = 3
    center_element_name = "G"

    _chi1 = CrystalStructure.from_file(struct_path1)#.struct_scale()
    _chi2 = CrystalStructure.from_file(struct_path2)#.struct_scale()

    sigma_rtrans = sigma_soap
    chi1  = _chi1.rtransform(sigma=sigma_rtrans, multi_species=True)
    chi2  = _chi2.rtransform(sigma=sigma_rtrans, multi_species=True)

    rcut = chi1.rcut
    soap = Soap(nbasis, nlmax, rcut, sigma_soap, nrad, ltype, center_element_name)

    print("rcut=", rcut)
##
    print("normal SOAP")
    print(soap.overlapInteg(chi1.chis[0],chi2.chis[0]))
    print(soap.similarity(chi1.chis[0],chi2.chis[0]))
    ##
    #print("dL/dt")
    #dL_dt = soap.doverlapInteg(chi1.chis[0], chi2.chis[0], variation = "dt")
    #for xi in range(3):
    #    for eta in range(3):
    #        print("   xi=",xi,"eta=",eta,": ",dL_dt[xi,eta])
    #        dL_dt[xi,eta]
    ##
    #print("dL/ds")
    #dL_ds = soap.doverlapInteg(chi1.chis[0], chi2.chis[0], variation = "ds")
    #for ii, dL_dsi in enumerate(dL_ds):
    #    print("   ii=", ii, ": ", dL_ds[ii])


    ### soap distance ###
    distance_ps=0.0
    distance_soap=0.0
    for i, (c1, c2) in enumerate(zip(chi1.chis,chi2.chis)):
        l11 = soap.overlapInteg(c1,c1)
        l12 = soap.overlapInteg(c1,c2)
        l22 = soap.overlapInteg(c2,c2)

        distance_ps2 =  (l11+l22-2.0*l12)/math.sqrt(l11*l22)
        distance_ps += math.sqrt(max(distance_ps2,0.0))
    print(" --- struc1, struc2, Distance = ", struct_path1, ',', struct_path2, ',', distance_ps)

    print("Normal ends.")
