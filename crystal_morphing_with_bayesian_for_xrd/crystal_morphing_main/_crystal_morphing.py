import sys
import os
import copy
import numpy as np
import math
import cmath
import itertools
from scipy.optimize import minimize, check_grad, approx_fprime
from crystal_structure import CrystalStructure, _CrystalStructure

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

# @stop_watch
def steepest_descent(self, itenum, chi1, chi2, step, d_des):
    """
    Update coordinate of atoms and reciprocal lattice vector by SD method.

    Paramters
    ----------
    itenum : int
        iteration number
    chi1 : CrystalStructure object
        target structure
    chi2 : CrystalStructure object
        updated structure
    step : str
        "t": update reciprocal lattice vector 
        "r": update coordinate of atoms
    d_des : float
        target distance

    Returns
    --------
    _chi : CrystalStructure object
        updated structure
    """
    self.d_des=d_des
    _check_reciprocal(chi1, chi2)

    if self.verbose:
        if step=="t":
            print("bb = \n", chi2.chis[0].lattice.matrix)
            # print("aa = \n", chi2.chis[0].lattice.reciprocal_lattice.matrix)
        if step=="r":
            print("rj = \n", chi2.rjs)

    if step=="t":
        bb = chi2.chis[0].lattice.matrix
        _dcost = self.jac(bb.flatten(), chi1, chi2, step, sd=True)
        _dcost = np.reshape(_dcost, (len(_dcost)//3, 3))

        # Update bb
        dbb   = _bb_prodcut_mat(bb, _dcost, itenum, self.mode, delta_t=self.delta_t, sd=True)
        _bb   = bb - dbb
        _chi  = self.update_chi(chi2, _bb.flatten(), step)

    elif step=="r":
        rjs = chi2.rjs
        _dcost = self.jac(rjs.flatten(), chi1, chi2, step, sd=True)
        _dcost = np.reshape(_dcost, (len(_dcost)//3, 3))

        ## Update rjs
        _rjs = np.zeros_like(rjs)
        for j, rj in enumerate(rjs):
            # _rjs[j]  = rj - (self.delta_r/itenum**0.5) * (_dcost[j] - _dcost[0])
            _rjs[j]  = rj - self.delta_r * (_dcost[j] - _dcost[0])
        _chi  = self.update_chi(chi2, _rjs.flatten(), step)

    else:
        print("Error: step should be \"t\" or \"r\". Currently using ", step)
        exit()

    return _chi


def quasi_newton(self, itenum, chi1, chi2, step, d_des,
                 maxiter=1, method='L-BFGS-B'):
    """
    Update coordinate of atoms and reciprocal lattice vector by Quasi-Newton method.

    Paramters
    ----------
    itenum : int
        iteration number
    chi1 : CrystalStructure object
        target structure
    chi2 : CrystalStructure object
        updated structure
    step : str
        "t": update reciprocal lattice vector
        "r": update coordinate of atoms
    d_des : float
        target distance
    maxiter : int
        maximum number of iterations for BFGS method
    method : str
        optimization method


    Returns
    --------
    _chi2 : CrystalStructure object
        updated structure
    """
    def printx(x):
        if self.verbose:
            print(" ---------- Callback: update parameters ----------\n")
            # print(np.reshape(x, (len(x)//3, 3)))

    self.d_des=d_des
    _check_reciprocal(chi1, chi2)

    if step=="t":
        bb  = chi2.chis[0].lattice.matrix
        bb  = bb.flatten()
        res = minimize(self.fun, bb, args=(chi1,chi2,step),
                        method=method, jac=self.jac,
                        options={'gtol':1e-5,'maxls':10,'maxiter':maxiter}, callback=printx)
    elif step=="r":
        rjs = chi2.rjs
        r0_z = rjs[0,2]
        rjs  = rjs.flatten()
        res = minimize(self.fun, rjs, args=(chi1,chi2,step),
                        method=method, jac=self.jac,
                        options={'ftol':1e-20,'gtol':1e-10,'maxls':10,'maxiter':maxiter}, callback=printx)
        rjs = np.reshape(res.x, (len(res.x)//3, 3))
        # rjs[:,:2] -= rjs[0,:2]
        # rjs[:,2] -= rjs[0,2] - r0_z
        rjs[:,:3] -= rjs[0,:3] # for 3d
        res.x = rjs.flatten()
    else:
        print("Error: step should be \"t\" or \"r\". Currently using ", step)
        exit()

    _chi2 = update_chi(self, chi2, res.x, step)

    return _chi2

# @stop_watch
def fun(self, x, chi1=None, chi2=None, step=None):
    """
    caclulate target function value
    Parameters
    -----------
    x : 1D array
        Input matrix. The reciprocal lattice vector or coordinate of atom are descripted as 1D according to the step.
    chi1 : CrystalStructure object
        target structure
    chi2 : CrystalStructure object
        updated structure
    step : str
        "t": update reciprocal lattice vector
        "r": update coordinate of atoms

    Returns
    --------
    res : float
        target function 
    """
    _chi2 = self.update_chi(chi2, x, step)
    if self.verbose:
        if step=="t":
            print("bb = \n", _chi2.chis[0].lattice.matrix)
        if step=="r":
            print("rj = \n", _chi2.rjs)

    res=0.0
    for i, (c1, c2) in enumerate(zip(chi1.chis, _chi2.chis)):
        if self.dis_type=="soap":
            distance = self.soap.distance(c1,c2)
        else:
            # l11 = self.soap.overlapInteg(chi,chi)
            l11 = self.l11[i]
            l12 = self.soap.overlapInteg(c1,c2)
            l22 = self.soap.overlapInteg(c2,c2)
            # distance_ps_square = (l11+l22-2.0*l12)/l11
            distance_ps_square = (l11+l22-2.0*l12)/math.sqrt(l11*l22)
            if distance_ps_square<0.0:
                distance_ps_square=0.0
            distance = math.sqrt(distance_ps_square)
        res += distance

    res = abs(res - self.d_des)**2
    if self.verbose:
        print("\nf(x) =", res)
    return res

# @stop_watch
def jac(self, x, chi1=None, chi2=None, step=None, sd=False):
    """
    Calculate differentiation of target function

    Parameters
    -----------
    x : 1D array
        Input matrix. The reciprocal lattice vector or coordinate of atom are descripted as 1D according to the step.
    chi1 : CrystalStructure object
        target structure
    chi2 : CrystalStructure object
        updated structure
    step : str
        "t": update reciprocal lattice vector
        "r": update coordinate of atoms
    sd : bool
       True : SD method

    Returns
    --------
    res : float
        Differentiated value of target function
    """
    cost_type="square"
    if sd:
        _chi2 = chi2
    else:
        _chi2 = self.update_chi(chi2, x, step)

    if self.verbose:
        print("\n ----Calculate derivatives-----")
    if step=="t":
        bb=_chi2.chis[0].lattice.matrix
        _dcost=np.zeros_like(bb)
        distance=0.0
        for i, (c1, c2) in enumerate(zip(chi1.chis, _chi2.chis)):
            if self.multi_species and self.verbose:
                print("i =", i, ": dD_dt for pair", c1.pair)
            if self.dis_type=='soap':
                dD_dt, dis = self.soap.ddistance_soap(c1, c2, variation = "dt",
                                                    cost_type=cost_type)
            else:
                dD_dt, dis = self.soap.ddistance_ps(c1, c2, variation = "dt",
                                                    cost_type=cost_type, l11=self.l11[i])
            _dcost += dD_dt
            distance += dis
        if self.verbose:
            print(" --> sum all {} components".format(len(_chi2.chis)))
        c=self.d_des
        if distance!=0.0:
            _dcost = (1.0 - c/distance) * _dcost

        if sd:
            denom = np.where(_dcost==0.0, 1.0, _dcost)
            _dcost = _dcost/abs(denom)**0.5
        else:
            _dcost = _bb_prodcut_mat(bb, _dcost, 1, self.mode)

    elif step=="r":
        _dcost=np.zeros_like(_chi2.rjs)
        distance=0.0
        for i, (c1, c2) in enumerate(zip(chi1.chis, _chi2.chis)):
            if self.multi_species and self.verbose:
                print("i =", i, ": dD_dr for pair", c1.pair)
            if self.dis_type=='soap':
                dD_dr, dis = self.soap.ddistance_soap(c1, c2, variation = "dr",
                                                    cost_type=cost_type, rjs=_chi2.rjs)
            else:
                dD_dr, dis = self.soap.ddistance_ps(c1, c2, variation = "dr",
                                                    cost_type=cost_type, l11=self.l11[i], rjs=_chi2.rjs)
            _dcost += dD_dr
            distance += dis
        if self.verbose:
            print(" --> sum all {} components".format(len(_chi2.chis)))
        c=self.d_des
        if distance!=0.0:
            _dcost = (1.0 - c/distance) * _dcost

        if sd:
            denom = np.where(_dcost==0.0, 1.0, _dcost)
            _dcost = _dcost/abs(denom)**0.5

        if self.mode=="xy":
            _dcost = _dcost * np.array([1,1,0])

    else:
        print("Error: step should be \"t\" or \"r\". Currently using ", step)
        exit()

    if self.verbose:
        print("jac = \n", _dcost)
        print(" -------- \n")
    res = _dcost.flatten()
    return res

# @stop_watch
def update_chi(self, chi2, x, step):
    """
    Update information of crystal structure

    Parameters
    -----------
    chi2 : CrystalStructure object
        updated structure
    x : 1D array
        Input matrix. The reciprocal lattice vector or coordinate of atom are descripted as 1D according to the step.
    step : str
        "t": update reciprocal lattice vector
        "r": update coordinate of atoms

    Returns
    --------
    _chi2 : CrystalStructure object
        Updated CrystalStructure object。
    """
    from _operator_local import _get_ref_gvecs, _get_ref_chigs, _index_in_rcut
    from _soap import rfunction
    ## local functions
    def _setup(chi2,x,step):
        xx = np.reshape(x, (len(x)//3, 3))
        if step=="t":
            bb=xx
            rjs=chi2.rjs
        elif step=="r":
            bb = chi2.chis[0].lattice.matrix
            rjs = xx
            chi2.rjs = rjs
        else:
            print("Error: step should be \"t\" or \"r\". Currently using ", step)
            exit()
        return  chi2, bb, rjs
    
    def _update_cs(bb, r_elems, mvecs, chi, gvecs, si, chigs, ref_chigs ):
        _cs = _CrystalStructure(bb, r_elems, mvecs)
        _cs.r_elems = r_elems
        _cs.rcut    = chi.rcut
        _cs.gvecs   = gvecs
        _cs.si      = si
        _cs.chigs   = chigs
        _cs.sigma   = chi.sigma
        _cs.ref_mvecs = chi.ref_mvecs
        _cs.ref_chigs = ref_chigs
        return _cs

    def _t_loop(chi, bb, iis):
        mvecs=[]
        gvecs=[]
        chigs=[]
        for ii in iis:
            mm=chi.ref_mvecs[ii]
            gvec = mm@bb
            chig = chi.ref_chigs[ii]
            mvecs.append(mm)
            gvecs.append(gvec)
            chigs.append(chig)
        return mvecs, gvecs, chigs
    
    ## Multiple species treatment
    if self.multi_species:
        chi2, bb, rjs = _setup(chi2,x,step)
        sp_names=chi2.chis[0].sp_names
        _chi2_list=[]
        for chi in chi2.chis:
            exp=math.e
            sigma=chi.sigma
            vu=chi.lattice.reciprocal_lattice.volume
            mvecs=[]
            gvecs=[]
            chigs=[]
            ref_chigs=[]

            pair=chi.pair
            bool_list=[True if s in pair else False for s in sp_names]
            comp_sp_names=list(itertools.compress(sp_names, bool_list))

            rcut = chi.rcut
            ref_gvecs = _get_ref_gvecs(chi.ref_mvecs, bb)
            iis = _index_in_rcut(ref_gvecs, rcut)
            if step=="t":
                mvecs, gvecs, chigs = _t_loop(chi,bb, iis)                        
                ref_chigs = chi.ref_chigs
            elif step=="r":
                Ajs=[]
                for i, rj in enumerate(rjs[bool_list]):
                    if comp_sp_names[i]==pair[0]:
                        Ajs.append(-1.0)
                        continue
                    Ajs.append(1.0)
                 
                ref_chigs = _get_ref_chigs(ref_gvecs, rjs[bool_list], Ajs, sigma, vu)
                
                for ii in iis:
                    mm=chi.ref_mvecs[ii]
                    gvec = ref_gvecs[ii]
                    chig = rfunction.chig(rjs[bool_list], gvec, Ajs, sigma, vu)
                    mvecs.append(mm)
                    gvecs.append(gvec)
                    chigs.append(chig)
            else:
                print("Error: step should be \"t\" or \"r\". Currently using ", step)
                exit()
            si = np.abs(chigs)                
            r_elems=[]
            for gvec in gvecs:
                if np.linalg.norm(gvec) == 0.0:
                    r_elems.append("G")
                else:
                    r_elems.append("X")
            
            ### update
            _cs = _update_cs(bb, r_elems, mvecs, chi, gvecs, si, chigs, ref_chigs )
            _cs.sp_names= chi.sp_names
            _cs.pair    = chi.pair

            _chi2_list.append(_cs)

        chi2.chis = _chi2_list
        
        return chi2

    ## Single species treatment
    else:
        chi2, bb, rjs = _setup(chi2,x,step)

        chi = chi2.chis[0]

        exp=math.e
        sigma=chi.sigma
        vu=chi.lattice.reciprocal_lattice.volume
        mvecs=[]
        gvecs=[]
        chigs=[]
        ref_chigs=[]
        
        rcut = chi.rcut
        ref_gvecs = _get_ref_gvecs(chi.ref_mvecs, bb)
        iis = _index_in_rcut(ref_gvecs, rcut)
        if step=="t":    
            mvecs, gvecs, chigs = _t_loop(chi,bb, iis)                        
            ref_chigs = chi.ref_chigs
        elif step=="r":
            Ajs = [1.0 for elem in rjs]
            ref_chigs = _get_ref_chigs(ref_gvecs, rjs, Ajs, sigma, vu)
                
            for ii in iis:
                mm=chi.ref_mvecs[ii]
                gvec = ref_gvecs[ii]
                chig = rfunction.chig(rjs, gvec, Ajs, sigma, vu)
                mvecs.append(mm)
                gvecs.append(gvec)
                chigs.append(chig)        
        else:
            print("Error: step should be \"t\" or \"r\". Currently using ", step)
            exit()
        si = np.abs(chigs)
        r_elems=[]
        for gvec in gvecs:
            if np.linalg.norm(gvec) == 0.0:
                r_elems.append("G")
            else:
                r_elems.append("X")

        ### update
        _cs = _update_cs(bb, r_elems, mvecs, chi, gvecs, si, chigs, ref_chigs )

        chi2.chis = [_cs]

        return chi2



def _bb_prodcut_mat(bb, mat, itenum, mode, delta_t=1.0, sd=False):
    """
    Calculate matrix product between reciprocal lattice vector and transformation matrix

    Parameters
    ----------
    bb : ndarray
        Reciprocal lattice vector
    mat : ndarray
        transformation matrix
    itenum : int
        number of iterations
    mode : str
        How to multiply a transform matrix
    delta_t : float
        step size

    Returns
    --------
    _bb : ndarray
        Matrix product of reciprocal lattice vector and transform matrix.
    """
    if mode == "all":
        coef = 1.0
        # coef_mat = (delta_t/itenum**0.5) * np.array([[1.0,coef,coef],[coef,1.0,coef],[coef, coef, 1.0]])
        coef_mat = delta_t * np.array([[1.0,coef,coef],[coef,1.0,coef],[coef, coef, 1.0]])
    elif mode == "diagonal":
        coef = 0.0
        coef_mat = delta_t*np.array([[1.0,coef,coef],[coef,1.0,coef],[coef, coef, 1.0]])
    elif mode == "non-diagonal":
        coef = 1.0
        coef_mat = delta_t*np.array([[0.0,coef,coef],[coef,0.0,coef],[coef, coef, 0.0]])
    elif mode == "xy":
        diag, non = 1.0, 1.0
        if sd:
            diag = 0.1
        # coef_mat = (delta_t/itenum**0.5)*np.array([[diag,non,0.0],[non,diag,0.0],[0.0, 0.0, 0.0]])
        coef_mat = delta_t*np.array([[diag,non,0.0],[non,diag,0.0],[0.0, 0.0, 0.0]])
    else:
        print("Error: mode can take all, diagonal, non-diagonal. Currently using ", mode)

    _bb = np.zeros((3,3))
    for xi in range(3):
        _bxi = np.zeros(3)
        for xid in range(3):
            _bxi += bb[xid]*mat[xid,xi] *coef_mat[xid,xi]
        _bb[xi] = _bxi

    return _bb


def _check_reciprocal(chi1, chi2):
    """
    Check CrystalStructure object is in reciprocal space

    Parameters
    ----------
    chi1 : CrystalStructure object
        Target structure.
    chi2 : CrystalStructure object
        Updated structure.
    """
    if chi1.is_reciprocal() == False:
        print("Error: CrystalInterpolate routines are not allowed for real space chi.")
        exit()
    if chi2.is_reciprocal() == False:
        print("Error: CrystalInterpolate routines are not allowed for real space chi.")
        exit()
    return
