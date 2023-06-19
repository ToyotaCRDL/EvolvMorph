import sys
import os
import copy
import numpy as np
import math
import cmath
import itertools
from scipy.optimize import minimize, check_grad, approx_fprime

##### S.Kajita and J.Oba construced crystal morphing code.    Ref : J. Oba, S. Kajita, Phys. Rev. Mater. 6, 023801 (2022).

homepath = os.path.dirname(__file__)
sys.path.append(os.path.join(homepath, "crystal_structure"))
sys.path.append(os.path.join(homepath, "soap"))

class CrystalInterpolate():
    """
    To interpolate between different crystal structures.
    """

    from _crystal_morphing import steepest_descent, quasi_newton, update_chi, fun, jac

    def __init__(self, soap, delta_t, delta_r, dis_type, mode, sub_steps, chi1, verbose=True):
        """
        Parameters
        ----------
        soap : Soap object
            Do overlap integration of SOAP, distance, power spectrum, and so on.
        delta_t : float
            Arbitrary coefficient for SD method to update reciprocal lattice vector (step size)
        delta_r : float
            Arbitrary coefficient for SD method to update coordinate of atoms (step size)
        dis_type : str
            type of distance
        mode : str
            mode for updating reciprocal lattice vector(_bb_prodcut_mat in _crystal_morphing.py)。
        sub_steps : int
            t, r step for one iteration for SD method
        chi1 : CrystalStructure object
            CrystalStructure object of target structure
        """
        self.soap          = soap
        self.delta_t       = delta_t
        self.delta_r       = delta_r
        self.dis_type      = dis_type
        self.mode          = mode
        self.sub_steps     = sub_steps
        self.multi_species = chi1.multi_species
        self.verbose = verbose

        l11_list=[]
        for chi in chi1.chis:
            l11 = self.soap.overlapInteg(chi,chi)
            l11_list.append(l11)
        self.l11 = l11_list

    def optimize(self, itenum, to_chi, from_chi, d_des, do_SD):
        """
        Update reciprocal lattice vector and coordinates of atoms

        Paramters
        ----------
        itenum : int
            number of iterations
        to_chi : CrystalStructure object
            target structure
        from_chi : CrystalStructure object
            updated structure
        d_des : float
            target distance
        do_SD : bool
            True : using SD method

        Returns
        --------
        from_chi : CrystalStructure object
            updated structure
        """
        import time
        start = time.time()
        if self.verbose:
            print("\n     Destination : ", d_des)
        if do_SD:
            # steepest decent
            if self.verbose:
                print("\n                --- t step "+self.mode, "---")
            for _ in range(self.sub_steps):
                from_chi = self.steepest_descent(itenum, to_chi, from_chi, "t", d_des=d_des)
            if self.verbose:
                d_soap, d_ps = self.report_metrics(to_chi, from_chi)
            from_chi = self.update_chi(from_chi, from_chi.rjs.flatten(), "r")

            if self.verbose:
                print("\n                --- r step ---")
            for _ in range(self.sub_steps):
                from_chi = self.steepest_descent(itenum, to_chi, from_chi, "r", d_des=d_des)
        else:
            # BFGS
            maxiter=int(math.sqrt(itenum))*5

            if self.verbose:
                print("\n                    t step "+self.mode)
            from_chi = self.quasi_newton(itenum, to_chi, from_chi, "t", d_des=d_des, maxiter=maxiter)
            if self.verbose:
                d_soap, d_ps = self.report_metrics(to_chi, from_chi)

            from_chi = self.update_chi(from_chi, from_chi.rjs.flatten(), "r")

            if self.verbose:
                print("\n                    r step ")
            from_chi = self.quasi_newton(itenum, to_chi, from_chi, "r", d_des=d_des, maxiter=maxiter)
        elapsed_time =  time.time() - start
#        print("")
#        print(f"----------- Iteration: {round(float(elapsed_time),3)} sec")
        return from_chi


    def report_metrics(self, chi1, chi2, bool_result=False):
        """
        report during transformation

        Parameters
        ----------
        chi1 : CrystalStructure object
            target structure
        chi2 : CrystalStructure object
            updated structure

        Returns
        -------
        distance_soap : float
            SOAP distance
        distance_ps : float
            power spectrum distance
        """        
        if self.verbose or bool_result:
            print("  ------------------ Diff report ------------------")
            print("    chi2  =>  chi1")
            lengths1 = chi1.chis[0].lattice.lengths
            angles1  = chi1.chis[0].lattice.angles
            lengths2 = chi2.chis[0].lattice.lengths
            angles2  = chi2.chis[0].lattice.angles
            print("    Bases")
            print("      lengths=", lengths2, "=>", lengths1)
            print("      angles =", angles2,  "=>", angles1)

        distance_ps=0.0
        distance_soap=0.0
        for i, (c1, c2) in enumerate(zip(chi1.chis,chi2.chis)):
            # l11 = self.soap.overlapInteg(chi1,chi1)
            l11 = self.l11[i]
            l12 = self.soap.overlapInteg(c1,c2)
            l22 = self.soap.overlapInteg(c2,c2)

            distance_soap2 = 2.0 - 2.0*l12/(l11*l22)**0.5
            distance_soap += math.sqrt(max(distance_soap2,0.0))

            # distance_ps =  (l11+l22-2.0*l12)/l11
            distance_ps2 =  (l11+l22-2.0*l12)/math.sqrt(l11*l22)
            distance_ps += math.sqrt(max(distance_ps2,0.0))
            # print(math.sqrt(max(distance_ps2,0.0)), math.sqrt(max(distance_soap2,0.0)))
        if self.verbose or bool_result:
            print("    Distance")
            print("      Power Spectra = ", distance_ps)
            print("               SOAP = ",distance_soap )
            print("                L11 = ", l11)
            print("                L12 = ", l12)
            print("                L22 = ", l22)

            print()

        return distance_soap, distance_ps

    def coords_report(self, chi1, chi2):
        """
        Output point intensity and coordinates in the reciprocal space

        Parameters
        ----------
        chi1 : CrystalStructure object
            target structure
        chi2 : CrystalStructure object
            updated structure
        """
        def _coords_report(chi):
            if hasattr(chi, 'pair'):
                print("      For pair :", chi.pair)
            fracs = chi.frac_coords
            gvecs = chi.gvecs
            elems = []
            for site in chi.sites:
                elems.append( site.specie.symbol )
            si = chi.get_si()
            for elem, frac, gvec, _si in zip(elems, fracs, gvecs, si):
                gvec = ['{:.5f}'.format(gv) for gv in gvec]
                _si = '{:.5f}'.format(_si)
                print("      ", elem,  _si, "(", "  ".join(gvec), ")", frac)

        print("    Coordinates")
        print("      chi1=")
        for c1 in chi1.chis:
            _coords_report(c1)
        print("      chi2=")
        for c2 in chi2.chis:
            _coords_report(c2)



if __name__ == "__main__":
    import time
    start = time.time()

    from crystal_structure import CrystalStructure
    from soap import Soap

    ##--- input structure file ---##
    # struct_path1 = './examination/CONTCAR_square_test'
    # struct_path2 = './examination/CONTCAR_bigsquare_test'
    # struct_path1 = './examination/C_graphen.cif'
    # struct_path2 = './examination/C_hexagonal.cif'
    # multi_species = False


    #--- input info loading ---#    220309
    import pandas as pd
    import numpy as np

    inputfile = 'INPUT.csv'

    inputdataraw = open(inputfile,'r').readlines()

    inputdata = []
    for i in range(len(inputdataraw)):
      l = inputdataraw[i].replace('\n','').replace(' ','').split(',')
      inputdata.append(l)

    inputdata = pd.DataFrame(inputdata)

    try:
        if inputdata[inputdata.iloc[:,0] == 'multi_species'].iloc[0,1] == 'True': 
           multi_species = True
        else:
           multi_species = False
    except IndexError:
       multi_species = False

    struct_path1 = inputdata[inputdata.iloc[:,0] == 'struct1'].iloc[0,1] #'./examination/MgO_1.cif'
    struct_path2 = inputdata[inputdata.iloc[:,0] == 'struct2'].iloc[0,1] #'./examination/MgO_2.cif'
#    multi_species = True # True: multi species, False: single species

    #--- read structure file ---
    chi1_real = CrystalStructure.from_file(struct_path1)
    chi2_real = CrystalStructure.from_file(struct_path2)

    ##--- set SOAP parameters ---##
    nbasis = 10
    nlmax  = 6
    sigma_soap  = 0.5
    nrad   = 60
    ltype  = 3
    center_element_name = "G"

    #--- transform to the reciprocal space ---
    sigma_rtrans = sigma_soap
    print("\n chi1: real ---> reciprocal")
    chi1  = chi1_real.rtransform(sigma_rtrans, multi_species)
    print("\n chi2: real ---> reciprocal")
    chi2  = chi2_real.rtransform(sigma_rtrans, multi_species)

    #--- create Soap object ---
    rcut = chi1.rcut
    soap = Soap(nbasis, nlmax, rcut, sigma_soap, nrad, ltype, center_element_name)
    print("\n rcut=", rcut)

    ##--- set CrystalInterpolate parameters ---##
    try:
      max_steps = int(inputdata[inputdata.iloc[:,0] == 'soap_max_steps'].iloc[0,1]) #20  # 3
    except IndexError:max_steps = 20

#    max_steps = 20
    sub_steps = 3
    delta_t = 2.*10**(-2)
    delta_r = 2.*10**(-2)
    sw_SD = 0.2
    dis_type = "ps"
    mode = "all" # set "xy" for 2d crystal morphing
    verbose = False #True

    #--- create CrystalInterpolate object ---
    ci = CrystalInterpolate(soap, delta_t, delta_r, dis_type, mode, sub_steps, chi1, verbose)
    if verbose:
        ci.coords_report(chi1, chi2)

    ##--- set distance parameters ---##
    d_soap, d_ps = ci.report_metrics(chi1, chi2, bool_result=True)
    d_ini = d_ps if dis_type=="ps" else d_soap
    d_des = 0.0*d_ini
    tol_d = 1.e-2


    #--- main loop ---
    start_main = time.time()
    #------------------------------------------------------------------#
    print("\n Crystal Interpolation starts: Max iterations =", max_steps)
    for _itenum in range(max_steps):
        itenum = _itenum + 1
        print("\n <-- Number of iterations =", itenum)
        do_SD = True if d_ps>sw_SD else False

        #--- optimization ---
        chi2 = ci.optimize(itenum, chi1, chi2, d_des, do_SD)
        d_soap, d_ps = ci.report_metrics(chi1, chi2, bool_result=True)

        #--- output CIF file ---
        cif_name="examination/tmp/test_"+str(itenum)+".cif"
        chi2.write_CIF(cif_name) # normal output
        # chi2.write_CIF("scaled_"+cif_name, scaling=True, scale_to_chi=chi1_real, scaling_method="area") # scaling output
        
        #--- check convergence ---
        d_x = d_ps if dis_type=="ps" else d_soap
        if abs(d_des-d_x) < d_des*tol_d:
            break
    #------------------------------------------------------------------#
    elapsed_time = time.time() - start_main
    print("\n <-- End of the interations -->")

    ci.report_metrics(chi1, chi2, bool_result=True)
    if verbose:
        ci.coords_report(chi1, chi2)

#    print ("Time in the main iterations: {:.6f}".format(elapsed_time) + "[sec]")
    elapsed_time = time.time() - start
#    print ("Total time: {:.6f}".format(elapsed_time) + "[sec]")

    print("Normal ends.")
