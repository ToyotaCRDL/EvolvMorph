import sys
import os

homepath = os.path.join(os.path.dirname(__file__), '../..')
sys.path.append(os.path.join(homepath, "crystal_structure"))
sys.path.append(os.path.join(homepath, "soap"))

from crystal_structure import CrystalStructure
from soap import Soap


struct_path1 = './CONTCAR_square_test'
struct_path2 = './CONTCAR_smallsquare_test'

#
nbasis = 10
nlmax  = 4
sigma_soap  = 0.5
nrad   = 60
ltype  = 3
center_element_name = "G"

_chi1 = CrystalStructure.from_file(struct_path1)#.struct_scale()
_chi2 = CrystalStructure.from_file(struct_path2)#.struct_scale()

sigma_rtrans = sigma_soap
chi1  = _chi1.rtransform(sigma=sigma_rtrans)
chi2  = _chi2.rtransform(sigma=sigma_rtrans)

rcut = chi1.rcut
soap = Soap(nbasis, nlmax, rcut, sigma_soap, nrad, ltype, center_element_name)

print("rcut=", rcut)

##-------------------------------------------------------------

test_flag  = "test0"

'''
reciprocal soap のテスト
Test 0-3:
    "test0": execution test
    "test1": L(chi1, chi2)                 = L(xhi, R chi2)
    "test2": d L(chi1, chi2=chi1)/dchi2(t) = 0
    "test3": dL(chi1, chi2)/dchi2          = dL(chi1, R chi2)/dchi2
axis specifies the axis to rotate in test1,3
'''

axis = [1,0,0]

print(test_flag)

if test_flag == "test0":
    print(chi1)
    print(chi2)
    print("bases")
    print(chi1.lattice.matrix)
    print(chi2.lattice.matrix)

    ##
    print("normal SOAP")
    print("        L12:", soap.overlapInteg(chi1,chi2))
    print("  tilde L12: ", soap.similarity(chi1,chi2))
    ##
    print("dL11/dt")
    dL_dt = soap.doverlapInteg(chi1, chi1, variation = "dt")
    for xi in range(3):
        for eta in range(3):
            print("   xi=",xi,"eta=",eta,": ",dL_dt[xi,eta])
            dL_dt[xi,eta]
    print("dL12/dt")
    dL_dt_a = soap.doverlapInteg(chi1, chi2, variation = "dt")
    for xi in range(3):
        for eta in range(3):
            print("   xi=",xi,"eta=",eta,": ",dL_dt_a[xi,eta])
            dL_dt_a[xi,eta]
    print("dL22/dt")
    dL_dt_b = soap.doverlapInteg(chi2, chi2, variation = "dt")
    for xi in range(3):
        for eta in range(3):
            print("   xi=",xi,"eta=",eta,": ",dL_dt_b[xi,eta])
            dL_dt_b[xi,eta]
    print("dsimilarity")
    dS_dt = soap.dsimilarity(chi1, chi2, variation = "dt")
    print(dS_dt)
    ##
    print("dL/ds")
    dL_ds = soap.doverlapInteg(chi1, chi2, variation = "ds")
    for ii, dL_dsi in enumerate(dL_ds):
        print("   ii=", ii, ": ", dL_ds[ii])

elif test_flag == "test1":
    output1 = soap.overlapInteg(chi1, chi2 )
    for deg in range(0, 361, 40):
        rt = stf.RotationTransformation(axis,deg)
        chi2r = rt.apply_transformation(_chi2)
        chi2r = chi2r.rtransform(sigma=sigma_rtrans)

        output2 = soap.overlapInteg(chi1, chi2r)
        print("   ",deg,"deg: ",output1,output2,output2/output1)

elif test_flag == "test2":
    print("> chi1,chi1")
    dL_dt = soap.doverlapInteg(chi1, chi1, variation = "dt")
    for xi in range(3):
        for eta in range(3):
            print("   xi=",xi,"eta=",eta,": ",dL_dt[xi,eta])

    print("> chi1,chi2")
    dL_dt = soap.doverlapInteg(chi1, chi2, variation = "dt")
    for xi in range(3):
        for eta in range(3):
            print("   xi=",xi,"eta=",eta,": ",dL_dt[xi,eta])

elif test_flag == "test3":
    output1_dt = soap.doverlapInteg(chi1, chi2,  variation="dt")
    output1_ds = soap.doverlapInteg(chi1, chi2,  variation="ds")
    for deg in range(0, 361, 40):

        rt = stf.RotationTransformation(axis,deg)
        chi2r = rt.apply_transformation(_chi2)
        chi2r = chi2r.rtransform(sigma=sigma_rtrans)

        print("> doverlapInteg_dt")
        output2 = soap.doverlapInteg(chi1, chi2r, variation="dt")
        print("   ",deg,"deg: ",output1_dt,output2,output2/output1_dt)

        print("> doverlapInteg_ds")
        output2 = soap.doverlapInteg(chi1, chi2r, variation="ds")
        print("   ",deg,"deg: ",output1_ds,output2,output2/output1_ds)


print("Normal ends.")
