import sys,os,glob,copy,re,shutil
import math,cmath
import warnings
import numpy as np
import matplotlib.pyplot as plt
import GPyOpt

import pandas as pd
import itertools

##### S.Kajita and J.Oba construced crystal morphing code.    Ref : J. Oba, S. Kajita, Phys. Rev. Mater. 6, 023801 (2022).
##### J.Lee and J.Oba modified the code for XRD application.        Ref : J. Lee, J. Oba, N. Ohba, S. Kajita, submitted., preprint arXiv 2302.10464 (2023).

##--- set path to your code ---##
path1 = "./crystal_morphing_main"    # modifying required
path2 = "./spectra"   # modified required

inputfile = "./INPUT.csv"
#inputdata = pd.read_csv(inputfile,header=None)

alphabet_list = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'   # 220404

## import library
sys.path.append(path1)
sys.path.append(path2)
sys.path.append(os.path.join(path1, "crystal_structure"))
sys.path.append(os.path.join(path1, "soap"))

from crystal_morphing import CrystalInterpolate
from crystal_structure import CrystalStructure
from soap import Soap
from spectra_gsas import Spectra
import spectra_gsas

## 220921 all files check ##

all_filename = [] ; all_cossim = []


#--- input info loading ---#    220309

inputdataraw = open(inputfile,'r').readlines()

inputdata = []
for i in range(len(inputdataraw)):
    l = inputdataraw[i].replace('\n','').replace(' ','').split(',')
    inputdata.append(l)

inputdata = pd.DataFrame(inputdata)


##--- set SOAP parameters ---##
nbasis = 10
nlmax  = 6
sigma_soap  = 0.5
nrad   = 60
ltype  = 3
center_element_name = "G"
#multi_species = True

multi_species = inputdata[inputdata.iloc[:,0] == 'multi_species'].iloc[0,1]  #False #True

## FAST case ##
# nbasis = 5
# nlmax  = 3
# nrad   = 30

##--- set CrystalInterpolate parameters ---##
#max_steps = 15  # 3

try:
  max_steps = int(inputdata[inputdata.iloc[:,0] == 'soap_max_steps'].iloc[0,1]) #20  # 3
except IndexError:max_steps = 20
sub_steps = 3
SD_d = 0.2
SD_i = 16 # int(int(max_steps)*0.80)  #15 # 30
delta_t = 2.*10**(-2)
delta_r = 2.*10**(-2)
tol_d = 0.1    #0.01
dis_type = "ps"
mode = "all"

##--- (optional) set scaling parameters ---##
scaling = False # set True if you want to fix a volume of output structures.
if scaling==True:
    chi_fix = CrystalStructure.from_file("path/to/structure") # a volume of an output structure will be changed to this structure.
sm = "volume" # if you want to do 2D scaling, set sm = "area".

def ci(struct_path1, struct_path2, ratio):
    print ('at ci ratio --- : ', ratio,'----------') # 220309
    print ('---------------------------------') # 220721
    chi1_real = CrystalStructure.from_file(struct_path1)
    chi2_real = CrystalStructure.from_file(struct_path2)
    if scaling==True:
        chi1_real = chi1_real.scale_chi(chi_fix, scaling_method=sm)
        chi2_real = chi2_real.scale_chi(chi_fix, scaling_method=sm)
    struct_path1 = struct_path1.replace('struct/','')
    struct_path2 = struct_path2.replace('struct/','')
    m1 = re.match(r'([a-zA-Z0-9_]+)\.cif', struct_path1.replace(dir_name+'/',''))
    m2 = re.match(r'([a-zA-Z0-9_]+)\.cif', struct_path2.replace(dir_name+'/',''))
    cif_name = dir_name+"/"+m2.groups()[0]+"_to_"+m1.groups()[0]+"_"+str(int(ratio*100)).zfill(4)+"_"+str(bo_i)+".cif"

    all_filename.append(cif_name)   ### 220921 added 

    #--- transform to the reciprocal space ---
    sigma_rtrans = sigma_soap

    if multi_species == 'True':
      chi1  = chi1_real.rtransform(sigma_rtrans, multi_species=True)   # multi species
      chi2  = chi2_real.rtransform(sigma_rtrans, multi_species=True)   # multi species

    else:
      chi1  = chi1_real.rtransform(sigma_rtrans, multi_species=False)   # multi species
      chi2  = chi2_real.rtransform(sigma_rtrans, multi_species=False)   # multi species


    #--- create Soap object ---
    rcut = chi1.rcut
    soap = Soap(nbasis, nlmax, rcut, sigma_soap, nrad, ltype, center_element_name)

    #--- create CrystalInterpolate object ---
    ci = CrystalInterpolate(soap, delta_t, delta_r, dis_type, mode, sub_steps, chi1, verbose=False)
    d_soap, d_ps = ci.report_metrics(chi1, chi2)
    d_ini = d_ps if dis_type=="ps" else d_soap
    d_x = d_ps if dis_type=="ps" else d_soap
    d_des = d_ini*(1.0 - ratio/100.0)

    if d_des==0.0:
        _chi1_real = CrystalStructure(chi1.chis[0].lattice.reciprocal_lattice,
                                      chi1_real.species, chi1.rjs, coords_are_cartesian=True)
        if scaling==True:
            _chi1_real = _chi1_real.scale_chi(chi_fix, scaling_method=sm)
        _chi1_real.write_CIF(cif_name)
        return _chi1_real

    #--- main loop ---
    #------------------------------------------------------------------#
    for _itenum in range(max_steps):
        if abs(d_des-d_ini) < 0.5*tol_d*d_des:
            break
        itenum = _itenum + 1
        do_SD = True if (d_ps>SD_d or d_des>SD_d) and itenum<SD_i else False
        chi2 = ci.optimize(itenum, chi1, chi2, d_des, do_SD)
        d_soap, d_ps = ci.report_metrics(chi1, chi2)

        # no improvement, then stop. (assuming dis_type=="ps")
        if abs(d_ps-d_x) < 1.e-5:
            break

        d_x = d_ps if dis_type=="ps" else d_soap
        if abs(d_des-d_x) < tol_d*d_des:
            break
    #------------------------------------------------------------------#
    print("\n <-- End of the interations -->")

    _chi2_real = CrystalStructure(chi2.chis[0].lattice.reciprocal_lattice,
                                  chi2_real.species, chi2.rjs, coords_are_cartesian=True)
    if scaling==True:
        _chi2_real = _chi2_real.scale_chi(chi_fix, scaling_method=sm)
    _chi2_real.write_CIF(cif_name)

    return _chi2_real



### bring from uspex code 220929

try:
    structure_selection_pair = inputdata[inputdata.iloc[:,0] == 'structure_selection_pair'].iloc[0,1]
except IndexError:
    structure_selection_pair='False'

try:
    xray_wavelength = inputdata[inputdata.iloc[:,0] == 'xray_wavelength'].iloc[0,1]
except IndexError:
    xray_wavelength='CuKa'

try:
    xrd_sigma=float(inputdata[inputdata.iloc[:,0] == 'xrd_sigma'].iloc[0,1])
except IndexError:
    xrd_sigma=0.5
try:
  bool_gsas=str(inputdata[inputdata.iloc[:,0] == 'bool_gsas'].iloc[0,1])
except IndexError:
  bool_gsas = 'False'
try:
   volmatrix = np.array(inputdata[inputdata.iloc[:,0] == 'volume_list'].iloc[0,1:].dropna(),dtype=float).round(2)
except IndexError:
   volmatrix = np.array([0.85,0.90,0.95,0.97,1.00,1.03,1.05,1.10,1.15]).round(2)


#two_theta_range=(0,120)
#bool_gsas=True
bool_gsas=inputdata[inputdata.iloc[:,0] == 'bool_gsas'].iloc[0,1]    #True #False #True
#sigma=0.5

##--- set target structure ---##   220309

obj_struct_path = inputdata[inputdata.iloc[:,0] == 'target'].iloc[0,1]
#obj_struct_path = 'struct/MgO.cif'
#obj_chi = CrystalStructure.from_file(obj_struct_path)

#spectra_obj = Spectra(obj_chi, two_theta_range, d_two_theta, xray_wavelength, method, use_gsas=bool_gsas)
#spectra_obj.gauss_filter(sigma)
#df_obj = spectra_obj.to_df()
#df_obj.plot(x='2_theta', y='Intensity')
# plt.savefig('xrd_Ladybug.png')

chi1 = pd.read_csv(obj_struct_path)

two_theta_min = min(chi1.iloc[:,0]) ; two_theta_max = max(chi1.iloc[:,0]) # modified 210802
two_theta_range=(two_theta_min,two_theta_max)     # arbitrarily modified 210802
two_thetas = np.array(chi1.iloc[:,0])    #  modified 210802
sigma = xrd_sigma   # 0.5     # 220913 changed
d_two_theta = chi1.iloc[1,0] - chi1.iloc[0,0]       #0.1 ;         # modified 210802
method = "XRD"

spectra1 = chi1.iloc[:,1]     #Spectra(chi1, two_theta_range, d_two_theta, method)
spectra1_smeared = spectra_gsas.gauss_filter2(spectra1,d_two_theta, sigma)# .gauss_filter(sigma)



def get_cs(chi):
    # print ('at-get-cs',chi)
    # spectra_obj = Spectra(obj_chi, two_theta_range, d_two_theta, method, use_gsas=True)
    # spectra_obj.gauss_filter(sigma)
#    spectra = Spectra(chi, two_theta_range, d_two_theta, method)
#    spectra.gauss_filter(sigma)
#    cs = spectra.cos_similarity(spectra_obj)   original
    original_vol = chi.lattice.volume ; cs_matrix = []
#    vol_diff = np.arange(0.84, 1.16, 0.02)
    vol_diff = volmatrix

    for vol in vol_diff:
       chi.lattice = chi.lattice.scale(original_vol*vol)
       if bool_gsas == 'True':
          spectra = Spectra(chi, two_theta_range, d_two_theta, xray_wavelength, method, use_gsas=bool_gsas)
       else:
          spectra = Spectra(chi, two_theta_range, d_two_theta, xray_wavelength, method)
       spectra.gauss_filter(sigma)
#       cs = spectra.cos_similarity(spectra_obj)
       cs = spectra_gsas.cos_similarity2(spectra1_smeared,spectra.f_result)
       cs_matrix.append(cs)
       print ('volume, volratio, cs : ',chi.lattice.volume,vol,cs)
#    return cs
    print ('---cs-matrix-max-and-vol---',np.max(cs_matrix),vol_diff[np.argmax(cs_matrix)])
    return np.max(cs_matrix),vol_diff[np.argmax(cs_matrix)]

def fun(x):
    chi = ci(to_chi, from_chi, x[0,0])
    cs, arb = get_cs(chi)
    all_cossim.append(cs)

    return cs.ravel()

def get_chi2():
     cs_list=[]
     for s in struct_path:
         if s in prev_chi:    # 220309
             cs = 0.0         # 220309
         else:   # 220309
             chi = CrystalStructure.from_file(s)
             if scaling==True:
                 chi = chi.scale_chi(chi_fix, scaling_method=sm)
             cs, arb = get_cs(chi)
         cs_list.append(cs)     ## space inserted  220309
     print('-------------------------------------')
     print('-----cs_list of input structures : ',cs_list)
     print('-------------------------------------')
     sort_id = np.argsort(np.array(cs_list))
     from_chi = struct_path[sort_id[-1]] #一番似てる
     to_chi   = struct_path[sort_id[-2]] #二番目に似てる

     print("from", struct_path[sort_id[-1]])
     print("to", struct_path[sort_id[-2]])

     return to_chi, from_chi

## 220718
def get_chi():
    from_chi, to_chi = c_struct_path[bo_i]

    print("from", from_chi)
    print("to", to_chi)

    return to_chi, from_chi


def struct_path_append(x, i):
#    o = open('materlist.csv','a')
    from pymatgen.io import cif
    path = dir_name+'/*.cif'
    flist = glob.glob(path)
    for file in flist:
        if '_' + str(int(x*100)).zfill(4) + '_' + str(i) + '.cif' in file:
            print ('bestchi-filename-----',file)
            bestchi = cif.CifParser(file).get_structures(primitive=False)[0]
            matername = file.split('/')[-1].split('_')[0]
            firstname = file.split('.cif')[0]
            print ('----matername---',matername)
            print ('----firstname---',firstname)
#            o = open('output/'+matername+'_search/all_cossim_output.csv','a')
#            print ('bestchi----')
#            print (bestchi)
            original_vol = bestchi.lattice.volume
            print ('----original_vol---',original_vol)
            cossim_val, vol = get_cs(bestchi)
            print ('vol : ',vol,' cossim_val : ',cossim_val)
            bestchi.lattice = bestchi.lattice.scale(original_vol*float(vol))
            print ('----original_vol*vol---',original_vol*float(vol))
            print ('bestchi-vol---')
            print (bestchi)

#            alphabetnum = alphabet_list[i]
#            o.write(str(alphabetnum)+','+firstname+'\n')
#            o.write(file+','+str(cossim_val)+'\n')

#            file2 = 'output/'+matername+'_search/'+matername+'_'+alphabetnum+'_'+str(i)+'_vol.cif'
            file2 = firstname+'_vol.cif' #'output/'+matername+'_search/'+matername+'_'+alphabetnum+'_'+str(i)+'_vol.cif'
            w = cif.CifWriter(bestchi)
            w.write_file(file2)
#            o.write(file2+','+str(cossim_val)+'\n')

            print ('before append---',struct_path)
#            struct_path.append(file)
            struct_path.append(file2)
            print ('after append---',struct_path)




def struct_path_append2(x, i):
#    o = open('materlist.csv','a')
    from pymatgen.io import cif
    path = dir_name+'/*.cif'
    flist = glob.glob(path)
    for file in flist:
        if '_' + str(int(x*100)).zfill(4) + '_' + str(i) + '.cif' in file:
            print ('bestchi-filename-----',file)
            bestchi = cif.CifParser(file).get_structures(primitive=False)[0]
            matername = file.split('/')[-1].split('_')[0]
            firstname = file.split('.cif')[0]
            print ('----matername---',matername)
            print ('----firstname---',firstname)
            o = open('output/'+matername+'_search/output_list.csv','a')
            print ('bestchi----')
            print (bestchi)
            original_vol = bestchi.lattice.volume
            print ('----original_vol---',original_vol)
            cossim_val, vol = get_cs(bestchi)
            print ('vol : ',vol,' cossim_val :',cossim_val)
            bestchi.lattice = bestchi.lattice.scale(original_vol*float(vol))
            print ('----original_vol*vol---',original_vol*float(vol))
            print ('bestchi-vol---')
            print (bestchi)

            alphabetnum = alphabet_list[i]
            o.write(str(alphabetnum)+','+firstname+'\n')
#            o.write(file+','+str(cossim_val)+'\n')

            file2 = 'output/'+matername+'_search/'+matername+'_'+alphabetnum+'_'+str(i)+'_vol.cif'
#            file2 = firstname+'_vol.cif' #'output/'+matername+'_search/'+matername+'_'+alphabetnum+'_'+str(i)+'_vol.cif'
            w = cif.CifWriter(bestchi)
            w.write_file(file2)
#            o.write(file2+','+str(cossim_val)+'\n')

#            print ('before append---',struct_path)
#            struct_path.append(file)
            struct_path.append(file2)
#            print ('after append---',struct_path)




##--- set input structure ---##   220309
#%%capture cap --no-stderr
#%%time



if structure_selection_pair == 'Limited':
  input_struct_list1 =  inputdata[inputdata.iloc[:,0] == 'input_struct1'].iloc[0,1:].dropna()
  input_struct_list2 =  inputdata[inputdata.iloc[:,0] == 'input_struct2'].iloc[0,1:].dropna()
  struct_path1 = []; struct_path2 = []
  for k in range(len(input_struct_list1)):
    if input_struct_list1.iloc[k]==None:
        continue
    struct_path1.append(input_struct_list1.iloc[k])
  for k in range(len(input_struct_list2)):
      if input_struct_list2.iloc[k]==None:
          continue
      struct_path2.append(input_struct_list2.iloc[k])
  c_struct_path = [] ; struct_path = struct_path1 + struct_path2
  for st1 in input_struct_list1:
      for st2 in input_struct_list2:
          c_struct_path.append((st1,st2))
          c_struct_path.append((st2,st1))


else:
  input_struct_list = inputdata[inputdata.iloc[:,0] == 'input_struct'].iloc[0,1:].dropna()
  print ('input_struct_list------',input_struct_list) # = inputdata[inputdata.iloc[:,0] == 'input_struct'].iloc[0,1:]
  struct_path = []
  for j in range(len(input_struct_list)):
    if input_struct_list.iloc[j]==None:
        continue
    struct_path.append(input_struct_list.iloc[j])

## 220718
  if structure_selection_pair == 'True':
#  c_struct_path = list(itertools.combinations(struct_path, 2))
    c_struct_path = list(itertools.permutations(struct_path, 2))

    print ('c_struct_path is for pair-----',c_struct_path)
  else:
    c_struct_path = struct_path
    print ('c_struct_path is for optimization----',c_struct_path)





#input_struct_list = inputdata[inputdata.iloc[:,0] == 'input_struct'].iloc[0,1]
#struct_path.append(inputdata[inputdata.iloc[:,0] == 'input_struct'].iloc[0,1]) #'struct/MgO_3.cif')
#struct_path.append('struct/MgO_3.cif')
#struct_path.append(inputdata[inputdata.iloc[:,0] == 'struc2'].iloc[0,1]) #'struct/MgO_8.cif')  # 220309---
#struct_path.append(inputdata[inputdata.iloc[:,0] == 'struc3'].iloc[0,1]) #'struct/MgO_8.cif')  # 220309---
#struct_path.append(inputdata[inputdata.iloc[:,0] == 'struc4'].iloc[0,1]) #'struct/MgO_8.cif')  # 220309---
#struct_path.append('struct/MgO_8.cif')


##--- make output directory ---##
obj_struct_path = obj_struct_path.replace('struct/','')
#m = re.match(r'([a-zA-Z0-9_]+)\.cif', obj_struct_path)
m = re.match(r'([a-zA-Z0-9_]+)\.csv', obj_struct_path)   ### 220929 cif->csv ver
dir_name = 'output/'+m.groups()[0]+'_search'
os.makedirs(dir_name, exist_ok=True)
os.makedirs(dir_name+'/found_structure', exist_ok=True)

##--- set Bayesian optimization parameters ---##
try: 
   num_ini_data = int(inputdata[inputdata.iloc[:,0] == 'num_init_data'].iloc[0,1]) #4 # number of initial random data samples
except IndexError:num_ini_data = 4
try:
    num_core = int(inputdata[inputdata.iloc[:,0] == 'num_core'].iloc[0,1]) #4 # if >1, do parallel
except IndexError:num_core = 4
try:
    num_iter = int(inputdata[inputdata.iloc[:,0] == 'num_iter'].iloc[0,1]) #4 # max iteration number
except IndexError:num_iter = 4
try:
    initial_step_boolean = str(inputdata[inputdata.iloc[:,0] == 'initial_step_boolean'].iloc[0,1]) #4 # max iteration number
except IndexError:initial_step_boolean = 'False'
try:
    random_seed = int(inputdata[inputdata.iloc[:,0] == 'random_seed'].iloc[0,1])
except IndexError:random_seed = 99


##--- Main ---##
## Searching crystal structure with corresponding XRD spectrum by using bayesian optimization
prev_chi=[]

if structure_selection_pair == 'True' or structure_selection_pair == 'Limited':
    n_cycles = len(c_struct_path)
else:
    n_cycles = len(c_struct_path)-1

o2 = open(dir_name+'/all_cossim_output.csv','a')

for bo_i in range(n_cycles):  ## 220718
    print ('--- total # of cycles : ',n_cycles,'---')
    print ('--- selection cycle number: ', bo_i,'---')

    if structure_selection_pair == 'True' or structure_selection_pair == 'Limited':
      to_chi, from_chi = get_chi()
      prev_chi += [to_chi, from_chi]
    else:
      to_chi, from_chi = get_chi2()
      prev_chi += [to_chi, from_chi]

    from numpy.random import seed

    seed(random_seed)
    bounds = [{'name': 'sigma', 'type': 'continuous', 'domain': (0,100)}]

    print ('---start bo loop ---',bo_i)

    if initial_step_boolean == 'True' or initial_step_boolean == 'T':

       bo = GPyOpt.methods.bayesian_optimization.BayesianOptimization(f=fun, domain=bounds, X=np.array([[0],[25],[50],[75],[100]]), maximize=True, acquisition_type='LCB',
                                                                   initial_design_numdata=num_ini_data, evaluator_type='local_penalization',
                                                                   batch_size=num_core, de_duplication=True)      ### 220309 insert input structures

    else:
       bo = GPyOpt.methods.bayesian_optimization.BayesianOptimization(f=fun, domain=bounds, maximize=True, acquisition_type='LCB',
               initial_design_numdata=num_ini_data, evaluator_type='local_penalization',batch_size=num_core, de_duplication=True)      ### 220309 insert input structures

    bo.run_optimization(max_iter=num_iter)
#    bo.plot_acquisition(filename=dir_name+'/acq_'+str(bo_i)+'.png')   # if patch is applied, turn on
#    bo.plot_convergence(filename=dir_name+'/con_'+str(bo_i)+'.png')   # if patch is applied, turn on

    x_opt = bo.x_opt[0]
    print ('---x_opt : ',x_opt,'---')
    if structure_selection_pair == 'True' or structure_selection_pair == 'Limited':
       struct_path_append(x_opt, bo_i)
    else:
       struct_path_append2(x_opt, bo_i)

    tmp_chi = CrystalStructure.from_file(struct_path[-1])
    print ('struct_path[-1]---',struct_path[-1])
    shutil.copy(struct_path[-1], dir_name+'/found_structure/')
    final_cossim, final_vol = get_cs(tmp_chi)    # 220309 added
    print('Found : cos_similarity =', final_cossim, final_vol, 'at ratio = ', bo.x_opt)
    print('--------')
    o2.write(str(struct_path[-1])+','+str(final_cossim)+'\n')
print('struct_path')
print(struct_path)
print('\n===== Ending :', dir_name, '=====\n')


## print all raw cossim data ###

#print (all_filename, len(all_filename))
#print (all_cossim, len(all_cossim))
#with open('output/output.txt', 'w') as f:
#    f.write(cap.stdout)
#all_data = pd.DataFrame({'all_filename':all_filename,'all_cossim':all_cossim})
#print (all_data)
#all_data.to_csv('./all_cossim_data.csv',sep=',',index=None)
#all_data.to_csv('output/'+matername+'_search/all_data.csv',sep=',',index=None)

#cap.show()
