import os
import sys
import numpy as np
import pandas as pd
from pymatgen.io import cif
from pymatgen.io.vasp import Poscar
#############################
### Written by Seiji Kajita ###
### Modified for USPEX by Joohwi Lee ###
########################################

import GSASIIscriptable as G2sc


### input info loading #  220912 added
inputfile = 'additional_INPUT.csv'
inputdataraw = open(inputfile,'r').readlines()

inputdata = []
for i in range(len(inputdataraw)):
    l = inputdataraw[i].replace('\n','').replace(' ','').split(',')
    inputdata.append(l)

inputdata = pd.DataFrame(inputdata)
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

### 220914 input_xy_gsas ###

if bool_gsas == 'True':
    import os.path
    import shutil
    if os.path.exists('../INST_XRY.PRM'):
        shutil.copy('../INST_XRY.PRM','./INST_XRY.PRM')

# 210419 #

class Spectra():
    '''
    pymatgenのStructureを入れるとXRDのsimulation結果がかえってくる処理群。
    二つの構造類似度計算もできる。
    '''
#    def __init__(self,chi_real, two_theta_range, d_two_theta, method="XRD"):     ### 220914 modified
    def __init__(self,chi_real, two_thetas, method="XRD", use_gsas=False):     ### 2200914 modified
        two_theta_range = (min(two_thetas),max(two_thetas)) ; d_two_theta = two_thetas[1] - two_thetas[0]      # 220914 modified
#        two_theta_range = (min(two_thetas),max(two_thetas)) #  220914 modified; d_two_theta = two_thetas[1] - two_thetas[0]      # 220914 modified

        '''
        Parameters:
        ----------
        chi_real: pymatgen Structure object
        two_theta_range: 1d list
            Range of 2 theta diffraction angles.
        d_two_theta: float
            Width of the two_theta_range
        method: string
            "XRD" X-ray diffraction
            "neutronD" neutron diffraction
        '''
        if use_gsas == False:
          from pymatgen.analysis.diffraction.xrd import XRDCalculator
          from pymatgen.analysis.diffraction.neutron import NDCalculator
#         from pymatgen.analysis.diffraction.tem import TEMCalculator
          if method == "XRD":
              try:
                  xrd = XRDCalculator(wavelength=xray_wavelength, symprec=0, debye_waller_factors=None)
              except KeyError:
                  xrd = XRDCalculator(wavelength=float(xray_wavelength), symprec=0, debye_waller_factors=None)
              spectrum = xrd.get_pattern(chi_real, scaled=True, two_theta_range=two_theta_range)
          elif method == "neutronD":
              nd  = NDCalculator(wavelength=1.54184, symprec=0, debye_waller_factors=None) # 1.54184 corresponds to CuKa
              spectrum = nd.get_pattern(chi_real, scaled=True, two_theta_range=two_theta_range)
              pass
          # elif method == "TEM":
          #     tem  = TEMCalculator()
          #     spectrum = tem.get_pattern(chi_real, scaled=True, two_theta_range=two_theta_range)
          else:
              print("Error. undefined method parameter:", method)
              exit()

          '''
          d_two_theta度刻みの2thetaとintensityデータを作成
          '''
#          two_thetas = np.arange( two_theta_range[0], two_theta_range[1], d_two_theta)
          f_result      = np.zeros(len(two_thetas))
          for peak_two_theta, intensity in zip(spectrum.x, spectrum.y):
              two_theta, k_theta = self._getNearestValueIndex(two_thetas, peak_two_theta)
              f_result[k_theta]    = intensity
          ##
          self.d_two_theta = d_two_theta
          self.two_thetas  = two_thetas
          self.f_result    = f_result
          return

        elif use_gsas == True:
            ### 220602 added ###
            ## Use GSAS-II to simulate XRD pattern
            TMP_DIR = './'  ## temporary folder to store sim.gpx, sim.cif, etc.
            DATA_DIR = './' #gsas_folder  ## path where 'INST_XRY.PRM' file exists.
            gpx = G2sc.G2Project(newgpx=TMP_DIR+'sim.gpx')

            w = cif.CifWriter(chi_real)
            w.write_file(TMP_DIR+'sim.cif')
            TARGET_CIF = TMP_DIR + 'sim.cif'
            FORMULA    = chi_real.formula.replace(" ", "")
            phase0 = gpx.add_phase(os.path.join(TARGET_CIF),phasename=FORMULA,fmthint='CIF')

            two_theta_range = (min(two_thetas),max(two_thetas))
            d_two_theta = two_thetas[1] - two_thetas[0]
            hist1 = gpx.add_simulated_powder_histogram("simulation",
                                      os.path.join(DATA_DIR,'INST_XRY.PRM'),
                                      min(two_thetas),max(two_thetas),d_two_theta,
                                      phases=gpx.phases())
            hist1.SampleParameters['Scale'][0] = 1000.
            gpx.data['Controls']['data']['max cyc'] = 0 # refinement not needed
            gpx.do_refinements([{}])
            gpx.save()
            x = gpx.histogram(0).getdata('x')
            y = gpx.histogram(0).getdata('ycalc')
            ##
            self.d_two_theta = d_two_theta
            self.two_thetas  = x
            self.f_result    = y
            os.remove(TMP_DIR+'sim.cif')
        return





    def gauss_filter(self, sigma=0.5):
        '''
        gaussian畳み込みでピークを太らせる
        '''
        def gaussian(theta, sigma):
            gauss = np.exp(-(theta**2)/(2*sigma**2))
            return gauss

        g = gaussian(np.arange(-40*self.d_two_theta, 40*self.d_two_theta, self.d_two_theta), sigma )
        self.f_result = np.convolve(self.f_result, g, 'same')
        return

    def log10_operate(self, delta=1.0):
        self.f_result = np.log10(self.f_result+delta)
        return

    def cos_similarity(self, _xrd):
        denom = np.linalg.norm(self.f_result) * np.linalg.norm(_xrd.f_result)
        return np.dot(self.f_result, _xrd.f_result)/denom

    def cos_log10_similarity(self, _xrd, delta=1.0):
        f_result1 = np.log10(self.f_result+delta)
        f_result2 = np.log10(_xrd.f_result+delta)
        denom = np.linalg.norm(f_result1) * np.linalg.norm(f_result2)
        return np.dot(f_result1, f_result2)/denom

    def to_csv(self, csvname):
        df = self.to_df()
        df.to_csv(csvname)
        return

    def to_df(self):
        df = pd.DataFrame(data=[self.two_thetas, self.f_result]).T
        df.columns = ["2_theta","Intensity"]
        return df

    def _getNearestValueIndex(self, list, num):
        """
        概要: リストからある値に最も近い値を返却する関数
        @param list: データ配列
        @param num: 対象値
        @return 対象値に最も近い値
        """
        # リスト要素と対象値の差分を計算し最小値のインデックスを取得
        idx = np.abs(np.asarray(list) - num).argmin()
        return list[idx], idx



### 210801 added ### 220914 added

def gauss_filter2(spectra, d_two_theta, sigma):
        '''
        gaussian畳み込みでピークを太らせる
        '''
        def gaussian(theta, sigma):
            gauss = np.exp(-(theta**2)/(2*sigma**2))
            return gauss

        g = gaussian(np.arange(-40*d_two_theta, 40*d_two_theta, d_two_theta), sigma)
        spectra2 = np.convolve(spectra, g, 'same')
        return spectra2

def cos_similarity2(spectra1, spectra2):
        denom = np.linalg.norm(spectra1) * np.linalg.norm(spectra2)
        return np.dot(spectra1, spectra2)/denom


def getNearstValueIndex2(list,num):
        idx = np.abs(np.asarray(list) - num).argmin()
        return list[idx], idx


if __name__ == "__main__":
    from pymatgen.io.vasp import Poscar
    from pymatgen.io import cif

    chi1 = pd.read_csv(sys.argv[1],header=None)

    two_theta_min = min(chi1.iloc[:,0]) ; two_theta_max = max(chi1.iloc[:,0]) # modified 210802
    two_theta_range=(two_theta_min,two_theta_max)     # arbitrarily modified 210802
    two_thetas = np.array(chi1.iloc[:,0])    #  modified 210802
    sigma = xrd_sigma   # 0.5     # 220913 changed
    d_two_theta = chi1.iloc[1,0] - chi1.iloc[0,0]       #0.1 ;         # modified 210802
    method = "XRD"

    spectra1 = chi1.iloc[:,1]     #Spectra(chi1, two_theta_range, d_two_theta, method)
    spectra1_smeared = gauss_filter2(spectra1,d_two_theta, sigma)# .gauss_filter(sigma)

    chi2 = cif.CifParser(sys.argv[2]).get_structures(primitive=False)[0]
#    chi2 = Poscar.from_file(sys.argv[2], read_velocities=False).structure
    if bool_gsas == 'True':
        spectra2 = Spectra(chi2, two_thetas, method='XRD', use_gsas=True)
    else: #if bool_gsas == 'False':
#    spectra2 = Spectra(chi2, two_theta_range, d_two_theta, method)
        spectra2 = Spectra(chi2, two_thetas, method='XRD')
#        print (dir(spectra2))
#    spectra2.gauss_filter(sigma)    
#        print (spectra2.f_result)


    original_vol = chi2.lattice.volume ; cs_matrix = []

    vol_def = float(open('vol_default','r').readlines()[0])

    if original_vol >= 4*vol_def:       ### added 210806
       chi2.lattice = chi2.lattice.scale(2*vol_def)     ### added 210806
       original_vol = chi2.lattice.volume     ### added 210806

    vol_diff = volmatrix          # 220913 changed np.arange(0.75, 1.25, 0.01) #[0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
    for vol in vol_diff:
       chi2.lattice = chi2.lattice.scale(original_vol*vol)

       if bool_gsas == 'True':
#          spectra2 = Spectra(chi2, two_theta_range, method='XRD', use_gsas='True')
          spectra2 = Spectra(chi2, two_thetas, method='XRD', use_gsas=True)
       elif bool_gsas == 'False':
           spectra2 = Spectra(chi2, two_thetas, method='XRD')
       spectra2_smeared = gauss_filter2(spectra2.f_result,d_two_theta,sigma)
       cs = cos_similarity2(spectra1_smeared, spectra2_smeared)
#       print (chi2.lattice.volume,vol,cs)    # for uspex 220916
       cs_matrix.append(cs)

    vol_diff_arb = abs(vol_diff-1)
    originalindex = list(vol_diff_arb).index(np.min(vol_diff_arb)) #list(vol_diff).index(np.min) ##np.min(abs(vol_diff-1)))
    print ('cos-similarity, original volume',cs_matrix[originalindex], vol_diff[originalindex])

    maxindex = cs_matrix.index(np.max(cs_matrix))
    print ('max-cos-similarity, changed volume ratio',cs_matrix[maxindex], vol_diff[maxindex]) #np.max(cs_matrix))
    chi2.lattice = chi2.lattice.scale(original_vol*vol_diff[maxindex])

    w = cif.CifWriter(chi2)   ### 220914 modified
    w.write_file('CONTCAR-volchanged.cif')    ### 220914 modified

    if bool_gsas == 'True':
          spectra2 = Spectra(chi2, two_thetas, method='XRD', use_gsas=True)
    elif bool_gsas == 'False':
           spectra2 = Spectra(chi2, two_thetas, method='XRD')
    spectra2_smeared = gauss_filter2(spectra2.f_result,d_two_theta,sigma)

#    data2 = np.column_stack((spectra2.two_thetas,spectra2.f_result))    # for uspex 220916
#    np.savetxt(sys.argv[2].split('.cif')[0]+'_volchanged_xrd.csv',data2,delimiter=',')   # for uspex 220916
#    data2 = np.column_stack((spectra2.two_thetas,spectra2_smeared))   # for uspex 220916
#    np.savetxt(sys.argv[2].split('.cif')[0]+'_volchanged_smeared_xrd.csv',data2,delimiter=',')   # for uspex 220916

    o = open('exp_cos_similarity_xrd_volchanged.csv','w')
    o.write(str(cs_matrix[maxindex]))
