import sys, os
from pymatgen.io import cif
import numpy as np
import pandas as pd

import GSASIIscriptable as G2sc

### Joohwi Lee ummarized the code.

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
    xrd_2theta_min=float(inputdata[inputdata.iloc[:,0] == 'xrd_2theta_min'].iloc[0,1])
except IndexError:
    xrd_2theta_min=0
try:
    xrd_2theta_max=float(inputdata[inputdata.iloc[:,0] == 'xrd_2theta_max'].iloc[0,1])
except IndexError:
    xrd_2theta_max=180
try:
    d_2theta=float(inputdata[inputdata.iloc[:,0] == 'd_2theta'].iloc[0,1])
except IndexError:
    d_2theta=0.1
try:
    xrd_sigma=float(inputdata[inputdata.iloc[:,0] == 'xrd_sigma'].iloc[0,1])
except IndexError:
    xrd_sigma=0.5
try:
  bool_gsas=str(inputdata[inputdata.iloc[:,0] == 'bool_gsas'].iloc[0,1])
except IndexError:
  bool_gsas = 'False'

filename = sys.argv[1]
xrdoutname = sys.argv[2]
two_thetas = np.arange(xrd_2theta_min,xrd_2theta_max+0.5*d_2theta,d_2theta)
d_two_theta = d_2theta
sigma = xrd_sigma

def gsas_cif_to_xrd(filename):
  TMP_DIR = './'  ## temporary folder to store sim.gpx, sim.cif, etc.
  DATA_DIR = './' #gsas_folder  ## path where 'INST_XRY.PRM' file exists.
  gpx = G2sc.G2Project(newgpx=TMP_DIR+'sim.gpx')

  chi_real =  cif.CifParser(filename).get_structures(primitive=False)[0]
  w = cif.CifWriter(chi_real)
  w.write_file(TMP_DIR+'sim.cif')
  TARGET_CIF = TMP_DIR + 'sim.cif'
  FORMULA    = chi_real.formula.replace(" ", "")

  phase0 = gpx.add_phase(os.path.join(TARGET_CIF),phasename=FORMULA,fmthint='CIF')

  d_two_theta = two_thetas[1] - two_thetas[0]
  hist1 = gpx.add_simulated_powder_histogram("simulation",os.path.join(DATA_DIR,'INST_XRY.PRM'),
                                           min(two_thetas),max(two_thetas),d_two_theta,
                                           phases=gpx.phases())

  hist1.SampleParameters['Scale'][0] = 1000.
  gpx.data['Controls']['data']['max cyc'] = 0 # refinement not needed

  gpx.do_refinements([{}])
  gpx.save()

  x = gpx.histogram(0).getdata('x')
  y = gpx.histogram(0).getdata('ycalc')

  return x,y

def gauss_filter2(spectra, d_two_theta, sigma):
        def gaussian(theta, sigma):
            gauss = np.exp(-(theta**2)/(2*sigma**2))
            return gauss

        g = gaussian(np.arange(-40*d_two_theta, 40*d_two_theta, d_two_theta), sigma)
        spectra2 = np.convolve(spectra, g, 'same')
        return spectra2

def pmg_cif_to_xrd(filename,two_thetas):
   from pymatgen.analysis.diffraction.xrd import XRDCalculator
   chi_real =  cif.CifParser(filename).get_structures(primitive=False)[0]
   try:
                xrd = XRDCalculator(wavelength=xray_wavelength, symprec=0, debye_waller_factors=None)
   except KeyError:
                xrd = XRDCalculator(wavelength=float(xray_wavelength), symprec=0, debye_waller_factors=None)
   d_two_theta = two_thetas[1] - two_thetas[0]
   spectrum = xrd.get_pattern(chi_real, scaled=True, two_theta_range=(xrd_2theta_min,xrd_2theta_max+0.5*d_two_theta))
   two_thetas = np.arange(xrd_2theta_min, xrd_2theta_max+0.5*d_two_theta, d_two_theta)
   f_result      = np.zeros(len(two_thetas))
   def getNearestValueIndex(list, num):
        idx = np.abs(np.asarray(list) - num).argmin()
        return list[idx], idx

   for peak_two_theta, intensity in zip(spectrum.x, spectrum.y):
                two_theta, k_theta = getNearestValueIndex(two_thetas, peak_two_theta)
                f_result[k_theta]    = intensity
   return two_thetas,f_result


if bool_gsas == 'True':
  x,y = gsas_cif_to_xrd(filename)
else:
  x,y = pmg_cif_to_xrd(filename,two_thetas)

y_smeared = gauss_filter2(y, d_two_theta, sigma)
data = np.column_stack((x,y))
print (data)
np.savetxt(xrdoutname,data,delimiter=',')

data2 = np.column_stack((x,y_smeared))
np.savetxt(xrdoutname.split('.csv')[0]+'_smeared.csv',data2,delimiter=',')
