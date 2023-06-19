import os
import sys
import numpy as np
import pandas as pd

#############################
### Written by Seiji Kajita ###
### Modified for USPEX XRD application by Joohwi Lee ###
########################################

class Spectra():
    def __init__(self,chi_real, two_theta_range, d_two_theta, method="XRD"):
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
        from pymatgen.analysis.diffraction.xrd import XRDCalculator
        from pymatgen.analysis.diffraction.neutron import NDCalculator

        if method == "XRD":
            xrd = XRDCalculator(wavelength="CuKa", symprec=0, debye_waller_factors=None)
            spectrum = xrd.get_pattern(chi_real, scaled=True, two_theta_range=two_theta_range)
        elif method == "neutronD":
            nd  = NDCalculator(wavelength=1.54184, symprec=0, debye_waller_factors=None) # 1.54184 corresponds to CuKa
            spectrum = nd.get_pattern(chi_real, scaled=True, two_theta_range=two_theta_range)
            pass
        else:
            print("Error. undefined method parameter:", method)
            exit()

        two_thetas = np.arange( two_theta_range[0], two_theta_range[1], d_two_theta)
        f_result      = np.zeros(len(two_thetas))
        for peak_two_theta, intensity in zip(spectrum.x, spectrum.y):
            two_theta, k_theta = self._getNearestValueIndex(two_thetas, peak_two_theta)
            f_result[k_theta]    = intensity
        self.d_two_theta = d_two_theta
        self.two_thetas  = two_thetas
        self.f_result    = f_result
        return

    def gauss_filter(self, sigma=0.5):
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
        idx = np.abs(np.asarray(list) - num).argmin()
        return list[idx], idx

