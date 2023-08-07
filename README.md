# Evolv&Morph
(Last update : August 07, 2023)
- Project researchers in Toyota Central R&D Labs., Inc. : Joohwi Lee, Junpei Oba, Nobuko Ohba, and  Seiji Kajita.
 
Evolv&Morph is a useful inverse design scheme. <br>

This is applied to create a crystal structure reproducing a given XRD pattern. <br>

This program involves:
- Evolutionary algorithm : modified parts of USPEX for XRD application <br>

  See original code of USPEX 9.4.4 developed by Prof. Oganov et al., https://uspex-team.org)

- Crystal morphing (Oba and Kajita, Phys. Rev. Mater. 6, 023801 (2022)) and bayesian optimization for XRD application 

- Other XRD-related tools


## Installation

It is recommended to use ubuntu. (Mac might be fine.) <br>

The result in the paper is obtained using ubuntu18.04, python3.7. <br>

The last test was done (June 19, 2023) using ubuntu18.04 and python3.9 created in docker.

1. Download Evolv&Morph code from github.

2. Update apt-get and install common/required packages. (sudo apt-get could be needed sometimes.) <br>
   Note that some of them are essential, but some of them are not so.

>  apt-get update --yes <br>

>  apt-get install --yes --no-install-recommends wget ca-certificates sudo locales fonts-liberation <br>

>  apt-get install --yes --no-install-recommends file git less nano patch tzdata unzip vim <br>

>  apt-get install --yes --no-install-recommends gnuplot ghostscript epstool epstool fig2dev pstoedit <br>

>  apt-get install --yes --no-install-recommends octave <br>

>  apt-get install --yes --no-install-recommends openssh-client <br>

>  apt-get install --yes --no-install-recommends firefox <br>

3. install required programs via mambaforge or pip (other conda such as miniconda  might be fine.)

> wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh <br>

> bash ./Mambaforge-Linux-x86_64.sh <br>

> mamba install --quiet --yes ipywidgets widgetsnbextension <br>

> mamba install --quiet --yes numpy=1.23.1 scipy pandas matplotlib ipympl zstandard zstd openssl  <br>

> mamba install --quiet --yes ase networkx pymatgen==2020.12.31 sympy   <br>

> mamba install --quiet --yes blas lapack libblas  <br>

> mamba install --quiet --yes gpy==1.10.0 gpyopt==1.2.6 paramz==0.9.5    <br>

> mamba install --quiet --yes gcc gfortran make nglview   <br>

> mamba install --quiet --yes notebook nb_conda_kernels  <br>

> mamba clean --all -f -y    <br>

> pip install rise <br> 

  Note that rise is installed via pip at python3.9. Mamba make error. <br>

4. Install GSAS-II for simulation of XRD from cif file.

  For the detail, see installation pages of GSAS-II. <br>

  https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/LinuxSingleStepInstaller   <br>

  https://gsas-ii.readthedocs.io/en/latest/packages.html   <br>
<br>
> mamba install --quiet --yes wxpython <br>

> mamba install --quiet --yes pyopengl pillow h5py imageio subversion <br>

> svn co https://subversion.xray.aps.anl.gov/pyGSAS/trunk ~/GSASII <br> 

> python ~/GSASII/GSASIIscriptable.py <br>

  Note that wxpython is successfully installed at python3.9, but failed at python3.7. <br>

  After installation, the installed GSASII folder should be added at ~/.bashrc. <br>

> export PYTHONPATH=/home/[user name]/GSASII:$PYTHONPATH


## Installation of evolutionary algorithm (USPEX)

USPEX 9.4.4 with matlab/octave was modified for XRD application. <br>

In addition, some grammars of old octave and python2.7 was modified to work well in the newer version. <br>

See https://uspex-team.org/en/uspex/downloads for original code. <br>

1. According to the original code installation, USPEX 9.4.4 should be installed.

2. Overwrite all the files and folders of modified part of USPEX in this code (./uspex_modified_for_xrd) into the installed USPEX folder.<br>
> cp -a -r ./uspex_modified_for_xrd/* [the folder USPEX is installed] 

   Here, -a option preserve the permission. This keeps -x permision for USPEX and other python files in the subfolders.

3. At ~/.bashrc, the followings should be added. 
> export PATH=[the folder USPEX is installed]:$PATH  <br>

> export USPEXPATH=[the folder USPEX is installed]/src


## Additional installation of crystal morphing

Crystal morphing and its application can be directly performed at the folder. <br>

However, the following procedure is essential prior to using it. <br>

1. Compile SOAP.

> cd ./crystal_morphing_with_bayesian_for_xrd/crystal_morphing_main/soap/_soap <br>

> bash f2pyCompile.sh 1> stdout.log 2> >(tee stderr.log >&2) <br>

> cd ../../../..

2. Patch matplot and GpyOpt installed by mambaforge.

   Do patch for avoiding errors for matplotlib and gpyopt related to version issue. <br>

   It is needed to know where matplotlib and gpyopt were installed by mambaforge (or conda). <br>

> find [mamba-basic-folder]/envs/[mamba-env-name] -name batch_local_penalization.py -exec cp ./patches/batch_local_penalization.py {} \\; <br>

> patch -b \`find [mamba-basic-folder]/envs/[mamba-env-name] -name plot_definitions.py |grep matplot\` ./patches/plot_definitions.py.patch

  For detail, see https://github.com/SheffieldML/GPy/pull/960 <br>

  If the following error occurs, see grammar related to -exec at https://stackoverflow.com/questions/2961673/find-missing-argument-to-exec <br>
```
find: missing argument to `-exec 
```

  

3. If the following error occurs from Bayesian optimization, numpy>1.24 should be downgrade into numpy1.23.1.

```
AttributeError: module 'numpy' has no attribute 'bool'. <br>
`np.bool` was a deprecated alias for the builtin `bool`. ...
```
> mamba uninstall numpy ; mamba install numpy=1.23.1


## Run evolutionary algorithm (./uspex_example_for_xrd)

How to run is written in jupyter notebook, evolutionary_algorithm_vs_target_xrd_csv_without_DFT.ipynb

If DFT calculation (VASP) is used as a supporting structural optimizer, additional setting is necessary.

Detailed information is written in jupyter notebook, evolutionary_algorithm_vs_target_xrd_csv_with_DFT.ipynb.
  

## Run crystal_morphing (./crystal_morphing_with_bayesian_for_xrd/crystal_morphing_main)

Crystal morphing makes intermediate structures between two input structures.

Detailed information is written in jupyter notebook, crystal_morphing.ipynb.

Also see README.md in the folder. <br>


## Run crystal_morphing with Bayesian for XRD application (./crystal_morphing_with_bayesian_for_xrd)

This process tries to find intermediate structures with higher cosine similarity of XRD pattern with respect to the target XRD.

Detailed information is written in jupyter notebook, crystal_morphing_bayesian_vs_xrd_csv.ipynb.


## Utils for XRD (./utils_for_xrd & ./utils_for_xrd/cossim_vs_xrdcsv/)

Some useful scripts are included.

Detailed information is written in jupyter notebook.


## Citation

If you find it useful to use this program in your research, please cite the following paper.

- Evolv&Morph 
```
Joohwi Lee, Junpei Oba, Nobuko Ohba, and Seiji Kajita, npj Computational Materials, 9, 135 (2023). 
```
When accepted, the reference will be changed into an accepted journal.

- Crystal Morphing
```
Junpei Oba and Seiji Kajita, Physical Review Materials 6, 023801 (2022).
```

## NOTICE

Copyright (C) 2023 TOYOTA CENTRAL R&D LABS., INC. All Rights Reserved.

For non-commercial research purposes only, under this license, permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), the rights to use, copy, modify, merge of the Software, subject to the following conditions:

The above copyright notice and this license shall be included in all copies or substantial portions of the Software.

Except as expressly stated above, no rights or licenses from any copyright holder is granted under this license, whether expressly, by implication, estoppel or otherwise.
Under this license, the rights to distribute the Software, including a derivative work of the Software, are NOT granted.

The name and trademarks of copyright holder(s) may NOT be used in advertising or publicity pertaining to the Software or portions thereof, including modifications or derivatives, without specific, written prior permission. Title to copyright in the Software will at all times remain with the copyright holders.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.






