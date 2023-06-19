# Crystal morphing

(Last update: June 19th, 2023)  
- Project researchers in Toyota Central R&D Lab. Inc.: Junpei Oba and Seiji Kajita.
<br>
Interpolating two different crystal structures. <br> 
SOAP distance and its derivative are used to update lattice vectors and atom positions.


## Required work
Detailed installation is written in Evolv&Morph. <br>

Note that crystal morphing requires compile of soap as following.

> cd ./soap/_soap <br>
> bash f2pyCompile.sh 1> stdout.log 2> >(tee stderr.log >&2) <br>
> cd ../..

## Run crystal morphing

Crystal morphing makes intermediate structures between two input structures.

Detailed information is also written in jupyter notebook, crystal_morphing.ipynb

#### 1. Check INPUT.csv and crystal_morphing.py
You may change the following parameters in INPUT.csv and crystal_morphing.py 

- **struct1**: path to the objective structure's file

- **struct2**: path to the initial structure's file

- **soap_max_steps**: number of iterarion

- **multi_species**: True for multi species case (e.g. MgO), or False for single species case (e.g. C)

- **sub_steps**: number of step to update lattice vectors and atom positions in each iterarion

- **delta_t**: step size for updating lattice vectors

- **delta_r**: step size for updating atom positions

- **sw_SD**: parameter to switch stochastic descent and BFGS methods to update structure parameters

- **verbose**: if True, output detailed information while calculation

#### 2. Execute crystal_morphing.py
> python crystal_morphing.py |tee log

#### 3. Check result
Intermediate structures will be outputted in the "examination/tmp/" directory by default.
Output file's name is specified in crystal_interpolation.py as follows.
> cif_name="examination/tmp/test_"+str(itenum)+".cif"

<br>



## Citation

If you find it useful to use this program in your research, please cite the following paper.

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



