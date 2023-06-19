# Crystal morphing with Bayesian for reproducing XRD

(Last update: June 19th, 2023)  
- Project researchers in Toyota Central R&D Lab. Inc.: Joohwi Lee, Junpei Oba, Nobuko Ohba, and Seiji Kajita. <br>

Create intermediate structures by crystal morphing to increase their XRD similarity of vs. target based on Bayesian optimiation. <br>


## Required work
Crystal morphing should work well. <br>

See README.md and jupyter notebook at ./crystal_morphing_main. <br>

## Run crystal morphing

Basically, how to run is written in jupyter notebook.


#### 1. Check INPUT.csv and morphing_with_bo_vs_xrd_csv.py
You may change the following parameters in INPUT.csv and morphing_with_bo_vs_xrd_csv.py.<br>

The parameter written in morphing_with_bo_vs_xrd_csv.py takes precedence over those in morphing main code. 

#### Prepare INPUT.csv
- **target** : insert path of target csv file including XRD (2theta, XRD intensity)

- **input_struct** : insert path of structures
-  (optional : structure_selection_pair is "Limited" case) <br>
  **input_struct1** & **input_struct2** : insert path of structures seperated into two groups

- **num_init_data** : number of initial random data samples for Bayesian optimization <br>

- **num_core** : number of data for validation for Bayesian optimization, if >1, do parallel

- **num_iter** : number of iteration number of Bayesian optimization

- **initial_step_boolean** : generate initial points for Bayesian (0%, 25%, 50%, 75%, 100%)

- **random_seed** : random seed for initial points (if initial_step_boolean exists, second step)
<br>

- **soap_max_steps** : maximum iteration of soap (default optimizer is steepest descent, if >15, it change to L-BFGS)

- **multi_species** : if chemical elements >=2, set True
<br>

- **structure_selection_pair** : selection of input structures for exploration path <br>
  True : All-pairs investgation, False : greedy algorithm, Limited : Limited-pairs investigation (see input_struct part)
<br>

- **xray_wavelength** : xray_wavelength : XRD source wavelength for pymatgen option

- **volume_list** : insert volume list for cossine similarity, such as 0.90,0.95, 1.00, 1.05, 1.10

- **xrd_sigma** : smearing factor (default : 0.5)
<br>

- **bool_gsas** : XRD is generated from GSAS-II (True). Otherwise, it is generated from pymatgen (False)
<br>

- (Optional) Prepare INST_XRD.PRM (necessary when bool_gsas is True)


#### 2. Execute crystal_morphing.py
> python ./morphing_with_bo_vs_xrd_csv.py |tee log

#### 3. Check result
Intermediate structures will be outputted in the "output/[materialname]_search" directory by default.

See output_list.csv and all_cossim_output.csv


## Citation

If you find it useful to use this program in your research, please cite the following paper.

- Evolv&Morph
```
Joohwi Lee, Junpei Oba, Nobuko Ohba, and Seiji Kajita, submitted., preprint arXiv 2302.10464 (2023).
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



