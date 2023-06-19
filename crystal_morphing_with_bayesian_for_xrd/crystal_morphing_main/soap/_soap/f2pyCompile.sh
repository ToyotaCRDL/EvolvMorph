#!/bin/#!/bin/sh

python -m numpy.f2py -c -lblas rfunction.f90 -m rfunction

## Don't use 
## --opt='-Ofast' 
## because it changes the resutls, though it is fast. So dangerous.

