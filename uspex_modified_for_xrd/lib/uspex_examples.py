#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
USPEX 9.4.4 release
2015 Oganov's Lab. All rights reserved.
"""

import os
import sys
import distutils.core


# Process directories
run_dir = os.getcwd()
bin_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
os.chdir(bin_dir)
os.chdir('examples')
examples_dir = os.getcwd()
readme_file = examples_dir + '/README.txt'
os.chdir(run_dir)


# -------------------------------------------------------------------------------
def uspex_examples(example):
    if not example:
        return

    f = open(readme_file, 'rb')
    readme_content = f.readlines()
    f.close()

    rows_numbers = []
    for i, row in enumerate(readme_content):
#        if row.find('LIST OF EXAMPLES:') >= 0:
        if row.find(b'LIST OF EXAMPLES:') >= 0:       # modified on 210422
            for j in range(i, len(readme_content)):
#                if readme_content[j].find('EX') == 0:
                if readme_content[j].find(b'EX') == 0:   # modified on 210422
                    rows_numbers.append(j)
            break

    examples_list = []
    for i in range(len(rows_numbers)):
        examples_list.append(readme_content[rows_numbers[i]])
    num_examples = len(examples_list)

    example_name = None
#    if example <> 'all':
    if example != 'all':     # modified on 210414
        try:
            example = int(example)
            example -= 1
            if example >= 0 and example < num_examples:
                examples_list = [examples_list[example]]
#                example_name = examples_list[0].split(':')[0]
                example_name = examples_list[0].split(b':')[0].decode()   # modified on 210422
            else:
                examples_list = []
                print (('Example does not exist! Possible values: %i - %i.') % (1, num_examples))   # modified on 210414
        except ValueError:
            examples_list = []
            print ('Example value should be integer!')   # modified on 210414

    print ('')   # modified on 210414
    for row in examples_list:
#        print (''.join(row))
        print (''.join(row.decode()))   # modified on 210414

    return example_name


# -------------------------------------------------------------------------------

uspex_var = 'USPEXPATH'
uspexpath = os.environ[uspex_var]


# -------------------------------------------------------------------------------
def uspex_copy(example):
    if not example:
        return

    example2copy = uspex_examples(example)
    if example2copy:
        from_dir = examples_dir + '/' + example2copy
        if os.path.exists(from_dir):
            distutils.dir_util.copy_tree(from_dir, run_dir)
            print (example2copy + ' example copied to the current directory.')   # modified on 210414


def copy_submission():
    from_dir = uspexpath + '/Submission'
    to_dir = run_dir + '/Submission'
    if os.path.exists(from_dir) and not os.path.exists(to_dir):
        distutils.dir_util.copy_tree(from_dir, to_dir)
        print ('Submission dir copied to the current directory.')   # modified on 210414

# -------------------------------------------------------------------------------
