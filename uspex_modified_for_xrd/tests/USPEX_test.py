#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
USPEX 9.4.4 release
2015 Oganov's Lab. All rights reserved.
"""

import os
import sys
import distutils.core
import subprocess
import time
import datetime
import shutil
import filecmp
from optparse import OptionParser

from lib.parse_OUTPUT import parse_OUTPUT
from lib.parse_Individuals import get_individuals_content
from lib.read_POSCAR import read_POSCAR, count_POSCARS
from lib.valid_tests import *



# -------------------------------- Options --------------------------------------
parser = OptionParser()
parser.add_option("-t", "--test", dest="testid",
                  help="specify a test or a comma-separated sequence of tests to execute. '-t all' executes all tests",
                  metavar="TEST")
parser.add_option("-l", "--list", dest="tests_list", action="store_true",
                  help="show list of all available tests")
parser.add_option("-o", "--octave", dest="uspex_octave", action="store_true",
                  help="run USPEX calculation with Octave instead of Matlab", metavar="RUN")

(options, args) = parser.parse_args()

if not options.testid or options.tests_list:
    print ('Valid test id:')   # modified on 210421
    for key in sorted(valid_tests.keys()):
        print ('%5s\t%-30s' % (key, valid_tests[key]))  # modified on 210421

    sys.exit(1)

# For calc dir name:
init_option = options.testid

try:
    testid_tmp = int(options.testid)
except ValueError:
    testid_tmp = options.testid

if testid_tmp == 'all':
    testid = testid_tmp
elif testid_tmp in valid_tests.keys():
    testid = [testid_tmp]
else:
    if type(testid_tmp) is str:
        testid = []
        for i in testid_tmp.split(','):
            try:
                if int(i) in valid_tests.keys():
                    testid.append(i)
            except ValueError:
                error_str = 'Test with id ' + i + ' does not exist.'
        if testid == []:
            error_str = 'No valid test id found.'
            parser.error(error_str)
    else:
        error_str = 'Test with id ' + options.testid + ' does not exist.'
        parser.error(error_str)

uspex_octave = ''
if options.uspex_octave:
    uspex_octave = '-o'

# parser.destroy()
# -------------------------------------------------------------------------------

start_time = time.time()


# -------------------------------------------------------------------------------
def duration(seconds):
    minutes = 0
    if seconds >= 0 and seconds < 60:
        duration = '%.1f seconds' % seconds
    elif seconds >= 60:
        minutes = int(seconds / 60)
        seconds = seconds - minutes * 60
        units = 'minutes'
        if minutes == 1:
            units = 'minute'
        duration = '%i %s %.1f seconds' % (minutes, units, seconds)
    return duration


# -------------------------------------------------------------------------------
start_time = time.time()


# -------------------------------------------------------------------------------
def run_test(testid_local, uspex_octave=''):
    test_start_time = time.time()
    st = datetime.datetime.fromtimestamp(test_start_time).strftime('%Y-%m-%d %H:%M:%S')

    testid_str = 'T%02i' % (int(testid_local))

    # ---------------------------------------------------------------------------
    # Functions part:

    # Function to check if all the files from the specified list exist:
    def check_file(uspex_file, empty_dirs):
        stat_info = 0
        return_value = 0
        stat = 'FAIL'

        if os.path.exists(uspex_file):
            stat_info = os.stat(uspex_file).st_size
            if stat_info > 0:
                return_value = 1
                stat = 'OK'
            else:
                if uspex_file in empty_dirs:
                    return_value = 1
                    stat = 'OK'

        print_str = '\t%-8s %-50s %8i' % (stat, uspex_file, stat_info)
        log.write(print_str + '\n')
        print (print_str)   # modified on 210421

        return return_value

    # Function for comparison of files in the resultS1 folder and in the reference folder:
    def ref_compare(file1, file2):
        if os.path.exists(file1) and os.path.exists(file2):
            result = filecmp.cmp(file1, file2)
        else:
            result = False
        return result

    # Copy input files/folders:
    def copy_folder(from_dir, to_dir):
        if os.path.exists(from_dir) and not os.path.exists(to_dir):
            distutils.dir_util.copy_tree(from_dir, to_dir)
            print_str = 'Specific dir and INPUT.txt file copied to the current directory.'
            log.write(print_str + '\n')
            print (print_str)   # modified on 210421

    # ----------------------------------------------------------------------------

    logfile = testid_str + '_' + example_name + '.log'
    log = open(logfile, 'wb')

    print_str = r'Date and time of execution: %s' % st   # modified on 210421
    log.write(print_str + r'\n')   # modified on 210421
    print (print_str)   # modified on 210421

    print_str = 'The following test is executing: %s\n' % example_name
    log.write(print_str + '\n\n')
    print (print_str)   # modified on 210421

    # ---------------------------------------------------------------------------

    # Initialize values:
    matlab_err = ''
    matlab_out = []
    num_files = 0
    total_passed = 0
    total_failed = 0
    poscars_passed = 0
    poscars_failed = 0
    poscars_number = 0
    bad_poscars = []
    status = 'Failed'
    error = ''
    duration = 0

    return_dict = {
        'name': example_name,
        'matlab_err': matlab_err,
        'matlab_out': matlab_out,
        'passed_files': num_files,
        'total_files': len(files_list),
        'total_passed': total_passed,
        'total_failed': total_failed,
        'poscars_passed': poscars_passed,
        'poscars_failed': poscars_failed,
        'poscars_number': poscars_number,
        'bad_poscars': bad_poscars,
        'status': status,
        'error': error,
        'duration': duration,
    }

    # ---------------------------------------------------------------------------
    # Run a short USPEX job:
    init_dir = os.getcwd()
    templates_dir = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])) + '/templates')
    calc_dir = os.path.abspath(init_dir + '/' + testid_str + '_' + example_name)
    results_dir = os.path.abspath(calc_dir + '/results1')

    copy_folder(templates_dir + '/' + example_name + '/', calc_dir)

    if not os.path.exists(calc_dir):
        error = '%s directory doesn\'t exist after copying. Exit.' % calc_dir
        log.write(error + '\n')
        print (error)   # modified on 210421
        return_dict['error'] = error
        return return_dict

    os.chdir(calc_dir)
    print_str = 'Changing to %s and execute USPEX job...' % calc_dir
    log.write(print_str + '\n')
    print (print_str)   # modified on 210421

    # ---------------------------------------------------------------------------
    # Insert '1 : fixRandSeed' to INPUT.txt to avoid usage of fixed random
    # scheme by users:
    # with open('INPUT.txt', 'a') as input_file:
    #    input_file.write('1 : fixRandSeed')
    # ---------------------------------------------------------------------------

    exec_list = ['USPEX', '-r']
    if uspex_octave != '':
        exec_list.append(uspex_octave)

    proc = subprocess.Popen(exec_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Poll process for new output until finished
    while True:
        outline = proc.stdout.readline()
        if outline == '' and proc.poll() != None:
            break
        row2print = datetime.datetime.fromtimestamp(time.time()).strftime(
            '%Y-%m-%d %H:%M:%S') + ' [T' + testid_local + ' ' + example_name + ']: ' + outline
        sys.stdout.write(row2print)
        sys.stdout.flush()
        matlab_out.append(row2print)

    (out, matlab_err) = proc.communicate()

    exitCode = proc.returncode
    exit_str = '\nExit code: %s\n' % exitCode
    print (exit_str)   # modified on 210421

    if example_name.find('restart_GULP') >= 0:
        files2delete = [
            'ANTISEEDS.mat',
            'ANTISEEDS.mat.backup',
            'Current_ORG.mat',
            'Current_ORG.mat.backup',
            'Current_POP.mat',
            'Current_POP.mat.backup',
            'USPEX_IS_DONE',
        ]

        for file_i in files2delete:
            if os.path.exists(file_i):
                os.remove(file_i)

        shutil.copy('INPUT_restart.txt', 'INPUT.txt')
        proc = subprocess.Popen(exec_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # shutil.copy('INPUT_restart.txt', 'INPUT.txt')

        # Poll process for new output until finished
        # matlab_out = []
        while True:
            outline = proc.stdout.readline()
            if outline == '' and proc.poll() != None:
                break
            row2print = datetime.datetime.fromtimestamp(time.time()).strftime(
                '%Y-%m-%d %H:%M:%S') + ' [T' + testid_local + ' ' + example_name + ']: ' + outline
            sys.stdout.write(row2print)
            sys.stdout.flush()
            matlab_out.append(row2print)

        (out, matlab_err) = proc.communicate()

        exitCode = proc.returncode
        exit_str = '\nExit code: %s\n' % exitCode
        print (exit_str)   # modified on 210421

    print_str = r''   # modified on 210421
    print_str += r'======================================================\n'   # modified on 210421
    print_str += r'=================   O U T P U T   ====================\n'   # modified on 210421
    print_str += r'Matlab output:\n%s\n' % ''.join(matlab_out)   # modified on 210421
    print_str += r'======================================================\n\n'   # modified on 210421
    log.writelines(print_str)
    log.writelines(exit_str)

    if matlab_err:
        print_str = r''   # modified on 210421
        print_str += r'======================================================\n'   # modified on 210421
        print_str += r'==================   E R R O R   =====================\n'   # modified on 210421
        print_str += r'Matlab error:\n%s' % matlab_err   # modified on 210421
        print_str += r'======================================================\n\n'   # modified on 210421
        log.writelines(print_str)
        print (print_str)   # modified on 210421

    print_str = 'USPEX execution finished. Start checking...'
    log.write(print_str + '\n')
    print (print_str)   # modified on 210421

    # ---------------------------------------------------------------------------


    # ---------------------------------------------------------------------------
    # Check directories and files exist:

    # num_files   = 0
    failed_list = []
    print_str = 'Checking files/directories in ' + os.getcwd() + ':'
    log.write(print_str + '\n')
    print (print_str)   # modified on 210421

    print_str = '\t%-8s %-50s %-8s' % ('Status', 'File/dir', 'Size (bytes)')
    log.write(print_str + '\n')
    print (print_str)   # modified on 210421

    print_str = '\t%-8s %-50s %-8s' % ('------', '--------', '------------')
    log.write(print_str + '\n')
    print (print_str)   # modified on 210421

    for cur_file in files_list:
        cur_num = check_file(cur_file, empty_dirs)
        num_files += cur_num
        if cur_num == 0:
            failed_list.append(cur_file)

    print_str = 'Files/directories exist: %i/%i' % (num_files, len(files_list))
    log.write(print_str + '\n')
    print (print_str)   # modified on 210421

#    if failed_list <> []:
    if failed_list <> []:   # modified on 210421
        print_str = '\nThe following files are not found or empty:\n\t' + '\n\t'.join(failed_list)
        log.write(print_str + '\n')
        print (print_str)   # modified on 210421
    # ---------------------------------------------------------------------------


    # ---------------------------------------------------------------------------
    # Compare files in referece and results1 folders:
    failed_refs = []
    total_refs = 0
    if os.path.exists(calc_dir + '/reference'):
        ref_list = os.listdir(calc_dir + '/reference')
        ref_list.sort()
        ref_list.remove('.svn')

        print_str = '\nCompare with files in the reference folder:\n'
        log.write(print_str + '\n')
        print (print_str)   # modified on 210421

        total_refs = len(ref_list)
        for i in ref_list:
            if ref_compare(calc_dir + '/reference/' + i, calc_dir + '/results1/' + i):
                ref_stat = 'OK'
            else:
                ref_stat = 'FAILED'
                failed_refs.append(i)

            print_str = '\t%-8s %-50s\t' % (ref_stat, i)
            log.write(print_str + '\n')
            print (print_str)   # modified on 210421
    # ---------------------------------------------------------------------------


    # ---------------------------------------------------------------------------
    # Compare results in OUTPUT.txt and Individuals files in results1/:
    if not os.path.exists(results_dir):
        error = '%s directory doesn\'t exist after USPEX execution. Exit.' % results_dir
        log.write(error + '\n')
        print (error)   # modified on 210421
        return_dict['error'] = error
        return return_dict

    os.chdir(results_dir)

    # Improved Mahdi's code for parsing of OUTPUT.txt:
    if not os.path.exists('OUTPUT.txt'):
        error = 'OUTPUT.txt file doesn\'t exist after USPEX execution. Exit.'
        log.write(error + '\n')
        print (error)   # modified on 210421
        return_dict['error'] = error
        return return_dict

    if example_name.find('metadynamics_GULP') >= 0:
        from_OUTPUT = parse_OUTPUT('META')
    elif example_name.find('Tinker-protein') >= 0:
        from_OUTPUT = parse_OUTPUT('PROTEIN')
    elif example_name.find('vcNEB_111') >= 0:
        from_OUTPUT = []
    else:
        from_OUTPUT = parse_OUTPUT()

    # Maxim's code for parsing of Individuals:
    if not os.path.exists('Individuals') and example_name.find('vcNEB_111') < 0:
        error = 'Individuals file doesn\'t exist after USPEX execution. Exit.'
        log.write(error + '\n')
        print (error)   # modified on 210421
        return_dict['error'] = error
        return return_dict

    struct_format = 'POSCARS'
    if example_name.find('metadynamics_GULP') >= 0:
        from_Individuals = get_individuals_content('Individuals_relaxed', 'META')
    elif example_name.find('Tinker-protein') >= 0:
        from_Individuals = get_individuals_content('Individuals', 'PROTEIN')
        struct_format = 'PDB'
    elif example_name.find('vcNEB_111') >= 0:
        from_Individuals = []
    else:
        from_Individuals = get_individuals_content('Individuals')

    # Comparison of the results:
    if example_name.find('vcNEB_111') < 0:
        print_str = '\n\nCompare contents of OUTPUT.txt and Individuals:'
        log.write(print_str + '\n')
        print (print_str)   # modified on 210421

        for col in columns2check:
            print_str = '\nChecking %s...' % col
            log.write(print_str + '\n')
            print (print_str)   # modified on 210421

            failed = 0
            passed = 0
            if len(from_OUTPUT[col]) == len(from_Individuals[col]['values']):
                for i in range(len(from_OUTPUT[col])):
#                    if from_OUTPUT[col][i] <> from_Individuals[col]['values'][i]:
                    if from_OUTPUT[col][i] != from_Individuals[col]['values'][i]:   # modified on 210421
                        failed += 1
                        total_failed += 1
                        print_str = 'OUTPUT: %s <> %s :Individuals' % (
                            from_OUTPUT[col][i], from_Individuals[col]['values'][i])
                        log.write(print_str + '\n')
                        print (print_str)   # modified on 210421
                    else:
                        passed += 1
                        total_passed += 1
                print_str = 'Failed: %s\t Passed: %s' % (failed, passed)
                log.write(print_str + '\n')
                print (print_str)    # modified on 210421

            else:
                total_failed += len(from_OUTPUT[col]) + len(from_Individuals[col]['values'])
                print_str = 'Number of data rows in OUTPUT.txt and in Individuals differs: %i <> %i' % (
                    len(from_OUTPUT[col]), len(from_Individuals[col]['values']))
                log.write(print_str + '\n')
                print (print_str)   # modified on 210421
                break
    # ---------------------------------------------------------------------------

    # ---------------------------------------------------------------------------
    # Compare structures in gatheredPOSCARS/PDB and BESTgatheredPOSCARS/PDB:
    if example_name.find('vcNEB_111') < 0:
        print_str = '\n\nCompare structures in gathered' + struct_format + ' and BESTgathered' + struct_format + ':'
        log.write(print_str + '\n')
        print (print_str)   # modified on 210421

        if not os.path.exists('gathered' + struct_format) or not os.path.exists('BESTgathered' + struct_format):
            error = 'gathered' + struct_format + ' or BESTgathered' + struct_format + ' don\'t exist after USPEX execution. Exit.'
            log.write(error + '\n')
            print (error)   # modified on 210421
            return_dict['error'] = error
            return return_dict

        id_list = count_POSCARS('BESTgathered' + struct_format, struct_format)

        poscars_number = len(id_list)
        for id in id_list:
            a = read_POSCAR('BESTgathered' + struct_format, id, struct_format)
            b = read_POSCAR('gathered' + struct_format, id, struct_format)

            if len(a) == len(b):
                for i in range(len(a)):
                    failed_flag = 0
                    if a[i] != b[i]:
                        print_str = 'Different structures: %i' % id
                        log.write(print_str + '\n')
                        print (print_str)   # modified on 210421

                        bad_poscars.append(id)
                        poscars_failed += 1
                        failed_flag = 1
                    break
                if failed_flag == 0:
                    poscars_passed += 1
            else:
                poscars_failed += 1

        print_str = 'Number of structures checked: %i' % poscars_number
        log.write(print_str + '\n')
        print (print_str)   # modified on 210421
    # ---------------------------------------------------------------------------

    # ---------------------------------------------------------------------------
    # Summary on the whole test:
    print_str = '\n======================================================================='
    log.write(print_str + '\n')
    print (print_str)   # modified on 210421

    print_str = 'Summary:\n'
    print_str += 'Test                   : %s\n' % ('T' + str(testid_local) + ' ' + example_name)
    print_str += 'Data comparison        : %3i/%3i\n' % (total_passed, (total_passed + total_failed))
    print_str += 'Structures matches     : %3i/%3i\n' % (poscars_passed, poscars_number)
    print_str += 'Files/directories exist: %3i/%3i\n' % (num_files, len(files_list))
    print_str += 'Reference files        : %3i/%3i\n' % (total_refs - len(failed_refs), total_refs)
    log.write(print_str + '\n')
    print (print_str)   # modified on 210421

    if total_failed == 0 and total_passed > 0 and ((total_refs > 0 and len(failed_refs) == 0) or total_refs == 0) and (
                num_files == len(files_list)) and poscars_failed == 0 and poscars_passed == poscars_number:
        status = 'Passed'
        print_str = 'Test passed!'
    else:
        status = 'Failed'
        print_str = 'Test failed!'

    if example_name.find('vcNEB_111') >= 0:
        if num_files == len(files_list):
            status = 'Passed'
            print_str = 'Test passed!'
        else:
            status = 'Failed'
            print_str = 'Test failed!'

    log.write(print_str + '\n')
    print (print_str)   # modified on 210421
    # ---------------------------------------------------------------------------

    # ---------------------------------------------------------------------------
    # End of test. Change to initial directory:
    os.chdir(init_dir)
    # ---------------------------------------------------------------------------
    log.close()

    test_end_time = time.time()
    duration = test_end_time - test_start_time

    return_dict = {
        'id': testid_local,
        'name': example_name,
        'matlab_err': matlab_err,
        'matlab_out': matlab_out,
        'passed_files': num_files,
        'total_files': len(files_list),
        'failed_refs': len(failed_refs),
        'total_refs': total_refs,
        'total_passed': total_passed,
        'total_failed': total_failed,
        'poscars_passed': poscars_passed,
        'poscars_failed': poscars_failed,
        'poscars_number': poscars_number,
        'bad_poscars': bad_poscars,
        'status': status,
        'error': error,
        'duration': duration,
    }
    return return_dict


# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
if testid == 'all':
    testid = valid_tests.keys()
    testid = sorted(testid)

for i in range(len(testid)):
    testid[i] = str(testid[i])
print ('The following tests will be executed: %s\n' % ', '.join(testid))   # modified on 210421

# exit(9999)

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
init_dir = os.getcwd()
calc_dir = os.path.abspath(init_dir + '/' + st + '___T' + init_option)

if not os.path.exists(calc_dir):
    os.mkdir(calc_dir)

results_list = []
for i in testid:
    os.chdir(calc_dir)
    # Import tests configurations:
    module = 'tests.' + valid_tests[int(i)]
    module = __import__(module, globals(), locals(), ['*'])

    for k in dir(module):
        locals()[k] = getattr(module, k)

    # Run the test:
    result = run_test(i, uspex_octave)
    results_list.append(result)
    os.chdir(init_dir)

end_time = time.time()
duration_sec = end_time - start_time

logfile = calc_dir + '/summary.log'
log = open(logfile, 'wb')

print ('\n\n==================================================================\n')   # modified on 210421
print_str = 'Total statistics:'
for res in results_list:
    print_str += '\n' + '\n'
    print_str += '\n' + 'Test ID                 : %s' % res['id']
    print_str += '\n' + 'Test name               : %s' % res['name']
    print_str += '\n' + 'Test duration           : %s' % duration(res['duration'])
    print_str += '\n' + 'Status                  : %s' % res['status'].upper()
    print_str += '\n' + 'Data comparison         : %3i/%3i' % (
        res['total_passed'], res['total_passed'] + res['total_failed'])
    print_str += '\n' + 'Structures matches      : %3i/%3i' % (res['poscars_passed'], res['poscars_number'])
    print_str += '\n' + 'Files/directories exist : %3i/%3i' % (res['passed_files'], res['total_files'])
    print_str += '\n' + 'Reference files         : %3i/%3i' % (
        res['total_refs'] - res['failed_refs'], res['total_refs'])
    print_str += '\n' + 'Matlab error(s):\n%s\n' % res['matlab_err']

print_str += '\n\n' + 'Duration: %s' % duration(duration_sec)
print (print_str)   # modified on 210421
log.write(print_str)
log.close()
# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# End of testing
sys.exit(0)
