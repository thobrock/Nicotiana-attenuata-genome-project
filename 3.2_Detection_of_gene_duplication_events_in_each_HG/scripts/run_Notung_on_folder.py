#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
 Title:     run_Notung_on_folder.py
 Author:    Thomas Brockmöller

 Description:
     runs Notung on all phylogenetic trees in a folder

 Example:
     run_Notung_on_folder.py -i trees/ -s speciesTree.nwk -o output/ -p .nt_ali_cleaned.phy_phyml_tree.txt

'''

import os
import glob
import argparse
from random import shuffle
from datetime import datetime
import sys

path_NOTUNG = 'software/Notung-2.6/Notung-2.6.jar'

parser = argparse.ArgumentParser(
    description='If you find bugs please contact: tbrockmoeller@ice.mpg.de')

parser.add_argument('-i', '--input', metavar='file', type=str,
                    help="input folder",
                    required=True)
parser.add_argument('-s', '--species', metavar='file', type=str,
                    help="species tree",
                    required=True)
parser.add_argument('-p', '--prefix', metavar='string', type=str,
                    help="gene tree file prefix (optional) (default: .nwk)",
                    required=False)
parser.add_argument('-o', '--output', metavar='file', type=str,
                    help="output folder (optional) (default: here)",
                    required=False)
parser.add_argument('--check', action="store_true", help="check if file was processed before (optional)",
                    required=False)
parser.add_argument('--shuffle', action="store_true", help="shuffle files for processing (optional)",
                    required=False)
parser.add_argument('--sorted', action="store_true", help="sort files for processing (optional)", required=False)
parser.add_argument('--count', action="store_true", help="count processed (optional)", required=False)
args = parser.parse_args()

# Example Notung command: java -jar /data4/Programs/Notung-2.6/Notung-2.6.jar -b all_trees.txt
# --root --rootscores --treeoutput newick --stpruned --log --savepng
# --speciestag prefix --progressbar --silent --exact-losses --info --treestats
default_PARAMETER_FOR_NOTUNG = {'--root': '',
                                '--rootscores': '',
                                # '--treeoutput': 'newick',
                                '--stpruned': '',
                                '--log': '',
                                '--savepng': '',
                                '--speciestag': 'prefix',
                                '--silent': '',
                                '--exact-losses': '',
                                '--info': '',
                                '--treestats': '',
                                '--outputdir': '.'}

file_SPECIES_TREE = args.species
folder_INPUT = args.input

var_CHECK_FILES = False
if args.check:
    var_CHECK_FILES = True

var_COUNT = False
if args.count:
    var_COUNT= True

var_SORTING = True
if args.sorted:
    var_SORTING = False

var_SHUFFLE = False
if args.shuffle:
    var_SHUFFLE = True

var_FILE_PREFIX = '.nwk'

if args.prefix:
    var_FILE_PREFIX = args.prefix

if args.output:
    default_PARAMETER_FOR_NOTUNG['--outputdir'] = args.output


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'


def notung(speciestree, genetreefile, params={}):
    global path_NOTUNG
    parameterString = ''
    for p in params:
        parameterString += ' ' + p + ' ' + params[p]
    os.system('java -jar ' + path_NOTUNG + ' -s ' + speciestree + ' -g ' + genetreefile + ' ' + parameterString)
    return 'java -jar ' + path_NOTUNG + ' -s ' + speciestree + ' -g ' + genetreefile + ' ' + parameterString



if var_COUNT:
    i = 0
    processing_files = glob.glob(os.path.join(folder_INPUT, '*' + var_FILE_PREFIX))
    for infile in processing_files:
        if os.path.isfile('%s%s.rooting.0' % (args.output, infile.split('/')[-1])):
            print(infile.split('/')[-1].split('-')[0])
            i += 1
    print('Processed files: %s' % i)

    sys.exit()

if not var_SORTING:
    processing_files = sorted(glob.glob(os.path.join(folder_INPUT, '*' + var_FILE_PREFIX)))
elif var_SHUFFLE:
    processing_files = glob.glob(os.path.join(folder_INPUT, '*' + var_FILE_PREFIX))
    shuffle(processing_files)
else:
    processing_files = glob.glob(os.path.join(folder_INPUT, '*' + var_FILE_PREFIX))

for infile in processing_files:
    print(bcolors.OKGREEN + 'Processing file: ' + bcolors.FAIL + infile + bcolors.ENDC)
    if var_CHECK_FILES:
        if os.path.isfile('%s%s.rooting.0' % (args.output, infile.split('/')[-1])):
            print(bcolors.OKBLUE + 'File already processed: ' + bcolors.FAIL + infile + bcolors.ENDC)
            continue

    time_START = datetime.now()
    command = notung(file_SPECIES_TREE, infile, default_PARAMETER_FOR_NOTUNG)
    time_END = datetime.now()
    time_DIFF = (time_END - time_START).total_seconds()
    print('\t', time_DIFF, 'seconds')
    print(command)
