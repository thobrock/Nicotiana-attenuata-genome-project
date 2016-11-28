#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
 Title:     predicting_phylogenetic_trees.py
 Author:    Thomas Brockmöller

 Description:
     predict phylogenetic trees for groups
     this script used third party software
'''


import os
import argparse
import shutil

path_TRANSLATORX = 'third_party_software/translatorx_vLocal.pl'
path_CONVERTFASTATOPHYLIP = 'third_party_software/convertFastaToPhylip.py'
path_JMODELTEST2 = 'third_party_software/jmodeltest/jModelTest.jar'
path_TRIMAL = 'third_party_software/trimal'
path_PHYML = 'third_party_software/phyml/src/phyml'
path_PHYML_MPI = 'phyml'

##################################
###         USER INPUT         ###
##################################

parser = argparse.ArgumentParser(
    description='This script builds phylogenetic trees of a corresponding orthologous group (defined with groupID or geneID). This script uses MUSCLE, TRIMAL, jMODELTEST and PHYML. If you find bugs please contact: tbrockmoeller@ice.mpg.de')

parser.add_argument('-i', '--input', metavar='file', type=str,
                    help="input groups file",
                    required=True)
parser.add_argument('-t', '--tmp', metavar='file', type=str,
                    help="folder for tmp files (optional) (default: tmp/)",
                    required=False)
parser.add_argument('-s', '--spp', metavar='string', type=str,
                    help="select species (optional) (default: ara,can,cca,csa,mgu,nat,nis,nit,nio,ptr,sme,sly,stu,vvi)",
                    required=False)
parser.add_argument('-st', '--steps', metavar='string', type=str,
                    help="steps (optional) (default: 1,1,1,1,1,1) (0:sequence file, 1:multiple alignment, 2:best model, 3:extract best model, 4:build tree, 5:color file)",
                    required=False)
parser.add_argument('-b', '--boot', metavar='N', type=int, help="bootstrap value (optional) (default: -5)", required=False)
parser.add_argument('-n', metavar='N', type=int, help="number of substitution schemes (optional) (default: 3)",
                    required=False)
parser.add_argument('-o', '--output', metavar='file', type=str,
                    help="output folder (optional) (default: tmp/)", required=False)
parser.add_argument('-m', '--threads', metavar='N', type=int, help="threads for jModelTest and phyml (optional) (default: 1)",
                    required=False)
parser.add_argument('-log', '--log', metavar='file', type=str,
                    help="log files (optional)",
                    required=False)
parser.add_argument('-log2', '--log2', metavar='file', type=str,
                    help="log files 2 (optional)",
                    required=False)
parser.add_argument('-sp', '--specific', metavar='string', type=str,
                    help="predict phylogeny only for this specific clusters (optional) (default: cluster1,cluster2)",
                    required=False)
parser.add_argument('-min', '--min', metavar='N', type=int,
                    help="smallest group size (optional) (default: 3)",
                    required=False)
parser.add_argument('-max', '--max', metavar='N', type=int,
                    help="largest group size (optional) (default: 30)",
                    required=False)
parser.add_argument('-x', '--extra', metavar='string', type=str,
                    help="extra parameter [a: all, n: only groups with NIATT and NIOBT, x: exclude groups with NIATT] (optional) (default: a)",
                    required=False)
parser.add_argument('-mpi', '--mpi', action="store_true", help="do not use mpi for phyml (optional)", required=False)
parser.add_argument('-g', '--groups', metavar='file', type=str,
                    help="load file with groupIDs to process (optional)",
                    required=False)
parser.add_argument('-p', '--points', action="store_true", help="replace . in geneID with ___ (optional)", required=False)
parser.add_argument("-v", action="count", default=0, help='show more information (optional)')
args = parser.parse_args()

##################################
###         CONSTANTS          ###
##################################

# Steps:
# 0 -> write sequence file
# 1 -> build alignment
# 2 -> best model
# 3 -> extract model
# 4 -> build tree
# 5 -> color file
const_DO_STEPS = [1, 1, 1, 1, 1, 1]

const_USE_MPI = True
const_SMALLEST_GROUP_SIZE = 3
const_LARGEST_GROUP_SIZE = 10000000
const_EXTRA_PARAMETER = 'a'
const_ALLOWED_SPECIES = ['ara', 'can', 'cca', 'csa', 'mgu', 'nat', 'nis', 'nit', 'nio', 'ptr', 'sme', 'sly', 'stu', 'vvi']
const_NOT_ALLOWED_SPECIES_IN_TREE = []
const_ALLOWED_GROUPS = []
const_UPLOAD_INFO_CODE = {'ara': 'range	#C6C6FF	Arabidopsis thaliana',
                          'can': 'range	#4CCD1E	Capsicum annuum',
                          'cca': 'range	#E995B4	Coffee canephora',
                          'csa': 'range	#6755E3	Cucumis sativus',
                          'gma': 'range	#2F74D0	Glycine max',
                          'mgu': 'range	#996666	Mimulus guttatus',
                          'mtr': 'range	#44B4D5	Medicago truncatula',
                          'nat': 'range	#FFCC00	Nicotiana attenuata',
                          'nio': 'range	#FF9900	Nicotiana obtusifolia',
                          'nis': 'range	#FF6600	Nicotiana sylvestris',
                          'nit': 'range	#FF3300	Nicotiana tomentosiformis',
                          'ptr': 'range	#8ED6EA	Populus trichocarpa',
                          'sme': 'range	#6BCCA6	Solanum melongena',
                          'sly': 'range	#99CC66	Solanum lycopersicum',
                          'stu': 'range	#999966	Solanum tuberosum',
                          'vvi': 'range	#5EAE9E	Vitis vinifera'}

const_BOOTSTRAP = -5
const_NUMBER_OF_SCHEMES = 3
const_NUMBER_OF_THREADS = 1

const_REPLACE_POINTS = False

##################################
###           FOLDER           ###
##################################

folder_OUTPUT = 'tmp/'
folder_TMP_FILES = 'tmp/'

##################################
###           FILES            ###
##################################

file_FASTA_CDS_ARA = 'FASTA/Athaliana_167_cds_primaryTranscriptOnly.fa'
file_FASTA_CDS_CAN = 'FASTA/Capsicum.annuum.L_Zunla-1_v2.0_CDS.fa'
file_FASTA_CDS_CCA = 'FASTA/coffea_cds.fna'
file_FASTA_CDS_CSA = 'FASTA/Csativus_122_cds_primaryTranscriptOnly.fa'
file_FASTA_CDS_MGU = 'FASTA/Mguttatus_v2.0_256_cds_primaryTranscriptOnly.fa'
file_FASTA_CDS_NAT = 'FASTA/NIATTr2.AN5.clean.cds'
file_FASTA_CDS_NIO = 'FASTA/NIOBTv3.AN7.cds'
file_FASTA_CDS_PTR = 'FASTA/Ptrichocarpa_210_cds_primaryTranscriptOnly.fa'
file_FASTA_CDS_SME = 'FASTA/SME_r2.5.1_cds_ip.fa'
file_FASTA_CDS_SLY = 'FASTA/Slycopersicum_225_cds_primaryTranscriptOnly.fa'
file_FASTA_CDS_STU = 'FASTA/Stuberosum_206_cds_primaryTranscriptOnly.fa'
file_FASTA_CDS_VVI = 'FASTA/Vvinifera_145_cds_primaryTranscriptOnly.fa'

file_POTATO_ANNOTATION_FILE = 'FASTA/Stuberosum_206_annotation_info.txt'

##################################
###      GLOBAL VARIABLES      ###
##################################

var_GROUPS = {}
var_SEQUENCES = {}
tmpPotatoIDs = {}

file_LOG = None
file_LOG2 = None

if args.points:
    const_REPLACE_POINTS = True

if args.log:
    file_LOG = args.log

if args.log2:
    file_LOG2 = args.log2

if args.output:
    folder_OUTPUT = args.output

if args.input:
    file_GROUPS = args.input

if args.n:
    const_NUMBER_OF_SCHEMES = args.n

if args.tmp:
    folder_TMP_FILES = args.tmp

if args.spp:
    const_ALLOWED_SPECIES = args.spp.split(',')

if args.boot:
    const_BOOTSTRAP = args.boot

if args.threads:
    const_NUMBER_OF_THREADS = args.threads

if args.steps:
    const_DO_STEPS = map(int, args.steps.split(','))

if args.specific:
    const_ALLOWED_GROUPS = args.specific.split(',')

if args.min:
    const_SMALLEST_GROUP_SIZE = int(args.min)

if args.max:
    const_LARGEST_GROUP_SIZE = int(args.max)

if args.extra:
    const_EXTRA_PARAMETER = args.extra

if args.mpi:
    const_USE_MPI = False

##################################
###         FUNCTIONS          ###
##################################

def loadGroupIdsFile(filename):
    groupIDs = []
    f = open(filename, 'r')
    for line in f:
        groupIDs.append(line.strip())
    f.close()
    return groupIDs


def readFastaFile(filename):
    global var_SEQUENCES
    global const_REPLACE_POINTS
    name = ''
    seq = ''
    f = open(filename, 'r')
    for line in f:
        if line[0] == '>':
            var_SEQUENCES[name] = seq
            seq = ''
            tmp = line.strip().split()[0].split('|')[0]
            name = tmp[1:]
            if const_REPLACE_POINTS:
                name = name.replace('.', '___')
        else:
            seq += line.strip()
    var_SEQUENCES[name] = seq
    f.close()
    return


def readGroupsFile(filename):
    global var_GROUPS
    global const_REPLACE_POINTS
    gF = open(filename, 'r')
    for line in gF:
        if const_REPLACE_POINTS:
            line = line.replace('.', '___')
        groupName = line.split(':')[0]
        groupGenes = line.split(':')[1].split()
        var_GROUPS[groupName] = groupGenes  # ['cluster0001', ['gene1', 'gene2']]
    gF.close()


def cleanGroup(groupGenes):
    global const_ALLOWED_SPECIES
    newGroup = []
    for gene in groupGenes:
        if gene.split('|')[0] in const_ALLOWED_SPECIES:
            newGroup.append(gene)
    return newGroup


def checkGroup(groupGenes):
    global const_NOT_ALLOWED_SPECIES_IN_TREE
    for gene in groupGenes:
        if gene.split('|')[0] in const_NOT_ALLOWED_SPECIES_IN_TREE:
            return False
    return True


def multipleAlignment(inputFile, outputFile):
    global path_TRANSLATORX
    os.system(
        'perl ' + path_TRANSLATORX + ' -i ' + inputFile + ' -o ' + outputFile + ' -p M')


def trimal(infile, outfile, params={'gt': '0.8'}):
    global path_TRIMAL
    parameterString = ''
    for p in params:
        parameterString += ' -' + p + ' ' + params[p]
    os.system(path_TRIMAL + ' -in ' + infile + ' -out ' + outfile + parameterString)


def convertToPhylip(infile, outfile):
    global path_CONVERTFASTATOPHYLIP
    os.system('python ' + path_CONVERTFASTATOPHYLIP + ' ' + infile + ' ' + outfile)


def jmodeltest(infile, outfile, params={}):
    global path_JMODELTEST2
    parameterString = ''
    for p in params:
        parameterString += ' -' + p + ' ' + params[p]
    commandForModeltest = 'java -jar ' + path_JMODELTEST2 + ' -d ' + infile + parameterString
    print commandForModeltest
    os.system(commandForModeltest + ' > ' + outfile)


def findBestModel(infile):
    f = open(infile)
    lines = f.readlines()
    tmp = lines.index("* AIC MODEL SELECTION : Best Model's command line\n")
    commandIndex = tmp + 2
    command = lines[commandIndex].split()
    modelIndex = command.index("-m")
    modelDigits = command[modelIndex + 1]
    [method, model, fa, fc, fg, ft, kappa, titv, Ra, Rb, Rc, Rd, Re, Rf, pInv, gamma] = lines[-1].strip().split()
    f.close()
    return {'method': method, 'model': model, 'modelDigits': modelDigits, 'fa': fa, 'fc': fc, 'fg': fg, 'ft': ft,
            'kappa': kappa, 'titv': titv,
            'Ra': Ra, 'Rb': Rb, 'Rc': Rc, 'Rd': Rd, 'Re': Re, 'Rf': Rf, 'pInv': pInv, 'gamma': gamma}


def phyml(infile, modelParams, params={}, threads=1):
    global path_PHYML
    global path_PHYML_MPI
    parameterString = ''
    for p in params:
        parameterString += ' -' + p + ' ' + params[p]
    modelString = '-m ' + modelParams['modelDigits'] + ' -f ' + modelParams['fa'] + ',' + modelParams['fc'] + ',' + \
                  modelParams['fg'] + ',' + modelParams['ft'] + ' -t ' + modelParams['titv']
    if not modelParams['pInv'] == 'N/A':
        modelString += ' -v ' + modelParams['pInv']
    if not modelParams['gamma'] == 'N/A':
        modelString += ' -a ' + modelParams['gamma']
    if const_USE_MPI:
        print('mpirun -np ' + str(threads) + ' ' + path_PHYML_MPI + ' --no_memory_check -i ' + infile + parameterString + ' ' + modelString)
        os.system('mpirun -np ' + str(threads) + ' ' + path_PHYML_MPI + ' --no_memory_check -i ' + infile + parameterString + ' ' + modelString)
    else:
        os.system(path_PHYML + ' --no_memory_check -i ' + infile + parameterString + ' ' + modelString)


def buildColorFile(filename, groupGenes):
    f = open(filename, 'w')
    for gene in sorted(groupGenes):
        f.write(gene + '\t' + const_UPLOAD_INFO_CODE[gene.split('|')[0]] + '\n')
    f.write('\n')
    f.close()

##################################
###            MAIN            ###
##################################

if args.v >= 1: print 'READ GROUPS FILES'
readGroupsFile(file_GROUPS)

if args.v >= 1: print 'READ FASTA FILES'
readFastaFile(file_FASTA_CDS_ARA)
readFastaFile(file_FASTA_CDS_CAN)
readFastaFile(file_FASTA_CDS_CCA)
readFastaFile(file_FASTA_CDS_CSA)
readFastaFile(file_FASTA_CDS_MGU)
readFastaFile(file_FASTA_CDS_NAT)
readFastaFile(file_FASTA_CDS_NIO)
readFastaFile(file_FASTA_CDS_PTR)
readFastaFile(file_FASTA_CDS_SLY)
readFastaFile(file_FASTA_CDS_STU)
readFastaFile(file_FASTA_CDS_SME)
readFastaFile(file_FASTA_CDS_VVI)

f = open(file_POTATO_ANNOTATION_FILE, 'r')
for line in f:
    values = line.strip().split()
    tmpPotatoIDs[values[3]] = values[2]
f.close()

var_GROUPS_ID_PROCESS = None
if args.groups:
    const_ALLOWED_GROUPS += loadGroupIdsFile(args.groups)

for groupName in sorted(var_GROUPS, reverse=True):

    if const_ALLOWED_GROUPS:
        if not groupName in const_ALLOWED_GROUPS:
            continue

    if file_LOG2:
        os.system('echo %s >> %s' % (groupName, file_LOG2))

    groupGenes = var_GROUPS[groupName]
    groupGenes = cleanGroup(groupGenes)

    if len(groupGenes) < const_SMALLEST_GROUP_SIZE:
        print 'WARNING: group is too small (min=' + str(const_SMALLEST_GROUP_SIZE) + ') !'
        continue

    if len(groupGenes) > const_LARGEST_GROUP_SIZE:
        print 'WARNING: group is too large (max=' + str(const_LARGEST_GROUP_SIZE) + ') !'
        continue

    if not checkGroup(groupGenes):
        print 'WARNING: group is excluded (' + str(const_LARGEST_GROUP_SIZE) + ') !'
        continue

    if const_EXTRA_PARAMETER == 'n':
        found_NIATT = False
        for gene in groupGenes:
            spp = gene.split('|')[0]
            if spp == 'nat' or spp == 'nio':
                found_NIATT = True
        if not found_NIATT:
            print 'WARNING: group does not contain NIATT or NIOBT genes - but are forced to do so!'
            continue

    if const_EXTRA_PARAMETER == 'x':
        found_NIATT = False
        for gene in groupGenes:
            spp = gene.split('|')[0]
            if spp == 'nat' or spp == 'nio':
                found_NIATT = True
        if found_NIATT:
            print 'WARNING: group does contain NIATT or NIOBT genes - but are forced not to do so!'
            continue

    ##### STEP 0: create sequence file #####
    if const_DO_STEPS[0] == 1:
        if args.v >= 1: print 'WRITE SEQUENCE FILE'
        f = open(folder_TMP_FILES + groupName + '-sequences.fasta', 'w')
        for gene in sorted(groupGenes):
            # for STU genes do not have 'points' in the geneID
            if gene.split('|')[0] == 'stu':
                f.write('>' + gene + '\n' + var_SEQUENCES[tmpPotatoIDs[gene.split('|')[1]]] + '\n')
            elif gene.split('|')[0] == 'mgu':
                # for MGU genes have 'points' in the geneID and need to be treated differently
                if const_REPLACE_POINTS:
                    f.write('>' + gene[:-4] + '\n' + var_SEQUENCES[gene.split('|')[1][:-4]] + '\n')
                else:
                    f.write('>' + gene[:-2] + '\n' + var_SEQUENCES[gene.split('|')[1][:-2]] + '\n')
            else:
                f.write('>' + gene + '\n' + var_SEQUENCES[gene.split('|')[1]] + '\n')
        f.close()

        if not os.path.isfile(folder_TMP_FILES + groupName + '-sequences.fasta'):
            if file_LOG:
                with open(file_LOG, "a") as myfile:
                    myfile.write('ERROR' + '\t' + 'WRITE SEQUENCE FILE' + '\t' + groupName + '\t' + folder_TMP_FILES + groupName + '-sequences.fasta' + '\t' + 'not exists' + '\n')
                continue

        if not folder_OUTPUT == folder_TMP_FILES:
            shutil.copy(folder_TMP_FILES + groupName + '-sequences.fasta',
                        folder_OUTPUT + groupName + '-sequences.fasta')

    ##### STEP 1: build multiple alignment, trim and convert to PHYLIP format #####
    if const_DO_STEPS[1]:
        if args.v >= 1: print 'BUILT MULTIPLE ALIGNMENT'
        multipleAlignment(folder_TMP_FILES + groupName + '-sequences.fasta',
                          folder_TMP_FILES + groupName + '-sequences_align')

        if not os.path.isfile(folder_TMP_FILES + groupName + '-sequences_align.nt_ali.fasta'):
            if file_LOG:
                with open(file_LOG, "a") as myfile:
                    myfile.write('ERROR' + '\t' + 'BUILD ALIGNMENT FILE' + '\t' + groupName + '\t' + folder_TMP_FILES + groupName + '-sequences_align.nt_ali.fasta' + '\t' + 'not exists' + '\n')
                continue

        if args.v >= 1: print 'CLEAN SEQUENCES (trimal)'
        trimal(folder_TMP_FILES + groupName + '-sequences_align.nt_ali.fasta',
               folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.fasta')  # -cons 60

        if not os.path.isfile(folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.fasta'):
            if file_LOG:
                with open(file_LOG, "a") as myfile:
                    myfile.write('ERROR' + '\t' + 'CLEAN ALIGNMENT FILE' + '\t' + groupName + '\t' + folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.fasta' + '\t' + 'not exists' + '\n')
                continue

        if not folder_OUTPUT == folder_TMP_FILES:
            shutil.copy(folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.fasta',
                        folder_OUTPUT + groupName + '-sequences_align.nt_ali_cleaned.fasta')

        if args.v >= 1: print 'CONVER TO PHYLIP FORMAT'
        convertToPhylip(folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.fasta',
                        folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.phy')

        if not os.path.isfile(folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.phy'):
            if file_LOG:
                with open(file_LOG, "a") as myfile:
                    myfile.write('ERROR' + '\t' + 'CONVERT TO PHYLIP FORMAT' + '\t' + groupName + '\t' + folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.phy' + '\t' + 'not exists' + '\n')
                continue

    ##### STEP 2: calculate best substitution model ####
    if const_DO_STEPS[2]:
        if args.v >= 1: print 'FIND BEST MODEL'
        jmodeltest(folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.fasta',
                   folder_TMP_FILES + groupName + '-bestModel.out',
                   {'f': '', 'i': '', 'g': '4', 's': str(const_NUMBER_OF_SCHEMES), 'AIC': '', 'a': '',
                    'tr': str(const_NUMBER_OF_THREADS)})

        if not os.path.isfile(folder_TMP_FILES + groupName + '-bestModel.out'):
            if file_LOG:
                with open(file_LOG, "a") as myfile:
                    myfile.write('ERROR' + '\t' + 'CALCULATE SUBSTITUTION MODEL' + '\t' + groupName + '\t' + folder_TMP_FILES + groupName + '-bestModel.out' + '\t' + 'not exists' + '\n')
                continue

    ##### STEP 3: extract best substitution model
    if const_DO_STEPS[3]:
        bestModelParameters = findBestModel(folder_TMP_FILES + groupName + '-bestModel.out')

        if not folder_OUTPUT == folder_TMP_FILES:
            shutil.copy(folder_TMP_FILES + groupName + '-bestModel.out', folder_OUTPUT + groupName + '-bestModel.out')

    ##### STEP 4: build tree #####
    if const_DO_STEPS[4]:
        if args.v: print 'BUILD TREE'
        phyml(folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.phy', bestModelParameters,
              {'d': 'nt', 'b': str(const_BOOTSTRAP)}, const_NUMBER_OF_THREADS)

        # '-quiet': '',

        if not os.path.isfile(folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.phy_phyml_tree.txt'):
            if file_LOG:
                with open(file_LOG, "a") as myfile:
                    myfile.write('ERROR' + '\t' + 'BUILD TREE FILE' + '\t' + groupName + '\t' + folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.phy_phyml_tree.txt' + '\t' + 'not exists' + '\n')
                continue

        if not folder_OUTPUT == folder_TMP_FILES:
            shutil.copy(folder_TMP_FILES + groupName + '-sequences_align.nt_ali_cleaned.phy_phyml_tree.txt',
                        folder_OUTPUT + groupName + '-sequences_align.nt_ali_cleaned.phy_phyml_tree.txt')

    ##### STEP 5: build color file #####
    if const_DO_STEPS[5]:
        if args.v >= 1: print 'BUILD COLOR FILE'

        buildColorFile(folder_TMP_FILES + groupName + '-sequences.color', groupGenes)

        if not os.path.isfile(folder_TMP_FILES + groupName + '-sequences.color'):
            if file_LOG:
                with open(file_LOG, "a") as myfile:
                    myfile.write('ERROR' + '\t' + 'CREATE COLOR FILE' + '\t' + groupName + '\t' + folder_TMP_FILES + groupName + '-sequences.color' + '\t' + 'not exists' + '\n')
                continue

        if not folder_OUTPUT == folder_TMP_FILES:
            shutil.copy(folder_TMP_FILES + groupName + '-sequences.color',
                        folder_OUTPUT + groupName + '-sequences.color')
