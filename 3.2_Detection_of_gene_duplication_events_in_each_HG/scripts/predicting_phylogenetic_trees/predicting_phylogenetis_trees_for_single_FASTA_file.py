#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
 Title:     predicting_phylogenetic_trees_for_single_FASTA_file.py
 Author:    Thomas Brockmöller

 Description:
     predict phylogenetic trees for genes in a FASTA file
     this script used third party software
'''

import os
import argparse
import shutil
import math

path_TRANSLATORX = 'third_party_software/translatorx_vLocal.pl'
path_CONVERTFASTATOPHYLIP = 'third_party_software/convertFastaToPhylip.py'
path_JMODELTEST2 = 'third_party_software/jmodeltest/jModelTest.jar'
path_TRIMAL = 'third_party_software/trimal'
path_PHYML = 'third_party_software/phyml/src/phyml'
path_PHYML_MPI = 'phyml'

file_NIAT_EXPRESSION = 'additional_files/gene_tpm_matrix.txt'

var_NIAT_EXPRESSION_TISSUES = ['ROT', 'FLB', 'OFL', 'LEC', 'LET', 'STT', 'OVA', 'STI', 'NEC', 'PED', 'POL', 'SNP', 'STO', 'STS',
    'COE', 'COL', 'ANT', 'SED', 'SES', 'SEW']

##################################
###         USER INPUT         ###
##################################

parser = argparse.ArgumentParser(
    description='This script builds phylogenetic trees of a corresponding orthologous group (defined with groupID or geneID) and uploads them to http://itol.embl.de. This script uses MUSCLE, TRIMAL, jMODELTEST and PHYML. If you find bugs please contact: tbrockmoeller@ice.mpg.de')

parser.add_argument('-i', '--input', metavar='file', type=str,
                    help="input FASTA file",
                    required=True)
parser.add_argument('-t', '--tmp', metavar='file', type=str,
                    help="folder for tmp files (optional) (default: tmp/)",
                    required=False)
parser.add_argument('-u', '--upload', action="store_true", help="upload tree file (optional) (default: 0)",
                    required=False)
parser.add_argument('-st', '--steps', metavar='string', type=str,
                    help="steps (optional) (default: 1,1,1,1,1,1,1) (0:sequence file, 1:multiple alignment, 2:best model, 3:extract best model, 4:build tree, 5: write color file, 6: expression file)",
                    required=False)
parser.add_argument('-b', '--boot', metavar='N', type=int, help="bootstrap value (optional) (default: -5)",
                    required=False)
parser.add_argument('-n', metavar='N', type=int, help="number of substitution schemes (optional) (default: 3)",
                    required=False)
parser.add_argument('-o', '--output', metavar='file', type=str,
                    help="output folder (optional) (default: tmp/)", required=False)
parser.add_argument('-m', '--threads', metavar='N', type=int,
                    help="threads for jModelTest and phyml (optional) (default: 1)",
                    required=False)
parser.add_argument('-log', '--log', metavar='file', type=str,
                    help="log files (optional)",
                    required=False)
parser.add_argument('-log2', '--log2', metavar='file', type=str,
                    help="log files 2 (optional)",
                    required=False)
parser.add_argument('-mpi', '--mpi', action="store_true", help="do not use mpi for phyml (optional)", required=False)
parser.add_argument('-p', '--points', action="store_true", help="replace . in geneID with ___ (optional)",
                    required=False)
parser.add_argument('-tr', '--trim', action="store_true", help="do not clean the alignment (optional)",
                    required=False)
parser.add_argument('--extra_trimal', metavar='string', type=str, help="extra parameter for trimal (optional)",
                    required=False)
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
const_DO_STEPS = [1, 1, 1, 1, 1, 1, 1]

const_USE_MPI = True
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
                          'NIT': 'range	#FF3300	Nicotiana tomentosiformis',
                          'NIS': 'range	#FF6600	Nicotiana sylvestris',
                          'nit': 'range	#FF3300	Nicotiana tomentosiformis',
                          'ptr': 'range	#8ED6EA	Populus trichocarpa',
                          'sme': 'range	#6BCCA6	Solanum melongena',
                          'sly': 'range	#99CC66	Solanum lycopersicum',
                          'stu': 'range	#999966	Solanum tuberosum',
                          'pin': 'range	#6E5A79	Petunia inflata',
                          'pax': 'range	#744B68	Petunia axillaris',
                          'vvi': 'range	#5EAE9E	Vitis vinifera'}

const_BOOTSTRAP = -5
const_NUMBER_OF_SCHEMES = 3
const_NUMBER_OF_THREADS = 1

const_UPLOAD_TREE = 0

const_REPLACE_POINTS = False
const_TRIM_ALIGNMENT = True

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
file_FASTA_CDS_NIS = 'FASTA/OLD_NISYL.AN7.cds'
file_FASTA_CDS_NIT = 'FASTA/OLD_NITOM.AN7.cds'
file_FASTA_CDS_PTR = 'FASTA/Ptrichocarpa_210_cds_primaryTranscriptOnly.fa'
file_FASTA_CDS_SME = 'FASTA/SME_r2.5.1_cds_ip.fa'
file_FASTA_CDS_SLY = 'FASTA/Slycopersicum_225_cds_primaryTranscriptOnly.fa'
file_FASTA_CDS_STU = 'FASTA/Stuberosum_206_cds_primaryTranscriptOnly.fa'
file_FASTA_CDS_VVI = 'FASTA/Vvinifera_145_cds_primaryTranscriptOnly.fa'

file_POTATO_ANNOTATION_FILE = 'FASTA/Stuberosum_206_annotation_info.txt'

##################################
###      GLOBAL VARIABLES      ###
##################################

var_SEQUENCES = {}

file_LOG = None
file_LOG2 = None

var_ADDITIONAL_TRIMAL_PARAMETER = ''

if args.points:
    const_REPLACE_POINTS = True

if args.log:
    file_LOG = args.log

if args.log2:
    file_LOG2 = args.log2

if args.output:
    folder_OUTPUT = args.output

if args.input:
    file_INPUT_FASTA = args.input
    file_INPUT_FASTA_SHORT = args.input.split('/')[-1]

if args.n:
    const_NUMBER_OF_SCHEMES = args.n

if args.tmp:
    folder_TMP_FILES = args.tmp

if args.boot:
    const_BOOTSTRAP = args.boot

if args.upload:
    const_UPLOAD_TREE = int(args.upload)

if args.threads:
    const_NUMBER_OF_THREADS = args.threads

if args.steps:
    const_DO_STEPS = map(int, args.steps.split(','))

if args.mpi:
    const_USE_MPI = False

if args.trim:
    const_TRIM_ALIGNMENT = False

if args.extra_trimal:
    var_ADDITIONAL_TRIMAL_PARAMETER = args.extra_trimal

##################################
###         FUNCTIONS          ###
##################################

def readFastaFile(filename):
    global var_SEQUENCES
    global const_REPLACE_POINTS
    name = ''
    seq = ''
    f = open(filename, 'r')
    for line in f:
        if line[0] == '>':
            if name:
                var_SEQUENCES[name] = seq
            seq = ''
            tmp = line.strip()
            name = tmp[1:]
            if const_REPLACE_POINTS:
                name = name.replace('.', '___')
        else:
            seq += line.strip()
    if name:
        var_SEQUENCES[name] = seq
    f.close()
    return


def multipleAlignment(inputFile, outputFile):
    global path_TRANSLATORX
    os.system(
        'perl ' + path_TRANSLATORX + ' -i ' + inputFile + ' -o ' + outputFile + ' -p M')


def trimal(infile, outfile, params={'gt': '0.8'}, extraparams=''):
    global path_TRIMAL
    parameterString = ''
    for p in params:
        parameterString += ' -' + p + ' ' + params[p]
    os.system(path_TRIMAL + ' -in ' + infile + ' -out ' + outfile + parameterString + ' ' + extraparams)


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
        print('mpirun -np ' + str(
            threads) + ' ' + path_PHYML_MPI + ' --no_memory_check -i ' + infile + parameterString + ' ' + modelString)
        os.system('mpirun -np ' + str(
            threads) + ' ' + path_PHYML_MPI + ' --no_memory_check -i ' + infile + parameterString + ' ' + modelString)
    else:
        os.system(path_PHYML + ' --no_memory_check -i ' + infile + parameterString + ' ' + modelString)


def buildColorFile(filename, groupGenes):
    # natGeneList = [] # TODO: for gene expression
    # slyGeneList = []
    f = open(filename, 'w')
    f.write('TREE_COLORS\nSEPARATOR TAB\n\nDATA\n')
    for gene in sorted(groupGenes):
        f.write(gene + '\t' + const_UPLOAD_INFO_CODE[gene.split('|')[0][:3]] + '\n')
    f.write('\n')
    f.close()


def buildNatExpressionFile(filename, natGenes):
    global file_NIAT_EXPRESSION
    global var_NIAT_EXPRESSION_TISSUES
    expressionGenes = {}  # geneID <-> [cuffgene1, ...]
    f = open(file_NIAT_EXPRESSION, 'r')
    header = f.next().strip().split('\t')[1:]
    for line in f:
        values = line.strip().split()

        geneID = values[0]

        if not geneID in natGenes:
            continue
        geneID_original = natGenes[geneID]

        expressionGenes[geneID_original] = {}

        for i_expressionTissue in range(0, len(values)-1):
            if header[i_expressionTissue] in var_NIAT_EXPRESSION_TISSUES:
                expressionGenes[geneID_original][header[i_expressionTissue]] = math.log(float(values[i_expressionTissue+1]) + 1)
    f.close()

    f = open(filename, 'w')
    f.write('DATASET_HEATMAP\n\nSEPARATOR TAB\n\nDATASET_LABEL\tExpression: N. attenuata\n\nCOLOR\t#ff0000\n\nCOLOR_MIN\t#ffffff\nCOLOR_MAX\t#ff0000\n\n')
    f.write('FIELD_LABELS\t' + '\t'.join(var_NIAT_EXPRESSION_TISSUES))
    f.write('\n')
    f.write('\nDATA\n')
    for g in expressionGenes:
            tmp_EXP = []
            for tissue in var_NIAT_EXPRESSION_TISSUES:
                tmp_EXP.append(str(expressionGenes[g][tissue]))
            f.write(g + '\t' + '\t'.join(tmp_EXP) + '\n')
    f.close()


##################################
###            MAIN            ###
##################################

##### STEP 0: create sequence file #####
if const_DO_STEPS[0] == 1:
    if args.v >= 1: print 'WRITE SEQUENCE FILE'

    readFastaFile(file_INPUT_FASTA)

    f = open(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.fasta', 'w')
    for gene in sorted(var_SEQUENCES):
        f.write('>' + gene + '\n' + var_SEQUENCES[gene] + '\n')
    f.close()

    if not os.path.isfile(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.fasta'):
        if file_LOG:
            with open(file_LOG, "a") as myfile:
                myfile.write(
                    'ERROR' + '\t' + 'WRITE SEQUENCE FILE' + '\t' + folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.fasta' + '\t' + 'not exists' + '\n')

    if not folder_OUTPUT == folder_TMP_FILES:
        shutil.copy(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.fasta',
                    folder_OUTPUT + file_INPUT_FASTA_SHORT + '-sequences.fasta')

##### STEP 1: build multiple alignment, trim and convert to PHYLIP format #####
if const_DO_STEPS[1]:
    if args.v >= 1: print 'BUILT MULTIPLE ALIGNMENT'
    multipleAlignment(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.fasta',
                      folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align')

    if not os.path.isfile(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali.fasta'):
        if file_LOG:
            with open(file_LOG, "a") as myfile:
                myfile.write(
                    'ERROR' + '\t' + 'BUILD ALIGNMENT FILE' + '\t' + folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali.fasta' + '\t' + 'not exists' + '\n')

    if const_TRIM_ALIGNMENT:
        if args.v >= 1: print 'CLEAN SEQUENCES (trimal)'
        trimal(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali.fasta',
           folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.fasta', extraparams=var_ADDITIONAL_TRIMAL_PARAMETER)  # -cons 60
    else:
        shutil.copy(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali.fasta', folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.fasta')

    if not os.path.isfile(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.fasta'):
        if file_LOG:
            with open(file_LOG, "a") as myfile:
                myfile.write(
                    'ERROR' + '\t' + 'CLEAN ALIGNMENT FILE' + '\t' + folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.fasta' + '\t' + 'not exists' + '\n')

    if not folder_OUTPUT == folder_TMP_FILES:
        shutil.copy(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.fasta',
                    folder_OUTPUT + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.fasta')

    if args.v >= 1: print 'CONVER TO PHYLIP FORMAT'
    convertToPhylip(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.fasta',
                    folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.phy')

    if not os.path.isfile(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.phy'):
        if file_LOG:
            with open(file_LOG, "a") as myfile:
                myfile.write(
                    'ERROR' + '\t' + 'CONVERT TO PHYLIP FORMAT' + '\t' + folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.phy' + '\t' + 'not exists' + '\n')

##### STEP 2: calculate best substitution model ####
if const_DO_STEPS[2]:
    if args.v >= 1: print 'FIND BEST MODEL'
    jmodeltest(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.fasta',
               folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-bestModel.out',
               {'f': '', 'i': '', 'g': '4', 's': str(const_NUMBER_OF_SCHEMES), 'AIC': '', 'a': '',
                'tr': str(const_NUMBER_OF_THREADS)})

    if not os.path.isfile(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-bestModel.out'):
        if file_LOG:
            with open(file_LOG, "a") as myfile:
                myfile.write(
                    'ERROR' + '\t' + 'CALCULATE SUBSTITUTION MODEL' + '\t' + folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-bestModel.out' + '\t' + 'not exists' + '\n')

##### STEP 3: extract best substitution model
if const_DO_STEPS[3]:
    bestModelParameters = findBestModel(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-bestModel.out')

    if not folder_OUTPUT == folder_TMP_FILES:
        shutil.copy(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-bestModel.out', folder_OUTPUT + file_INPUT_FASTA_SHORT + '-bestModel.out')

##### STEP 4: build tree #####
if const_DO_STEPS[4]:
    if args.v: print 'BUILD TREE'
    phyml(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.phy', bestModelParameters,
          {'d': 'nt', 'b': str(const_BOOTSTRAP)}, const_NUMBER_OF_THREADS)

    # '-quiet': '',

    if not os.path.isfile(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.phy_phyml_tree.txt'):
        if file_LOG:
            with open(file_LOG, "a") as myfile:
                myfile.write(
                    'ERROR' + '\t' + 'BUILD TREE FILE' + '\t' + folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.phy_phyml_tree.txt' + '\t' + 'not exists' + '\n')

    if not folder_OUTPUT == folder_TMP_FILES:
        shutil.copy(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.phy_phyml_tree.txt',
                    folder_OUTPUT + file_INPUT_FASTA_SHORT + '-sequences_align.nt_ali_cleaned.phy_phyml_tree.txt')
        # tree file: -sequences_align.nt_ali_cleaned.phy_phyml_tree.txt

#### STEP 5: build color file #####
if const_DO_STEPS[5]:



    if args.v: print 'BUILD TREE'
    readFastaFile(file_INPUT_FASTA)

    buildColorFile(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.color.txt', var_SEQUENCES)

    if not os.path.isfile(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.color.txt'):
        if file_LOG:
            with open(file_LOG, "a") as myfile:
                myfile.write(
                    'ERROR' + '\t' + 'CREATE COLOR FILE' + '\t' + folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.color.txt' + '\t' + 'not exists' + '\n')

    if not folder_OUTPUT == folder_TMP_FILES:
        shutil.copy(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.color.txt',
                    folder_OUTPUT + file_INPUT_FASTA_SHORT + '-sequences.color.txt')

#### STEP 6: build expression file #####
if const_DO_STEPS[6]:
    if args.v: print 'BUILD TREE'
    readFastaFile(file_INPUT_FASTA)

    genes = {}

    for gene in var_SEQUENCES:
        if gene[:3] == 'nat':
            geneID_clean = gene.split('|')[1].split('.')[0].split('___')[0]

            genes[geneID_clean] = gene

    buildNatExpressionFile(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.nat-expression.txt', genes)

    if not os.path.isfile(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.nat-expression.txt'):
        if file_LOG:
            with open(file_LOG, "a") as myfile:
                myfile.write(
                    'ERROR' + '\t' + 'CREATE NAT EXPRESSION FILE' + '\t' + folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.nat-expression.txt' + '\t' + 'not exists' + '\n')

    if not folder_OUTPUT == folder_TMP_FILES:
        shutil.copy(folder_TMP_FILES + file_INPUT_FASTA_SHORT + '-sequences.nat-expression.txt',
                    folder_OUTPUT + file_INPUT_FASTA_SHORT + '-sequences.nat-expression.txt')





