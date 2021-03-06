#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
 Title:     detect_duplication_events_-_for_all_species.py
 Author:    Thomas Brockmöller

 Description:
     reads all gene duplication events from a folder containing Notung (v2.6) output files

 Limitations:
     developed for a specific set of species

 Example:
     detect_duplication_events_-_for_all_species.py -i folder_with_processed_NOTUNG_files/ -o recent_duplications.100pb.bootstrap.0.9.tsv -b 0.9 -p '.0' --exclude clusternames_with_wrong_alignments.txt

'''


from ete2 import PhyloTree
import re
import argparse
import os
import glob

parser = argparse.ArgumentParser(
    description='If you find bugs please contact: tbrockmoeller@ice.mpg.de')

parser.add_argument('-i', '--input', metavar='file', type=str,
                    help="input folder with Notung files",
                    required=True)
parser.add_argument('-o', '--output', metavar='file', type=str,
                    help="output file",
                    required=True)
parser.add_argument('-p', '--prefix', metavar='string', type=str,
                    help="gene tree file prefix (optional) (default: .nwk)",
                    required=False)
parser.add_argument('-b', '--bootstrap', metavar='float', type=str,
                    help="bootstrap value (default: 0.9)",
                    required=False)
parser.add_argument('-ds', '--duplication_score', metavar='float', type=float,
                    help="cutoff for species overlapping (optional) (default: 0)",
                    required=False)
parser.add_argument('-e', '--exclude', metavar='file', type=str,
                    help="exclude groups file",
                    required=False)
args = parser.parse_args()

folder_INPUT = args.input

file_OUTPUT = args.output

exclude_groups = []

if args.exclude:
    f = open(args.exclude)
    for line in f:
        exclude_groups.append(line.strip().split()[0])
    f.close()

var_FILE_PREFIX = '.nwk'

if args.prefix:
    var_FILE_PREFIX = args.prefix

var_DUPLICATION_SCORE = 0

if args.duplication_score:
    var_DUPLICATION_SCORE = float(args.duplication_score)

var_BOOTSTRAP = 0.9

if args.bootstrap:
    var_BOOTSTRAP = float(args.bootstrap)

var_SPECIES_LIST = ['ara', 'can', 'cca', 'csa', 'mgu', 'nat', 'nis', 'nit', 'nio', 'ptr', 'sme', 'sly', 'stu', 'vvi']
var_SELECTED_SPECIES_LIST = var_SPECIES_LIST


def remove_interal_nodes(tree):
    tree = re.sub(r'\)[n|r]\d+\:', '):', tree)
    tree = re.sub(r'\)[n|r]\d+\[', ')[', tree)
    return tree


def removeLost(genes):
    genes_new = set([])

    for gene in genes:
        if not gene[-4:] == 'LOST':
            genes_new.add(gene)
    return genes_new


def getSpecies(genes):
    species = set([])
    for gene in genes:
        if not gene[-4:] == 'LOST':
            spp = gene.split('|')[0]
            species.add(spp)
    return species


def returnDuplicationEventNewSpeciesTree(child1_spp, child2_spp):
    # Union... where is the duplication event (Vereinigung)
    speciesUnion = child1_spp.union(child2_spp)

    if {'ara'} == speciesUnion:
        return 'specific ARA'
    if {'can'} == speciesUnion:
        return 'specific CAN'
    if {'cca'} == speciesUnion:
        return 'specific CCA'
    if {'csa'} == speciesUnion:
        return 'specific CSA'
    if {'mgu'} == speciesUnion:
        return 'specific MGU'
    if {'nat'} == speciesUnion:
        return 'specific NAT'
    if {'nis'} == speciesUnion:
        return 'specific NIS'
    if {'nit'} == speciesUnion:
        return 'specific NIT'
    if {'nio'} == speciesUnion:
        return 'specific NIO'
    if {'ptr'} == speciesUnion:
        return 'specific PTR'
    if {'sme'} == speciesUnion:
        return 'specific SME'
    if {'sly'} == speciesUnion:
        return 'specific SLY'
    if {'stu'} == speciesUnion:
        return 'specific STU'
    if {'vvi'} == speciesUnion:
        return 'specific VVI'

    if {'ara', 'ptr'} == speciesUnion:
        return 'ARA-PTR'

    if ('csa' in speciesUnion and ('ptr' in speciesUnion or 'ara' in speciesUnion) and len(
            {'can', 'cca', 'mgu', 'nat', 'nio', 'nis', 'nit', 'sme', 'stu', 'sly', 'vvi'}.intersection(
                speciesUnion)) == 0):
        return 'CSA'

    if ('vvi' in speciesUnion and ('csa' in speciesUnion or 'ara' in speciesUnion or 'ptr' in speciesUnion) and len(
            {'can', 'cca', 'mgu', 'nat', 'nio', 'nis', 'nit', 'sme', 'stu', 'sly'}.intersection(
                speciesUnion)) == 0):
        return 'VVI'

    if len({'can', 'cca', 'mgu', 'nat', 'nio', 'nis', 'nit', 'sme', 'stu', 'sly'}.intersection(speciesUnion)) and len(
            {'vvi', 'ptr', 'csa', 'ara'}.intersection(speciesUnion)):
        return 'OLD'

    if ('mgu' in speciesUnion and len({'can', 'cca', 'nat', 'nio', 'nis', 'nit', 'sme', 'stu', 'sly'}.intersection(speciesUnion))):
        return 'MGU'

    if not ('nat' in speciesUnion or 'nio' in speciesUnion or 'nis' in speciesUnion or 'nit' in speciesUnion):
        if 'can' in speciesUnion:
            return 'CAN'
        if 'sme' in speciesUnion:
            return 'SME'
        if 'sly' in speciesUnion and 'stu' in speciesUnion:
            return 'SLY-STU'
        else:
            print('%s\t%s' % (','.join(sorted(child1_spp)), ','.join(sorted(child2_spp))))
            return 'ERROR!!!!!!!'

    if not ('can' in speciesUnion or 'sme' in speciesUnion or 'sly' in speciesUnion or 'stu' in speciesUnion):
        if not ('nio' in speciesUnion or 'nit' in speciesUnion):
            if 'nat' in speciesUnion and 'nis' in speciesUnion:
                return 'NAT-NIS'
            else:
                return 'ERROR!!!!!!!'
        elif not ('nat' in speciesUnion or 'nis' in speciesUnion):
            if 'nio' in speciesUnion and 'nit' in speciesUnion:
                return 'NIO-NIT'
            else:
                print('%s\t%s' % (','.join(sorted(child1_spp)), ','.join(sorted(child2_spp))))
                return 'ERROR!!!!!!!'
        else:
            return 'NICOTIANA'

    return 'SOLANACEAE'


######################################
### MAIN
######################################

fo = open(file_OUTPUT, 'w')
fo.write(
    'FILE\tDUPLICATION_ID\tSPECIES_1\tSPECIES_2\tDUPLICATION\tBOOTSTRAP\tBOOTSTRAP_CHILD1\tBOOTSTRAP_CHILD2\tOVERLAP\tOVERLAP_NUMBERS\tGENES_1\tGENES_2\tSURE\n')

i_duplication = 0

for file_INPUT in sorted(glob.glob(os.path.join(folder_INPUT, '*' + var_FILE_PREFIX))):
    file_INPUT_NAME = file_INPUT.split('/')[-1].split('-')[0]

    clustername = file_INPUT_NAME.split('-')[0]

    if clustername in exclude_groups:
        continue

    f = open(file_INPUT)
    var_NEXUS_TREE = f.next().strip()
    f.close()

    var_NEXUS_TREE_processed = remove_interal_nodes(var_NEXUS_TREE)

    t = PhyloTree(var_NEXUS_TREE_processed)

    duplications_list = {}
    duplications_dict = {}

    node_i = 0

    for node in t.traverse("preorder"):
        node_i += 1

        duplication_type_lost = 0  # 1: copy is lost, 2: both copies retain
        duplication_score = -1
        duplication_score_string = ''
        useThisNode = False
        nodeIsSure = False

        if len(node.get_children()) == 2:

            if hasattr(node, "D") and node.D == 'Y':
                child1 = node.get_children()[0]
                child2 = node.get_children()[1]

                CORRECT_CHILD_NODE = False

                # EXCLUDE THE NODE THAT IS NOT THE LOSS NODE
                while not CORRECT_CHILD_NODE:
                    if len(child1.get_children()) == 2:
                        child1_1 = child1.get_children()[0]
                        child1_1_genes = child1_1.get_leaf_names()
                        child1_1_genes = removeLost(child1_1_genes)

                        child1_2 = child1.get_children()[1]
                        child1_2_genes = child1_2.get_leaf_names()
                        child1_2_genes = removeLost(child1_2_genes)

                        if len(child1_1_genes) == 0:
                            child1 = child1_2
                            continue
                        if len(child1_2_genes) == 0:
                            child1 = child1_1
                            continue
                        CORRECT_CHILD_NODE = True
                    else:
                        CORRECT_CHILD_NODE = True

                CORRECT_CHILD_NODE = False

                # EXCLUDE THE NODE THAT IS NOT THE LOSS NODE
                while not CORRECT_CHILD_NODE:
                    if len(child2.get_children()) == 2:
                        child2_1 = child2.get_children()[0]
                        child2_1_genes = child2_1.get_leaf_names()
                        child2_1_genes = removeLost(child2_1_genes)

                        child2_2 = child2.get_children()[1]
                        child2_2_genes = child2_2.get_leaf_names()
                        child2_2_genes = removeLost(child2_2_genes)

                        if len(child2_1_genes) == 0:
                            child2 = child2_2
                            continue
                        if len(child2_2_genes) == 0:
                            child2 = child2_1
                            continue
                        CORRECT_CHILD_NODE = True
                    else:
                        CORRECT_CHILD_NODE = True

                child1_genes = child1.get_leaf_names()
                child2_genes = child2.get_leaf_names()

                child1_spp = getSpecies(child1_genes).intersection(var_SPECIES_LIST)
                child2_spp = getSpecies(child2_genes).intersection(var_SPECIES_LIST)

                child1_genes = removeLost(child1_genes)
                child2_genes = removeLost(child2_genes)

                bootstrap = 0

                speciesIntersection = child1_spp.intersection(child2_spp)
                speciesUnion = child1_spp.union(child2_spp)

                duplication_score = float(len(speciesIntersection)) / float(len(speciesUnion))
                duplication_score_string = str(len(speciesIntersection)) + ' / ' + str(len(speciesUnion))

                if speciesIntersection:
                    duplication_type_lost = 2
                else:
                    duplication_type_lost = 1

                # check bootstrap, if we can trust this node
                if hasattr(node, "B"):
                    bootstrap = node.B
                else:
                    bootstrap = node.support

                if hasattr(child1, "B"):
                    bootstrap_child1 = child1.B
                else:
                    bootstrap_child1 = child1.support

                if hasattr(child2, "B"):
                    bootstrap_child2 = child2.B
                else:
                    bootstrap_child2 = child2.support

                duplication = returnDuplicationEventNewSpeciesTree(child1_spp, child2_spp)

                if float(bootstrap) < var_BOOTSTRAP or float(bootstrap_child1) < var_BOOTSTRAP or float(
                        bootstrap_child2) < var_BOOTSTRAP:
                    if duplication_score < var_DUPLICATION_SCORE:
                        useThisNode = True
                        nodeIsSure = False
                    else:
                        useThisNode = True
                        nodeIsSure = True

                else:
                    useThisNode = True
                    nodeIsSure = True

                # list of gene of child 1 and child 2
                gene_new1 = []
                gene_new2 = []
                for gene in sorted(child1_genes):
                    gene_name = gene.split('|')[1].replace('___', '.')
                    if gene.split('|')[0] in var_SELECTED_SPECIES_LIST:
                        gene_new1.append(gene_name)
                for gene in sorted(child2_genes):
                    gene_name = gene.split('|')[1].replace('___', '.')
                    if gene.split('|')[0] in var_SELECTED_SPECIES_LIST:
                        gene_new2.append(gene_name)

                # write this duplication to the output file
                fo.write(file_INPUT_NAME + '\t' + str(i_duplication) + '\t')

                fo.write(','.join(sorted(child1_spp)))
                fo.write('\t')

                fo.write(','.join(sorted(child2_spp)))
                fo.write('\t')

                fo.write(duplication)
                fo.write('\t')
                fo.write(str(bootstrap))
                fo.write('\t')

                fo.write(str(bootstrap_child1))
                fo.write('\t')
                fo.write(str(bootstrap_child2))
                fo.write('\t')

                fo.write(str(duplication_score))
                fo.write('\t')
                fo.write(str(duplication_score_string))
                fo.write('\t')

                fo.write(','.join(sorted(gene_new1)))
                fo.write('\t')

                fo.write(','.join(sorted(gene_new2)))
                fo.write('\t')
                fo.write(str(nodeIsSure))

                fo.write('\n')

                # index for duplication
                i_duplication += 1

fo.close()
