#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
 Title:     filter_blast_file.py
 Author:    Thomas Brockmöller
 Description:
     reduces a BLAST output file by a geneList
'''
import argparse

parser = argparse.ArgumentParser(
    description='Extracts blastHits according to cutoffs. If you find bugs please contact: tbrockmoeller@ice.mpg.de')
parser.add_argument('-bi', metavar='string', type=str, help="BLAST infile", required=True)
parser.add_argument('-bo', metavar='string', type=str, help="BLAST outfile (optional)", required=False)
parser.add_argument('-cl', metavar='string', type=str, help="cutoff percentage of length (optional) (default: 60)", required=False)
parser.add_argument('-ci', metavar='string', type=str, help="cutoff percentage of identity (optional) (default: 50)", required=False)
parser.add_argument('-cc', metavar='string', type=str, help="cutoff percentage of coverage (optional) (default: 60)", required=False)
parser.add_argument('-ce', metavar='string', type=str, help="cutoff e-value (optional) (default: 1e-20)",
                    required=False)
args = parser.parse_args()

if args.bi:
    blastInFile = args.bi

if args.bo:
    blastOutFile = args.bo
else:
    blastOutFile = None

# CUTOFFS
if args.cl:
    cutoff_length = float(args.cl)
else:
    cutoff_length = 60

if args.ci:
    cutoff_identity = float(args.ci)
else:
    cutoff_identity = 50

if args.cc:
    cutoff_coverage = float(args.cc)
else:
    cutoff_coverage = 60

if args.ce:
    cutoff_eVal = float(args.ce)
else:
    cutoff_eVal = float(1e-20)

fastaCdsARA = 'FASTA/ara.fasta'
fastaCdsCAN = 'FASTA/can.fasta'
fastaCdsCSA = 'FASTA/csa.fasta'
fastaCdsMGU = 'FASTA/mgu.fasta'
fastaCdsNAT = 'FASTA/nat.fasta'
fastaCdsNIO = 'FASTA/nio.fasta'
fastaCdsPTR = 'FASTA/ptr.fasta'
fastaCdsSLY = 'FASTA/sly.fasta'
fastaCdsSME = 'FASTA/sme.fasta'
fastaCdsSTU = 'FASTA/stu.fasta'
fastaCdsVVI = 'FASTA/vvi.fasta'
fastaCdsGMA = 'FASTA/gma.fasta'
fastaCdsALY = 'FASTA/aly.fasta'
fastaCdsCRU = 'FASTA/cru.fasta'

var_SEQUENCES = {}

def readFasta(filename, addition=''):
    global var_SEQUENCES
    name = ''
    seq = ''
    f = open(filename, 'r')
    for line in f:
        if line[0] == '>':
            var_SEQUENCES[name] = seq + addition
            seq = ''
            tmp = line.strip().split()[0]
            name = tmp[1:]
        else:
            seq += line.strip()
        var_SEQUENCES[name] = seq
    f.close()


readFasta(fastaCdsARA)
readFasta(fastaCdsCAN)
readFasta(fastaCdsCSA)
readFasta(fastaCdsMGU)
readFasta(fastaCdsNAT)
readFasta(fastaCdsNIO)
readFasta(fastaCdsPTR)
readFasta(fastaCdsSLY)
readFasta(fastaCdsSME)
readFasta(fastaCdsSTU)
readFasta(fastaCdsVVI)
readFasta(fastaCdsGMA)
readFasta(fastaCdsALY)
readFasta(fastaCdsCRU)

geneLength = {}
for i in sorted(var_SEQUENCES):
    geneLength[i] = len(var_SEQUENCES[i])

f = open(blastInFile, 'r')
if blastOutFile:
    fo = open(blastOutFile, 'w')
    fro = open(blastOutFile + '.rejected', 'w')
for line in f:
    (queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart,
     subjectEnd, eVal, bitScore) = line.strip().split('\t')

    eVal = float(eVal)
    percIdentity = float(percIdentity)
    alnLength = float(alnLength)
    queryEnd = float(queryEnd)
    queryStart = float(queryStart)
    subjectEnd = float(subjectEnd)
    subjectStart = float(subjectStart)
    queryCoverage = 100 * (queryEnd - queryStart) / geneLength[queryId]
    subjectCoverage = 100 * (subjectEnd - subjectStart) / geneLength[subjectId]

    if eVal < cutoff_eVal and percIdentity > cutoff_identity and alnLength > cutoff_length and queryCoverage > cutoff_coverage and subjectCoverage > cutoff_coverage:
        if blastOutFile:
            fo.write(line)
        else:
            print(line.strip())
    else:
        if blastOutFile:
            fro.write(line)
if blastOutFile:
    fo.close()
    fro.close()

f.close()
