#!/usr/bin/env python

from __future__ import print_function
import sys
from subprocess import call
import csv
import Bio
import re
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
import argparse
from argparse import RawTextHelpFormatter

# silence Biopython warnings
import warnings
from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonWarning)

# Make file write mode compatible for Python 2 and Python 3
writemode = 'wb' if sys.version_info < (3,) else 'w'

# script help and usage
parser=argparse.ArgumentParser(
    description='Given a tab-delimited file providing genomic coordinates of transposon insertions, this script\nparses a Genbank file and tabulates the total number of transposon insertions and insertion sites\nper locus tag. It also tabulates intergenic insertions and locus tags without any insertions. \n\nNOTE: Organisms with multiple Genbank records (e.g. those with multiple chromosomes or plasmids)\nshould be concatenated into a single .gbk file before executing this script. For example:\n% cat NC_000001.gbk NC_000002.gbk [...NC_00000n.gbk] > concatenated.gbk\n\nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)',
    epilog='Output columns for coding hits:\n1-SeqRecord ID\n2-feature type\n3-locus tag\n4-strand\n5-total transposon hits\n6-hits in TA sites\n7-hits in non-TA sites\n8-total # of hittable TA sites\n9-percentage of TA sites hit\n10-reads per kbp\n11-product\n\nAuthor: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/tnseq\nJanuary 2016\nupdated February 2016\n ', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file containing, at a minimum, locus tags and corresponding genomic\ncoordinates')
parser.add_argument('[DATA FILE]', help='insertion data as three-column tab-delimited file containing chromosome,\ngenomic coordinate, and total number of insertions at that coordinate, e.g.\nchr1\t327491\t1639')
parser.add_argument('[prefix]', help='output file names will be prepended with this prefix')
parser.add_argument('[N-terminal trim]', type=float, help='percent (expressed as a decimal) of the total gene length to trim from the\nN-terminus of the translated protein')
parser.add_argument('[C-terminal trim]', type=float, help='percent (expressed as a decimal) of the total gene length to trim from the\nC-terminus of the translated protein')
args=parser.parse_args()

insertionEvents=[]

print("reading insertion data...")

# read insertion data into list
with open(sys.argv[2], 'r') as f:
	insertionPointsFile = csv.reader(f,delimiter='\t')
	for chromosome,coordinate,insertion in insertionPointsFile:
		insertionEvents.append([chromosome, coordinate, insertion])

print("done")
print("tabulating insertions by locus tag...")

outputFileName = "%s.tabulated_insertions.coding.txt" % sys.argv[3]
outputFile = open(outputFileName, writemode)

noHitsFileName = "%s.nohits.txt" % sys.argv[3]
noHitsFile = open(noHitsFileName, writemode)

# define list of CDS coordinates for calculating intergenic regions
CDS_list = [[0,0]]

NtermTrim = float(sys.argv[4])
CtermTrim = float(sys.argv[5])

# parse input genbank file with SeqIO and loop over all sequence records
genbankFile = open(sys.argv[1], 'r')
for sequenceRecord in SeqIO.parse(genbankFile, "genbank"):
	genome = sequenceRecord.seq
	genomeSequence = str(sequenceRecord.seq)
	print("working on %s %s (%i bp, %2.2f%% GC)...\n" % (''.join(sequenceRecord.id), ''.join(sequenceRecord.description).rstrip('.').replace(', complete genome', ''), len(genome), GC(genome)))

	for feature in sequenceRecord.features:

		if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA' or feature.type == 'ncRNA':

                        locusTag = ''.join(feature.qualifiers["locus_tag"])
                        product = ''.join(feature.qualifiers["product"])
                        strand = int(feature.location.strand)
                        hits = []
                        if 'old_locus_tag' in feature.qualifiers:
                                oldlocusTag = ''.join(feature.qualifiers["old_locus_tag"])


                                # deal with strand-specific coordinate definitions
                                if strand == 1:
                                        strandSign = '+'
                                        startCoord = int(feature.location.start.position)
                                        endCoord = int(feature.location.end.position)
                                        # correct coordinates based on positional arguments
                                        realStartCoord = int(round(startCoord + ((endCoord - startCoord) * NtermTrim)))
                                        realEndCoord = int(round(endCoord - ((endCoord - startCoord) * CtermTrim)))
                                        correctedGeneLen = abs(realEndCoord - realStartCoord)
                                        geneSequence = str(genome[realStartCoord:realEndCoord])
                                        # collect uncorrected coordinates in CDS_list for determining intergenic regions later
                                        CDS_list.append([startCoord, endCoord])
                                        for insertionEvent in insertionEvents:
                                                chromosome = insertionEvent[0]
                                                # essential test to match coordinates of seqRecord with the same chromosome in the input hits file
                                                if chromosome == ''.join(sequenceRecord.id):
                                                        if realStartCoord <= int(insertionEvent[1]) <= realEndCoord:
                                                                hits.append(int(insertionEvent[2]))
                                                                coordinate = int(insertionEvent[1])

                                if strand == -1:
                                        strandSign = '-'
                                        startCoord = int(feature.location.end.position) #this is not a typo - biopython defines start and end coordinates from left to right regardless of strandedness
                                        endCoord = int(feature.location.start.position) #this is not a typo
                                        # correct coordinates based on positional arguments
                                        realStartCoord = int(round(startCoord - ((startCoord - endCoord) * NtermTrim)))
                                        realEndCoord = int(round(endCoord + ((startCoord - endCoord) * CtermTrim)))
                                        correctedGeneLen = abs(realEndCoord - realStartCoord)
                                        geneSequence = str(genome[realEndCoord:realStartCoord])
                                        	# collect uncorrected coordinates in CDS_list for determining intergenic regions later
                                        CDS_list.append([endCoord, startCoord])
                                        for insertionEvent in insertionEvents:
                                                chromosome = insertionEvent[0]
                                                # essential test to match coordinates of seqRecord with the same chromosome in the input hits file
                                                if chromosome == ''.join(sequenceRecord.id):
                                                        if realEndCoord <= int(insertionEvent[1]) <= realStartCoord:
                                                                hits.append(int(insertionEvent[2]))
                                                                coordinate = int(insertionEvent[1])
                                                                
                                # write output statistics if hits found
                                if hits:
                                        # calculate RpK (reads per kilobase of gene) based on corrected gene length
                                        readsPerKb = sum(hits) / (float(correctedGeneLen) / 1000)

                                        # classify output into either TA site hits or non-TA site hits
                                        if len(hits) > 0:
                                                codingHits = "%s\t%s\t%s\t%s\t%s\t%i\t%2.1f\t%s" % (''.join(sequenceRecord.id), ''.join(feature.type), locusTag, oldlocusTag, strandSign, sum(hits), readsPerKb, product)
                                        else:
                                                codingHits = "%s\t%s\t%s\t%s\t%s\t%i\t%2.1f\t%s" % (''.join(sequenceRecord.id), ''.join(feature.type), locusTag, oldlocusTag, strandSign, sum(hits), readsPerKb, product)
                                        outputFile.write(codingHits+"\n")


print("...done")
print("insertions in coding regions written to %s" % outputFileName)

noHitsFile.close()
outputFile.close()
genbankFile.close()
