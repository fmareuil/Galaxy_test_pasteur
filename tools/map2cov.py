#! /usr/bin/env python
# -*- coding: utf-8 -*-

#  This file is a part of COV2HTML software
#  COV2HTML <http://sourceforge.net/projects/cov2html/> 
#  A visualization and analysis tool of Bacterial NGS data for biologists.
# 
#  ---------------------------------------------------------------------- 
#  Copyright (C) 2012 Marc Monot
#  
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version. 
#   <http://www.gnu.org/licenses/gpl-3.0.html>
# 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#  ---------------------------------------------------------------------- 
#  
#  @author Marc Monot <marc.monot@pasteur.fr>
#  @author Mickael Orgeur


### PYTHON Module
# Base
import sys
import re
from optparse import OptionParser, OptionGroup

# Project
from _modules.functions import *

# Add from Python Package Index
from _modules.ordereddict import *


### MAIN
print
usage = """\
Usage: %prog [option] reference_file [option] alignment_file -t 'DNA|RNA|TSS|CHIP' {options}

Program: MAP2COV - Convert mapping data into coverage format readable by the COV2HTML interface
Version: 4.4
Contact: Marc Monot <marc.monot@pasteur.fr>

Reference format supported:
    GenBank, EMBL, FASTA/GFF3

Alignment format supported:
    SAM/BAM, ELAND, WIG

Data type supported:
    DNA, RNA, TSS, CHIP
        
Optional Parameters
    Anonymize, Read-length, Strand-specific, Coverage Style\
"""

getopt = OptionParser(usage=usage)

annotation = OptionGroup(getopt, 'Annotated reference sequence')
annotation.add_option('-g', '--genbank', dest='gbk', metavar='INPUT_FILE', help='GenBank annotation file')
annotation.add_option('-e', '--embl', dest='embl', metavar='INPUT_FILE', help='EMBL annotation file')
getopt.add_option_group(annotation)

prediction = OptionGroup(getopt, 'Reference Sequence + CDS Annotation', 'Both parameters are inseparable')
prediction.add_option('-f', '--fasta', dest='fasta', metavar='INPUT_FILE', help='Fasta reference sequence')
prediction.add_option('-c', '--cds', dest='cds', metavar='INPUT_FILE', help='Fasta/GFF3 CDS regions')
getopt.add_option_group(prediction)

alignment = OptionGroup(getopt, 'Alignment map')
alignment.add_option('-s', '--sam', dest='sam', metavar='INPUT_FILE', help='SAM/BAM alignment file')
alignment.add_option('-i', '--eland', dest='eland', metavar='INPUT_FILE', help='ELAND alignment file')
alignment.add_option('-w', '--wig', dest='wig', metavar='INPUT_FILE', help='WIG coverage file')
getopt.add_option_group(alignment)

datatype = OptionGroup(getopt, 'Data type')
datatype.add_option('-t', '--type', dest='datatype', metavar='DNA|RNA|TSS|CHIP', help='DNA, RNA, TSS or CHIP data')
getopt.add_option_group(datatype)

output = OptionGroup(getopt, 'Optional parameters')
output.add_option('-a', '--confidential', dest='secret', metavar='YES|NO', default='no', \
help='Export data confidentially (default: NO)')
output.add_option('-b', '--stranded', dest='strand', metavar='YES|NO', default='no', \
help='Perform a strand-specific analysis of data (default: NO)')
output.add_option('-l', '--read_length', dest='length', metavar='INT', default='20', \
help='Ignore reads shorter than fixed threshold (default: 20; 0 to disable)')
output.add_option('-v', '--cov_style', dest='style', metavar='MEAN|COUNTS', default='MEAN', \
help='Genetic Elements Coverage style: MEAN coverage or Raw COUNTS (default: MEAN)')
output.add_option('-o', '--output', dest='outname', metavar='OUTPUT_PREFIX', default='export_map2cov', \
help='Prefix for the output text file compressed in bz2 format (default: export_map2cov). The extension \'.txt.bz2\' is automatically added')
getopt.add_option_group(output)

options, arguments = getopt.parse_args()
# No input file provided
if options.gbk == None and options.embl == None and options.cds == None and options.fasta == None \
and options.sam == None and options.eland == None and options.wig == None:
    getopt.error('\tNo input file provided\n\t\t\tUse -h or --help for more details\n')
# Options -c and -f are inseparable
if (options.cds == None and options.fasta) or (options.cds and options.fasta == None):
    getopt.error('\tOptions \'-c\' and \'-f\' cannot be used seperately\n\t\t\tUse -h or --help for more details\n')
# Reference sequence test
if options.gbk:
    pass
elif options.embl:
    pass
elif options.fasta:
    pass
else:
    getopt.error('\tNo reference sequence provided\n\t\t\tUse -h or --help for more details\n')
# Alignment file test
if options.sam:
    pass
elif options.eland:
    pass
elif options.wig:
    pass
else:
    getopt.error('\tNo alignment file provided\n\t\t\tUse -h or --help for more details\n')
# Datatype test
try:
    type = options.datatype
    if type.lower() != 'dna' and type.lower() != 'rna' and type.lower() != 'tss' and type.lower() != 'chip':
        getopt.error("""\tDatatype provided is not recognized \'-t\': %s
        \t\tChoose between \'DNA\', \'RNA\', \'TSS\' and \'CHIP\'\n\t\t\tUse -h or --help for more details\n""" % type)
except AttributeError:
    getopt.error('\tNo type of data provided \'-t\'\n\t\t\tUse -h or --help for more details\n')
# Read-length threshold test
try:
    length = int(options.length)
    if length < 0:
        getopt.error('\tThe read-length threshold \'-l\' must be a positive integer\n\t\t\tUse -h or --help for more details\n')
except ValueError:
        getopt.error('\tThe read-length threshold \'-l\' must be a positive integer\n\t\t\tUse -h or --help for more details\n')
# WIG tests
if options.wig != None and options.style == 'COUNTS':
    getopt.error('\tThe Coverage Style \'Raw COUNTS\' cannot be used with WIG input file\n\t\t\tUse -h or --help for more details\n')
if options.strand.lower() != 'no' and options.wig != None:
    print 'WARNING: The option \'-b\' is ignored with WIG input file\n'
    options.strand = 'no'
if length != 20 and length != 0 and options.wig != None:
    print 'WARNING: The option \'l\' is ignored with WIG input file\n'
    length = 0
# Output test
if options.strand.lower() == 'no':
    if options.outname[-8:] != '.txt.bz2':
        outname = options.outname + '.txt.bz2'
    else:
        outname = options.outname
# Two files are generated for option '-b': one for the strand plus and one for the strand minus
else:
    if options.outname[-8:] != '.txt.bz2':
        outname_stranded = options.outname + '_stranded.txt.bz2'
    else:
        outname_stranded = options.outname[:-8] + '_stranded.txt.bz2'

### GENOMIC FILES: GBK, EMBL or FASTA
# GBK FORMAT
if options.gbk:
    try:
        filin = open(options.gbk, 'rU')
    except IOError:
        sys.exit('ERROR: No such file %s' % options.gbk)
    
    # Check GenBank format
    if checkGenbankFile(filin):
        sys.exit('ERROR: Annotation file is not in GenBank format')
    
    # CDS infos: name, CDS product, CDS first, CDS end, CDS way, CDS coverage (0)
    print('START Importing CDS Information')
    dict_cds, exon = cdsGenbankRecovery(filin, options.secret.lower())
    print('%i CDS regions extracted\n' % len(dict_cds))
    if exon > 0:
        print('Number of excluded CDS: %i\n' % exon)

    # Genome sequence forward
    print('START Importing Reference Sequence File')
    genome_forward = genomeGenbankRecovery(filin)
    if genome_forward == []:
        sys.exit('ERROR: No genome sequence in GenBank file')
    print('Length of the sequence: %i bp\n' % len(genome_forward))

    filin.close()

# EMBL FORMAT
elif options.embl:
    try:
        filin = open(options.embl, 'rU')
    except IOError:
        sys.exit('ERROR: No such file %s' % options.embl)
    
    # Check EMBL format
    if checkEmblFile(filin):
        sys.exit('ERROR: Annotation file is not in EMBL format')
    
    print('START Importing CDS Information')
    dict_cds, exon = cdsEmblRecovery(filin, options.secret.lower())
    print('%i CDS regions extracted\n' % len(dict_cds))
    if exon > 0:
        print('Number of excluded CDS: %i\n' % exon)
    
    # Genome sequence forward
    print('START Importing Reference Sequence')
    genome_forward = genomeEmblRecovery(filin)
    if genome_forward == []:
        sys.exit('ERROR: No genome sequence in EMBL file')
    print('Length of the sequence: %i bp\n' % len(genome_forward))
    
    filin.close()

# FASTA FORMAT
elif options.fasta:
    # Check genome and CDS files format (fasta)
    if checkFastaFile(options.fasta):
        sys.exit('ERROR: Reference file is not in Fasta format')
    filin = open(options.cds, 'rU')
    fileformat = checkFastaOrGFF(filin)
    if fileformat == 0:
        sys.exit('ERROR: CDS file is not in Fasta or GFF3 format')
    filin.close()
    
    # Genome sequence forward and reverse
    print('START Importing Genome File')
    genome_forward, genome_reverse = genomeFastaRecovery(options.fasta)
    print('Length of the sequence: %i bp\n' % len(genome_forward))
    
    # CDS infos: name, CDS product, CDS first, CDS end, CDS way, CDS coverage (plus/minus)
    print('START Importing CDS File')
    if fileformat == 'F':
        dict_cds = cdsFastaRecovery(options.cds, genome_forward, genome_reverse, options.secret.lower())
    if fileformat == 'G':
        dict_cds = cdsGFFRecovery(options.cds, genome_forward, genome_reverse, options.secret.lower())
    print('%i CDS regions extracted\n' % len(dict_cds))
    
### READING FILES: SAM/BAM, ELAND or WIG
# SAM/BAM FORMAT
if options.sam:
    # Check SAM/BAM alignment format
    format = checkSamFile(options.sam)
    if format == 1:
        sys.exit('ERROR: Alignment file is not in SAM/BAM format')
    
    print('START Reading SAM/BAM File')
    genome_coverage, reads, tss_coverage = readSamFile(options.sam, format, length, genome_forward, type.lower())
    if genome_coverage == 1:
        sys.exit('ERROR: Several reference-sequence names were detected in SAM/BAM file')
    elif genome_coverage == 2:
        sys.exit('ERROR: Length of reference sequence provided does not correspond to the mapping reference')
    elif genome_coverage == 3:
        sys.exit('ERROR: Mapping position out of the reference sequence')
    print('Number of mapped reads: %i\n' % reads)

# ELAND FORMAT
elif options.eland:
    # Check ELAND alignment format
    format = checkElandFile(options.eland)
    if format == 1:
        sys.exit('ERROR: Alignment file is not in ELAND format')

    print('START Reading ELAND File')
    genome_coverage, reads, tss_coverage = readElandFile(options.eland, format, length, genome_forward, type.lower())
    if genome_coverage == 1:
        sys.exit('ERROR: Mapping position out of the reference sequence')
    print('Number of mapped reads: %i\n' % reads)

# WIG FORMAT
elif options.wig:
    # Check WIG coverage format
    if checkWigFile(options.wig):
        sys.exit('ERROR: Coverage file is not in WIG format')
    
    print('START Reading WIG File')
    genome_coverage, positions, tss_coverage = readWigFile(options.wig, genome_forward)
    if genome_coverage == 1:
        sys.exit('ERROR: Coverage position out of the reference sequence')
    print('Number of coverage positions: %i\n' % positions)

### RESULTS FILE
# CDS coverage
print('START Creating CDS Coverage')
if options.style.lower() == 'mean':
    dict_cds = cdsCoverage(genome_coverage, dict_cds, type.lower(), options.style.lower())
if options.style.lower() == 'counts':
    dict_cds = cdsCoverage(tss_coverage, dict_cds, type.lower(), options.style.lower())
print('%i CDS coverages calculated\n' % len(dict_cds))

# Intergenic region coverage
print('START Creating Intergenic-Region Coverage')
if options.style.lower() == 'mean':
    dict_ig = intergenicCoverage(genome_coverage, dict_cds, options.style.lower())
if options.style.lower() == 'counts':
    dict_ig = intergenicCoverage(tss_coverage, dict_cds, options.style.lower())
print('%i intergenic-region coverages calculated\n' % len(dict_ig))

# Output editing and compression
print('START Writing Result File')
output = writeResults(genome_forward, genome_coverage, tss_coverage, dict_cds, dict_ig, type.lower(), options.strand.lower())
if len(output) == 1:
    compressResults(output[0], outname)
    print('New file created: %s\n' % outname)
else:
    compressResults(output[0] + "\n" + output[1], outname_stranded)
    print('New file created: %s\n' % outname_stranded)
