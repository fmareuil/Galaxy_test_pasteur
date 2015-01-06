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
import os
import shlex
import subprocess
import re
import sys
import struct
import gzip
import bz2

# Add from Python Package Index
from ordereddict import *


### SHELL function    
def shellCmd(cmd):
    """Launch a command line to shell, waiting for execution."""
    return subprocess.Popen(shlex.split(cmd)).wait()


### FILE format control
def checkFastaFile(filin):
    """Check sequence Fasta file given by user"""
    first_line = 1
    try:
        for line in open(filin, 'rU'):
            # Test '>'
            if first_line :
                if line[0] != '>':
                    return 1
                else:
                    first_line = 0
                    continue
            # Test 1er base par ligne : '>ATGCN'
            base = changeCase(line[0])
            if base != 'A' and base != 'T' and base != 'C' and base != 'G' and base != 'N' and base != '>':
                return 1
        return 0
    except IOError:
        sys.exit('ERROR: No such file %s' % filin)

def checkFastaOrGFF(filin):
    """Check CDS Fasta or GFF3 file given by user"""
    line = filin.readline()
    # Test '##gff-version 3'
    if line[0:13] == '##gff-version':
        return 'G'
    # Test Fasta file
    elif line[0] == '>':
        return 'F'
    else:
        return 0

def checkGenbankFile(filin):
    """Check GenBank annotation file given by user"""
    line = filin.readline()
    # Test 'LOCUS'
    if line[0:5] != 'LOCUS':
        return 1
    else:
        return 0

def checkEmblFile(filin):
    """Check EMBL annotation file given by user"""
    line = filin.readline()
    # TEST 'ID'
    if line[0:2] != 'ID':
        return 1
    else:
        return 0

def checkSamFile(filin):
    """Check SAM/BAM alignment file given by user"""
    try:
        # BAM format
        try:
            bam = gzip.open(filin, 'rb')
            magic = bam.read(4)           
            if magic[0:3] == 'BAM':
                bam.close()
                return 'bam'
            else:
                return 1
        # SAM format
        except IOError:
            for line in open(filin, 'rU'):
                # Headers
                if line[0] == '@':
                    continue
                # Alignment line
                else:
                    if re.search('^[ !-?A-~]+\t[0-9]+\t[!-~]+\t[0-9]+\t[0-9]+', line):
                        return 'sam'
                    else:
                        return 1
    except IOError:
        sys.exit('ERROR: No such file %s' % filin)

def checkElandFile(filin):
    """Check ELAND alignment file given by user"""
    try:
        for line in open(filin ,'rU'):
            # Old ELAND format
            if re.search('^[!-~]+\t[ATGCNatgcn.]+', line):
                return 'old'
            # New ELAND format
            elif re.search('^[!-~]+\t([0-9]+\t){5}', line):
                return 'new'
            else:
                return 1
    except IOError:
        sys.exit('ERROR: No such file %s' % filin)

def checkWigFile(filin):
    """Check WIG coverage file given by user"""
    try:
        i = 0
        for line in open(filin, 'rU'):
            i += 1
            if line[0:5] == 'track':
                return 0
            elif i == 100000:
                return 1
    except IOError:
        sys.exit('ERROR: No such file %s' % filin)


### SEQUENCE parsing
def genomeGenbankRecovery(filin):
    """Get genome sequence from Genbank file"""
    genome_forward = []
    for line in filin:
        if re.search('^\s+[0-9]+\s[atgcnATGCN\s]+\n', line):
            line = re.sub(' ', '', line[10:-1])
            genome_forward += changeCase(line)
    return genome_forward

def cdsGenbankRecovery(filin, secret):
    """Get a dictionnary from GenBank file that contains:
    CDS name, CDS sequence, CDS product, CDS first, CDS end, CDS way, CDS coverage, CDS type."""
    cds = OrderedDict()
    line = filin.readline()
    cds_type = ''
    id = exon = 0
    while line != '':
        # CDS parsing
        if line[5:8] == 'CDS' or line[5:9] == 'rRNA' or line[5:9] == 'tRNA' or line[5:10] == 'ncRNA':
            cds_id = gene = ''
            cds_product = []
            inc = 1
            # Multi-exon gene detected
            if re.search('join', line) or re.search(',', line):
                exon += 1
                line = filin.readline()
                continue
            # Wrong characters identifed in GenBank file
            line = re.sub('[>|<]', '', line)
            # CDS, rRNA or tRNA on strand +
            if re.search('\s{5}CDS\s{13}[0-9]+\.\.[0-9]+\n', line):
                position = line[21:-1].split('..')
                strand = 1
                cds_type = 'CDS'
            elif re.search('\s{5}CDS\s{13}complement\([0-9]+\.\.[0-9]+\)\n', line):
                position = line[32:-2].split('..')
                strand = -1
                cds_type = 'CDS'
            elif re.search('\s{5}rRNA\s{12}[0-9]+\.\.[0-9]+\n', line):
                position = line[21:-1].split('..')
                strand = 1
                cds_type = 'rRNA'
            elif re.search('\s{5}rRNA\s{12}complement\([0-9]+\.\.[0-9]+\)\n', line):
                position = line[32:-2].split('..')
                strand = -1
                cds_type = 'rRNA'
            elif re.search('\s{5}tRNA\s{12}[0-9]+\.\.[0-9]+\n', line):
                position = line[21:-1].split('..')
                strand = 1
                cds_type = 'tRNA'
            elif re.search('\s{5}tRNA\s{12}complement\([0-9]+\.\.[0-9]+\)\n', line):
                position = line[32:-2].split('..')
                strand = -1
                cds_type = 'tRNA'
            elif re.search('\s{5}ncRNA\s{11}[0-9]+\.\.[0-9]+\n', line):
                position = line[21:-1].split('..')
                strand = 1
                cds_type = 'ncRNA'
            elif re.search('\s{5}ncRNA\s{11}complement\([0-9]+\.\.[0-9]+\)\n', line):
                position = line[32:-2].split('..')
                strand = -1
                cds_type = 'ncRNA'                
            else:
                inc = 0
            # Features parsing
            if inc == 1:
                line = filin.readline()
                while line[0:21] == ' ' * 21:
                    # Gene name
                    if line[21:26] == '/gene':
                        gene = line.split('"')[1]
                    # Locus tag
                    if line[21:31] == '/locus_tag':
                        cds_id = line.split('"')[1]
                    # CDS product
                    if line[21:29] == '/product':
                        # Closing quote on the first line
                        if line[-2] == '"':
                            cds_product += line[31:-2].split(' ')
                        # Research or cloding quote
                        else:
                            cds_product += line[31:-1].split(' ')
                            line = filin.readline()
                            while line[-2] != '"':
                                cds_product += line[21:-1].split(' ')
                                line = filin.readline()
                            cds_product += line[21:-2].split(' ')
                    line = filin.readline()
                # Data storing
                if cds_id == '':
                    cds_id = gene
                elif gene != '':
                    cds_id = cds_id + '|' + gene
                if secret != 'no':
                    id += 1
                    cds_id = 'CDS' + str(id)
                    cds_product = ''
                cds[cds_id] = [cds_product, '', int(position[0]), int(position[1]), strand, [0, 0], cds_type]
                gene = cds_product = ''
        # Sequence
        elif line[0:6] == 'ORIGIN':
            break
        else:
            line = filin.readline()
    return cds, exon

def genomeEmblRecovery(filin):
    """Get genome sequence from EMBL file"""
    genome_forward = []
    for line in filin:
        if re.search('^\s{5}[atgcnATGCN\s]+\s+[0-9]+\n', line):
            line = re.sub(' ', '', line[5:-10])
            genome_forward += changeCase(line)
    return genome_forward

def cdsEmblRecovery(filin, secret):
    """Get a dictionnary from EMBL file that contains:
    CDS name, CDS sequence, CDS product, CDS first, CDS end, CDS way, CDS coverage, CDS type."""
    cds = OrderedDict()
    line = filin.readline()
    cds_type = ''
    id = exon = 0
    while line != '':
        # CDS parsing
        if line[5:8] == 'CDS' or line[5:9] == 'rRNA' or line[5:9] == 'tRNA' or line[5:10] == 'ncRNA':
            cds_id = gene = ''
            cds_product = []
            inc = 1
            # Multi-exon gene detected
            if re.search('join', line) or re.search(',', line):
                exon += 1
                line = filin.readline()
                continue
            # Wrong characters identifed in EMBL file
            line = re.sub('[>|<]', '', line)
            # CDS, rRNA or tRNA on strand +
            if re.search('FT\s{3}CDS\s{13}[0-9]+\.\.[0-9]+\n', line):
                position = line[21:-1].split('..')
                strand = 1
                cds_type = 'CDS'
            elif re.search('FT\s{3}CDS\s{13}complement\([0-9]+\.\.[0-9]+\)\n', line):
                position = line[32:-2].split('..')
                strand = -1
                cds_type = 'CDS'
            elif re.search('FT\s{3}rRNA\s{12}[0-9]+\.\.[0-9]+\n', line):
                position = line[21:-1].split('..')
                strand = 1
                cds_type = 'rRNA'
            elif re.search('FT\s{3}rRNA\s{12}complement\([0-9]+\.\.[0-9]+\)\n', line):
                position = line[32:-2].split('..')
                strand = -1
                cds_type = 'rRNA'
            elif re.search('FT\s{3}tRNA\s{12}[0-9]+\.\.[0-9]+\n', line):
                position = line[21:-1].split('..')
                strand = 1
                cds_type = 'tRNA'
            elif re.search('FT\s{3}tRNA\s{12}complement\([0-9]+\.\.[0-9]+\)\n', line):
                position = line[32:-2].split('..')
                strand = -1
                cds_type = 'tRNA'
            elif re.search('FT\s{3}ncRNA\s{11}[0-9]+\.\.[0-9]+\n', line):
                position = line[21:-1].split('..')
                strand = 1
                cds_type = 'ncRNA'
            elif re.search('FT\s{3}ncRNA\s{11}complement\([0-9]+\.\.[0-9]+\)\n', line):
                position = line[32:-2].split('..')
                strand = -1
                cds_type = 'ncRNA'
            else:
                inc = 0
            # Features parsing
            if inc == 1:
                line = filin.readline()
                while line[2:21] == ' ' * 19 and line[0:2] != 'SQ':
                    # Gene name
                    if line[21:26] == '/gene':
                        gene = line.split('"')[1]
                    # Locus tag
                    if line[21:31] == '/locus_tag':
                        cds_id = line.split('"')[1]
                    # CDS product
                    if line[21:29] == '/product':
                        # Closing quote on the first line
                        if line[-2] == '"':
                            cds_product += line[31:-2].split(' ')
                        # Research or cloding quote
                        else:
                            cds_product += line[31:-1].split(' ')
                            line = filin.readline()
                            while line[-2] != '"':
                                cds_product += line[21:-1].split(' ')
                                line = filin.readline()
                            cds_product += line[21:-2].split(' ')
                    line = filin.readline()
                # Data storing
                if cds_id == '':
                    cds_id == gene
                elif gene != '':
                    cds_id = cds_id + '|' + gene
                if secret != 'no':
                    id += 1
                    cds_id = 'CDS' + str(id)
                    cds_product = ''
                cds[cds_id] = [cds_product, '', int(position[0]), int(position[1]), strand, [0, 0], cds_type]
                gene = cds_product = ''
        # Sequence
        elif line[0:2] == 'SQ':
            break
        else:
            line = filin.readline()
    return cds, exon

def genomeFastaRecovery(filin):
    """Get genome sequence from fasta file"""
    genome_forward = ''.join(map(str,[changeCase(line).replace(' ', '').split('\r')[0].split('\n')[0] for line in open(filin) if line[0] != '>']))
    return genome_forward, reverseComplement(genome_forward)

def cdsFastaRecovery(filin, genome_forward, genome_reverse, secret):
    """Get a dictionnary from CDS fasta file that contains :
    CDS name, CDS sequence, CDS product, CDS first, CDS end, CDS way, CDS coverage, CDS type."""
    cds = OrderedDict()
    id = 0
    # Stock in CDS : Product, sequence
    for line in open(filin, 'rU'):
        if line[0] == '>':
            cds_info = line.split('>')[1].split()
            cds_id = cds_info[0]
            if len(cds_info) > 1: 
                cds_product = cds_info[1:]
            else:
                cds_product = ''
            if secret != 'no':
                id += 1
                cds_id = 'CDS' + str(id)
                cds_product = ''
            cds[cds_id] = [cds_product, '', 0, 0, 1, [0, 0], 'CDS']
        else:
            cds[cds_id][1] += changeCase(line.split('\r')[0].split('\n')[0])            
    # Stock in CDS : first, end, sens and remove sequence
    for cds_id in cds:
        # Find in forward
        position = genome_forward.find(cds[cds_id][1])
        # If in forward
        if position != -1:
            # Duplication Mask : change 10 first bases of match by N
            genome_forward = genome_forward[:position] + 'N'*10 + genome_forward[position+10:]
        # If not in forward, find in reverse
        else:
            position_reverse = genome_reverse.find(cds[cds_id][1])
            # If no match pass the CDS
            if position_reverse == -1:
                # Remove no match gene from dictionnary
                cds.pop(cds_id)
                continue
            # If in reverse
            else:
                position = len(genome_reverse) - len(cds[cds_id][1]) - position_reverse
                cds[cds_id][4] = -1
                # Duplication Mask : change 10 last bases of match by N                
                genome_reverse = genome_reverse[:position_reverse+len(cds[cds_id][1])-10] + 'N'*10 + genome_reverse[position_reverse+len(cds[cds_id][1]):]
        # Ajout de 1 (+1 Biologique)        
        cds[cds_id][2] = position + 1
        cds[cds_id][3] = position + len(cds[cds_id][1])
    return cds

def cdsGFFRecovery(filin, genome_forward, genome_reverse, secret):
    """Get a dictionnary from CDS GFF3 file that contains :
    CDS name, CDS sequence, CDS product, CDS first, CDS end, CDS way, CDS coverage, CDS type."""
    cds = OrderedDict()
    id = 0
    # Stock in CDS : Product, sequence
    for line in open(filin, 'rU'):
        if line[0] != '#':
            cds_info = line.split()
            if cds_info[2] == 'CDS' or cds_info[2] == 'rRNA' or cds_info[2] == 'tRNA' or cds_info[2] == 'ncRNA':
                if len(cds_info[8].split(';')) == 1:
                    cds_id = cds_info[8].split(';')[0].split('=')[1]
                    cds_product = ''
                else:
                    cds_id = cds_info[8].split(';')[0].split('=')[1]
                    cds_product = cds_info[8].split(';')[1:] + cds_info[9:]
                if cds_info[6] == '-':
                    strand = -1
                else:
                    strand = 1
                cds[cds_id] = [cds_product, '', int(cds_info[3]), int(cds_info[4]), strand, [0, 0], cds_info[2]]
        if line[0:7] == '##FASTA':
            return cds
    return cds


### SEQUENCE functions
def changeCase(seq):
    """Change case to capital of a sequence string."""
    return seq.replace('a','A').replace('t','T').replace('c','C').replace('g','G').replace('n','N')

def reverseComplement(seq):
    """Reverse Complement a fonction, must be in capital letter."""
    return seq.replace('A','t').replace('T','a').replace('C','g').replace('G','c').replace('a','A').replace('t','T').replace('c','C').replace('g','G')[::-1]


### ALIGNMENT functions
def readSamFile(filin, format, read_threshold, genome_forward, datatype):
    """Convert SAM / BAM information into coverage, do not use SNP information"""
    genome_coverage = ['plus', 'minus']
    genome_coverage[0] = [0]*len(genome_forward)
    genome_coverage[1] = [0]*len(genome_forward)
    tss_coverage = ['plus', 'minus']
    tss_coverage[0] = [0]*len(genome_forward)
    tss_coverage[1] = [0]*len(genome_forward)
    reads = 0
    
    # BAM format
    refID_multi_select = 0
    refID_select = 0
    m_SQ_l = []
    m_SQ_n = []
    if format == 'bam':
        bam = gzip.open(filin, 'rb')
        # Headers
        magic = bam.read(4)
        ltext = bam.read(4)
        ltext = struct.unpack('i', ltext)       
        text = bam.read(ltext[0])
        n_ref = bam.read(4)
        n_ref = struct.unpack('i', n_ref)        
        for i in range(n_ref[0]):
            # Reference sequence
            l_name = bam.read(4)
            l_name = struct.unpack('i', l_name)
            name = bam.read(l_name[0])
            m_SQ_n.append(name)
            lref = bam.read(4)       
            lref = struct.unpack('i', lref)
            m_SQ_l.append(lref[0])
        # Multi-BAM : automatic use of the right mapping reference sequence  
        for i in range(len(m_SQ_l)):
            if m_SQ_l[i] == len(genome_forward):
                name = m_SQ_n[i]
                lref = m_SQ_l[i]
                refID_select = i
                refID_multi_select += 1
        if lref != len(genome_forward) or refID_multi_select > 1:  # Length of mapping sequence != Length of reference sequence
            return 2, 0, 0
        # Alignments
        block_size = bam.read(4)
        n_out = 0
        while block_size != '':
            # Data extraction
            block_size = struct.unpack('i', block_size)
            refID = bam.read(4)
            refID = struct.unpack('i',refID)            
            pos = bam.read(4)
            pos = struct.unpack('i', pos)[0]
            skip = bam.read(4)
            # Parsing of bitwise flag
            flag_nc = bam.read(4)
            flag_nc = struct.unpack('I', flag_nc)
            flag = (flag_nc[0] & 0xffff0000) >> 16
            flag = str(bin(flag))[2:].zfill(11)
            l_seq = bam.read(4)
            l_seq = struct.unpack('i', l_seq)[0]
            next_block = bam.read(block_size[0]-(5*4))
            block_size = bam.read(4)
            # Test RefID
            if refID_select != refID[0]:
                continue
            # Threshold of read length fixed with options '-i'
            if l_seq < read_threshold:
                continue
            # Ignore unmapped reads (position is set at -1, flag is set at 1)
            if pos == -1 or flag[-3] == '1':
                continue
            last = pos + l_seq
            # Test position out of the sequence
            if pos > len(genome_forward):
                pos = len(genome_forward)
                return 3, 0, 0
            if last > len(genome_forward):
                last = len(genome_forward)
            reads += 1
            # Genome coverage (Matching positions are set from 0)
            for i in range(pos, last):
                # Strand plus
                if flag[-5] == '0':
                    genome_coverage[0][i] += 1
                # Strand minus: flag is set at 1
                else:
                    genome_coverage[1][i] += 1
            # TSS coverage (Matching positions are set from 0)
            # Strand plus
            if flag[-5] == '0':
                tss_coverage[0][pos] += 1
            # Strand minus: flag is set at 1
            else:
                tss_coverage[1][last-1] += 1
        return genome_coverage, reads, tss_coverage
    
    # SAM format
    if format == 'sam':
        refID_multi_select = 0
        inc = 0
        multi = 0
        m_SQ_n = []
        m_SQ_l = []
        for line in open(filin, 'rU'):
            # Headers
            if line[0] == '@':
                # Header of mapping reference sequence
                if line[0:3] == '@SQ':
                    inc += 1
                    header = line[:-1].split('\t')
                    for ele in header:
                        if ele[0:3] == 'SN:':
                            n_seq = ele[3:]
                            m_SQ_n.append(n_seq)
                        if ele[0:3] == 'LN:':
                            l_seq = int(ele[3:])
                            m_SQ_l.append(l_seq)
                    # Multi-SAM : automatic use of the right mapping reference sequence         
                    if inc > 1:
                        for i in range(len(m_SQ_l)):
                            if m_SQ_l[i] == len(genome_forward):
                               n_seq = m_SQ_n[i]
                               l_seq = m_SQ_l[i]
                               multi = 1
                               refID_multi_select += 1
                               continue
                        # Length of mapping sequence != Length of reference sequence
                        if l_seq != len(genome_forward) or refID_multi_select > 1:
                            return 2, 0, 0
                            continue
            else:
                line = str(line).split('\t')
                flag = str(bin(int(line[1])))[2:].zfill(11)
                name_match = str(line[2])
                if multi and name_match != n_seq:
                    continue
                match_position = int(line[3])
                read_length = len(line[9])
                # Threshold of read length fixed with options '-i'
                if read_length < read_threshold:
                    continue
                first = match_position
                last = match_position + read_length
                # Ignore unmapped reads (position is set at 0, flag is set at 1)
                if first == 0 or flag[-3] == '1':
                    continue
                # Test position out of the sequence
                if pos > len(genome_forward):
                    pos = len(genome_forward)
                    return 3, 0, 0
                if (last -1) > len(genome_forward):
                    last = len(genome_forward) - 1
                reads += 1
                # Genome coverage (Matching position are set from 1)
                for i in range(first, last):
                    # Strand plus
                    if flag[-5] == '0':
                        genome_coverage[0][i-1] += 1
                    # Strand minus: flag is set at 1
                    else:
                        genome_coverage[1][i-1] += 1
                # TSS coverage (Matching position are set from 1)
                # Strand plus
                if flag[-5] == '0':
                    tss_coverage[0][first-1] += 1
                # Strand minus: flag is set at 1
                else:
                    tss_coverage[1][last-2] += 1

        return genome_coverage, reads, tss_coverage

def readElandFile(filin, format, read_threshold, genome_forward, datatype):
    """Convert ELAND information into coverage, do not use SNP information"""
    genome_coverage = ['plus', 'minus']
    genome_coverage[0] = [0]*len(genome_forward)
    genome_coverage[1] = [0]*len(genome_forward)
    tss_coverage = ['plus', 'minus']
    tss_coverage[0] = [0]*len(genome_forward)
    tss_coverage[1] = [0]*len(genome_forward)
    reads = 0

    for line in open(filin, 'rU'):
        line = line.split('\t')
        # Old ELAND format
        if format == 'old':
            read_length = len(line[1])
            try:
                match_position = int(line[7])
                strand = line[8]
            # Ignore unmapped reads
            except IndexError:
                continue
        # New ELAND format
        else:
            read_length = len(line[8])
            if line[12] != '' and line[13] != '':
                match_position = int(line[12])
                strand = line[13]
            # Ignore unmapped reads
            else:
                continue
        # Threshold of read length fixed with options '-i'
        if read_length < read_threshold:
            continue
        first = match_position
        last = match_position + read_length        
        # Test positon out of the genome
        if first < 1:
            return 1, 0, 0
        if last > len(genome_forward):
            last = len(genome_forward)
        reads += 1
        # Genome coverage (Matching position are set from 1)
        for i in range(first, last):
            # Strand plus
            if strand == 'F':
                genome_coverage[0][i-1] += 1
            # Strand minus: line[13] is set at R
            else:
                genome_coverage[1][i-1] += 1
        # TSS coverage (Matching position are set from 1)
        # Strand plus
        if strand == 'F':
            tss_coverage[0][first-1] += 1
        # Strand minus: line[13] is set at R
        else:
            tss_coverage[1][last-2] += 1
    
    return genome_coverage, reads, tss_coverage

def readWigFile(filin, genome_forward):
    """Convert WIG information into coverage"""
    genome_coverage = ['plus', 'minus']
    genome_coverage[0] = [0]*len(genome_forward)
    genome_coverage[1] = [0]*len(genome_forward)
    positions = 0
    tss_coverage = ['plus', 'minus']

    for line in open(filin, 'rU'):
        if re.search('^[0-9]+\t[0-9]+', line):
            positions += 1
            data = line[:-1].split('\t')[0:2]
            position = int(data[0])
            coverage = int(data[1])
            # Test positon out of the genome
            if position < 1 or position > len(genome_forward):
                return 1, 0, 0
            genome_coverage[0][position-1] = coverage
    return genome_coverage, positions, tss_coverage


### COVERAGE functions
def cdsCoverage(genome_coverage, dict_cds, datatype, coverage):
    """Return Mean Coverage or Raw Counts for each CDS, or their promotor regions for tss and chip"""
    genome_coverage = [map(int, genome_coverage[0]), map(int, genome_coverage[1])]
    # CDS coverage is calculated from genome coverage on the entire gene
    if datatype != 'tss' and datatype != 'chip':
        for cds_id in dict_cds:
            # Strand plus
            plus = sum(genome_coverage[0][int(dict_cds[cds_id][2]-1):int(dict_cds[cds_id][3])])
            if coverage == 'mean':
                dict_cds[cds_id][5][0] = float(plus) / len(genome_coverage[0][int(dict_cds[cds_id][2]-1):int(dict_cds[cds_id][3])])
            elif coverage == 'counts':
                dict_cds[cds_id][5][0] = float(plus)
            # Strand minus
            minus = sum(genome_coverage[1][int(dict_cds[cds_id][2]-1):int(dict_cds[cds_id][3])])
            if coverage == 'mean':
                dict_cds[cds_id][5][1] = float(minus) / len(genome_coverage[1][int(dict_cds[cds_id][2]-1):int(dict_cds[cds_id][3])])
            elif coverage == 'counts':
                dict_cds[cds_id][5][1] = float(minus)
        return dict_cds
    # CDS coverage is calculated from genome coverage on the region [-250:ATG:+100]
    else:
        for cds_id in dict_cds:
            # Strand plus
            if int(dict_cds[cds_id][4]) == 1:
                start = int(dict_cds[cds_id][2]) - 250
                # Test position out of the first base
                if start < 1:
                    start = 1
                stop = int(dict_cds[cds_id][2]) + 2 + 100
                # Test position out of the last base
                if stop > len(genome_coverage[0]):
                    stop = len(genome_coverage[0])
                plus = sum(genome_coverage[0][start-1:stop])
                if coverage == 'mean':
                    dict_cds[cds_id][5][0] = float(plus) / len(genome_coverage[0][start-1:stop])
                elif coverage == 'counts':
                    dict_cds[cds_id][5][0] = float(plus)
                minus = sum(genome_coverage[1][start-1:stop])
                if coverage == 'mean':
                    dict_cds[cds_id][5][1] = float(minus) / len(genome_coverage[1][start-1:stop])
                elif coverage == 'counts':
                    dict_cds[cds_id][5][1] = float(minus)
            # Strand minus: strand is set at -1
            else:
                start = int(dict_cds[cds_id][3]) + 250
                # Test position out of the last base
                if start > len(genome_coverage[0]):
                    start = len(genome_coverage[0])
                stop = int(dict_cds[cds_id][3]) - 2 - 100
                # Test position out of the first base
                if stop < 1:
                    stop = 1
                plus = sum(genome_coverage[0][stop-1:start])
                if coverage == 'mean':
                    dict_cds[cds_id][5][0] = float(plus) / len(genome_coverage[0][stop-1:start])
                elif coverage == 'counts':
                    dict_cds[cds_id][5][0] = float(plus)
                minus = sum(genome_coverage[1][stop-1:start])
                if coverage == 'mean':
                    dict_cds[cds_id][5][1] = float(minus) / len(genome_coverage[1][stop-1:start])
                elif coverage == 'counts':
                    dict_cds[cds_id][5][1] = float(minus)
        return dict_cds

def intergenicCoverage(genome_coverage, dict_cds, coverage):
    """Return Mean Coverage or Raw Counts for each intergenic region"""
    genome_coverage = [map(int, genome_coverage[0]), map(int, genome_coverage[1])]
    dict_ig = OrderedDict()
    pos = 1
    id = 1
    cds_id = dict_cds.keys()
    for i in range(len(cds_id)):
        # CDS start position
        start = int(dict_cds[cds_id[i]][2]) - 1
        if start >= pos:
            cov = ['plus', 'minus']
            # Strand plus
            plus = sum(genome_coverage[0][pos-1:start])
            if coverage == 'mean':
                cov[0] = float(plus) / len(genome_coverage[0][pos-1:start])
            elif coverage == 'counts':
                cov[0] = float(plus)
            # Strand minus
            if genome_coverage[1] != []:
                minus = sum(genome_coverage[1][pos-1:start])
                if coverage == 'mean':
                    cov[1] = float(minus) / len(genome_coverage[1][pos-1:start])
                elif coverage == 'counts':
                    cov[1] = float(minus)
            else:
                cov[1] = 0
            # Intergenic region at the beginning of the sequence
            if i == 0:
                dict_ig['IG%i'%id] = ['begin', cds_id[i], pos, start, cov]
            # Intergenic regions
            else:
                dict_ig['IG%i'%id] = [cds_id[i-1], cds_id[i], pos, start, cov]
            id += 1
        pos = int(dict_cds[cds_id[i]][3]) + 1
        if i == len(cds_id)-1:
            # Intergenic region at the end of the sequence
            if len(genome_coverage[0])-1 >= pos:
                cov = ['plus', 'minus']
                # Strand plus
                plus = sum(genome_coverage[0][pos-1:len(genome_coverage[0])])
                if coverage == 'mean':
                    cov[0] = float(plus) / len(genome_coverage[0][pos-1:len(genome_coverage[0])])
                elif coverage == 'counts':
                    cov[0] = float(plus)
                # Strand minus
                minus = sum(genome_coverage[1][pos-1:len(genome_coverage[1])])
                if coverage == 'mean':
                    cov[1] = float(minus) / len(genome_coverage[1][pos-1:len(genome_coverage[1])])
                elif coverage == 'counts':
                    cov[1] = float(minus)
                dict_ig['IG%i'%id] = [cds_id[i], 'end', pos, len(genome_coverage[0])-1, cov]
    return dict_ig


### EXPORT functions
def writeResults(genome_forward, strand_coverage, tss_coverage, dict_cds, dict_ig, datatype, strand):
    """Export analysed data in output file : 
    (All) Coverage split in 10,000 fasta fragment
    (All) CDS information : name, product, first, end, way, coverage, type"""

    # Strand-specific analysis
    if strand != 'no':
        output = ['plus', 'minus']
        output[0] = ('@Type=' + str(datatype.upper()) + '\n' + '@Strand=+\n')
        output[1] = ('@Type=' + str(datatype.upper()) + '\n' + '@Strand=-\n')
        # Genome coverage
        output[0] += exportDataSplit('Coverage', strand_coverage[0], 10000)
        output[1] += exportDataSplit('Coverage', strand_coverage[1], 10000)
        # TSS coverage
        if datatype == 'tss':
            if tss_coverage[0] == 'plus':
                tss_coverage[0] = tss_coverage[1] = [0] * len(genome_forward)
            output[0] += exportDataSplit('Start', tss_coverage[0], 10000)
            output[1] += exportDataSplit('Start', tss_coverage[1], 10000)
        # CDS coverage
        output[0] += '@CDS\n'
        output[1] += '@CDS\n'
        for cds_id in dict_cds:
            cov = int(round(dict_cds[cds_id][5][0]))
            if cov == 0: cov = 1    # Uncovered genes are set at 1
            output[0] += ('>' + str(cds_id) + ' ' + str(dict_cds[cds_id][2]) + ' ' + str(dict_cds[cds_id][3]) + ' ' + str(dict_cds[cds_id][4])  + ' ' + str(cov) + ' ' + str(dict_cds[cds_id][6]) + '\n' + ' '.join(map(str, dict_cds[cds_id][0])) + '\n')
            cov = int(round(dict_cds[cds_id][5][1]))
            if cov == 0: cov = 1    # Uncovered genes are set at 1
            output[1] += ('>' + str(cds_id) + ' ' + str(dict_cds[cds_id][2]) + ' ' + str(dict_cds[cds_id][3]) + ' ' + str(dict_cds[cds_id][4]) + ' ' + str(cov) + ' ' + str(dict_cds[cds_id][6]) + '\n' + ' '.join(map(str, dict_cds[cds_id][0])) + '\n')
        # Intergenic-region coverage
        output[0] += '@IG\n'
        output[1] += '@IG\n'
        for ig_id in dict_ig:
            cov = int(round(dict_ig[ig_id][4][0]))
            if cov == 0: cov = 1    # Uncovered intergenic regions are set at 1
            output[0] += ('>' + str(ig_id) + ' ' + str(dict_ig[ig_id][2]) + ' ' + str(dict_ig[ig_id][3]) + ' ' + str(cov) + '\n' + str(dict_ig[ig_id][0]) + ' ' + str(dict_ig[ig_id][1]) + '\n')
            cov = int(round(dict_ig[ig_id][4][1]))
            if cov == 0: cov = 1    # Uncovered intergenic regions are set at 1
            output[1] += ('>' + str(ig_id) + ' ' + str(dict_ig[ig_id][2]) + ' ' + str(dict_ig[ig_id][3]) + ' ' + str(cov) + '\n' + str(dict_ig[ig_id][0]) + ' ' + str(dict_ig[ig_id][1]) + '\n')

    # Not strand-specific analysis
    else:
        output = ['both']
        output[0] = ('@Type=' + str(datatype.upper()) + '\n')
        # Genome coverage
        genome_coverage = [0]*len(strand_coverage[0])
        for i in range(len(strand_coverage[0])):
            genome_coverage[i] = strand_coverage[0][i] + strand_coverage[1][i]
        output[0] += exportDataSplit('Coverage', genome_coverage, 10000)
        # TSS coverage
        if datatype == 'tss':
            if tss_coverage[0] == 'plus':
                tss_coverage[0] = tss_coverage[1] = [0] * len(genome_forward)        
            genome_coverage = [0]*len(tss_coverage[0])
            for i in range(len(tss_coverage[0])):
                genome_coverage[i] = tss_coverage[0][i] + tss_coverage[1][i]
            output[0] += exportDataSplit('Start', genome_coverage, 10000)
        # CDS coverage
        output[0] += '@CDS\n'
        for cds_id in dict_cds:
            if int(round(dict_cds[cds_id][5][0])) == 0 and int(round(dict_cds[cds_id][5][1])) == 0:
                cov = 1    # Uncovered genes are set at 1
            else:
                cov = int(round(dict_cds[cds_id][5][0] + dict_cds[cds_id][5][1]))
            output[0] += ('>' + str(cds_id) + ' ' + str(dict_cds[cds_id][2]) + ' ' + str(dict_cds[cds_id][3]) + ' ' + str(dict_cds[cds_id][4]) + ' ' + str(cov) + ' ' + str(dict_cds[cds_id][6]) + '\n' + ' '.join(map(str, dict_cds[cds_id][0])) + '\n')
        # Intergenic-region coverage
        output[0] += '@IG\n'
        for ig_id in dict_ig:
            if int(round(dict_ig[ig_id][4][0])) == 0 and int(round(dict_ig[ig_id][4][1])) == 0:
                cov = 1    # Uncovered intergenic regions are set at 1
            else:
                cov = int(round(dict_ig[ig_id][4][0] + dict_ig[ig_id][4][1]))
            output[0] += ('>' + str(ig_id) + ' ' + str(dict_ig[ig_id][2]) + ' ' + str(dict_ig[ig_id][3]) + ' ' + str(cov) + '\n' + str(dict_ig[ig_id][0]) + ' ' + str(dict_ig[ig_id][1]) + '\n')
        
    return output

def exportDataSplit(datatype, data, split):
    """Return data in export format :
    @datatype
    >data1
    >data2 
    cut in split length."""
    sep = ' '
    text = '@' + str(datatype) + '\n'
    for i in range((len(data)/split)+1):
        text += sep.join(map(str,data[i*split:(i+1)*split])) + '\n'
    return text

def compressResults(source, outname):
    """Compress the output file using a bz2 compressor object."""
    data = bz2.compress(source)
    filout = file(outname, 'wb')
    filout.write(data)
    filout.close()
