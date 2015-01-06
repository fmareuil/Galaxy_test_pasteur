PROGRAM
=======

map2cov.py - run as command line in a shell
map2cov_tk.py - run as graphical user interface in Tcl/Tk


VERSION
=======

Version 4.3


INTRODUCTION
============

The MAP2COV program converts mapping data into coverage file readable by the 
COV2HTML interface. From an alignment map of NGS reads and an annotation file 
associated, MAP2COV calculates the coverage per base of the reference sequence, 
as well as the mean coverage for each CDS and intergenic region.
MAP2COV is able to parse mapping alignments in SAM/BAM, ELAND and WIG formats. 
The CDS list can be provided from an GenBank or EMBL annotation file, or as a 
FASTA or GFF3 file along with a reference genome file in FASTA format.
The output text file generated is compressed in bzip2 (bz2) format in order to 
be easily load on the COV2HTML interface for a futher analysis and visualization 
of data.
The COV2HTML interface is available at http://mmonot.eu/COV2HTML/connexion.php


PREREQUISITES
=============

Unix/Linux
- None (Python 2.X including Tk/Tcl are pre-integrated)

Python and Tk/Tcl source code
http://www.python.org/getit/releases/2.7.5/
http://www.tcl.tk/software/tcltk/download.html


COMMAND LINE
============

    ./map2cov.py [option] reference_file [option] alignment_file -t 'DNA|RNA' {options} args
    
    Help:   
    
        ./map2cov.py -h
        ./map2cov.py --help
    
    Options:
        
        * Reference information
            Annotated reference sequence:
            -g|--genbank    INPUT_FILE: GenBank annotation file
            -e|--embl       INPUT_FILE: EMBL annotation file

            Reference Sequence + CDS Annotation: (both parameters are inseparable)
            -f|--fasta      INPUT_FILE: Fasta reference sequence
            -c|--cds        INPUT_FILE: Fasta or GFF3 CDS regions

        * Alignment map:
            -s|--sam        INPUT_FILE: SAM/BAM alignment file
            -i|--eland      INPUT_FILE: ELAND alignment file
            -w|--wig        INPUT_FILE: WIG coverage file

        * Data type:
            -t|--type       STRING: 'DNA' (DNA-, ChIP-seq) or 
                                    'RNA' (RNA-seq) or
                                    'TSS' (TSS-seq) or
				    'CHIP' (ChIP-seq) data

        * Optional parameters:
            -a|--confidential   STRING: 'YES' or 'NO'
                                        Export data confidentially (default: 'NO')
                                        CDS names are replaced by anonymous identifiers
            -b|--stranded       STRING: 'YES' or 'NO'
                                        Perform a strand-specific processing (default: 'NO')
                                        Two output files are generated: one for the sense strand and 
                                        one for the antisens strand (extension '_plus.txt.bz2' and 
                                        '_minus.txt.bz2' respectively)
            -l|--read_length    INTEGER: Ignore reads shorter than fixed threshold 
                                         (default: 20; 0 to disable; accept only positive integer)
	    -v|--cov_style     STRING: 'MEAN' or 'COUNTS'
	                       		Genetic Elements Coverage style: MEAN coverage or Raw COUNTS
               	 	       		(default: MEAN)
            -o|--output         OUTPUT_PREFIX: Prefix for the output text file compressed in 
                                               bz2 format (default: export_map2cov.txt.bz2)


TK INTERFACE
============

    ./map2cov_tk.py &

    (i)   Reference information fields:
          Browse GenBank, EMBL or FASTA file of the reference sequence used for the mapping 
          assembly. In case of FASTA sequence, a FASTA file containing the CDS regions is 
          required.
    
    (ii)  Alignment file fields:
          Browse SAM/BAM, ELAND or WIG file corresponding to the mapping of reads against 
          the reference sequence provided.
    
    (iii) Data type field:
          Choose 'DNA' for DNA- or ChIP-sequencing data.
          Choose 'RNA' for RNA-sequencing data.
          Choose 'TSS' for TSS-sequencing data (not available for WIG files, no sequence reads).
          Choose 'CHIP' for ChIP-sequencing data.

    (iv)  Output directory field:
          Browse a directory where the output file will be generated. By default, the output 
          file is saved in the current directory, where the program is executed.
          The output text file is named 'export_map2cov.txt.bz2' and bzip2 compressed.
          
    (v)   Optional parameters:
          (a) Anonymize data field:
              Choose 'Yes' to replace CDS names by anonymous identifiers in the output file.

          (b) Read-length threshold field:
              Provide a minimum size for filtering reads too short. By default, the threshold is 
              set at 20 (0 to disable; accept only positive integer). This parameter is ignored 
              for WIG files.

          (c) Strand-specific data field:
              Choose 'Yes' in case of strand-specific sequencing data. This parameter ignored 
              for WIG files). Two output files will be then generated:
                  'export_map2cov_plus.txt.bz2' for sense strand;
                  'export_map2cov_minus.txt.bz2 for antisense strand.

          (d) Genetic Elements Coverage Style:
              Choose 'MEAN' : Coverage of genes and intergenic regions is calculated with the mean coverage
              Choose 'COUNTS : Coverage of genes and intergenic regions is the number of read that match on them (Raw Counts).
          
    (vi)  Launch analysis
          Execute the analysis. During the processing a pop-up window is displayed in order 
          to follow the files parsing.


SUPPLEMENTARY INFORMATION
=========================

* The program does not handle multi-exon genes in GenBank and EMBL files.

* The TSS datatype, the stranded analysis and the read-length filtering are not 
  available for WIG files.

* The extension '.txt.bz2' is added to the output file if missing.


TEST FILES
==========

* MAP2COV is provided with input test files which contain:

    (i) Reference information from 20 kb of the Clostridium difficile 630 genome (NC_009089):
        (a) Annotated sequence
            GENBANK 'Annotation_GenBank.gbk'
            EMBL 'Annotation_EMBL.gbk'

        (b) Sequence and annotation
            FASTA 'Reference_FASTA.fasta'
            FASTA 'CDS_FASTA.fasta' or GFF3 'CDS_GFF3.gff'

    (ii) Alignment map of 300 in reverse strand from position 2001 to 2045, 300 in forward strand from position 17956 to 18000:
        (a) SAM 'Format_SAM.sam' or BAM 'Format_BAM.bam'

        (b) ELAND 'Format_ELAND.txt'

        (c) WIG 'Format_WIG.txt'


CONTACT
=======

Marc Monot <marc.monot@pasteur.fr>

