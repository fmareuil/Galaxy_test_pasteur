import sys
import getopt

def usage():
    print '\n\n########################################################################'
    print '# GemSIM - Generic Error Model based SIMulator of N.G. sequencing data #'
    print '########################################################################\n'
    print '\nGemReads.py:\n'
    print 'Takes a reference genome, an empirical error model, and a haplotype file'
    print 'listing SNP locations and frequencies, and creates a simulated data set'
    print 'of random reads, as would be produced by a next-gen sequencing run.'
    print 'Output is in fastq format, suitable for input into popular alignment'
    print 'software.' 
    print '\nOptions:'
    print '      -h prints these instructions.'
    print '      -r reference genome, in fasta format.' 
    print '      -R Only for metagenome projects. Directory containing references.' 
##MAREUIL
    print '      -L Only for metagenome projects. Liste of references.'
##MAREUIL
    print '      -a Only for metagenome projects. Species-abundance file.'
    print '      -n number of reads to produce. For paired end reads, number of pairs.'
    print '      -g haplotype file, specifying location and frequency of snps.'
    print '      -G directory of haplotype files (metagenomics mode only).'
    print '      -l length of reads. Integer value, or -l d for empirical distribution.'
    print '      -m error model file *_single.gzip or *_paired.gzip.'
    print '      -c use this flag if you wish to draw reads from a circular genome.'
    print '      -q quality score offset. Usually 33 or 64 (see manual).' 
    print '      -o output file name prefix.'
    print '      -u Mean fragment length for paired end reads. -u d for empirical.'
    print '      -s standard deviation for fragment length. Use only with -u and -p.' 
    print '      -p use only to create paired end reads.\n\n'
    
def main(argv):
    Liste=''   
    direct=''
    refFile=''
    try:
        opts, args = getopt.getopt(argv, "hR:L:r:")
        print opts, args
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts: 
        print opt, 'opt', arg, 'arg'      
        if opt == '-h':
            usage()
            sys.exit()
        if opt == '-R':
            print arg
	elif opt == '-L':           
            print arg
        elif opt == '-r':
            print arg   

if __name__=="__main__":
    main(sys.argv[1:])
