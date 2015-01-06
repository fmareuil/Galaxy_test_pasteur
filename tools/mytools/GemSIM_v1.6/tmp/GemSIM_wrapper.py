#!/usr/bin/python -W ignore

import sys
import optparse
import subprocess
import tarfile
from os.path import abspath, realpath, dirname, join as joinpath

resolved = lambda x: realpath(abspath(x))

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()
 
def untar():
     pass
     return
def launcher(cmd):
    pass
    return
    
def parsetar():
    pass
    return

def _badpath(path, base):
    # joinpath will ignore base if path is absolute
    print resolved(joinpath(base,path)).startswith(base)
    return not resolved(joinpath(base,path)).startswith(base)

def _badlink(info, base):
    # Links are interpreted relative to the directory containing the link
    tip = resolved(joinpath(base, dirname(info.name)))
    return _badpath(info.linkname, base=tip)

def safemembers(members):
    base = resolved(".")

    for finfo in members:
        print finfo
        if _badpath(finfo.name, base):
            print >>sys.stderr, finfo.name, "is blocked (illegal path)"
        elif finfo.issym() and _badlink(finfo,base):
            print >>sys.stderr, finfo.name, "is blocked: Hard link to", finfo.linkname
        elif finfo.islnk() and _badlink(finfo,base):
            print >>sys.stderr, finfo.name, "is blocked: Symlink to", finfo.linkname
        else:
            yield finfo
    
def __main__():
    parser = optparse.OptionParser()
    parser.add_option('-r','', dest='r', help='reference genome, in fasta format.' )
    parser.add_option('-R','', dest='R', help='Only for metagenome projects. Archive tar containing references.' )
    parser.add_option('-a','', dest='a', help='Only for metagenome projects. Species-abundance file.' )
    parser.add_option('-n','', dest='n', help='number of reads to produce. For paired end reads, number of pairs.' )
    parser.add_option('-g','', dest='g', help='haplotype file, specifying location and frequency of snps.')
    parser.add_option('-G','', dest='G', help='Archive tar of haplotype files (metagenomics mode only.')
    parser.add_option('-l','', dest='l', help='length of reads. Integer value, or -l d for empirical distribution.' )
    parser.add_option('-m','', dest='m', help='error model file *_single.gzip or *_paired.gzip.' )
    parser.add_option('-c','', dest='c', help='use this flag if you wish to draw reads from a circular genome.' )
    parser.add_option('-q','', dest='q', help='quality score offset. Usually 33 or 64 (see manual).' )
    parser.add_option('-o','',  dest='o', help='output file name prefix.' )
    parser.add_option('-u','', dest='u', help='Mean fragment length for paired end reads. -u d for empirical.' )
    parser.add_option('-s','', dest='s', help='standard deviation for fragment length. Use only with -u and -p.' )
    parser.add_option('-p','', dest='p', help='use only to create paired end reads.' )
    (options, args) = parser.parse_args()
    
    command = {}
    if options.r:
        command["-r"] = options.r
    if options.R:
        tarmembers = tarfile.open(options.R)
        tarmembers.extractall(path="temporaire", members=safemembers(tarmembers))
	tarmembers.close()
        command["-R"] = options.R
    if options.a:
        command["-a"] = options.a    
    if options.n:
        command["-n"] = options.n   
    if options.g:
        command["-g"] = options.g	
    if options.G:
        command["-G"] = options.G	
    if options.l:
        command["-l"] = options.l		
    if options.m:
        command["-m"] = options.m	
    if options.c:
        command["-c"] = options.c		
    if options.q:
        command["-q"] = options.q
    if options.o:
        command["-o"] = options.o
    if options.u:
        command["-u"] = options.u	
    if options.s:
        command["-s"] = options.s		
    if options.p:
        command["-p"] = options.p
    cmd = "python "
    for key in command.keys():
        cmd = cmd + key +' '+ command[key] + ' '
    print cmd
    
				
if __name__=="__main__": __main__()
