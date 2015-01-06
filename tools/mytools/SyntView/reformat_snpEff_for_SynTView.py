#! /local/gensoft2/adm/bin/python2.6 

import os
import pwd
import re
from sys import argv, exit, stdout, stderr, path

from optparse import OptionParser

import SynTView.Formatter as Formatter
import PTT.PTTReader as PTTReader
import SnpEffUtils.Reader as SnpEffReader

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-p", dest="ptt_file", 
        help=".ptt file", default=None)
    parser.add_option("-s", dest="snp_file", 
        help="snpEff output file", default=None)
    parser.add_option("-l", dest="log_file", 
        help="Optional log file", default=None)
    parser.add_option("-o", dest="out_prefix", 
        help="Prefix for output files", default=None)
    
    # grab options and arguments in an usable form
    if not argv[1:] :
       argv.append('-h')

    prog_name = os.path.basename(argv[0])
    (options, args) = parser.parse_args()

                
    if not options.snp_file :
        parser.error("\nsnpEff output file missing.")
    else :
        if not os.access(options.snp_file, os.R_OK):
            parser.error("Unable to read snpEff output file %s " %options.snp_file)
                
    if not options.out_prefix :
        parser.error("\nPrefix for output file is missing.")
    else :
        out_snp   = options.out_prefix+".snp"
        out_indel = options.out_prefix+".indel"
        if os.access(out_snp, os.F_OK) and not os.access(out_snp, os.W_OK):
            parser.error("\nSNP output file %s exists but is not writable." %out_snp)
        if os.access(out_indel, os.F_OK) and not os.access(out_indel, os.W_OK):
            parser.error("\nINDEL output file %s exists but is not writable." %out_indel)

    log_fh = stderr
    if options.log_file:
        if os.access(options.log_file, os.F_OK) and not os.access(options.log_file, os.W_OK):
            parser.error("\n log file %s exists but is not writable." %options.log_file)
        else:
            log_fh = open(options.log_file,"w")

    if not options.ptt_file :
        parser.error("\nPTT file missing.\n")
    else :
        if not os.access(options.ptt_file, os.R_OK):
            parser.error("Unable to read PTT file %s " %options.ptt_file)

    ptt_reader      = PTTReader.PTTReader()
    snp_reader      = SnpEffReader.SnpEffReader()
    snp_reformatter = Formatter.SynTViewSnpFormatter()

    with open(options.ptt_file, "r") as ptt_fh:
        ptt_dict    = ptt_reader.read_and_check(ptt_fh)

    with open(options.snp_file, "r") as snp_fh:
        snp_list = snp_reader.read_and_check(snp_fh)
        
    with open(out_snp, "w") as out_snp_fh:
        with open(out_indel, "w") as out_indel_fh:
            for snp in snp_list:
                try:
                    type, reformatted_line = snp_reformatter.format(snp, ptt_dict)
                    if type == "SNP":
                        print >>out_snp_fh, reformatted_line
                    if type == "INS" or type == "DEL":
                        print >>out_indel_fh, reformatted_line
                except Formatter.SynTViewSnpFormatterGeneError, msg:
                    print >>log_fh, str(msg)
                except Formatter.SynTViewSnpFormatterIntragenicError, msg:
                    print >>log_fh, str(msg)+" ignored"
