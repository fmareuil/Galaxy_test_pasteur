#!/usr/bin/env python

import optparse, string, os

def foo_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split())
    
def test_files(pathfile):
    namefile = os.path.basename(pathfile)
    path = os.path.dirname(pathfile)
    prefix_file = "%s.dat" % string.split(namefile,".")[0]
    strresponse = ""
    if os.path.exists(os.path.join(path, "%s.bam" % prefix_file)):
        strresponse = ''.join([strresponse, "bam"])
    if os.path.exists(os.path.join(path, "%s.bai" % prefix_file)):
        strresponse = ''.join([strresponse, "bai"])
    return strresponse

if __name__ == '__main__':
    op = optparse.OptionParser()
    op.add_option('-i', '--input', type='string')
    op.add_option('-o', '--output', default=None)
    opts, args = op.parse_args()
    flagexist = test_files(opts.input)
    print flagexist == "bam"
    if flagexist == "bambai":
        out = open(opts.output,"w")
	out.write("il y a les 2 fichiers bam et bai")
	out.close()
    elif flagexist == "bam":
        out = open(opts.output,"w")
	out.write("il n'y a que le fichier bam")
	out.close()
    elif flagexist == "bai":
        out = open(opts.output,"w")
	out.write("il n'y a que le fichier bai")
	out.close()
    else:
        out = open(opts.output,"w")
        out.write("il y a ni le bam ni le bai")
	out.close()
