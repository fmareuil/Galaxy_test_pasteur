
"""
------------------------------------------
Author: Olivia Doppelt-Azeroual
Organism: CIB, Institut Pasteur
Version: 0.2

------------------------------------------

formatdb and blast wrapper:
Runs the either:  
- formatdb 
- mblastall


usage: blast_wrapper.py --program [options]
-program: formatdb or blastall 



"""


import string, sys, os, optparse, subprocess, glob, tarfile,shutil

def getPathOnly(filePath):
    try:
        return string.join(string.split(filePath,"/")[0:len(string.split(filePath,"/"))-1],"/")
    except OSError:
        return "."

def getFileNameOnly(filePath):
    return string.join(string.split(filePath,'/')[len(string.split(filePath,"/"))-1],"")

def main():
    
    # Command line parser creation
    parser = optparse.OptionParser()

    parser.add_option("--progName",action="store",
                      help="mblastall or formatdb")
    
    # formatdb options
    parser.add_option("--i", action="store",
                      help="Fasta Sequence to enter in the new DB")
    parser.add_option("--p", action="store",
                      help="SequenceType")
    parser.add_option("--n", action="store", 
                      help="output prefix name")


    # blastall options
    parser.add_option( '--threads', action="store", 
                      help='The number of threads to use' )
    parser.add_option("--blast_type",action="store",
                      help="Blast type blastn or blastp")
    parser.add_option("--query",action="store",
                      help="Nucleotide query sequence(s)")
    parser.add_option("--db",action="store",
                      help="available DB")
    parser.add_option("--evalue",action="store",
                      help="Set expectation value cutoff")
    parser.add_option("--out_format",action="store",
                      help="Output format")
    parser.add_option("--filter_query",action="store",
                      help="Filter out low complexity regions (with DUST)")
    parser.add_option("--max_hits",action="store",
                      help="Maximum hits to show regarless of the hit species")
    parser.add_option("--best_hits",action="store",
                      help="Keeps the best hit only ")
    parser.add_option("--word_size",action="store",
                      help="Word size for wordfinder algorithm")
    # parser.add_option("--word_sizeAdv",action="store",
    #                   help="Word size for wordfinder algorithm")
    parser.add_option("--output",action="store",
                      help="Blast output")
    
    opts,args=parser.parse_args()



    dico=vars(opts)
    argumentDict = {}
    # Assiging the values to the right variable names
    for k,v in dico.iteritems():
        if v!= None and v!= False :
            argumentDict[k] = str(v)


    # FormatDB launching
    if opts.progName == "formatdb":
        dbName = argumentDict['n']
#        sys.stderr.write("dbName:"+ dbName + "\n")
        cmd = " formatdb -i " + argumentDict['i'] + " -p " + argumentDict['p'] + " -n " + dbName
        print cmd
        proc = subprocess.Popen( cmd, shell=True, stderr=sys.stderr, stdout=sys.stdout)
        proc.wait()

        # compression of the blast DB created
        pathDB=dbName+'.*'
        listFile=glob.glob(pathDB)
        path=os.getcwd()
        newPath=getPathOnly(dbName)
        tarDB=tarfile.open(dbName+'.tar.gz','w|gz')
        os.chdir(newPath)
        for i in listFile:
            iname=getFileNameOnly(i)
            tarDB.add(iname, recursive=False)
            os.remove(iname)
            
        tarDB.close()
        os.chdir(path)
        
        # renaming of the output for Galaxy to copy the result
#        sys.stderr.write("tarDBNAME"+tarDB.name+"\n")
        shutil.move(tarDB.name,dbName)

    # blastall launching
    else:
        dbName=argumentDict["db"]
      
        try:
            tarfile.is_tarfile(dbName)
            localDB=tarfile.open(dbName,'r')
            localDB.extractall()
            
            # Adjust the DB name according to Galaxy Working Dir Name.
            # VERY UGLY ....
            dbNames = localDB.getnames()
            nsq = string.split(dbNames[0],'.')
            argumentDict["db"]=nsq[0]+"."+nsq[1]
            
        except IOError:
            ()

        arguments=[""]*20
        for k,v in argumentDict.iteritems():
            
            if k == "progName":
                arguments[0] = ('%s ' % (v))
            if k == "threads":
                arguments[1] =  (' -a %s ' % (v))
            if k == "blast_type":
                arguments[2] = (' -p %s ' % (v))
            if k == "query":
                arguments[6] = ('-i %s ' % (v))
            if k == "db":
                arguments[8] = ('-d %s ' % (v))
            if k == "evalue":
                arguments[10] = (' -e %s ' % (v))
            if k =="out_format":
                arguments[12] = (' -m %s ' % (v))
            if k =="filter_query":
                arguments[14] = (' -F %s ' % (v))
            if k =="max_hits":
                arguments[14] = (' -v %s ' % (v))
            if k =="best_hits":
                arguments[17] = (' -b %s ' % (v))
            if k =="word_size":
                arguments[18] = (' -W %s ' % (v))
            if k =="output":
                arguments[19] = (' -o %s ' % (v))

            
        cmd = ''.join(arguments)
        print cmd
        proc = subprocess.Popen( cmd, shell=True, stderr=sys.stderr, stdout=sys.stdout)
        proc.wait()    
        
        
        



if __name__=="__main__": 
    main()
