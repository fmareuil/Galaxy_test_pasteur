"""
an utility function to run PROSA-II with os.system()

author: Jens Linge
"""
import os, shutil, string, time, sys, traceback
## from Parameters import *

def RunProsa(pdbList, tmpDir, reportFN,\
             trashDirectory,\
             prosaExe,\
             slowNetworks):
    """
    runs the PROSA-II program on a given list of pdb files
    pdbList should be a list of filenames in tmpDir
    tmpDir is the directory (string) to which all the big output files
      are written to
    reportFN is the absolute filename of the files with the raw results
    """
    os.chdir(tmpDir)
    ## Mareuil <
    #prosaExeDF = prosaExe + ' -d -f '  #no display, enable file input
    prosaExeDF = prosaExe + ' '
    ## Mareuil >
    inputFN = 'prosa_input.txt'      #no need to change this temp name

    print 'starting PROSA-II checks, using temporary directory:'
    prosaOutputHandle = open(reportFN, 'w')
    for file in pdbList:
        print '  working with file:', file
        # the input string for PROSA-II
        whatToDo = 'read pdb ' + file + ' ' + file + ' \n' +\
                   'analyse energy ' + file + ' \n' +\
                   'print energy ' + file + ' ' + file + '\n' +\
                   'delete ' + file + ' \nquit\n'

        commandFileName = os.path.join(tmpDir,inputFN)
        commandFile = open(commandFileName, 'w')
        commandFile.write(whatToDo)
        commandFile.close()
        #starting prosa:
        ## Mareuil <
        print 'starting prosa with command:'
        #print prosaExeDF + commandFileName
        print prosaExeDF + inputFN
        #os.system(prosaExeDF + commandFileName)
        os.system(prosaExeDF + inputFN)
        ## Mareuil >
        #appending the results:
        anaFN = file + '.ana'
        try:
            anaHandle = open(anaFN,'r')
        except:
            print 'WARNING: .ana file could not be found. skipping Prosa check.'
            print '-'*60
            traceback.print_exc(file=sys.stdout)
            print '-'*60
            return
            
        anaString = string.join(anaHandle.readlines())
        prosaOutputHandle.write(anaString)
        anaHandle.close()

    prosaOutputHandle.close()
