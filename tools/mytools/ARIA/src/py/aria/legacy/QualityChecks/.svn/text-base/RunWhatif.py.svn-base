"""
an utility function to run WHATIF with os.system()

author: Jens Linge
"""
import os, shutil, string, time
## from Parameters import *

## TODO: code is still in QualityChecks.py  ...
## def PreparePdb():
##     """
##     use PreparePdb() to select a certain range of your pdb files,
##     e.g residues 13-74 instead of taking all residues for the checks
##     """
##     pass
    
def RunWhatif(pdbList, tmpDir, fulchkReportFN, checkdbFN,\
              trashDirectory,\
              whatIfExe,\
              slowNetworks):
    """
    runs the WHATIF WHATCHECK program on a given list of pdb files
    pdbList should be a list of filenames with absolute paths (strings)
    tmpDir is the directory (string) to which all the big output files
      are written to
    reportFN is the filename of the files with the raw results
    """
    print 'starting WHATIF FULCHK'
    os.chdir(tmpDir)
    fulchkTmpName = 'FULCHK_' + str(time.time())[:-5] + str(time.time())[-3:]
    checkdbTmpName = 'CHECKDB_' + str(time.time())[:-5] + str(time.time())[-3:]
    open(tmpDir + '/' + fulchkTmpName, 'w').close()
    open(tmpDir + '/' + checkdbTmpName, 'w').close()

    for file in pdbList:
        print 'starting whatif with file:', file
        time.sleep(2)
        if os.path.exists('check.db'):
            os.remove('check.db')
        whatToDo = ' %FULCHK \n ' + file + ' \n ' + file + '\n exit\n y\n'
        startUpFil = open(tmpDir + '/STARTUP.FIL', 'w')
        startUpFil.write(whatToDo)
        startUpFil.close()

	## 161106: modified in order to use WhatCheck instead of WhatIF

        if os.path.split(whatIfExe)[1] in ['whatcheck','DO_WHATCHECK.COM']:
          os.system('cd ' + tmpDir + ';' + whatIfExe + ' ' + file)
        else:
          os.system('cd ' + tmpDir + ';' + whatIfExe)
        os.system('cd ' + tmpDir + ';' + 'cat pdbout.txt >> ' + fulchkTmpName)
        os.system('cd ' + tmpDir + ';' + 'cat check.db >> ' + checkdbTmpName)
    shutil.copy(tmpDir + '/' + fulchkTmpName, fulchkReportFN)
    shutil.copy(tmpDir + '/' + checkdbTmpName, checkdbFN)
    print 'wrote file FULCHK_REPORT and CHECKDB_REPORT.'


