"""
an utility function to run PROCHECK with os.system()

author: Jens Linge
"""
import glob, os, shutil, string, time
## from Parameters import *

def RunProcheck(pdbList, tmpDir, reportFN, postscriptDir,\
                trashDirectory,\
                procheckExe,\
                cshExe,\
                pdbSuffix,\
                procheckResolution,\
                slowNetworks):
    """
    runs the PROCHECK program on a given list of pdb files
    pdbList should be a list of filenames with absolute paths (strings)
    tmpDir is the directory (string) to which all the big output files
      are written to
    reportFN is the filename of the file with the raw results
    postscriptDir is the directory where the .ps files are written to
    procheckResolution needed to run procheck
    """
    inputFN = 'procheck.csh' #no need to change this filename

    os.chdir(tmpDir)
    print 'starting PROCHECK with resolution: ' + procheckResolution
    reportHandle = open(os.path.join(tmpDir, reportFN), 'w')
    for file in pdbList:
        print 'working on file', file
        #writing temporary csh script:
        cshScript = """cd %s
    %s %s %s
    """ % (tmpDir, procheckExe, file, procheckResolution)
        cshFN = os.path.join(tmpDir, inputFN)

        cshHandle = open(cshFN, 'w')
        cshHandle.write(cshScript)
        cshHandle.close()

        os.system(cshExe + ' ' + cshFN)
        
        #appending the results:
        time.sleep(slowNetworks)
        upTo = -len(pdbSuffix)
        outputFN = file[:upTo] + '.sum'
        print outputFN, os.getcwd()

        outputHandle = open(outputFN,'r')
        outputString = string.join(outputHandle.readlines())
        print outputString #test
        reportHandle.write(outputString)
        outputHandle.close()
    reportHandle.close()

    #copying the postscript files:
    if postscriptDir and os.path.exists(postscriptDir):
        psFiles = glob.glob('*.ps')
        print 'copying postscript files'
        for fileN in psFiles:
            print 'copying', fileN, 'to', postscriptDir
            shutil.copyfile(fileN, os.path.join(postscriptDir,fileN))
            os.system('gzip -f ' + fileN) # -f for overwriting
        newpsFiles = glob.glob(postscriptDir + '/*ps')
        for fileN in newpsFiles:
            print 'running gzip on file: ' + fileN
            os.system('gzip -f ' + fileN) # -f for overwriting

