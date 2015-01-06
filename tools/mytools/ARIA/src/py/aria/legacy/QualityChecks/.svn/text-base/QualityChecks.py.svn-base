"""
QualityChecks.py

a module to validate an ensemble of structures with the programs:
  PROCHECK
  WHATIF
  PROSA

QualityChecks.py runs these programs and provides some statistics of the
results.

One tip for the setup of PROSA and PROCHECK. They both need some system variables set.
Thus, it's convenient to have a script like:

-----------

#!/bin/csh
setenv PROSA_BASE /some/directory/prosa/prosabase/
/some/directory/prosa/bin/prosaII $*

-----------

instead of the 'pure' binary /home/Bis/shared/rh71/prosa/bin/prosaII

For PROCHECK, use something like:

-----------

#!/bin/csh
setenv  PROCHECK /some/directory/procheck-3.5.4
setenv prodir /some/directory/procheck-3.5.4
$PROCHECK/procheck.scr $*       

-----------

author: Jens Linge, Pasteur Institute, Paris
        linge@pasteur.fr
"""
__author__   = "$Author: bardiaux $"
__revision__ = "$Revision: 1.2 $"

__help__ = """
starts the following checks on all pdb files in the specified directory:
    PROCHECK:     Ramachandran plot, bad contacts
    WHATIF:       FULCHK
    PROSA:        inverse Boltzmann energy

all the output files of these programs are written to a specified
temporary directory (since they waste a lot of disk space!)

INSTALLATION NOTES:
please edit the header of the runChecks() function below
according to the setup on your machine
trashDirectory       should be a trash disk
procheckOnOff        1 if you want to use procheck, 0 if not
prosaOnOff           1 if you want to use prosa, 0 if not
whatifOnOff          1 if you want to use whatcheck, 0 if not
whatIfExe            absolute path of whatif executable 
procheckExe          absolute path of procheck executable 
cshExe               absolute path of csh executable 
prosaExe             absolute path of prosa executable 
latexExe             absolute path of latex executable 
dvipsExe             absolute path of dvips executable 
aminoAcidRange       usually 'all' or something like '23 58'
pdbSuffix            usually '.pdb' (should have three letters)
procheckResolution   default is '2.0'
postscriptFilesOnOff for copying over the PROCHECK .ps files 
useFileNam           whether you want to use a file.nam file
howManyPdb           number of pdb structures you want to include
removeTrashOnOff     removes the trash files automatically
slowNetworks         waiting for slow networks, usually 5-20 sec
commentStringForOutput  appears in .tex and .ascii summary files

linge@pasteur.fr
"""

## modifications: WR

FILENAME_REPORT = 'quality_checks'
FILENAME_REPORT_TEXT = FILENAME_REPORT
FILENAME_REPORT_TEX = FILENAME_REPORT + '.tex'
FILENAME_REPORT_PROSA = FILENAME_REPORT + '.prosa'
FILENAME_REPORT_PROCHECK = FILENAME_REPORT + '.procheck'
FILENAME_REPORT_WHATIF = FILENAME_REPORT + '.whatif'
FILENAME_REPORT_WHATIF_FULL = FILENAME_REPORT + '.whatif_fulchk'
# BARDIAUX 2.3
FILENAME_REPORT_MP_CLASHSCORE = FILENAME_REPORT + '.clashscore'

import glob, os, re, shutil, string, sys, time, traceback, random
## BARDIAUX
from aria.ariabase import AriaBaseClass
ARIA_VERSION = AriaBaseClass().get_version_string()
## Mareuil <
def copypdb_for_prosa(filename, dest_directory, verbose):
    readfile = open(filename,'r')
    listfile = readfile.readlines()
    readfile.close()
    writefile = open(dest_directory +'/'+filename,'w')
    for i in listfile:
        if i[13:16] == 'OT1':
            writefile.write(i[0:13]+'O  '+i[16:])
	elif i[13:16] == 'OT2':
            if verbose: print 'pdb modification for prosa2003'
        else:
            writefile.write(i)
    writefile.close()
## Mareuil >
#FOR INSTALLATION PLEASE CHANGE THE NEXT FEW LINES:
def runChecks(workingDirectory = '/tmp/test_pdb',\
              trashDirectory = '/tmp/trash',\
              procheckOnOff = 1,\
              prosaOnOff = 1,\
              whatifOnOff = 1,\
              clashlistOnOff = 1,\
              whatIfExe = 'whatif',\
              procheckExe = 'procheck',\
              cshExe = 'csh',\
              prosaExe = "prosaII",\
              clashlistExe = "clashlist",\
              latexExe = 'latex',\
              dvipsExe = 'dvips',\
              aminoAcidRange = 'all',\
              pdbSuffix = '.pdb',\
              procheckResolution = '2.0',\
              postscriptFilesOnOff = 0,\
              useFileNam = 1,\
              howManyPdb = 20,\
              removeTrashOnOff = 0,\
              slowNetworks = 2,
              commentStringForOutput = '',
              verbose = 0,
              return_procheck_results = 0,
              fileList=None):

    from DelTrailingSlash import DelTrailingSlash
    from RunProcheck import RunProcheck
    from RunProsa import RunProsa
    from RunWhatif import RunWhatif
    from Descriptive import Descriptive
    from ReadCheckDB import readCheckDB

    from aria.Molprobity import MolprobityClashlist

    #before take-off, check whether procheck, whatcheck, prosa executables do exist:
    if not os.path.exists(procheckExe):
        if verbose: print 'can not find PROCHECK executable on disk => PROCHECK disabled.'
        procheckOnOff = 0  #overwriting the procheckOnOff switch
    if not os.path.exists(whatIfExe):
        if verbose: print 'can not find WHATIF executable on disk => WHATIF disabled.'
        whatifOnOff = 0  #overwriting the whatifOnOff switch
    if not os.path.exists(prosaExe):
        if verbose: print 'can not find PROSA executable on disk => PROSA disabled.'
        prosaOnOff = 0  #overwriting the prosaOnOff switch
    if not os.path.exists(clashlistExe):
        if verbose: print 'can not find CLASHLIST executable on disk => CLASHLIST disabled.'
        clashlistOnOff = 0  #overwriting the clashlistOnOff switch
        
    # type conversion of howManyPdb:
    if howManyPdb == None:
        howManyPdb = 0
    if type(howManyPdb) == type('a'):
        howManyPdb = string.atoi(howManyPdb)
        
        
    if (procheckOnOff == 0) and (whatifOnOff == 0) and (prosaOnOff == 0) and (clashlistOnOff == 0):
        if verbose: print 'no executables found for PROCHECK, PROSA, WHATIF, clashlist. giving up.'
        return
    

    ###############################################################################
    ## #0. specify the working directory and a list of files to work on:
    ## #if a working directory is specified on command line, take it:
    ## if len(sys.argv) == 2:
    ##     workingDirectory = sys.argv[1]
    ## elif len(sys.argv) == 1:
    ##     workingDirectory = os.getcwd()
    ##     print 'working on directory', workingDirectory
    ## else:
    ##     print 'TOO MANY COMMANDLINE ARGUMENTS GIVEN!'
    ##     print 'READ THE INSTRUCTIONS:', __help__
    ##     print 'PROGRAM ABORTED.'
    ##     sys.exit(1)

    #chdir to working directory:
    print workingDirectory
    os.chdir(workingDirectory)

    #get a list of pdb files:

    ## WR: try to read file.nam from [working_dir]/molmol

    molmol_path = os.path.join(workingDirectory, 'molmol')
    fileNam = os.path.join(molmol_path, 'file.nam')

    if useFileNam and os.path.exists(fileNam):
        if verbose: print 'USING THE BEST ' + str(howManyPdb) + ' STRUCTURES FROM file.nam' 
        fileHandle = open(fileNam)
        nameList = fileHandle.readlines()
        fileList = []
        qqq = 0
        for eachLine in nameList:
            qqq = qqq + 1
            if qqq > howManyPdb:
                break
            if string.strip(eachLine):
                fileList.append(string.strip(eachLine))
    else:

        if not fileList:
        
            fileList = glob.glob('*' + pdbSuffix)
            fileList.sort()

    try:
        fileList = fileList[:howManyPdb]
        
    except:
        # if the list is shorter than howManyPdb shorter,
        # no need to worry...
        pass

    ## WR: use basenames only

    fileList = [os.path.basename(fn) for fn in fileList]

    startText= '''QualityCheck.py
    a module to run WHATCHECK, PROCHECK, PROSA for an ensemble of files
    '''
    if verbose: print startText

    if verbose: print "PYTHONPATH set to:"
    if verbose: print sys.path

    if verbose: print '\nworking on the files:'
    for eachFN in fileList:
        if verbose: print eachFN

    ###############################################################################
    #1. copy the PDB files to the temporary directory:

    # delete some whitespace:
    aminoAcidRange = string.lower(string.strip(aminoAcidRange))
    pdbSuffix = string.strip(pdbSuffix)

    #mkdir in trashDirectory:

    NN = random.random() 
    tempDirectory = trashDirectory + '/TEMP_' + str(time.time()) + str(NN)[2:]
    os.mkdir(tempDirectory)

    #write a README file to the TEMP_ directory:
    readmeHandle = open(tempDirectory + '/README', 'w')
    readmeHandle.write('QualityChecks.py was started in:\n  ')
    readmeHandle.write(workingDirectory + '\n')
    readmeHandle.write('trash directory:\n  ')
    readmeHandle.write(tempDirectory + '\n')
    readmeHandle.write('-> please remove this temporary directory after the program has finished.\n')
    readmeHandle.write('working with the ' + str(len(fileList)) + ' files:\n')
    for eachFile in fileList:
        readmeHandle.write('  ' + eachFile + '\n')
    readmeHandle.close()

    #copy some files:
    for file in fileList:
        try:
            if verbose: print 'copying', file, 'to', tempDirectory
            ## Mareuil <
            copypdb_for_prosa(file, tempDirectory, verbose)
            #shutil.copyfile(file, tempDirectory + '/' + file)
            ## Mareuil >
        except:
            if verbose: print 'ERROR: could not copy the file:'
            if verbose: print '       ' + file
            if verbose: print '       to the directory:'
            if verbose: print '       ' + tempDirectory
            return

    ###############################################################################
    #2. prepare PDB files if necessary (using WHATIF):

    #use WHATIF to generate a new .pdb file which contains only the specified range:
    if aminoAcidRange != 'all' or (not aminoAcidRange):
        os.chdir(tempDirectory)
        if verbose: print 'waiting ' + str(slowNetworks) + ' sec for slow networks'
        time.sleep(slowNetworks) #waiting for slow networks 
        if verbose: print 'starting whatif to generate pdb files'
        aminoAcidRange = re.sub('-', ' ', aminoAcidRange)
        aminoAcidRange = re.sub(':', ' ', aminoAcidRange)
        aminoAcidRange = re.sub(';', ' ', aminoAcidRange)
        rangeList = string.split(aminoAcidRange)
        rangeList[0] = str(string.atoi(rangeList[0]) - 1)
        rangeList[1] = str(string.atoi(rangeList[1]) + 1)
        endNewRange = str(string.atoi(rangeList[1]) - 2)
        if len(rangeList) != 2:
            if verbose: print 'WARNING: specified amino acid range is not valid!'
            if verbose: print 'PROGRAM ABORTED'
            return
        for file in fileList:
            if verbose: print 'converting', file, 'with range:', aminoAcidRange
            whatToDo = ' getmol ' + file + ' 1\n %delete ' +\
                       rangeList[1] + '\n %delete ' + rangeList[0] + '\n %makmol ' +\
                       file + '\n ' + file + '\n y \n' + rangeList[0] + ' ' +\
                       endNewRange + ' 0 \n \n end y \n'
            startUpFil = open(tempDirectory + '/STARTUP.FIL', 'w')
            startUpFil.write(whatToDo)
            startUpFil.close()
            os.system(whatIfExe)
        os.chdir(workingDirectory)
    else:
        if verbose: print 'all residues chosen => I do not have to delete residues using WHATIF.'


    ###############################################################################
    #3. run PROCHECK:
    if procheckOnOff:
        if verbose: print 'starting PROCHECK:'
        reportFN = os.path.join(workingDirectory, FILENAME_REPORT_PROCHECK)
        if postscriptFilesOnOff:
            postscriptDir = workingDirectory
        else:
            postscriptDir = None
        try:
            RunProcheck(pdbList=fileList,\
                        tmpDir=tempDirectory,\
                        reportFN=reportFN,\
                        postscriptDir=postscriptDir,\
                        trashDirectory = trashDirectory,\
                        procheckExe = procheckExe,\
                        cshExe = cshExe,\
                        pdbSuffix = pdbSuffix,\
                        procheckResolution = procheckResolution,\
                        slowNetworks = slowNetworks)
        except:
            if verbose: print 'problem occurred while running procheck, skipping it.'
            if verbose: print '-'*60
            traceback.print_exc(file=sys.stdout)
            if verbose: print '-'*60
            procheckOnOff = 0 #overwriting procheckOnOff !!!

        os.chdir(workingDirectory) #just to make sure
        
    ###############################################################################
    #4. run PROSA:
    if prosaOnOff:
        if verbose: print 'starting PROSA-II:'
        os.chdir(workingDirectory) #just to make sure
        reportFN = os.path.join(workingDirectory, FILENAME_REPORT_PROSA)
        try:
            RunProsa(fileList, tempDirectory, reportFN,\
                     trashDirectory = trashDirectory,\
                     prosaExe = prosaExe,\
                     slowNetworks = slowNetworks)
        except:
            if verbose: print 'problem occurred while running prosa, skipping it.'
            if verbose: print '-'*60
            traceback.print_exc(file=sys.stdout)
            if verbose: print '-'*60
            prosaOnOff = 0 #overwriting prosaOnOff !!!

        os.chdir(workingDirectory) #just to make sure

    ###############################################################################
    #5. run WHATCHECK:
    if whatifOnOff:
        if verbose: print 'starting WHATCHECK:'
        reportFN = os.path.join(workingDirectory, FILENAME_REPORT_WHATIF_FULL)
        checkdbFN = os.path.join(workingDirectory, FILENAME_REPORT_WHATIF)
        try:
            RunWhatif(pdbList=fileList,\
                      tmpDir=tempDirectory,\
                      fulchkReportFN=reportFN,\
                      checkdbFN=checkdbFN,\
                      trashDirectory = trashDirectory,\
                      whatIfExe = whatIfExe,\
                      slowNetworks = slowNetworks)
        except:
            if verbose: print 'problem occurred while running whatcheck, skipping it.'
            if verbose: print '-'*60
            traceback.print_exc(file=sys.stdout)
            if verbose: print '-'*60
            whatifOnOff = 0 #overwriting whatifOnOff !!!

        os.chdir(workingDirectory) #just to make sure

    ###############################################################################
    # BARDIAUX 2.3
    #5. run Molprobity CLASHLIST:
    if clashlistOnOff:
        if verbose: print 'starting Molprobity Clashlist:'
        reportFN = os.path.join(workingDirectory, FILENAME_REPORT_MP_CLASHSCORE)
        
        mpcl = MolprobityClashlist(clashlistExe)

        os.chdir(tempDirectory)
        try:
            clashlist_scores = mpcl.runClashlist(fileList, reportFN)
        except:
            if verbose: print 'problem occurred while running clashlist, skipping it.'
            if verbose: print '-'*60
            traceback.print_exc(file=sys.stdout)
            if verbose: print '-'*60
            clashlistOnOff = 0 #overwriting clashlistOnOff !!!
            
        os.chdir(workingDirectory) #just to make sure

    ###############################################################################
    #6. parsing PROCHECK output:
    if procheckOnOff:
        reportFN = os.path.join(workingDirectory, FILENAME_REPORT_PROCHECK)
        reportHandle = open(reportFN)
        reportLines = reportHandle.readlines()
        reportString = string.join(reportLines)
        reportHandle.close()

        inputFileName = re.compile('>>>-----.*?\n.*?\n\s*\|\s*(\S+)\s+')
        residues = re.compile('(\d+)\s*residues\s\|')
        ramachandranPlot = re.compile('Ramachandran\splot:\s*(\d+\.\d+)%\s*core\s*(\d+\.\d+)%\s*allow\s*(\d+\.\d+)%\s*gener\s*(\d+\.\d+)%\s*disall')
        labelledAll = re.compile('Ramachandrans:\s*(\d+)\s*.*?out\sof\s*(\d+)')
        labelledChi = re.compile('Chi1-chi2\splots:\s*(\d+)\s*.*?out\sof\s*(\d+)')
        badContacts = re.compile('Bad\scontacts:\s*(\d+)')
        gFactors = re.compile('G-factors\s*Dihedrals:\s*([0-9-+.]+)\s*Covalent:\s*([0-9-+.]+)\s*Overall:\s*([0-9-+.]+)')
        beginSum = re.compile('----------<<<')

        procheckNameList = []
        residuesList = []
        ramCoreList = []
        ramAllowList = []
        ramGenerList = []
        ramDisallList = []
        labAllList = []
        labAllOutList = []
        labChiList = []
        labChiOutList = []
        badConList = []
        gFacDihList = []
        gFacCovList = []
        gFacOveList = []

        startPosition = 0
        while  beginSum.search(reportString, startPosition):
            startPosition = beginSum.search(reportString, startPosition).end()
            if inputFileName.search(reportString, startPosition):
                matchedInputFileName = inputFileName.search(reportString, startPosition)
                matchedResidues = residues.search(reportString, startPosition)
                matchedRamachandranPlot = ramachandranPlot.search(reportString, startPosition)
                matchedLabelledAll = labelledAll.search(reportString, startPosition)
                matchedLabelledChi = labelledChi.search(reportString, startPosition)
                matchedbadContacts = badContacts.search(reportString, startPosition)
                matchedGFactors = gFactors.search(reportString, startPosition)

                procheckNameList.append(matchedInputFileName.group(1))
                residuesList.append(string.atof(matchedResidues.group(1)))
                ramCoreList.append(string.atof(matchedRamachandranPlot.group(1)))
                ramAllowList.append(string.atof(matchedRamachandranPlot.group(2)))
                ramGenerList.append(string.atof(matchedRamachandranPlot.group(3)))
                ramDisallList.append(string.atof(matchedRamachandranPlot.group(4)))
                labAllList.append(string.atof(matchedLabelledAll.group(1)))
                labAllOutList.append(string.atof(matchedLabelledAll.group(2)))
                labChiList.append(string.atof(matchedLabelledChi.group(1)))
                labChiOutList.append(string.atof(matchedLabelledChi.group(2)))
                badConList.append(string.atof(matchedbadContacts.group(1)))
                gFacDihList.append(string.atof(matchedGFactors.group(1)))
                gFacCovList.append(string.atof(matchedGFactors.group(2)))
                gFacOveList.append(string.atof(matchedGFactors.group(3)))




    ###############################################################################
    # 7. parsing PROSA-II output:
    if prosaOnOff:
        os.chdir(workingDirectory)
        reportFN = os.path.join(workingDirectory, FILENAME_REPORT_PROSA)
        reportHandle = open(reportFN)

        proteinName = re.compile('#Protein:\s*(\S+)\s*')
        comment = re.compile('#')

        totalEnergy = []
        prosaFiles = []
        prosaEnergies = []


        for line in reportHandle.readlines():
            matchedProteinName = proteinName.match(line)
            if matchedProteinName:
                #for the first line, insert if-statement:
                if len(totalEnergy) == 0:
                    prosaFiles.append(matchedProteinName.group(1))
                    continue
                #otherwise, write out the data:
                desc = Descriptive()
                desc.addData(totalEnergy)
                prosaEnergies.append(desc.getMean())
                #initialization for the next file:
                prosaFiles.append(matchedProteinName.group(1))
                totalEnergy = []
            elif comment.search(line):
                continue
            elif string.strip(line) == '':
                continue
            else:
                totalEnergy.append(string.atof(string.split(line)[3]))

        reportHandle.close()


    ###############################################################################
    # 8. parsing WHATCHECK output:
    if whatifOnOff:
        #specify filenames:
        fullchkFN = os.path.join(workingDirectory, FILENAME_REPORT_WHATIF_FULL)
        fullchkHandle = open(fullchkFN)

        #regular expressions:
        reCompound = re.compile('====\s*Compound\s*code\s+(\S*)\s+===')
        reRamachandran = re.compile('Ramachandran\s*Z-score\s*:\s*([0-9.Ee-]+)')
        reChi12 = re.compile('chi-1\S*chi-2\s*correlation\s*Z-score\s*:\s*([0-9.Ee-]+)')
        reInsideOutside = re.compile('inside\S*outside\s*RMS\s*Z-score\s*:\s*([0-9.Ee-]+)')
        re39 = re.compile('#\s*39')
        reOnedirection = re.compile('one\s*direction.')
        re40 = re.compile('#\s*40\s*#')
        reBB = re.compile('Backbone\s*conformation\s*Z-score\s*:\s*([0-9.Ee-]+)')
        re1st = re.compile('1st\s*generation\s*packing\s*quality\s*:\s*([0-9.Ee-]+)')
        re2nd = re.compile('2nd\s*generation\s*packing\s*quality\s*:\s*([0-9.Ee-]+)')
        reAppearance = re.compile('Ramachandran\s*plot\s*appearance\s*:\s*([0-9.Ee-]+)')
        reRotamer = re.compile('chi-1\S*chi-2\s*rotamer\s*normality\s*:\s*([0-9.Ee-]+)')
        reBBconformation = re.compile('Backbone\s*conformation\s*:\s*([0-9.Ee-]+)')
        reBondLenghts = re.compile('Bond\s*lengths\s*:\s*([0-9.Ee-]+)')
        reBondAngles = re.compile('Bond\s*angles\s*:\s*([0-9.Ee-]+)')
        reOmegaAngle = re.compile('Omega\s*angle\s*restraints\s*:\s*([0-9.Ee-]+)')
        reSideChain = re.compile('Side\s*chain\s*planarity\s*:\s*([0-9.Ee-]+)')
        reImproper = re.compile('Improper\s*dihedral\s*distribution\s*:\s*([0-9.Ee-]+)')
        reInOut = re.compile('Inside\S*Outside\s*distribution\s*:\s*([0-9.Ee-]+)')

        #initialize the output lists:
        listCompound = []
        listRamachandran = []
        listChi12 = []
        listInsideOutside = []
        list39 = []
        listBB = []
        list1st = []
        list2nd = []
        listAppearance = []
        listRotamer = []
        listBBconformation = []
        listBondLenghts = []
        listBondAngles = []
        listOmegaAngle = []
        listSideChain = []
        listImproper = []
        listInOut = []

        ## #initialize the toggle for test 39 (1=inside test39, 0=outside test39):
        ## toggle39 = 0 
        ## startCount39 = 0
        ## counter39 = 0

        #loop over the lines of the FULLCHK file:
        for eachLine in fullchkHandle.readlines():
            seCompound = reCompound.search(eachLine)
            if seCompound:
                listCompound.append(seCompound.group(1))
                toggle39 = 0
            seRamachandran = reRamachandran.search(eachLine)
            if seRamachandran:
                listRamachandran.append(string.atof(seRamachandran.group(1)))
            seChi12 = reChi12.search(eachLine)
            if seChi12:
                listChi12.append(string.atof(seChi12.group(1)))
            seInsideOutside = reInsideOutside.search(eachLine)
            if seInsideOutside:
                listInsideOutside.append(string.atof(seInsideOutside.group(1)))

        ##     se39 = re39.search(eachLine)
        ##     #in test #39 count the linenumbers:
        ##     if se39:
        ##         toggle39 = 1
        ##         continue
        ##     if toggle39:
        ##         seOnedirection = reOnedirection.search(eachLine)
        ##         if seOnedirection:
        ##             startCount39 = 1
        ##             continue
        ##     if startCount39:
        ##         #the end is reached when:
        ##         se40 = re40.search(eachLine)
        ##         if se40:
        ##             #now append the list value:
        ##             list39.append(counter39)
        ##             counter39 = 0
        ##             startCount39 = 0
        ##             toggle39 = 0
        ##             continue
        ##         if string.split(eachLine) > 5:
        ##             counter39 = counter39 + 1

            seBB = reBB.search(eachLine)
            if seBB:
                listBB.append(string.atof(seBB.group(1)))
            se1st = re1st.search(eachLine)
            if se1st:
                list1st.append(string.atof(se1st.group(1)))
            se2nd = re2nd.search(eachLine)
            if se2nd:
                list2nd.append(string.atof(se2nd.group(1)))
            seAppearance = reAppearance.search(eachLine)
            if seAppearance:
                listAppearance.append(string.atof(seAppearance.group(1)))
            seRotamer = reRotamer.search(eachLine)
            if seRotamer:
                listRotamer.append(string.atof(seRotamer.group(1)))
            seBBconformation = reBBconformation.search(eachLine)
            if seBBconformation:
                listBBconformation.append(string.atof(seBBconformation.group(1)))
            seBondLenghts = reBondLenghts.search(eachLine)
            if seBondLenghts:
                listBondLenghts.append(string.atof(seBondLenghts.group(1)))
            seBondAngles = reBondAngles.search(eachLine)
            if seBondAngles:
                listBondAngles.append(string.atof(seBondAngles.group(1)))
            seOmegaAngle = reOmegaAngle.search(eachLine)
            if seOmegaAngle:
                listOmegaAngle.append(string.atof(seOmegaAngle.group(1)))
            seSideChain = reSideChain.search(eachLine)
            if seSideChain:
                listSideChain.append(string.atof(seSideChain.group(1)))
            seImproper = reImproper.search(eachLine)
            if seImproper:
                listImproper.append(string.atof(seImproper.group(1)))
            seInOut = reInOut.search(eachLine)
            if seInOut:
                listInOut.append(string.atof(seInOut.group(1)))


        checkdbFN = os.path.join(workingDirectory, FILENAME_REPORT_WHATIF)
        bmpchks, bh2chks, ba2chks = readCheckDB(checkdbFN)
        bmpchkList = []
        bh2chkList = []
        ba2chkList = []

        for eachV in bmpchks.values():
            bmpchkList.append(len(eachV))
        for eachV in bh2chks.values():
            bh2chkList.append(len(eachV))
        for eachV in ba2chks.values():
            ba2chkList.append(len(eachV))


        #some checks for consistency:
        allLists = [listCompound, listRamachandran, listChi12, listInsideOutside, listBB, list1st, list2nd, listAppearance, listRotamer, listBBconformation, listBondLenghts, listBondAngles, listOmegaAngle, listSideChain, listImproper, listInOut]

        testiii = -1
        for eachList in allLists:
            testiii = testiii + 1
            if len(eachList) != len(listCompound):
                if verbose: print 'WARNING:\n'
                if verbose: print str(len(listCompound)) + ' structures, but ' + str(len(eachList)) + ' elements (list index no ' + str(testiii) + ' )\n'
                if verbose: print eachList

    ###############################################################################
    # BARDIAUX 2.3
    #8. get clshscore statistics
    if clashlistOnOff:
        clashlist_Zscores, clashlist_percentile = mpcl.doStatistics(clashlist_scores)

        clashlistChecks = [[clashlist_scores, 'Clashscore'],\
                           [clashlist_percentile, 'Clashscore percentile'],\
                           [clashlist_Zscores, 'Clashscore Z-score']]

        #print clashlistChecks

    ###############################################################################
    #9. write the results as ASCII and LATEX files, give some comments:

    #file handle and header for ASCII
    statOutFile = os.path.join(workingDirectory, FILENAME_REPORT_TEXT)
    if verbose: print 'writing all the statistics to ' + statOutFile
    statOutHandle = open(statOutFile, 'w')
    statOutHandle.write('QualityChecks summary:\n\n')
    statOutHandle.write(commentStringForOutput + '\n\n')
    statOutHandle.write('Temporary directory:\n  ')
    statOutHandle.write(tempDirectory + '\n\n')
    statOutHandle.write('Working directory:\n  ')
    statOutHandle.write(workingDirectory + '\n\n')
    statOutHandle.write('The following files have been analysed:\n')
    for eachFile in fileList:
        statOutHandle.write('  ' + eachFile + '\n')
    statOutHandle.write('\n')



    #file handle and header for LaTeX:
    latexOutFile =  os.path.join(workingDirectory, FILENAME_REPORT_TEX)
    if verbose: print 'writing a LaTex file to ' + latexOutFile
    latexOutHandle = open(latexOutFile, 'w')

    startTexString= r"""\documentclass{article} 
    \begin{document}
    \sloppy
    \begin{table}
    \begin{center}
    \caption[]
    {Quality indices generated for your aria%s project """ % ARIA_VERSION

    startTexString = startTexString + commentStringForOutput + r""" }\label{quality_indices}
    {\footnotesize
    \begin{tabular}{l|cccc}
    % \hline
    performed checks   &value &error &min  &max \\

    \hline
    """ 
    latexOutHandle.write(startTexString)

    ###############################################################################
    #9a. PROCHECK:
    if procheckOnOff:
    #    procheckNames = [procheckNameList, 'filenames']
    #    procheckResidues = [residuesList, 'residues']
        procheckChecks = [[ramCoreList, 'most favoured regions'],\
                          [ramAllowList, 'allowed regions'],\
                          [ramGenerList, 'generously allowed regions'],\
                          [ramDisallList, 'disallowed regions'],\
                          [labAllList, 'labelled residues (all Ramachandrans)'],\
    #                      [labAllOutList, 'out of'],\
                          [labChiList, 'labelled residues (Chi1-chi2 plots)'],\
    #                      [labChiOutList, 'out of'],\
                          [badConList, 'bad contacts'],\
                          [gFacDihList, 'G-factor dihedrals'],\
                          [gFacCovList, 'G-factor covalent'],\
                          [gFacOveList, 'G-factor overall']]

        ## WR (11/03/04):

        if return_procheck_results:

            import numpy
            
            _d = {}

            for values, descr in procheckChecks:
                _d[descr] = numpy.array(values)
            
            return _d

    ###############################################################################
    #9b. WHATCHECK:
    if whatifOnOff:
        whatcheckZ = [[list1st, '1st generation packing quality Z-score (QUACHK)'],\
                      [list2nd, '2nd generation packing quality Z-score (NQACHK)'],\
                      [listAppearance, 'Ramachandran plot appearance Z-score (RAMCHK)'],\
                      [listRotamer, 'Chi-1 chi-2 rotamer normality Z-score (C12CHK)'],\
                      [listBBconformation, 'Backbone conformation Z-score (BBCCHK)']]
        whatcheckRMSZ= [[listBondLenghts, 'Bond lengths RMS Z-score (BNDCHK)'],\
                        [listBondAngles, 'Bond angles RMS Z-score (ANGCHK)'],\
                        [listOmegaAngle, 'Omega angle restraints RMS Z-score (OMECHK)'],\
                        [listSideChain, 'Side chain planarity RMS Z-score (PLNCHK)'],\
                        [listImproper, 'Improper dihedral distribution RMS Z-score (HNDCHK)'],\
                        [listInOut, 'Inside/outside distribution RMS Z-score (INOCHK)']]
        whatcheckN = [[bmpchkList, 'Inter-atomic bumps (BMPCHK)'],\
                      [bh2chkList, 'Unsatisfied hydrogen donors (BH2CHK)'],\
                      [ba2chkList, 'Unsatisfied hydrogen acceptors (BA2CHK)']]

    ###############################################################################
    #9c. PROSA:
    if prosaOnOff:
        prosaChecks = [[prosaEnergies, 'PROSA-II mean force energy [avg, sd, min, max]:']]
    else:
        prosaChecks = [None]


    ###############################################################################
    #9d. writing PROCHECK results:

    if procheckOnOff:
        statOutHandle.write('\n'+55*' ' + 'value     error       min       max\n')
        statOutHandle.write('PROCHECK:\n')
        latexOutHandle.write('PROCHECK: & & & ' + r"\\" + '\n')
        for eachL in procheckChecks:
            Desc = Descriptive()
            Desc.addData(eachL[0])
            statOutHandle.write(eachL[1] + (52-len(eachL[1]))* ' ')
            statOutHandle.write(Desc.ascii() + '\n')
            latexOutHandle.write(eachL[1] + ' & ')
            latexOutHandle.write(Desc.latex())
        statOutHandle.write('\n')

    ###############################################################################
    #9e. writing WHATCHECK results:
    if whatifOnOff:
        statOutHandle.write('WHATCHECK:\n')
        latexOutHandle.write(r'\hline' + '\n')
        latexOutHandle.write('WHATCHECK: & & & ' + r"\\" + '\n')
        for eachL in whatcheckZ:
            Desc = Descriptive()
            Desc.addData(eachL[0])
            statOutHandle.write(eachL[1] + (52-len(eachL[1]))* ' ')
            statOutHandle.write(Desc.ascii() + '\n')
            latexOutHandle.write(eachL[1] + ' & ')
            latexOutHandle.write(Desc.latex())
        #statOutHandle.write('\n')
        #statOutHandle.write('WHATCHECK RMS Z-scores:\n')
        #latexOutHandle.write('WHATCHECK RMS Z-scores: & & & ' + r"\\" + '\n')
        for eachL in whatcheckRMSZ:
            Desc = Descriptive()
            Desc.addData(eachL[0])
            statOutHandle.write(eachL[1] + (52-len(eachL[1]))* ' ')
            statOutHandle.write(Desc.ascii() + '\n')
            latexOutHandle.write(eachL[1] + ' & ')
            latexOutHandle.write(Desc.latex())
        #statOutHandle.write('\n')
        #latexOutHandle.write(' & & & ' + r"\\" + '\n')
        for eachL in whatcheckN:
            Desc = Descriptive()
            Desc.addData(eachL[0])
            statOutHandle.write(eachL[1] + (52-len(eachL[1]))* ' ')
            statOutHandle.write(Desc.ascii() + '\n')
            latexOutHandle.write(eachL[1] + ' & ')
            latexOutHandle.write(Desc.latex())
        statOutHandle.write('\n')



    ###############################################################################
    #9f. writing PROSA results:

    if prosaOnOff:
        latexOutHandle.write(r'\hline' + '\n')
        #latexOutHandle.write('PROSA-II & & & ' + r"\\" + '\n')
        statOutHandle.write('PROSA:\n')
        for eachL in prosaChecks:
            if not eachL: continue
            Desc = Descriptive()
            Desc.addData(eachL[0])
            statOutHandle.write(eachL[1] + (52-len(eachL[1]))* ' ')
            statOutHandle.write(Desc.ascii() + '\n')
            latexOutHandle.write(eachL[1] + ' & ')
            latexOutHandle.write(Desc.latex())
        statOutHandle.write('\n')


    ###############################################################################
    #9f. writing CLASHLIST results:

    if clashlistOnOff:
        latexOutHandle.write(r'\hline' + '\n')
        #latexOutHandle.write('PROSA-II & & & ' + r"\\" + '\n')
        statOutHandle.write('MOLPROBITY:\n')
        for eachL in clashlistChecks:
            if not eachL: continue
            Desc = Descriptive()
            Desc.addData(eachL[0])
            statOutHandle.write(eachL[1] + (52-len(eachL[1]))* ' ')
            statOutHandle.write(Desc.ascii() + '\n')
            latexOutHandle.write(eachL[1] + ' & ')
            latexOutHandle.write(Desc.latex())
        statOutHandle.write('\n')
        
    ###############################################################################
    #9g. closing output filehandles:
    statOutHandle.write('\n\n')
    statOutHandle.close()
    endTexString = r"""
    % \hline
    \end{tabular} \\ }
    \end{center}
    \end{table}
    \end{document}
    """
    latexOutHandle.write(endTexString)
    latexOutHandle.close()

    ###############################################################################
    #9h. running LaTeX & dvips:
    if os.path.exists(latexExe):
        os.chdir(workingDirectory)
        os.system(latexExe + ' ' + latexOutFile)
        if os.path.exists(dvipsExe):
            #wait for slow networks:
            time.sleep(slowNetworks) 
            #get .dvi filename:
            dviFN = latexOutFile[:-4] + '.dvi'
            psFN = latexOutFile[:-4] + '.ps'
            os.system(dvipsExe + ' '  + dviFN + ' -o ' + psFN)

    ###############################################################################
    #10. remove temporary files (if switched on):
    ## if removeTrashOnOff:
    ##     os.system('rm -rf ' + tempDirectory) #dangerous, handle with care!


    ###############################################################################
    #11. bye:
    endText = """QualityChecks.py finished."""
    if verbose: print endText



###############################################################################
if __name__ == '__main__':
    #if a working directory is specified on command line, take it:
    if len(sys.argv) == 2:
        workingDirectory = sys.argv[1]
    elif len(sys.argv) == 1:
        workingDirectory = os.getcwd()
        if verbose: print 'working on directory', workingDirectory
    else:
        if verbose: print 'TOO MANY COMMANDLINE ARGUMENTS GIVEN!'
        if verbose: print 'READ THE INSTRUCTIONS:', __help__
        if verbose: print 'PROGRAM ABORTED.'
    runChecks(workingDirectory = workingDirectory,\
              commentStringForOutput = '')
