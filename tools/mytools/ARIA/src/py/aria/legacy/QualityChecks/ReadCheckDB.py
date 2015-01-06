"""
a parser for checkdb output files from WHATIF

parses:
BMPCHK
BH2CHK
BA2CHK

linge@pasteur.fr
October 2001
"""

import re, string

def readCheckDB(checkdbFN):
    #checkdbFN = '/home/Bis/linge/il4_run3_water/water/CHECKDB_REPORT'
    checkdbHandle = open(checkdbFN)

    #compile some patterns:
    idPattern = re.compile('ID\s+:\s+(\S+)\s*\n')
    checkIDPattern = re.compile('CheckID\s+:\s+(\S+)\s*\n')
    namePattern = re.compile('Name\s+:\s+(\S+)\s*\n')
    valuePattern =  re.compile('Value\s+:\s+(\S+)\s*\n')

    withinBMPCHK = 0
    withinBH2CHK = 0
    withinBA2CHK = 0

    bmpchks = {}
    bh2chks = {}
    ba2chks = {}

    for line in checkdbHandle.readlines():
        #search all the patterns:
        idSearch = idPattern.match(line)
        checkIDSearch = checkIDPattern.search(line)
        nameSearch = namePattern.search(line)
        valueSearch = valuePattern.search(line)
        # get the filename -> currentID:
        if idSearch:
            currentID = idSearch.group(1)
    #        print currentID #test
            withinBMPCHK = 0
            withinBH2CHK = 0
            withinBA2CHK = 0
            continue
        # get the current check -> currentCheck:
        if checkIDSearch:
            currentCheck = checkIDSearch.group(1)
            if currentCheck == 'BMPCHK':
                withinBMPCHK = 1
                bmpchks[currentID] = []
            elif currentCheck == 'BH2CHK':
                withinBH2CHK = 1
                bh2chks[currentID] = []
            elif currentCheck == 'BA2CHK':
                withinBA2CHK = 1
                ba2chks[currentID] = []
        if withinBMPCHK:
            if valueSearch:
                bmpchks[currentID].append(valueSearch.group(1))
        elif withinBH2CHK:
            if valueSearch:
                bh2chks[currentID].append(valueSearch.group(1))
        elif withinBA2CHK:
            if valueSearch:
                ba2chks[currentID].append(valueSearch.group(1))

##     print bmpchkList
##     print bh2chkList
##     print ba2chkList

##     bmpchkList = []
##     bh2chkList = []
##     ba2chkList = []

##     for eachV in bmpchks.values():
##         bmpchkList.append(len(eachV))
##     for eachV in bh2chks.values():
##         bh2chkList.append(len(eachV))
##     for eachV in ba2chks.values():
##         ba2chkList.append(len(eachV))

##     print bmpchks
##     print bh2chks
##     print ba2chks 

    return bmpchks, bh2chks, ba2chks
