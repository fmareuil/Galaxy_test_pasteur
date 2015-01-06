"""
The general ARIA data format for the chemical shifts
together with the converting scripts to XEASY, NMRVIEW and AURELIA.
All the chemical shift data handling is done with this module.
It contains two classes:
-PpmList: the class represents an atom chemical shift list
-Atom:    represents one atom (with ppm, atomname, etc.)

Instantiate PpmList, then read and write the data with its methods:
#example code:
import PpmList
PPM = PpmList.PpmList('test', 3)  #name='test', dimension=3
PPM.ReadXeasyProt('/home/linge/aria/example/n15.prot')
PPM.WriteAureliUserInfo('/home/linge/writeppm.aui')
PPM.WriteChem('/home/linge/writeppm.chem')
PPM.WritePpm('/home/linge/writeppm.ppm')
PPM.WriteXeasyProt('/home/linge/writeppm.prot')
PPM.Stdout()
"""
__author__   = "$Author: bardiaux $"
__revision__ = "$Revision: 1.1.1.1 $"
__date__     = "$Date: 2010/03/23 15:27:24 $"

import math, os, re, string
import DeleteCnsComments, DeleteComments, SequenceList
import FortranFormat, TextFile
import AminoAcid, Nomenclature, PseudoAtom

###############################################################################
class PpmList:
    """
    Contains the chemical shifts of all atoms of one spectrum
    supplies methods for input and output in various formats
    
    public attributes:
        atomlist       list of atom objects
        atomdicfa      dictionary (key: (residuenumber, atomname),
                       value: atom object)
        comment        comment for the NOE list
        dimension      dimensionality of the spectrum (2, 3 or 4)
        name           name of the spectrum
        fileName       name of the file from which the data come from
        residuelist    list of residue names
        
    public methods:
        AddAtom               adds one atomobject 
        AddSequence           adds sequence information
        ConvertWildcards      converts wildcards to single atomnames
        ReadAnsig             reads an ANSIG crosspeaks export file
        ReadAureliaUserInfo   reads an user info file of AURELIA
        ReadChem              reads an ARIA .chem file
        ReadNmrView           reads a NMRView .out file
        ReadPipp              reads a PIPP .shifts file
        ReadPpm               reads an ARIA .ppm file
        ReadRegine            reads an Regine chemical shift file
        ReadSparky            reads a Sparky .list file
        ReadXeasyProt         reads a XEASY .prot file
        ReadXeasyProtSeq      reads a XEASY .prot file with a sequence file
        RemoveAtom            removes one atomobject
        RemoveDoubleQuotes    removes double quotation marks from all the atomnames
        Stdout                writes the data to stdout
        StdoutAll             writes all the data to stdout
        WriteAureliaUserInfo  writes an user info file of AURELIA
        WriteBioMagResBank    writes a file for BioMagResBank deposition
        WriteChem             writes an .chem file
        WriteId               writes an ARIA .id file
        WritePpm              writes an ARIA .ppm file
        WriteSparky           writes a Sparky .list file
        WriteXeasyProt        writes a XEASY .prot file
        WriteXML2String       returns a string containing the sequence in XML format
        WriteXML2File         writes the sequence to an XML file  
    NOTE:
        all fields are treated internally as strings, except:
        atomnames are always tuples of strings
        dimension is always an integer
    """
    def __init__(self, name = 'chemical shifts' , dimension = 0,\
                 comment = 'Aria Data Format'):
        self.atomlist = []
        self.atomdicfa = {}
        self.comment = comment
        self.dimension = dimension
        self.fileName = ''
        self.name = name
        self.residuelist = []
        

    def AddAtom(self, aobj):
        self.atomlist.append(aobj)
        self.atomdicfa[(aobj.residuenumber, aobj.atomname)] = aobj

    def AddSequence(self, residueList):
        """
        residueList contains the sequence, e.g. like
        ['ARG', 'GLY', 'LEU']
        in 3-letter code
        sequence number 1 is residueList[0]
        """
        self.residuelist = residueList
        for atom in self.atomlist:
            if atom.aminoacid:
                if len(residueList) < string.atoi(atom.residuenumber):
                    print 'WARNING: the provided residue list contains only ' +\
                          str(len(residueList)) + ' residues,'
                    print '         but the ppm list contains a residuenumber ' + atom.residuenumber
                else:
                    if atom.aminoacid != residueList[string.atoi(atom.residuenumber) - 1]:
                        print 'WARNING: ' + atom.aminoacid + ' in ' + atom.residuenumber +\
                              'is not the same as ' +\
                              residueList[string.atoi(atom.residuenumber) - 1] +\
                              ' in  the provided residue list'
            else:
#                print 'adding ' + residueList[string.atoi(atom.residuenumber) - 1] + ' ' + atom.residuenumber #test
                try:
                    atom.aminoacid = str(residueList[string.atoi(atom.residuenumber) - 1])
                except IndexError:
                    print 'WARNING:  the provided residue list contains only ' +\
                          str(len(residueList)) + ' residues,'
                    print '          but the ppm list contains a residuenumber ' + atom.residuenumber
                    
    def ConvertWildcards(self):
        for entry in self.atomlist:
            newList = []
            for atom in entry.atomname:
                if not entry.aminoacid:
                    newList.append(atom)
                else:
                    convertedTuple = PseudoAtom.Pseudo2IupacTuple(entry.aminoacid, atom)
                    for converted in convertedTuple:
                        newList.append(converted)
#                    print atom, newList #test
            entry.atomname = newList
            

    def ReadAnsig(self, fileName):
        """
        reads an ANSIG crosspeaks export file or ANSIG crosspeaks storage file
        just extracts the chemical shift assignments
        """
        if _DoesFileExist(fileName) == 0:
            return
        print 'reading the ANSIG crosspeaks export file:\n ', fileName
        
        #important - clean atomlist and atomdicfa:
        self.atomlist = []
        self.atomdicfa = {}
        self.fileName = fileName
        
        cpHandle = TextFile.TextFile(fileName)
        
        #there are two possible cases:
        #1. ANSIG v3.3 export crosspeaks file
        #   second line contains 'totalNumber dimensionality'
        #   spectrumName appears in every line
        #2. ANSIG v3.3 storage crosspeaks file
        #   second line contains 'spectrumName totalNumber dimensionality'
        #   spectrumName appears only in the header
        
        #get total number of crosspeaks and dimensionality from the first two lines:
        eachLine = cpHandle.readline()
        eachLine = cpHandle.readline()
        
        totAndDim = string.split(eachLine)
        totalNumber = totAndDim[-2]
        self.dimension = int(totAndDim[-1])
        
        if len(totAndDim) == 2:  #for the ANSIG crosspeaks files           
            format2 = FortranFormat.FortranFormat('3E13.6,A12,7I6,6A4')
            format3 = FortranFormat.FortranFormat('4E13.6,A12,9I6,9A4')
            format4 = FortranFormat.FortranFormat('5E13.6,A12,11I6,12A4')

            #read the assignments:
            ppmdic = {}
            if self.dimension == 2:
                for eachLine in cpHandle:
                    if len(eachLine) < 40:
                        continue
                    line = FortranFormat.FortranLine(eachLine, format2)
                    # {(residuenumber, 3-lettercode, atomname): ppm}
                    ppmdic[(string.strip(str(line[11])), string.strip(string.upper(line[13])),\
                            string.strip(line[15]))] = string.strip(str(line[0]))
                    ppmdic[(string.strip(str(line[12])), string.strip(string.upper(line[14])),\
                            string.strip(line[16]))] = string.strip(str(line[1]))

            elif self.dimension == 3:
                for eachLine in cpHandle:
                    if len(eachLine) < 40:
                        continue
                    line = FortranFormat.FortranLine(eachLine, format3)
##                    for eachBla in line: #test
##                        print eachBla, #test
##                    print '' #test
                    # {(residuenumber, 3-lettercode, atomname): ppm}
                    ppmdic[(string.strip(str(line[14])), string.strip(string.upper(line[17])),\
                            string.strip(line[20]))] = string.strip(str(line[0]))
                    ppmdic[(string.strip(str(line[15])), string.strip(string.upper(line[18])),\
                            string.strip(line[21]))] = string.strip(str(line[1]))
                    ppmdic[(string.strip(str(line[16])), string.strip(string.upper(line[19])),\
                            string.strip(line[22]))] = string.strip(str(line[2]))

            elif self.dimension == 4:
                for eachLine in cpHandle:
                    if len(eachLine) < 40:
                        continue
                    line = FortranFormat.FortranLine(eachLine, format4)
                    # {(residuenumber, 3-lettercode, atomname): ppm}
                    ppmdic[(string.strip(str(line[17])), string.strip(string.upper(line[21])),\
                            string.strip(line[25]))] = string.strip(str(line[0]))
                    ppmdic[(string.strip(str(line[18])), string.strip(string.upper(line[22])),\
                            string.strip(line[26]))] = string.strip(str(line[1]))
                    ppmdic[(string.strip(str(line[19])), string.strip(string.upper(line[23])),\
                            string.strip(line[27]))] = string.strip(str(line[2]))
                    ppmdic[(string.strip(str(line[20])), string.strip(string.upper(line[24])),\
                            string.strip(line[28]))] = string.strip(str(line[3]))

            #set some default values:
            segid = None
            shifterror = '0.0'

            #transfer the data from ppmdic to the Atom objects:
            #don't use residuenumber 0 => all the ppms of residue 0 are not used!
##            print ppmdic.keys() #test
            for keyword in ppmdic.keys():
                if (string.strip(keyword[0]) == '') or \
                   (string.strip(keyword[0]) == '0') or \
                   (len(string.strip(keyword[1])) != 3) or \
                   (string.strip(keyword[2]) == ''):
                    continue
                ATOM = Atom(keyword[0], keyword[1], segid, (keyword[2],),\
                            ppmdic[keyword], shifterror)
                self.AddAtom(ATOM)
          
        else:  #for the ANSIG storage files
            format2 = FortranFormat.FortranFormat('3E13.6,5I6,6A4')
            format3 = FortranFormat.FortranFormat('4E13.6,7I6,9A4')
            format4 = FortranFormat.FortranFormat('5E13.6,9I6,12A4')
            spectrumName = totAndDim[0]

            #read the assignments:
            ppmdic = {}
            if self.dimension == 2:
                for eachLine in cpHandle:
                    if len(eachLine) < 40:
                        continue
                    line = FortranFormat.FortranLine(eachLine, format2)
                    # {(residuenumber, 3-lettercode, atomname): ppm}
                    ppmdic[(string.strip(str(line[8])), string.strip(string.upper(line[10])),\
                            string.strip(line[12]))] = string.strip(str(line[0]))
                    ppmdic[(string.strip(str(line[9])), string.strip(string.upper(line[11])),\
                            string.strip(line[13]))] = string.strip(str(line[1]))

            elif self.dimension == 3:
                for eachLine in cpHandle:
                    if len(eachLine) < 40:
                        continue
                    line = FortranFormat.FortranLine(eachLine, format3)
                    # {(residuenumber, 3-lettercode, atomname): ppm}
                    ppmdic[(string.strip(str(line[10])), string.strip(string.upper(line[13])),\
                            string.strip(line[16]))] = string.strip(str(line[0]))
                    ppmdic[(string.strip(str(line[11])), string.strip(string.upper(line[14])),\
                            string.strip(line[17]))] = string.strip(str(line[1]))
                    ppmdic[(string.strip(str(line[12])), string.strip(string.upper(line[15])),\
                            string.strip(line[18]))] = string.strip(str(line[2]))

            elif self.dimension == 4:
                for eachLine in cpHandle:
                    if len(eachLine) < 40:
                        continue
                    line = FortranFormat.FortranLine(eachLine, format4)
                    # {(residuenumber, 3-lettercode, atomname): ppm}
                    ppmdic[(string.strip(str(line[14])), string.strip(string.upper(line[18])),\
                            string.strip(line[22]))] = string.strip(str(line[0]))
                    ppmdic[(string.strip(str(line[15])), string.strip(string.upper(line[19])),\
                            string.strip(line[23]))] = string.strip(str(line[1]))
                    ppmdic[(string.strip(str(line[16])), string.strip(string.upper(line[20])),\
                            string.strip(line[24]))] = string.strip(str(line[2]))
                    ppmdic[(string.strip(str(line[17])), string.strip(string.upper(line[21])),\
                            string.strip(line[25]))] = string.strip(str(line[3]))

            #set some default values:
            segid = None
            shifterror = '0.0'

            #transfer the data from ppmdic to the Atom objects:
            #don't use residuenumber 0 => all the ppms of residue 0 are not used!
            for keyword in ppmdic.keys():
                if (string.strip(keyword[0]) == '') or \
                   (string.strip(keyword[0]) == '0') or \
                   (len(string.strip(keyword[1])) != 3) or \
                   (string.strip(keyword[2]) == ''):
                    continue
                ATOM = Atom(keyword[0], keyword[1], segid, keyword[2],\
                            ppmdic[keyword], shifterror)
                self.AddAtom(ATOM)

            
    def ReadAureliaUserInfo(self, fileName):
        """
        reads an Aurelia User Info File
        """
        if _DoesFileExist(fileName) == 0:
            return
        print 'reading an Aurelia User Info File:\n  ', fileName
        print 'We always use the following format for the User Info Files:'
        print '  # 8.17 NH 7 2FMR'
        print '  # ppm atomname residuenumber segid'
        print '  segid should contain 4 letters or should be blank'
        print '  other formats can not be read in by this method!'
        #important - clean atomlist and atomdicfa:
        self.atomlist = []
        self.atomdicfa = {}
        self.fileName = fileName
        auihandle = TextFile.TextFile(fileName)
        for line in auihandle:
            linelist = string.split(line)
            if len(linelist) < 4:
                continue
            ATOM = Atom()
            ATOM.shift = linelist[1]
            ATOM.atomname = linelist[2]
            ATOM.residuenumber = linelist[3]
            try:
                ATOM.segid = linelist[4]
            except:
                ATOM.segid = '    '
            self.AddAtom(ATOM)
        auihandle.close()


    def ReadSparky(self, fileName):
        """
        reads a sparky file
        Sparky format

        E36 H1' 1H 6.029
        
        """
        if _DoesFileExist(fileName) == 0:
            return

        #important - clean atomlist and atomdicfa:
        self.atomlist = []
        self.atomdicfa = {}
        self.fileName = fileName
        
        completelist = open(fileName)
        for line in completelist.readlines():
            linelist = string.split(line)

            #in case some lines are not complete 
            if len(linelist) < 4:
                continue
            if linelist[0] == "Group":
                continue
            
            # create an istance of class Atom
            ATOM = Atom()
##             print line
##             print linelist   # ?

            ATOM.aminoacid = AminoAcid.AminoAcid(linelist[0][0])[1]
            ATOM.residuenumber = linelist[0][1:]
            ATOM.atomname = (linelist[1],)
            ATOM.atomtype = linelist[2]
            ATOM.shift = linelist[3]
            ATOM.shifterror = linelist [4]
            
            self.AddAtom(ATOM)
        completelist.close()

##     def ReadSparky(self, fileName):
##         """
##         reads a sparky file
##         """
##         if _DoesFileExist(fileName) == 0:
##             return

##         #important - clean atomlist and atomdicfa:
##         self.atomlist = []
##         self.atomdicfa = {}
##         self.fileName = fileName
##         auihandle = open(fileName)
##         for line in auihandle.readlines():
##             linelist = string.split(line)
##             if len(linelist) < 4:
##                 continue
##             ATOM = Atom()
##             print line
##             print linelist
            
##             ATOM.shift = linelist[1]
##             ATOM.atomname = linelist[2]
##             ATOM.residuenumber = linelist[3]
##             try:
##                 ATOM.segid = linelist[4]
##             except:
##                 ATOM.segid = '    '
##             self.AddAtom(ATOM)
##         auihandle.close()

    
    def ReadChem(self, fileName):
        """
        reads an ARIA .chem file

        example for one line:
        do ( store1 = 60.709 ) ( resid 1 and name CA )

        comments like '! bla until lineend' or '{bla}' are neglected
        use DeleteCnsComments to parse comments
        """
        if _DoesFileExist(fileName) == 0:
            return
        print 'reading an ARIA chemical shift file', fileName

        #important - clean atomlist and atomdicfa:
        self.atomlist = []
        self.atomdicfa = {}
        self.fileName = fileName

        #get the file without the comments:
        bigstring = DeleteCnsComments.GetString(fileName)
#        print bigstring #test
        #split the string in lines:
        lines = string.split(bigstring, '\n')

        ppmAssign = re.compile('do\s*\(\s*store1\s*=\s*([0-9-+.Ee]+)\s*\)\s*\(\s*resid\s*(\d+)\s*and\s*name\s*(\S+)\s*\)')
        
        for line in lines:
            #for wrong or empty lines:
            if len(line) < 20:
                continue
#            print line #test
            linelist = string.split(line)
            ATOM = Atom()
            ppmSearch = ppmAssign.search(line)

            # new for store5 * store6 -> skip if it's not store1
            # and pattern doesn't match:
            if not ppmSearch:
                continue

            ATOM.residuenumber = ppmSearch.group(2)
            ATOM.aminoacid = None
            ATOM.segid = None
            ATOM.atomname = (ppmSearch.group(3), )
            ATOM.shift = ppmSearch.group(1)
            ATOM.shifterror = '0.0'
            self.AddAtom(ATOM)
    
    def ReadNmrView(self, fileName):
        """
        reads a NMRView chemical shift file

        each line contains:
        residuenumber.atomname shift

        example (one line of an .out file):
        2.HA       4.166 0

        comments like '# bla until lineend' are neglected
        use DeleteComments to parse comments
        """
        if _DoesFileExist(fileName) == 0:
            return
        print 'reading a NMRView .out file', fileName

        #important - clean atomlist and atomdicfa:
        self.atomlist = []
        self.atomdicfa = {}
        self.fileName = fileName

        #get the file without the comments:
        bigstring = DeleteComments.GetString(fileName)

        #split the string in lines:
        lines = string.split(bigstring, '\n')

        for line in lines:
            linelist = string.split(line)
            #for wrong or empty lines:
            if len(linelist) < 3:
                continue
            ATOM = Atom()
            firstFieldList = string.split(linelist[0], '.')
            ATOM.residuenumber = firstFieldList[0]
            ATOM.aminoacid = None
            ATOM.segid = None
            ATOM.atomname = (PseudoAtom.Pseudo2Atom(firstFieldList[1]),)
            ATOM.shift = linelist[1]
            ATOM.shifterror = '0.0'
            self.AddAtom(ATOM)

    def ReadPipp(self, fileName): 
        """
        reads a PIPP .shifts file
        """
        if _DoesFileExist(fileName) == 0:
            return
        print 'reading a PIPP .shifts file', fileName
        
        #important - clean atomlist and atomdicfa:
        self.atomlist = []
        self.atomdicfa = {}
        self.fileName = fileName
        assignNow = 0

        shiftHandle = TextFile.TextFile(fileName)
        
        #split the string in lines:
        for line in shiftHandle:
            if line[0] == '#':
                continue
            linelist = string.split(line)
                
            if linelist == []:
                continue
            elif linelist[0] == 'RES_ID':
                residueID = linelist[1]
            elif linelist[0] == 'RES_TYPE':
                residueType = linelist[1]
            elif linelist[0] == 'SPIN_SYSTEM_ID':
                spinSystemID = linelist[1]
                assignNow = 1
                continue
            elif linelist[0] == 'END_RES_DEF':
                assignNow = 0
            elif linelist[0] == 'SHIFT_FL_FRMT':
                fileFormat = linelist[1]
            elif linelist[0] == 'FIRST_RES_IN_SEQ':
                firstResidue = linelist[1]

            if assignNow:
                ATOM = Atom()
                ATOM.residuenumber = residueID
                ATOM.aminoacid = residueType
                ATOM.atomname = tuple(string.split(linelist[0], '|'))
                ATOM.shift = linelist[1]
                ATOM.shifterror = '0.0'    #default
                ATOM.segid = None         #default
                self.AddAtom(ATOM)

        print fileName, 'used the format', fileFormat, 'with first residue', firstResidue
        
    
    
    def ReadPpm(self, fileName):
        """
        reads an ARIA .ppm file
        each line contains:
        residuenumber
        aminoacid
        segid
        atomname
        shift
        shifterror

        unused fields are set to None
        
        uses DeleteComments, the following comments are recognized:
        { bla }
        { bla { bla } bla }
        ! bla until end of line
        # bla until end of line
        
        """
        if _DoesFileExist(fileName) == 0:
            return
        print 'reading an ARIA. ppm file', fileName
        #important - clean atomlist and atomdicfa:
        self.atomlist = []
        self.atomdicfa = {}
        self.fileName = fileName
        #get the file without the comments:
        bigstring = DeleteComments.GetString(fileName)
        #split the string in lines:
        lines = string.split(bigstring, '\n')
        for line in lines:
            linelist = string.split(line)
            #for wrong or empty lines:
            if len(linelist) < 6:
                continue
            ATOM = Atom()
            ATOM.residuenumber = linelist[0]
            ATOM.aminoacid = linelist[1]
            ATOM.segid = linelist[2]
            ATOM.atomname = (linelist[3], )
            ATOM.shift = linelist[4]
            ATOM.shifterror = linelist[5]
            self.AddAtom(ATOM)
        
    
            
    def ReadRegine(self, fileName):
        """
        reads an Regine chemical shift file

        format:

            ALA 12 HA 3.5131

            3-lettercode residueNumber atomName ppm

        Please note that the atomnames are converedt from IUPAC to CNS
        and vice versa!
        """
        if _DoesFileExist(fileName) == 0:
            return
        print 'reading an Regine ppm file:\n  ', fileName
        #important - clean atomlist and atomdicfa:
        self.atomlist = []
        self.atomdicfa = {}
        self.fileName = fileName
        fileStream = open(fileName)
        for eachLine in fileStream.readlines():
            lineList = string.split(eachLine)
            if len(lineList) < 4:
                continue
            ATOM = Atom()
            ATOM.shifterror = "0.0"
            if len(lineList) == 5: 
                ATOM.shifterror = lineList[4]
            ATOM.shift = lineList[3]
            ATOM.aminoacid = string.upper(lineList[0])
            ATOM.atomname = (Nomenclature.ConvertCnsProtonNames(ATOM.aminoacid, lineList[2]),)
            ATOM.residuenumber = lineList[1]
            ATOM.segid = '    '
            self.AddAtom(ATOM)
        fileStream.close()

    
    def ReadXeasyProt(self, fileName):
        """
        reads a XEASY .prot file
        uses the ReadXeasy module
        """
        #for the XEASY
        import ReadXeasy
        if _DoesFileExist(fileName) == 0:
            return
        #important - clean atomlist and atomdicfa:
        self.atomlist = []
        self.atomdicfa = {}
        print 'reading the .prot file', fileName
        self.fileName = fileName
        XPROT = ReadXeasy.XeasyProt()
        XPROT.ReadProt(fileName)
        for EACH in XPROT.atomlist:
            ATOM = Atom()
            ATOM.residuenumber = EACH.fragmentnumber
            ATOM.atomname = EACH.ariaatomname
            if EACH.shift == '999.000':
                ATOM.shift = None
            else:
                ATOM.shift = EACH.shift
            ATOM.shifterror = EACH.shifterror
            ATOM.xeasyatomname = EACH.xeasyatomname
            ATOM.xeasyatomnumber = EACH.atomnumber
            self.AddAtom(ATOM)
        self.RemoveDoubleQuotes() #conversion of " into ''


    def ReadXeasyProtSeq(self, fileName, seqFileName):
        """
        reads a XEASY .prot file and a sequence file
        uses the ReadXeasy and SequenceList modules
        """
        #for the XEASY files:
        import ReadXeasy
        if _DoesFileExist(fileName) == 0:
            return
        if _DoesFileExist(seqFileName) == 0:
            return

        #read the sequence:
        print 'reading the .prot file', seqFileName
        SL = SequenceList.SequenceList()
        SL.ReadSeq(seqFileName)
        
        #important - clean atomlist and atomdicfa:
        self.atomlist = []
        self.atomdicfa = {}
        print 'reading the .prot file', fileName
        self.fileName = fileName
        XPROT = ReadXeasy.XeasyProt()
        XPROT.ReadProt(fileName)
        for EACH in XPROT.atomlist:
            ATOM = Atom()
            ATOM.residuenumber = EACH.fragmentnumber
            try:
                ATOM.aminoacid = SL.aalist[string.atoi(EACH.fragmentnumber)-1]
            except:
                pass
            ATOM.atomname = EACH.ariaatomname
            if EACH.shift == '999.000':
                ATOM.shift = None
            else:
                ATOM.shift = EACH.shift
            ATOM.shifterror = EACH.shifterror
            ATOM.xeasyatomname = EACH.xeasyatomname
            ATOM.xeasyatomnumber = EACH.atomnumber
            self.AddAtom(ATOM)
        self.RemoveDoubleQuotes() #conversion of " into ''


    def RemoveAtom(self, aobj):
        """
        removes a given atom object from atomlist and atomdic
        """
        if aobj in self.atomlist:
            self.atomlist.remove(aobj)
        if self.atomdicfa.has_key((aobj.residuenumber, aobj.atomname)):
            del(self.atomdicfa[(aobj.residuenumber, aobj.atomname)])

        
    def RemoveDoubleQuotes(self):
        """
        removes double quotation marks from all the atomnames
        this is useful e.g. for XEASY RNA atomnames
        all the double quotation marks are replaced by two single
        quotation marks
        """
        doubleQuote =  re.compile('"')
        for eachP in self.atomlist:
            tmpList = []
            for eachA in eachP.atomname:
                if eachA:
                    tmpList.append(doubleQuote.sub("''", eachA))
                else:
                    tmpList.append(eachA)
            eachP.atomname = tuple(tmpList)
                

    
    def Stdout(self):
        print self.__repr__()
        
    def __repr__(self):
        outString = 'residuenumber aminoacid segid atomname ppm dppm comment \n'
        for EACH in self.atomlist:
            outString = outString + ' ' + str(EACH.residuenumber) + ' ' + str(EACH.aminoacid) + ' ' +\
                        str(EACH.segid) + ' ' + str(EACH.atomname)+ ' ' + str(EACH.shift)+ ' ' +\
                        str(EACH.shifterror) + ' ' + str(EACH.comment) + '\n'
        return outString
        

    def StdoutAll(self):
        """prints out all the attributes + the xeasy stuff"""
        print 'residuenumber aminoacid segid atomname ppm dppm comment',\
              'xatomname xatomnumber'
        for EACH in self.atomlist:
            print EACH.residuenumber, EACH.aminoacid, EACH.segid,\
                  EACH.atomname, EACH.shift, EACH.shifterror, EACH.comment,\
                  EACH.xeasyatomname, EACH.xeasyatomnumber
        
    
    def WriteAureliaUserInfo(self, fileName):
        """
        writes an Aurelia User Info File
        Does not use the chemical shifts with values 999.000
        if there are more than one atomnames in the tuple, only the
        first will be used!
        """
        print 'writing an Aurelia User Info File:\n  ', fileName
        print 'We always use the following format for the User Info Files:'
        print '  # 8.17 NH 7 2FMR'
        print '  # ppm atomname residuenumber segid'
        print '  segid should contain 4 letters or 4 spaces'
        auihandle = TextFile.TextFile(fileName, 'w')
        for EACH in self.atomlist:
            #those with 999.000 don't have an assignment:
            if EACH.shift != '999.000':
                if EACH.segid == None:
                    outsegid = '    '  #4 spaces
                else:
                    outsegid = EACH.segid
                auihandle.write('# ' + EACH.shift + ' ' +\
                                EACH.atomname[0] + ' ' +\
                                EACH.residuenumber +\
                                outsegid + '\n')
        
    
    def WriteAureliaUserInfo(self, fileName):
        """
        writes an Aurelia User Info File
        Does not use the chemical shifts with values 999.000
        if there are more than one atomnames in the tuple, only the
        first will be used!
        """
        print 'writing an Aurelia User Info File:\n  ', fileName
        print 'We always use the following format for the User Info Files:'
        print '  # 8.17 NH 7 2FMR'
        print '  # ppm atomname residuenumber segid'
        print '  segid should contain 4 letters or 4 spaces'
        auihandle = TextFile.TextFile(fileName, 'w')
        for EACH in self.atomlist:
            #those with 999.000 don't have an assignment:
            if EACH.shift != '999.000':
                if EACH.segid == None:
                    outsegid = '    '  #4 spaces
                else:
                    outsegid = EACH.segid
                auihandle.write('# ' + EACH.shift + ' ' +\
                                EACH.atomname[0] + ' ' +\
                                EACH.residuenumber +\
                                outsegid + '\n')
        

    
    def _Write_BioMagResBank_Polymer_Saveframe( self, fileHandle ):
        """write out the amino acid sequence in the future"""
        
        ## Number of residues in residue loop per line
        ## Please specify as a real
        number_residues_per_line            =  5.0
        number_residues_per_line_1Letter    = 20.0
        
        ## Tuple of three letter codes        
        aaSequence          = ()
        aaSequence1Letter   = ""
        oldResidueNumber    = -999
        aaSequenceCount     = 0

        tmpStr = """
    ##############################
    #  Monomeric polymers        #
    ##############################

    
save_my_protein
   _Saveframe_category     monomeric_polymer

   _Mol_type               polymer
   _Mol_polymer_class      protein
"""
        fileHandle.write( tmpStr )

        for EACHR in self.residuelist:
            aaSequenceCount     = aaSequenceCount + 1
            aaSequence          = aaSequence + ( EACHR, )
            aaSequence1Letter   = aaSequence1Letter + AminoAcid.AminoAcid( EACHR )[0]
            oldResidueNumber    = aaSequenceCount
            if ( math.modf( ( aaSequenceCount )/
                number_residues_per_line_1Letter )[0] == 0.0 ):
                aaSequence1Letter = aaSequence1Letter + '\n'

        
        if ( aaSequence == '' ):
            print 'ERROR: aminoacid sequence not known from input data!'
        else:
            resLoop = """
   
        ##############################
        #  Polymer residue sequence  #
        ##############################


"""
            resLoop = resLoop + "   _Residue_count          %s\n" % aaSequenceCount
            resLoop = resLoop + "   _Mol_residue_sequence\n;\n"
            fileHandle.write( resLoop + aaSequence1Letter + "\n;\n\n" )
            resLoop = """
     loop_
        _Residue_seq_code
        _Residue_label

        """
            iii = 1
            for eachR in aaSequence:
                resLoop = resLoop + ' %3d  %3s     ' % ( iii, eachR )
                if math.modf( (iii)/number_residues_per_line )[0] == 0.0:
                    resLoop = resLoop + '\n        '
                iii = iii + 1
            resLoop = resLoop +  '\n\n   stop_\n'
            fileHandle.write( resLoop )
        fileHandle.write("\n\nsave_\n\n")


    def _WriteBioMagResBank_Molecular_System_Saveframe( self, fileHandle ):
        
        
        tmpStr = """

##################################
#  Molecular system description  #
##################################


save_my_system
   _Saveframe_category      molecular_system

   _Mol_system_name        'my system'

   loop_
      _Mol_system_component_name
      _Mol_label

       "my protein" $my_protein

   stop_

save_
"""
        fileHandle.write( tmpStr )

    
    
    def WriteBioMagResBank(self, fileName):
        """
        writes a file for BioMagResBank deposition
        
        Does not use the chemical shifts with values 999.000
        if there are more than one atomnames in the tuple, only the
        first will be used!
        """
        
        print 'Writing the BioMagResBank STAR formatted file', fileName
        bioHandle = open(fileName, 'w')

        # Header for STAR format
        bioHandle.write( "data_chemical_shift_set\n" )
        # Write a generic molecular system saveframe
        self._WriteBioMagResBank_Molecular_System_Saveframe( bioHandle )
        # Write the sequence in a polymer saveframe
        self._Write_BioMagResBank_Polymer_Saveframe( bioHandle )

        tmpStr = \
"""
    ###################################
    #  Assigned chemical shift lists  #
    ###################################

###################################################################
#       Chemical Shift Ambiguity Index Value Definitions          #
#                                                                 #
#   Index Value            Definition                             #
#                                                                 #
#      1             Unique                                       #
#      2             Ambiguity of geminal atoms or geminal methyl #
#                         proton groups                           #
#      3             Aromatic atoms on opposite sides of the ring #
#                        (e.g. Tyr HE1 and HE2 protons)           #
#      4             Intraresidue ambiguities (e.g. Lys HG and    #
#                         HD protons)                             #
#      5             Interresidue ambiguities (Lys 12 vs. Lys 27) #
#      9             Ambiguous, specific ambiguity not defined    #
#                                                                 #
###################################################################


save_chemical_shift_set
   _Saveframe_category               assigned_chemical_shifts

"""  
        # Write some comments
        tmpStr = tmpStr + \
"""
   _Details
;
Derived from the file: """
        tmpStr = tmpStr + self.fileName + """
;
   _Mol_system_component_name        "my protein"
"""
        bioHandle.write(tmpStr)
        
        
        ## Chemical shift table
        ## SegId left out for now
        commentsForPpm = """
loop_
     _Atom_shift_assign_ID
     _Residue_seq_code
     _Residue_label
     _Atom_name
     _Atom_type
     _Chem_shift_value
     _Chem_shift_value_error
     _Chem_shift_ambiguity_code

"""
        bioHandle.write(commentsForPpm)

        atomShiftAssignCounter = 1
        
        for EACH in self.atomlist:
            if EACH.shift != '999.000' and EACH.shift != None:
                # segid
                if EACH.segid:
                    outSEGID = EACH.segid
                else:
                    outSEGID = '.'  # .  instead of whitespace or None

                # residue type
                if EACH.aminoacid:
                    outAA = EACH.aminoacid
                else:
                    outAA = '.'  # .  instead of whitespace or None

                #for the ambiguity:  #TODO : ambiguity in STAR files!!!
##                 if len(EACH.atomname) > 1:
##                     ambiguity = 9
##                 else:
##                     ambiguity = 1
                ambiguity = "."

                if EACH.shifterror:
                    outError = EACH.shifterror
                else:
                    outError = '.' # .  instead of whitespace or None

                ## All upper case please
                tmpAtomName = EACH.atomname[0][:]
                tmpAtomName = string.upper( tmpAtomName )
                
                ## Convert CNS pseudo atoms to IUPAC, e.g. LEU HD1# -> MD1
                ## or leave untouched if they are not CNS pseudo atom names
                ## for the standard amino acids and nucleic acids
                tmpAtomName = Nomenclature.Convert_PseudoAtomName_CNS_2_IUPAC(
                                EACH.aminoacid, tmpAtomName )
                
                ## Convert IUPAC atoms of a methyl group it's constituing pseudo atom
                ## E.g. LEU HD11 -> MD1
                ## This does NOT convert anything else than methyl and
                ## amino groups since the atoms in these groups can not be assigned
                ## individually.
                ## The first argument is the type of pseudo atom as defined in the library
                ## file originating from the AQUA software.
                tmpAtomName = Nomenclature.Convert_IUPAC_AtomName_2_PseudoAtomName(
                                2, EACH.aminoacid, tmpAtomName )
                
                ## Convert IUPAC pseudo atoms to BMRB, e.g. LEU MD1 -> HD1
                ## This does NOT convert anything else than methyl and
                ## amino groups. I.e. it doesn't convert methylene and some
                ## other groups. Those need multiple rows in BMRB
                ## chemical shift table.
                tmpAtomName = Nomenclature.Convert_AtomName_IUPAC_2_BMRB_ChemShift(
                                EACH.aminoacid, tmpAtomName )
                
                # Take care of special characters in STAR
                if ( re.compile("#").match(tmpAtomName, 1) ):
                    tmpAtomName = "'%s'" % tmpAtomName
                elif ( re.compile("\"").match(tmpAtomName, 1) ):
                    tmpAtomName = "'%s'" % tmpAtomName
                if ( re.compile("'").match(tmpAtomName, 1) ):
                    tmpAtomName = "\"%s\"" % tmpAtomName
                                
                ## Note that taking the atom nucleus from the first char is
                ## not garanteed to always work!
                nucleus_type = string.upper( EACH.atomname[0][0] )
                
                outString = "    %6s %4s %3s  %-4s %-2s %7s %6s  %s\n" %\
                            (   atomShiftAssignCounter,
##                                outSEGID,
                                EACH.residuenumber,
                                outAA,
                                tmpAtomName,
                                nucleus_type,
                                EACH.shift,
                                outError,
                                ambiguity )
                bioHandle.write(outString)
                
                atomShiftAssignCounter = atomShiftAssignCounter + 1

        lessComments = """
stop_
"""

        ## Just leave it out for now if it's not used any way.
        moreComments = """
# The following loop is used to define sets of Atom-shift assignment IDs that
# represent related ambiguous assignments taken from the above list of
# assigned chemical shifts.  Each element in the set should be separated by a
# comma, as shown in the example below, and is the assignment ID for a chemical
# shift assignment that has been given as ambiguity code of 4 or 5.  Each set
# indicates that the observed chemical shifts are related to the defined 
# atoms, but have not been assigned uniquely to a specific atom in the set.

loop_
  _Atom_shift_assign_ID_ambiguity   

#
#    Sets of Atom-shift Assignment Ambiguities
               #              
#    ------------------------------------------
# Example:    5,4,7
#
                @
stop_
"""
        bioHandle.write(lessComments)

        # write out a tail for STAR format
        bioHandle.write("\nsave_\n\n")
        bioHandle.close()


    


    def WriteChem(self, fileName):
        """
        writes a .chem list which can be read by cns
        Does not use the chemical shifts with values 999.000
        if there are more than one atomnames in the tuple, only the
        first will be used!
        """
        print 'writing a .chem file', fileName
        chemhandle = TextFile.TextFile(fileName, 'w')
        
        chemhandle.write('! derived from the file:\n')
        chemhandle.write('! ' + self.fileName + '\n')
        for EACH in self.atomlist:
            #those with 999.000 don't have an assignment:
            if EACH.shift and EACH.shift != '999.000':
##                 chemhandle.write('do ( store1 = ' + EACH.shift +\
##                                  ' ) ( resid ' + EACH.residuenumber +\
##                                  ' and name ' + EACH.atomname[0] + ' )\n')


## SHALL WE USE STORE 5 and 6 on top of store1 for the errors???
                if  EACH.shifterror:
                    outShiftError = EACH.shifterror
                else:
                    outShiftError = '0.0'

                midshift = string.atof(EACH.shift)
                lowshift = string.atof(EACH.shift) - string.atof(outShiftError)
                upshift = string.atof(EACH.shift) + string.atof(outShiftError)
                chemhandle.write('do ( store1 = ' + str(midshift) +\
                                 ' ) ( resid ' + EACH.residuenumber +\
                                 ' and name ' + EACH.atomname[0] + ' )\n')
                chemhandle.write('do ( store5 = ' + str(lowshift) +\
                                 ' ) ( resid ' + EACH.residuenumber +\
                                 ' and name ' + EACH.atomname[0] + ' )\n')
                chemhandle.write('do ( store6 = ' + str(upshift) +\
                                 ' ) ( resid ' + EACH.residuenumber +\
                                 ' and name ' + EACH.atomname[0] + ' )\n')
    
        chemhandle.write('\n')
        chemhandle.close()

    def WriteSparky(self, fileName):
        """
        writes a .list list which can be read by cns
        Does not use the chemical shifts with values 999.000
        if there are more than one atomnames in the tuple, only the
        first will be used!
        """
        print 'writing a .list file', fileName
        chemhandle = TextFile.TextFile(fileName, 'w')
        
        chemhandle.write('%s   %s  %s    %s   %s  %s' %('Group','Atom','Nuc',\
                                                        'Shift','Sdev',\
                                                        'Assignments'))
        chemhandle.write('\n')
        chemhandle.write('\n')
        
        for EACH in self.atomlist:
            #those with 999.000 don't have an assignment:
            if EACH.shift and EACH.shift != '999.000':
                   a=AminoAcid.AminoAcid(EACH.aminoacid)[0]+str(EACH.residuenumber)
                   az = len(a)
                   ak = 5 - az
                   
                   b=EACH.atomname[0]
                   bz=len(b)
                   bk = 7 - bz
                   
                   c=EACH.atomtype
                   cz=len(c)
                   ck = 5 - cz
                   
                   d=str(EACH.shift)
                   dz=len(d)
                   dk = 9 - dz
                   
                   e=str(EACH.shifterror)
                   ez=len(e)
                   ek = 7 - ez
                   
                   f='1'
                   fz=len(f)
                   fk = 7 - fz

                   first = (ak) * ' '
                   second= (bk)  * ' '
                   third = (ck) * ' '
                   fourth = (dk) * ' '
                   fifth = (ek) * ' '
                   sixth = (fk) * ' '
                   
                   first=first+a
                   second=second+b
                   third=third+c
                   fourth=fourth+d
                   fifth=fifth+e
                   sixth=sixth+f
                   
                   chemhandle.write(first+second+third+fourth+fifth+sixth+'\n')


        chemhandle.write('\n')
        chemhandle.close()
        
    def WriteId(self, fileName):
        """
        writes an ARIA .id file which can be read by cns
        if there are more than one atomname in the tuple, only the
        first will be used!
        if there aren't any atomnumbers, they will be added (counting from 1)
        """
        print 'writing an .id file', fileName
        idhandle = TextFile.TextFile(fileName, 'w')
        idhandle.write('! derived from the file:\n')
        idhandle.write('! ' + self.fileName + '\n')
        atomCounter = 1
        for EACH in self.atomlist:
            if EACH.xeasyatomnumber == None:
                outAtomNumber = str(atomCounter)
            else:
                outAtomNumber = EACH.xeasyatomnumber
            idhandle.write('do ( store2 = ' + outAtomNumber +\
                           ' ) ( resid ' + EACH.residuenumber +\
                           ' and name ' + EACH.atomname[0] + ' )\n')
            atomCounter = atomCounter + 1

    
    def WritePpm(self, fileName):
        """
        writes a list in ARIA's .ppm format
        CNS will read the restraints in store1
        does not use the chemical shifts with values 999.000
        if there are more than one atomnames in the tuple, only the
        first will be used!
        """
        print 'writing a .ppm file', fileName
        ppmhandle = TextFile.TextFile(fileName, 'w')
        ppmhandle.write('! derived from the file:\n')
        ppmhandle.write('! ' + self.fileName + '\n')
        for EACH in self.atomlist:
            #those with 999.000 don't have an assignment:
            if EACH.shift != '999.000':
                if EACH.residuenumber:
                    outResidueNumber = EACH.residuenumber
                else:
                    outResidueNumber = '-'
                if EACH.aminoacid:
                    outAminoAcid = EACH.aminoacid
                else:
                    outAminoAcid = '-'
                if EACH.segid:
                    outSegid = EACH.segid
                else:
                    outSegid = '-'
                if EACH.atomname:
                    outAtomname = EACH.atomname[0]
                else:
                    outAtomname = '-'
                if EACH.shift:
                    outShift = EACH.shift
                else:
                    outShift = '-'
                if EACH.shifterror:
                    outShiftError = EACH.shifterror
                else:
                    outShiftError = '-'
##                print outResidueNumber + ' ' +\
##                      outAminoAcid + ' ' +\
##                      outSegid + ' ' +\
##                      outAtomname + ' ' +\
##                      outShift + ' ' +\
##                      outShiftError
                ppmhandle.write(outResidueNumber + ' ' +\
                                outAminoAcid + ' ' +\
                                outSegid + ' ' +\
                                outAtomname + ' ' +\
                                outShift + ' ' +\
                                outShiftError + '\n')
        


    def WriteXeasyProt(self, fileName, useall=1):
        """
        writes a XEASY .prot file

        useall=1:    writes out all the atoms (default)
        otherwise:   does not write atoms with shifts 999.000

        if there are more than one atomnames in the tuple, only the
        first will be used!
        """
        print 'writing a XEASY .prot file', fileName
        ppmhandle = TextFile.TextFile(fileName, 'w')
        for EACH in self.atomlist:
            if EACH.xeasyatomnumber == None:
                EACH.xeasyatomnumber = 0
            if EACH.shift == None:
                outshift = '0'
            else:
                outshift = EACH.shift
            if EACH.shifterror == None:
                outshifterror = '0.0'
            else:
                outshifterror = EACH.shifterror
            if EACH.residuenumber == None:
                EACH.residuenumber = '0'
            if (EACH.shift != '999.000') or (useall == 1):
                ppmhandle.write('%4i %7.3f %5.3f %-5s %3i \n' %\
                                (int(EACH.xeasyatomnumber),
                                 float(outshift),
                                 float(outshifterror),
                                 EACH.atomname[0],
                                 int(EACH.residuenumber)))
                                #(string.atoi(EACH.xeasyatomnumber),\
                                # string.atof(outshift),\
                                # string.atof(outshifterror),\
                                # EACH.atomname[0],\
                                # string.atoi(EACH.residuenumber)))

    

    def WriteXML2String(self, xmlVersion="1.0", encoding="UTF-8",\
                        dtdFileName="ppm1.1.dtd", dtdVersionName="1.1"):
        """
        returns a string containing the ppm data in XML format

        xmlVersion       appears in the XML declaration
        encoding         appears in the XML declaration
        dtdFileName      points to the corresponding DTD
        dtdVersionName   version name/number of the DTD
        """
        #beginning and end:
        beginSection = """<?xml version="%s" encoding="%s"?>
<!DOCTYPE aria_ppm:list SYSTEM "%s">
<aria_ppm:list>
  <aria_ppm:version>%s</aria_ppm:version>
""" % (xmlVersion, encoding, dtdFileName, dtdVersionName)

        endSection = """</aria_ppm:list>
"""
        #loop over all the atoms:
        outString = ''
        for eachObj in self.atomlist:
            beginAssignment = """  <aria_ppm:assignment>
"""
            endAssignment = """  </aria_ppm:assignment>
"""
            beginList = """    <aria_ppm:atom_list>
"""
            endList = """    </aria_ppm:atom_list>
"""
            outString = outString + beginAssignment + beginList
            for eachAtom in eachObj.atomname:
#                print eachAtom #test
                if eachObj.segid:
                    segmentName = eachObj.segid
                else:
                    segmentName = '' #default
                atomString = """      <aria_ppm:atom>
        <aria_ppm:segment_name>%s</aria_ppm:segment_name>
        <aria_ppm:residue_number>%s</aria_ppm:residue_number>
        <aria_ppm:atom_name>%s</aria_ppm:atom_name>
      </aria_ppm:atom>
""" % (segmentName, eachObj.residuenumber, eachAtom)
                outString = outString + atomString
            ppmString = """    <aria_ppm:ppm_list>
      <aria_ppm:ppm>%s</aria_ppm:ppm>
      <aria_ppm:ppm_error>%s</aria_ppm:ppm_error>
    </aria_ppm:ppm_list>
""" % (eachObj.shift, eachObj.shifterror)
            outString = outString + endList + ppmString + endAssignment
        outString = beginSection + outString + endSection
        return outString

    
    def WriteXML2File(self, fileName):
        """writes the ppm data to an XML file"""
        try:
            outhandle = TextFile.TextFile(fileName, 'w')
        except IOError:
            print 'could not open the file', fileName
            print 'Abort WriteXML2File method.'
            return
        print 'writing to the file:', fileName
        outString = self.WriteXML2String()
        outhandle.write(outString)
        outhandle.close()
        
###############################################################################
class Atom:
    """
    represents the chemical shift assignment of one atom

    ambiguities lead to tuples of atomnames for a single chemical shift
    
    This object contains all the data as attributes, there are no
    methods available.

    all the attributes contain strings, except atomname which contains a tupel
    of strings!
    
    attributes:
        general attributes:
        residuenumber        residuenumber
        aminoacid            3-letter code
        segid                segid
        atomname             tuple of atomnames (ARIA/CNS nomenclature)
        shift                chemical shift
        shifterror           chemical shift error
        comment
        
        other attributes:
        xeasyatomname        XEASY atomname
        xeasyatomnumber      XEASY atomnumber
    """
    def __init__(self, residuenumber = None, aminoacid = None, segid = None,\
                 atomname = None, shift = None, shifterror = None,\
                 comment = None, xeasyatomname = None, xeasyatomnumber = None):
        self.residuenumber = residuenumber
        self.aminoacid = aminoacid
        self.segid = segid
        self.atomname = atomname
        self.shift = shift
        self.shifterror = shifterror
        self.comment = comment

        self.xeasyatomname = xeasyatomname
        self.xeasyatomnumber = xeasyatomnumber
        
    
###############################################################################
def _DoesFileExist(fileName):
    if os.path.exists(fileName) == 0:
        print 'WARNING:', fileName, 'does not exist.'
        return 0
    return 1

###############################################################################
#test code:
if __name__ == "__main__":
    print 'testing module:'
    AT1 = Atom('123', 'LEU', 'ABCD', ('HA', 'HB1'), '3.321', '0.002', 'test', 'HA', '321')
    AT2 = Atom('124', 'ILE', 'ABCD', ('HA', 'HB2'), '3.123', '0.003', 'test', 'HA', '654')
    PL=PpmList()
    PL.atomlist = [AT1, AT2]
    PL.atomdicfa = {}
    PL.dimension = '2'
    PL.name = 'test'
    PL.fileName = 'test'
#    print '1. The test ppm list is:'
#    print PL.atomlist
    print '\nThe Stdout() method prints:'
    PL.Stdout()
    print '\nThe XML output is:'
    print PL.WriteXML2String()
    print '\ntest done. bye.'
