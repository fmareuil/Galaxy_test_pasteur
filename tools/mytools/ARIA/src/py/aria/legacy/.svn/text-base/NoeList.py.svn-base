"""
A module for representing NOE data in ARIA's general NOE data format.
Contains methods to read and write NOE data in various formats
(e.g. XEASY, ANSIG, PIPP)
All the NOE data handling is done with this module.
"""
__author__   = "$Author: bardiaux $"
__revision__ = "$Revision: 1.1.1.1 $"
__date__     = "$Date: 2010/03/23 15:27:24 $"

import copy, os, re, string

import AminoAcid, PseudoAtom
import SequenceList, Comparisons, HeteronucleusPlusProton
import DeleteCnsComments
import FortranFormat, TextFile

###############################################################################
class NoeList:
    """
    Contains all NOE data of one spectrum in ARIA's general NOE data format
    supplies methods for input and output in various formats
    Author: Jens Linge, EMBL Heidelberg, 1998-1999

    public attributes:
        name           name of the spectrum
        comment        comment for the NOE list
        dimension      dimensionality of the spectrum (2, 3 or 4)
        fileNames      list of fileNames from which the data come from
        peaksdic       dictionary (key: (spectrumName, peakNumber), value: peakobject)
        peakslist      list of NOE peaks (with their contributions)
        intraList      list of NOE peaks (with their contributions)
        seqList        list of NOE peaks (with their contributions)
        mediumList     list of NOE peaks (with their contributions)
        longList       list of NOE peaks (with their contributions)
        unassignedList list of NOE peaks (with their contributions)
        ambiguousList  list of NOE peaks (with their contributions)
    public methods:
        AddPeak               adds one peak (peakobject)
        AddSequence           adds sequence information
        ConvertWildcards      converts wildcards to IUPAC
        CountNontrivial       count the number of assignments of each peak
        IntraSeqMediumLong    splits the peak list into intra, medium and long
        ReadAnsigCpe          reads a crosspeaks export file
        ReadAnsigNoe          reads a general NOE restraints file
        ReadAssignPeaksProt   insert a spectrum, reading XEASY .assign, .peaks,
                              and .prot files
        ReadAureliaUserPeak   reads an AURELIA user peak list
        ReadList              insert a spectrum, reading ARIA's .list file
        ReadListPeaksProt     insert a spectrum, reading XEASY .peaks, .prot,
                              and ARIA's .list files
        ReadLolUpl            reads DYANA .lol and .upl files
        ReadOldListPeaksProt  insert a spectrum, reading XEASY .peaks, .prot,
                              and ARIA's old>
        ReadNMRViewPeaks      insert a spectrum, reading a NMRView extracted
                              peak file .list files
        ReadPeaksProt         insert a spectrum, reading XEASY .peaks and .prot
        ReadPippPck           insert a spectrum, reading PIPP's .PCK file
        ReadRegine            insert a spectrum, reading a Regine file
        ReadSparky            insert a spectrum, reding a Sparky's .list file
        ReadTbl               insert a spectrum, reading ARIA's .tbl files
        RemovePeak            removes one peak (specified by its peakNumber)
        RemoveDoubleQuotes    removes double quotation marks from the atomnames
        RestraintCombination  transforms long range restraints in ambiguous ones.
        SetSpectrumName       sets the spectrum name for all contributions
        Sort                  sorting after spectrumName and peakNumber
        Stdout                writes the data in the new .list format to stdout
        WriteAnsigNoe         writes an ANSIG general NOE restraints file
        WriteAureliaUserPeak  writes an AURELIA user peak list
        WriteList             writes a list in ARIA's .list format
        WriteNMRViewPeaks     writes a list in NMRView .xpk format
        WriteOldList          writes a list in ARIA's old .list format
        WriteSparky           writes a list in Sparky's .list format
        WriteTbl              writes .tbl in ARIA format
        WriteUplLol           writes .upl and .lol files for MolMol
        WriteXeasyAssignPeaks writes .peaks and .assign files
        WriteXML2String       returns a string containing the sequence in XML format
        WriteXML2File         writes the sequence to an XML file

    NOTE:
        -all fields are treated internally as strings, except:
          atomnames are always tuples of strings
          dimension is always an integer
        -empty fields are set to None, showns as '-' in stdout
        -there are internal objects for XEASY .assign, .peaks and .prot
    """
    def __init__(self, name = '', dimension = 0, comment = 'ARIA peak list'):
        self.fileNames = []
        self.name = name
        self.comment = comment
        self.dimension = dimension
        self.peakslist = []
        self.intraList = []
        self.seqList = []
        self.mediumList = []
        self.longList = []
        self.unassignedList = []
        self.ambiguousList = []
        self.peaksdic = {}
        #internal attributes for XEASY .peaks files:
        self._iname1 = '?'
        self._iname2 = '?'
        self._iname3 = '?'
        self._iname4 = '?'
        

    def __repr__(self):
        return self.name

    def _cleanUp(self):
        #delete all the existing peaks:
        self.peakslist = []
        self.intraList = []
        self.seqList = []
        self.mediumList = []
        self.longList = []
        self.unassignedList = []
        self.ambiguousList = []
    
    def AddPeak(self, peakobject):
        self.peakslist.append(peakobject)
        self.peaksdic[(peakobject.spectrumName, peakobject.peakNumber)] = peakobject

    def AddSequence(self, sequenceDictionary):
        """
        sequenceDictionary = {'16': 'ARG', '17': 'GLY'}
        is a dictionary containing the residue numbers as keys and the
        three letter code amino acids as values (both should be strings)
        """
        for eachP in self.peakslist:
            for eachC in eachP.contributions:
                if sequenceDictionary.has_key(eachC.residue1):
                    eachC.aa1 = sequenceDictionary[eachC.residue1]
                if sequenceDictionary.has_key(eachC.residue2):
                    eachC.aa2 = sequenceDictionary[eachC.residue2]
                    
    def AddSequenceFromFile(self, sequenceFilename, startNumber = 1):
        """
        read a sequence file and adds the information to the NOE peak list
        startNumber (integer) for the residue numbering
        sequence must be in three-letter code format
        """
        if _DoesFileExist(sequenceFilename) == 0:
            return
        SL=SequenceList.SequenceList()
        SL.ReadSeq(sequenceFilename)
        sequenceDictionary = {}
        for aa in SL.aalist:
            sequenceDictionary[str(startNumber)] = aa
            startNumber = startNumber + 1
        self.AddSequence(sequenceDictionary)
        
    def ConvertWildcards(self):
        for eachP in self.peakslist:
            for eachC in eachP.contributions:
                for eachA in (('atomnameh1', 'aa1'), ('atomnamep1', 'aa1'),\
                              ('atomnameh2', 'aa2'), ('atomnamep2', 'aa2')):
                    newList = []
                    for atom in getattr(eachC, eachA[0]):
#                        print getattr(eachC, eachA[1]), atom #test
                        if getattr(eachC, eachA[0]) != None and atom != None:
                            convertedTuple = PseudoAtom.Pseudo2IupacTuple(getattr(eachC, eachA[1]), atom)
                            for converted in convertedTuple:
                                newList.append(converted)
                            setattr(eachC, eachA[0], tuple(newList))
            
    def CountNontrivial(self):
        if len(self.peakslist) == 0:
            print 'peakslist is empty. CountNontrivial method aborted.'
            return
        print 'counting the contributions for each peak and setting the attribute nta for each contribution.'
        for EACHP in self.peakslist:
            for EACHCON in EACHP.contributions:
                #convert to string:
                EACHCON.nta = str(len(EACHP.contributions))
        

    def IntraSeqMediumLong(self):
        """
        splits self.peakslist into shorter lists:
          self.intraList          |i-j| = 0
          self.seqList            |i-j| = 1
          self.mediumList         1 < |i-j| < 5
          self.longList           |i-j| > 4
          self.unassignedList
          self.ambiguousList

        the self.peakslist attribute remains untouched
        """
        print 'sorting the peaks: intra, seq, medium, long, unassigned, ambiguous'
        #IMPORTANT:the methode has to start with empty lists
        self.intraList=[]
        self.seqList=[]
        self.mediumList=[]
        self.longList=[]
        self.unassignedList=[]
        self.ambiguousList=[]
        
        for eachPeak in self.peakslist:
            if len(eachPeak.contributions) > 1:
                self.ambiguousList.append(eachPeak)
            elif len(eachPeak.contributions) == 0:
                continue
            else:
                if eachPeak.contributions[0].residue1 != None and \
                   eachPeak.contributions[0].residue2 != None:
                    if abs(string.atoi(eachPeak.contributions[0].residue1) -\
                           string.atoi(eachPeak.contributions[0].residue2)) == 0:
                        self.intraList.append(eachPeak)
                    elif abs(string.atoi(eachPeak.contributions[0].residue1) -\
                           string.atoi(eachPeak.contributions[0].residue2)) == 1:
                        self.seqList.append(eachPeak)
                    elif abs(string.atoi(eachPeak.contributions[0].residue1) -\
                             string.atoi(eachPeak.contributions[0].residue2)) < 5:
                        self.mediumList.append(eachPeak)
                    else:
                        self.longList.append(eachPeak)
                else:
                    self.unassignedList.append(eachPeak)

        #some messages for stdout:
        print len(self.intraList), len(self.seqList), len(self.mediumList),\
              len(self.longList), len(self.unassignedList), len(self.ambiguousList)
              


    def ReadAnsigCpe(self, fileName, het1, pro1, het2, pro2):
        """
        reads an ANSIG crosspeaks export file or an ANSIG storage crosspeaks file
        """
        if _DoesFileExist(fileName) == 0:
            return
        print 'reading the ANSIG crosspeaks export file:\n ', fileName
        
        #important - delete all the existing peaks:
        self._cleanUp()

        
        cpHandle = TextFile.TextFile(fileName)
        
        bang = re.compile('#')

        #some default values:
        peakType = None
        segid1 = None
        segid2 = None
        volumeError = None
        distance = None
        distanceError = None
        lowerBound = None
        upperBound = None
        dh1ppm = None
        dp1ppm = None
        dh2ppm = None
        dp2ppm = None

        #get total number of crosspeaks and dimensionality from the first two lines:
        eachLine = cpHandle.readline()
        eachLine = cpHandle.readline()
        
        totAndDim = string.split(eachLine)
        totalNumber = totAndDim[-2]
        self.dimension = int(totAndDim[-1])

        #there are two possible cases:
        #1. ANSIG v3.3 export crosspeaks file
        #   second line contains 'totalNumber dimensionality'
        #   spectrumName appears in every line
        #2. ANSIG v3.3 storage crosspeaks file
        #   second line contains 'spectrumName totalNumber dimensionality'
        #   spectrumName appears only in the header

        if len(totAndDim) == 2:  #for the ANSIG crosspeaks files           
            format2 = FortranFormat.FortranFormat('3E13.6,A12,7I6,6A4')
            format3 = FortranFormat.FortranFormat('4E13.6,A12,9I6,9A4')
            format4 = FortranFormat.FortranFormat('5E13.6,A12,11I6,12A4')

            #use (lineNumber - 2) as peak number:
            lineNumber = 2

            for eachLine in cpHandle:
                lineNumber = lineNumber + 1
                if bang.match(eachLine) or len(eachLine) < 40:
                    continue

                #some more default values:
                atomnameh1 = (None,)
                atomnameh2 = (None,)
                atomnamep1 = (None,)
                atomnamep2 = (None,)
                h1ppm = None
                h2ppm = None
                p1ppm = None
                p2ppm = None


                exec 'lineList = FortranFormat.FortranLine(eachLine, format' + str(self.dimension) +')'
                outS = ''
                if het1 != 'N':
                    outS = outS + 'h1ppm = string.strip(str(lineList[' + het1 + '-1]));'
                    outS = outS + 'atomnameh1 = (string.strip(lineList[' + str(self.dimension * 5 + 4) + '+' + het1 + ']),);' 
                if het2 != 'N':
                    outS = outS + 'h2ppm = string.strip(str(lineList[' + het2 + '-1]));'
                    outS = outS + 'atomnameh2 = (string.strip(lineList[' + str(self.dimension * 5 + 4) + '+' + het2 + ']),);'
                if pro1 != 'N':
                    outS = outS + 'p1ppm = string.strip(str(lineList[' + pro1 + '-1]));'
                    outS = outS + 'atomnamep1 = (string.strip(lineList[' + str(self.dimension * 5 + 4) + '+' + pro1 + ']),);' 
                if pro2 != 'N':
                    outS = outS + 'p2ppm = string.strip(str(lineList[' + pro2 + '-1]));'
                    outS = outS + 'atomnamep2 = (string.strip(lineList[' + str(self.dimension * 5 + 4) + '+' + pro2 + ']),);' 
                #ANSIG intensity => ARIA volume:
                outS = outS + 'volume = string.strip(str(lineList[' + str(self.dimension) + ']));' 
                outS = outS + 'spectrumName = string.strip(lineList[' + str(self.dimension + 1) + ']);' 

                if pro1 != 'N':
                    outS = outS + 'residue1 = string.strip(str(lineList[' + str(self.dimension * 3 + 4) + '+' + pro1 + ']));'
                else:
                    outS = outS + 'residue1 = "";'
                    
                if pro2 != 'N':
                    outS = outS + 'residue2 = string.strip(str(lineList[' + str(self.dimension * 3 + 4) + '+' + pro2 + ']));' 
                else:
                    outS = outS + 'residue2 = "";'
                    
                if pro1 != 'N':
                    outS = outS + 'aa1 = string.upper(string.strip(lineList[' +\
                           str(self.dimension * 4 + 4) + '+' + pro1 + ']));'
                else:
                    outS = outS + 'aa1 = "";'
                    
                if pro2 != 'N':
                    outS = outS + 'aa2 = string.upper(string.strip(lineList[' +\
                           str(self.dimension * 4 + 4) + '+' + pro2 + ']));' 
                else:
                    outS = outS + 'aa2 = "";'
                    
                exec outS

                if h1ppm == '': h1ppm = None
                if p1ppm == '': p1ppm = None
                if h2ppm == '': h2ppm = None
                if p2ppm == '': p2ppm = None
                if volume == '': volume = None
                if spectrumName == '': spectrumName = None
                if residue1 == '': residue1 = None
                if residue2 == '': residue2 = None
                if aa1 == '': aa1 = None
                if aa2 == '': aa2 = None
                if atomnameh1 == ('',): atomnameh1 = (None,)
                if atomnameh2 == ('',): atomnameh2 = (None,)
                if atomnamep1 == ('',): atomnamep1 = (None,)
                if atomnamep2 == ('',): atomnamep2 = (None,)

                #for the strange residue numbers, e.g. '72?'
                if residue1:
                    if residue1[-1] == '?':
                        residue1 = residue1[:-1]
                if residue2:
                    if residue2[-1] == '?':
                        residue2 = residue2[:-1]
                
                #put it in a contribution instance:
                CONT = NoeContribution(peakNumber = str(lineNumber - 2),\
                                       peakType = peakType,\
                                       residue1 = residue1,\
                                       aa1 = aa1 ,\
                                       segid1 = segid1,\
                                       atomnameh1 = atomnameh1,\
                                       atomnamep1 = atomnamep1,\
                                       residue2 = residue2,\
                                       aa2 = aa2,\
                                       segid2 = segid2,\
                                       atomnameh2 = atomnameh2,\
                                       atomnamep2 = atomnamep2,\
                                       volume = volume,\
                                       volumeError = volumeError,\
                                       distanceAve = distance,\
                                       distanceStd = distanceError,\
                                       lowerBound = lowerBound,\
                                       upperBound = upperBound,\
                                       h1ppm = h1ppm,\
                                       dh1ppm = dh1ppm,\
                                       p1ppm = p1ppm,\
                                       dp1ppm = dp1ppm,\
                                       h2ppm = h2ppm,\
                                       dh2ppm = dh2ppm,\
                                       p2ppm = p2ppm,\
                                       dp2ppm = dp2ppm)

                #ANSIG specific attributes:
                CONT.spectrumName = spectrumName

                #create a peak instance and add contribution:
                NOE = NoePeak()
                NOE.AddContribution(CONT)

                #add peak to the NOE list:
                self.AddPeak(NOE)
            
        else:  #for the ANSIG storage files
            format2 = FortranFormat.FortranFormat('3E13.6,5I6,6A4')
            format3 = FortranFormat.FortranFormat('4E13.6,7I6,9A4')
            format4 = FortranFormat.FortranFormat('5E13.6,9I6,12A4')
            spectrumName = totAndDim[0]
            
            #use (lineNumber - 2) as peak number:
            lineNumber = 2

            for eachLine in cpHandle:
                lineNumber = lineNumber + 1
                if bang.match(eachLine) or len(eachLine) < 40:
                    continue

                #some more default values:
                atomnameh1 = (None,)
                atomnameh2 = (None,)
                atomnamep1 = (None,)
                atomnamep2 = (None,)
                h1ppm = None
                h2ppm = None
                p1ppm = None
                p2ppm = None

                exec 'lineList = FortranFormat.FortranLine(eachLine, format' + str(self.dimension) +')'
                outS = ''
                if het1 != 'N':
                    outS = outS + 'h1ppm = string.strip(str(lineList[' + het1 + '-1]));'
                    outS = outS + 'atomnameh1 = (string.strip(lineList[' + str(self.dimension * 5 + 1) + '+' + het1 + ']),);' 
                if het2 != 'N':
                    outS = outS + 'h2ppm = string.strip(str(lineList[' + het2 + '-1]));'
                    outS = outS + 'atomnameh2 = (string.strip(lineList[' + str(self.dimension * 5 + 1) + '+' + het1 + ']),);'
                if pro1 != 'N':
                    outS = outS + 'p1ppm = string.strip(str(lineList[' + pro1 + '-1]));'
                    outS = outS + 'atomnamep1 = (string.strip(lineList[' + str(self.dimension * 5 + 1) + '+' + pro1 + ']),);' 
                if pro2 != 'N':
                    outS = outS + 'p2ppm = string.strip(str(lineList[' + pro2 + '-1]));'
                    outS = outS + 'atomnamep2 = (string.strip(lineList[' + str(self.dimension * 5 + 1) + '+' + pro2 + ']),);' 

                #ANSIG intensity => ARIA volume:
                outS = outS + 'volume = string.strip(str(lineList[' + str(self.dimension) + ']));' 
                outS = outS + 'residue1 = string.strip(str(lineList[' + str(self.dimension * 3 + 1) + '+' + pro1 + ']));' 
                outS = outS + 'residue2 = string.strip(str(lineList[' + str(self.dimension * 3 + 1) + '+' + pro2 + ']));' 

                outS = outS + 'aa1 = string.upper(string.strip(lineList[' +\
                     str(self.dimension * 4 + 1) + '+' + pro1 + ']));' 
                outS = outS + 'aa2 = string.upper(string.strip(lineList[' +\
                     str(self.dimension * 4 + 1) + '+' + pro2 + ']));' 

                exec outS

                if h1ppm == '': h1ppm = None
                if p1ppm == '': p1ppm = None
                if h2ppm == '': h2ppm = None
                if p2ppm == '': p2ppm = None
                if volume == '': volume = None
                if residue1 == '': residue1 = None
                if residue2 == '': residue2 = None
                if aa1 == '': aa1 = None
                if aa2 == '': aa2 = None
                if atomnameh1 == ('',): atomnameh1 = (None,)
                if atomnameh2 == ('',): atomnameh2 = (None,)
                if atomnamep1 == ('',): atomnamep1 = (None,)
                if atomnamep2 == ('',): atomnamep2 = (None,)

                #put it in a contribution instance:
                CONT = NoeContribution(peakNumber = str(lineNumber - 2),\
                                       peakType = peakType,\
                                       residue1 = residue1,\
                                       aa1 = aa1 ,\
                                       segid1 = segid1,\
                                       atomnameh1 = atomnameh1,\
                                       atomnamep1 = atomnamep1,\
                                       residue2 = residue2,\
                                       aa2 = aa2,\
                                       segid2 = segid2,\
                                       atomnameh2 = atomnameh2,\
                                       atomnamep2 = atomnamep2,\
                                       volume = volume,\
                                       volumeError = volumeError,\
                                       distanceAve = distance,\
                                       distanceStd = distanceError,\
                                       lowerBound = lowerBound,\
                                       upperBound = upperBound,\
                                       h1ppm = h1ppm,\
                                       dh1ppm = dh1ppm,\
                                       p1ppm = p1ppm,\
                                       dp1ppm = dp1ppm,\
                                       h2ppm = h2ppm,\
                                       dh2ppm = dh2ppm,\
                                       p2ppm = p2ppm,\
                                       dp2ppm = dp2ppm)

                #ANSIG specific attributes:
                CONT.spectrumName = spectrumName

                #create a peak instance and add contribution:
                NOE = NoePeak()
                NOE.AddContribution(CONT)

                #add peak to the NOE list:
                self.AddPeak(NOE)


    def ReadAnsigNoe(self, noefile, het1, pro1, het2, pro2):
        """
        reads an ANSIG general NOE restraint file
        
        I use the same convention as in the .html file for ANSIG:
        het1 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro1 = '1', '2', '3', '4' or 'N' (column number or not used)
        het2 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro2 = '1', '2', '3', '4' or 'N' (column number or not used)

        I use the ANSIG 'absolute intensity' field for the ARIA volumes
        
        """
        #some messages:
        print 'reading the ANSIG general NOE file:'
        print '  ' + noefile

        #check if the files exist:
        _DoesFileExist(noefile)
        self._cleanUp()

        
        #important - delete all the existing peaks:
        self._cleanUp()
        
        #read the general NOE file:
        noeHandle = TextFile.TextFile(noefile)
        bang = re.compile('!')
        format2 = FortranFormat.FortranFormat('I2,3A4,F8.3,3A4,F8.3,I6,A12,F6.1,2E11.4')
        format3 = FortranFormat.FortranFormat('I2,3A4,F8.3,3A4,F8.3,3A4,F8.3,I6,A12,F6.1,2E11.4')
        format4 = FortranFormat.FortranFormat('I2,3A4,F8.3,3A4,F8.3,3A4,F8.3,3A4,F8.3,I6,A12,F6.1,2E11.4')
        
        #unavailable data set to default:
        peakType = None
        segid1 = None
        segid2 = None
        volumeError = None
        distance = None
        distanceError = None
        lowerBound = None
        upperBound = None
        dh1ppm = None
        dp1ppm = None
        dh2ppm = None
        dp2ppm = None

        for eachLine in noeHandle:
            if bang.match(eachLine) or len(eachLine) < 10:
                continue
##            #replace tabs with whitespace (circumvent TR1 problem):
##            print eachLine, #test
##            eachLine = re.sub('\t', ' ', eachLine)
            #get the right format for the corresponding dimensionality:
            exec 'lineList = FortranFormat.FortranLine(eachLine, format' + eachLine[0:1] +')'
##            print lineList, lineList[3] #test
            self.dimension = lineList[0]
            outS = ''
            #residue1:
            outS = outS + 'residue1 = lineList[' + str(4 * string.atoi(pro1) - 3) + '];'
            outS = outS + 'aa1 = lineList[' + str(4 * string.atoi(pro1) - 2) + '];'
            if het1 != 'N':
                outS = outS + 'atomnameh1 = lineList[' + str(4 * string.atoi(het1) - 1) + '];'
                outS = outS + 'h1ppm = lineList[' + str(4 * string.atoi(het1)) + '];'
            else:
                atomnameh1 = None
                h1ppm = None
            if pro1 != 'N':
                outS = outS + 'atomnamep1 = lineList[' + str(4 * string.atoi(pro1) - 1) + '];'
                outS = outS + 'p1ppm = lineList[' + str(4 * string.atoi(pro1)) + '];'
            else:
                atomnamep1 = None
                p1ppm = None
            #residue2:
            outS = outS + 'residue2 = lineList[' + str(4 * string.atoi(pro2) - 3) + '];'
            outS = outS + 'aa2 = lineList[' + str(4 * string.atoi(pro2) - 2) + '];'
            if het2 != 'N':
                outS = outS + 'atomnameh2 = lineList[' + str(4 * string.atoi(het2) - 1) + '];'
                outS = outS + 'h2ppm = lineList[' + str(4 * string.atoi(het2)) + '];'
            else:
                atomnameh2 = None
                h2ppm = None
            if pro2 != 'N':
                outS = outS + 'atomnamep2 = lineList[' + str(4 * string.atoi(pro2) - 1) + '];'
                outS = outS + 'p2ppm = lineList[' + str(4 * string.atoi(pro2)) + '];'
            else:
                atomnamep2 = None
                p2ppm = None
            outS = outS + 'peakNumber = lineList[' + str(self.dimension * 4 + 1) + '];'
            outS = outS + 'spectrumName = lineList[' + str(self.dimension * 4 + 2) + '];'
            outS = outS + 'mixtime = lineList[' + str(self.dimension * 4 + 3) + '];'
            #use absolute intensities for the volumes:
            outS = outS + 'volume = lineList[' + str(self.dimension * 4 + 4) + '];'
            outS = outS + 'relativeint = lineList[' + str(self.dimension * 4 + 5) + '];'

            exec outS
            
            #put it in a contribution instance:
            CONT = NoeContribution(peakNumber = string.strip(str(peakNumber)),\
                                   peakType = peakType,
                                   residue1 = string.strip(str(residue1)),\
                                   aa1 = string.upper(string.strip(aa1)),\
                                   segid1 = segid1,\
                                   atomnameh1 = (string.strip(atomnameh1),),\
                                   atomnamep1 = (string.strip(atomnamep1),),\
                                   residue2 = string.strip(str(residue2)),\
                                   aa2 = string.upper(string.strip(aa2)),\
                                   segid2 = segid2,\
                                   atomnameh2 = (string.strip(atomnameh2),),\
                                   atomnamep2 = (string.strip(atomnamep2),),\
                                   volume = string.strip(str(volume)),\
                                   volumeError = volumeError,\
                                   distanceAve = distance,\
                                   distanceStd = distanceError,\
                                   lowerBound = lowerBound,\
                                   upperBound = upperBound,\
                                   h1ppm = string.strip(str(h1ppm)),\
                                   dh1ppm = dh1ppm,\
                                   p1ppm = string.strip(str(p1ppm)),\
                                   dp1ppm = dp1ppm,\
                                   h2ppm = string.strip(str(h2ppm)),\
                                   dh2ppm = dh2ppm,\
                                   p2ppm = string.strip(str(p2ppm)),\
                                   dp2ppm = dp2ppm)
            
            #ANSIG specific attributes:
            CONT.spectrumName = string.strip(str(spectrumName))
            CONT.mixtime = string.strip(str(mixtime))
            CONT.relativeint = string.strip(str(relativeint))
            
            #create a peak instance and add contribution:
            NOE = NoePeak()
            NOE.AddContribution(CONT)

            #add peak to the NOE list:
            self.AddPeak(NOE)
    
    def ReadAssignPeaksProt(self, assignfile, peaksfile, protfile,\
                            het1, pro1, het2, pro2):
        """
        reads Xeasy .assign, .peaks and .prot files
        
        I use the same convention as in the .html file for XEASY:
        het1 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro1 = '1', '2', '3', '4' or 'N' (column number or not used)
        het2 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro2 = '1', '2', '3', '4' or 'N' (column number or not used)

        I use pro1 and pro2 for the matrix -> binary format conversion

        the strategy for combining the information is:
        1. take all the peaks in the .peaks files
        2. all the peaks, for which there is an assignment
           in the .assign file and no complete assignment in the
           .peaks file, will get the possible assignments from the
           .assign file
           that means that fully assigned peaks from .peaks
           stay untouched
        """
        #some messages:
        print 'reading the XEASY files:'
        print '  ' + assignfile
        print '  ' + peaksfile
        print '  ' + protfile

        #for XEASY I need:
        #from Aria.DataformatConversion import ReadXeasy
        import ReadXeasy
        
        #internal XEASY data stores:
        self._ASSIGN = ReadXeasy.XeasyAssign() #internal store of XEASY .assign data
        self._PEAKS = ReadXeasy.XeasyPeaks()   #internal store of XEASY .peaks data
        self._PROT = ReadXeasy.XeasyProt()     #internal store of XEASY .prot data
        
        #check if the files exist:
        _DoesFileExist(assignfile)
        _DoesFileExist(peaksfile)
        _DoesFileExist(protfile)
        
        #for the documentation:
        self.fileNames.append(assignfile)
        self.fileNames.append(peaksfile)
        self._cleanUp()

        #important - delete all the existing peaks:
        self._cleanUp()
        
        #read in .assign, use pro1 and pro2 for the matrix:
        #.first = pro1 and .second = pro2 !!!
        self._ASSIGN.ReadAssign(assignfile, pro1, pro2)
        
        #read in .peaks and .prot:
        self._PEAKS.ReadPeaks(peaksfile)
        self._PROT.ReadProt(protfile)
        
        #get the dimensionality from the .peaks file:
        self.dimension = self._PEAKS.dimension
        
        #get INAME entries from .peaks:
        self._iname1 = self._PEAKS.iname1
        self._iname2 = self._PEAKS.iname2
        self._iname3 = self._PEAKS.iname3
        self._iname4 = self._PEAKS.iname4
        
        #put all the information together (like described above):
        #start with the peaks from the .peaks file:
        for EACHP in self._PEAKS.peakslist:
            EACHP.spectrumName = None #for the spectrumName
            NEWP = NoeContribution()
            NEWP.peakNumber = EACHP.peakNumber

            #copy the user-defined spectrum type:
            getudst = copy.copy(EACHP.spectrumType)
            getinte = copy.copy(EACHP.integration)

##            print getudst, getinte #test
            
            #get the right columns for the atomnumbers (0 => None):
            if het1 != 'N':
                NEWP.atomnumberh1 = getattr(EACHP, 'atomnumber' + het1)
                if NEWP.atomnumberh1 == '0':
                    NEWP.atomnumberh1 = None
            if pro1 != 'N':
                NEWP.atomnumberp1 = getattr(EACHP, 'atomnumber' + pro1)
                if NEWP.atomnumberp1 == '0':
                    NEWP.atomnumberp1 = None
            if het2 != 'N':
                NEWP.atomnumberh2 = getattr(EACHP, 'atomnumber' + het2)
                if NEWP.atomnumberh2 == '0':
                    NEWP.atomnumberh2 = None
            if pro2 != 'N':
                NEWP.atomnumberp2 = getattr(EACHP, 'atomnumber' + pro2)
                if NEWP.atomnumberp2 == '0':
                    NEWP.atomnumberp2 = None
            
            NEWP.peakType = EACHP.peakType
            NEWP.volume = EACHP.volume
            NEWP.volumeError = EACHP.volumeError
            NEWP.xeasyudst = EACHP.spectrumType
            NEWP.xeasyinte = EACHP.integration
            
            #get the right columns for the frequencies:
            if het1 != 'N':
                NEWP.h1ppm = getattr(EACHP, 'w' + het1)
            if pro1 != 'N':
                NEWP.p1ppm = getattr(EACHP, 'w' + pro1)
            if het2 != 'N':
                NEWP.h2ppm = getattr(EACHP, 'w' + het2)
            if pro2 != 'N':
                NEWP.p2ppm = getattr(EACHP, 'w' + pro2)
            
            #now copy this contribution to the peak:
            WHOLEP = NoePeak()
            WHOLEP.AddContribution(NEWP)
            
            #now copy the whole peak to the self.peakslist attribute:
            self.AddPeak(WHOLEP)

            
        #now have a look at the .assign file:
        print 'ARIA uses all the peaks in the .peaks files.'
        print '  all the peaks, for which there is an assignment',\
              'in the .assign file'
        print '  and no complete assignment in the .peaks file, will get'
        print '  the possible assignments from the .assign file:'

        OLDPEAK = None
        for EACHP in self._ASSIGN.assignments:
            EACHP.spectrumName = None #trivial
            #if peak is unambiguously assigned in .peaks, keep this peak:
            if self.peaksdic.has_key((EACHP.spectrumName, EACHP.peakNumber)):
                #check if all the assignments are there:
                if len(self.peaksdic[(EACHP.spectrumName, EACHP.peakNumber)].contributions) == 1:
                    OLDPEAK = None
                    fullyAssigned = 1
                    if (het1 != 'N') and \
                       (self.peaksdic[(EACHP.spectrumName, EACHP.peakNumber)].contributions[0].atomnameh1 == (None, )):
                        fullyAssigned = 0 
                    if (pro1 != 'N') and \
                       (self.peaksdic[(EACHP.spectrumName, EACHP.peakNumber)].contributions[0].atomnamep1 == (None, )):
                        fullyAssigned = 0 
                    if (het2 != 'N') and \
                       (self.peaksdic[(EACHP.spectrumName, EACHP.peakNumber)].contributions[0].atomnameh2 == (None, )):
                        fullyAssigned = 0 
                    if (pro2 != 'N') and \
                       (self.peaksdic[(EACHP.spectrumName, EACHP.peakNumber)].contributions[0].atomnamep2 == (None, )):
                        fullyAssigned = 0

                    if fullyAssigned:
                        continue
                    
                #if ambiguous peak is in .peaks, it will be removed now:
                print '    replacing peak', EACHP.peakNumber, 'from .peaks with',\
                      'the corresponding peak from .assign'
                OLDPEAK = self.peaksdic[(EACHP.spectrumName, EACHP.peakNumber)]
                self.RemovePeak(EACHP.spectrumName, EACHP.peakNumber)
##             print EACHP.peakNumber, EACHP.dimension, EACHP.binlist, EACHP.first,\
##                   EACHP.second, EACHP.possible1, EACHP.possible2, EACHP.possible3,\
##                   EACHP.possible4   #test
            
##            print OLDPEAK.peakNumber #test
                
            #I have to know which column belongs to h1, p1, h2, p2
            #.first = pro1 and .second = pro2 !!! see above

            #initialize the lists:
            h1list = [None]
            p1list = [None]
            h2list = [None]
            p2list = [None]

            #find the right lists for h1, p1, h2, p2:
            if het1 != 'N':
                h1list = getattr(EACHP, 'possible' + het1)
            if pro1 != 'N':
                p1list = getattr(EACHP, 'possible' + pro1)
            if het2 != 'N':
                h2list = getattr(EACHP, 'possible' + het2)
            if pro2 != 'N':
                p2list = getattr(EACHP, 'possible' + pro2)

            #create a NoePeak instance and get the peak number:
            NOEP = NoePeak()
            NOEP.peakNumber = EACHP.peakNumber
                
            #set the attributes of NoeContribution,
            #loop through all possible assignments:

            #add some of the peak attributes from the removed peak of .peaks:
            def _addAtt(oldPeak, noeC):
                noeC.volume = oldPeak.contributions[0].volume
                noeC.volumeError = oldPeak.contributions[0].volumeError
                noeC.h1ppm = oldPeak.contributions[0].h1ppm
                noeC.dh1ppm = oldPeak.contributions[0].dh1ppm
                noeC.p1ppm = oldPeak.contributions[0].p1ppm
                noeC.dp1ppm = oldPeak.contributions[0].dp1ppm
                noeC.h2ppm = oldPeak.contributions[0].h2ppm
                noeC.dh2ppm = oldPeak.contributions[0].dh2ppm
                noeC.p2ppm = oldPeak.contributions[0].p2ppm
                noeC.dp2ppm = oldPeak.contributions[0].dp2ppm
                noeC.peakType = oldPeak.contributions[0].peakType


            #consider all the possible assignments, use the matrix, luke!
            #I have to go through each possible assignment:
            
            #2D case:
            if (het1 == 'N') and (pro1 != 'N') and (het2 == 'N') \
               and (pro2 != 'N'):
                p1Counter = -1
                for eachp1 in p1list:
                    p1Counter = p1Counter + 1
                    p2Counter = -1
                    for eachp2 in p2list:
                        p2Counter = p2Counter + 1
                        NOEC = NoeContribution()
                        NOEC.peakNumber = EACHP.peakNumber
                        NOEC.atomnumberp1 = eachp1
                        NOEC.atomnumberp2 = eachp2
                        NOEC.xeasyudst = getudst
                        NOEC.xeasyinte = getinte
                        #add old attributes:
                        if OLDPEAK:
                            _addAtt(OLDPEAK, NOEC)

                        #check if this contribution is switched on in the assignment window:
                        if string.atoi(pro1) < string.atoi(pro2):
                            matrixC = p1Counter
                            matrixR = p2Counter
                        else:
                            matrixC = p2Counter
                            matrixR = p1Counter
                            
                        #if the element is 1, it's used:
                        if EACHP.assimatrix[matrixC][matrixR]:
                            #pool the contributions in NoePeak:
                            NOEP.AddContribution(NOEC)
                            
                        #pool the contributions in NoePeak:
                        NOEP.AddContribution(NOEC)
                            
            #3D case, heteronucleus in het1:
            if (het1 != 'N') and (pro1 != 'N') and (het2 == 'N') \
               and (pro2 != 'N'):
                h1Counter = -1
                for eachh1 in h1list:
                    h1Counter = h1Counter + 1
                    p1Counter = -1
                    for eachp1 in p1list:
                        p1Counter = p1Counter + 1
                        p2Counter = -1
                        for eachp2 in p2list:
                            p2Counter = p2Counter + 1
                            NOEC = NoeContribution()
                            NOEC.peakNumber = EACHP.peakNumber
                            NOEC.atomnumberh1 = eachh1
                            NOEC.atomnumberp1 = eachp1
                            NOEC.atomnumberp2 = eachp2
                            NOEC.xeasyudst = getudst
                            NOEC.xeasyinte = getinte
                            #add old attributes:
                            if OLDPEAK:
                                _addAtt(OLDPEAK, NOEC)
                            #check if this contribution is switched on in the assignment window:
                            if string.atoi(pro1) < string.atoi(pro2):
                                matrixC = p1Counter
                                matrixR = p2Counter
                            else:
                                matrixC = p2Counter
                                matrixR = p1Counter
                            
                            #if the element is 1, it's used:
                            if EACHP.assimatrix[matrixC][matrixR]:
                                #pool the contributions in NoePeak:
                                NOEP.AddContribution(NOEC)
                            
                            
            #3D case, heteronucleus in het2:
            if (het1 == 'N') and (pro1 != 'N') and (het2 != 'N') \
               and (pro2 != 'N'):
                h2Counter = -1
                for eachh2 in h2list:
                    h2Counter = h2Counter + 1
                    p1Counter = -1
                    for eachp1 in p1list:
                        p1Counter = p1Counter + 1
                        p2Counter = -1
                        for eachp2 in p2list:
                            p2Counter = p2Counter + 1
                            NOEC = NoeContribution()
                            NOEC.peakNumber = EACHP.peakNumber
                            NOEC.atomnumberp1 = eachp1
                            NOEC.atomnumberh2 = eachh2
                            NOEC.atomnumberp2 = eachp2
                            NOEC.xeasyudst = getudst
                            NOEC.xeasyinte = getinte
                            #add old attributes:
                            if OLDPEAK:
                                _addAtt(OLDPEAK, NOEC)
                            #check if this contribution is switched on in the assignment window:
                            if string.atoi(pro1) < string.atoi(pro2):
                                matrixC = p1Counter
                                matrixR = p2Counter
                            else:
                                matrixC = p2Counter
                                matrixR = p1Counter
                            
                            #if the element is 1, it's used:
                            if EACHP.assimatrix[matrixC][matrixR]:
                                #pool the contributions in NoePeak:
                                NOEP.AddContribution(NOEC)

                            
            #4D case:
            if (het1 != 'N') and (pro1 != 'N') and (het2 != 'N') \
               and (pro2 != 'N'):
                for eachh1 in h1list:
                    h1Counter = h1Counter + 1
                    h2Counter = -1
                    for eachh2 in h2list:
                        h2Counter = h2Counter + 1
                        p1Counter = -1
                        for eachp1 in p1list:
                            p1Counter = p1Counter + 1
                            p2Counter = -1
                            for eachp2 in p2list:
                                p2Counter = p2Counter + 1
                                NOEC = NoeContribution()
                                NOEC.peakNumber = EACHP.peakNumber
                                NOEC.atomnumberh1 = eachh1
                                NOEC.atomnumberp1 = eachp1
                                NOEC.atomnumberh2 = eachh2
                                NOEC.atomnumberp2 = eachp2
                                NOEC.xeasyudst = getudst
                                NOEC.xeasyinte = getinte
                                #add old attributes:
                                if OLDPEAK:
                                    _addAtt(OLDPEAK, NOEC)
                                #check if this contribution is switched on in the assignment window:
                                if string.atoi(pro1) < string.atoi(pro2):
                                    matrixC = p1Counter
                                    matrixR = p2Counter
                                else:
                                    matrixC = p2Counter
                                    matrixR = p1Counter

                                #if the element is 1, it's used:
                                if EACHP.assimatrix[matrixC][matrixR]:
                                    #pool the contributions in NoePeak:
                                    NOEP.AddContribution(NOEC)


            #add the peak to NoeList:
            self.AddPeak(NOEP)
            
        #now get the assigned shifts, shifterrors, atomnames and residue numbers
        #from .prot, loop through the whole NoeList:
        for EPEAK in self.peakslist:
            delCon = []
            for ECON in EPEAK.contributions:
                if self._PROT.atomdican.has_key(ECON.atomnumberh1):
                    ECON.ah1ppm = self._PROT.atomdican[ECON.atomnumberh1].shift
                    ECON.dah1ppm = self._PROT.atomdican[ECON.atomnumberh1].shifterror
                    ECON.atomnameh1 = self._PROT.atomdican[ECON.atomnumberh1].ariaatomname
##                    ECON.residue1 = self._PROT.atomdican[ECON.atomnumberh1].fragmentnumber
                else:
                    if not ((ECON.atomnumberh1 == None) or (ECON.atomnumberh1 == '0')):
                        print 'WARNING: atomnumber', ECON.atomnumberh1, 'not found in .prot'
                if self._PROT.atomdican.has_key(ECON.atomnumberp1):
                    ECON.ap1ppm = self._PROT.atomdican[ECON.atomnumberp1].shift
                    ECON.dap1ppm = self._PROT.atomdican[ECON.atomnumberp1].shifterror
                    ECON.atomnamep1 = self._PROT.atomdican[ECON.atomnumberp1].ariaatomname
                    ECON.residue1 = self._PROT.atomdican[ECON.atomnumberp1].fragmentnumber
                    #if residuenumber het and pro don't fit together, remove contribution:
                    if het1 != 'N' and ECON.atomnameh1 != (None, ) and \
                       ECON.atomnamep1 != (None, ):
                        if self._PROT.atomdican[ECON.atomnumberh1].fragmentnumber != \
                           self._PROT.atomdican[ECON.atomnumberp1].fragmentnumber and \
                           len(EPEAK.contributions) > 1:
                            delCon.append(ECON)
                else:
                    if not ((ECON.atomnumberp1 == None) or (ECON.atomnumberp1 == '0')):
                        print 'WARNING: atomnumber', ECON.atomnumberp1, 'not found in .prot'
                if self._PROT.atomdican.has_key(ECON.atomnumberh2):
                    ECON.ah2ppm = self._PROT.atomdican[ECON.atomnumberh2].shift
                    ECON.dah2ppm = self._PROT.atomdican[ECON.atomnumberh2].shifterror
                    ECON.atomnameh2 = self._PROT.atomdican[ECON.atomnumberh2].ariaatomname
##                    ECON.residue2 = self._PROT.atomdican[ECON.atomnumberh2].fragmentnumber
                else:
                    if not ((ECON.atomnumberh2 == None) or (ECON.atomnumberh2 == '0')):
                        print 'WARNING: atomnumber', ECON.atomnumberh2, 'not found in .prot'
                if self._PROT.atomdican.has_key(ECON.atomnumberp2):
                    ECON.ap2ppm = self._PROT.atomdican[ECON.atomnumberp2].shift
                    ECON.dap2ppm = self._PROT.atomdican[ECON.atomnumberp2].shifterror
                    ECON.atomnamep2 = self._PROT.atomdican[ECON.atomnumberp2].ariaatomname
                    ECON.residue2 = self._PROT.atomdican[ECON.atomnumberp2].fragmentnumber
                    #if residuenumber het and pro don't fit together, remove contribution:
                    if het2 != 'N' and ECON.atomnameh2 != (None, ) and \
                       ECON.atomnamep2 != (None, ):
                        if self._PROT.atomdican[ECON.atomnumberh2].fragmentnumber != \
                           self._PROT.atomdican[ECON.atomnumberp2].fragmentnumber and \
                           len(EPEAK.contributions) > 1:
                            delCon.append(ECON)
                else:
                    if not ((ECON.atomnumberp2 == None) or (ECON.atomnumberp2 == '0')):
                        print 'WARNING: atomnumber', ECON.atomnumberp2, 'not found in .prot'
            for eachDel in delCon:
                EPEAK.contributions.remove(eachDel)
            del delCon
        self.RemoveDoubleQuotes()

        #finally, throw out the contributions (carefully, conservative),
        #in which the heteronucleus does not fit to the proton:
        for EACHP in self.peakslist:
            delCon = []
            for ECON in EPEAK.contributions:
                if het1 != 'N' and ECON.atomnameh1 != (None, ):
                    if not HeteronucleusPlusProton.checkProton(ECON.atomnameh1[0], ECON.atomnamep1[0]):
                        delCon.append(ECON)
                if het2 != 'N' and ECON.atomnameh2 != (None, ):
                    if not HeteronucleusPlusProton.checkProton(ECON.atomnameh2[0], ECON.atomnamep2[0]):
                        delCon.append(ECON)
            for eachDel in delCon:
                EACHP.contributions.remove(eachDel)
            del delCon


    def ReadList(self, fileName):
        """
        reads an ARIA .list file

        the only allowed comments start a line with a '#'
        read the manual for a description of the format
        """
        self.fileNames.append(fileName)

        #check, if file exists:
        _DoesFileExist(fileName)
        
        print 'Reading NOE list from file: %s' % fileName

        #important: delete all the existing peaks
        self._cleanUp()

        #fileHandle:
        fileHandle = TextFile.TextFile(fileName)

        #for the first peak (see if statment below):
        PEAK=None

        for eachLine in fileHandle:
            lineList = string.split(eachLine)
            #getting rid of the delimiters: '-' => None:
            while '-' in lineList:
                index = lineList.index('-')
                lineList.remove('-')
                lineList.insert(index, None)
            #throwing out nonsense and comments:
            if (len(lineList) < 4):
                continue
            if lineList[0] == 'p':
                #add the last peak to peakslist:
                if PEAK:
                    self.AddPeak(PEAK)
                #new peak instance:
                PEAK = NoePeak()
                #get the data:
                spectrumName = lineList[1]
                peakNumber = lineList[2]
                testSet = lineList[3]
                inWeight = lineList[4]
                volume = lineList[5]
                volumeError = lineList[6]
                intensity = lineList[7]
                intensityError = lineList[8]
                h1ppm = lineList[9]
                dh1ppm = lineList[10]
                p1ppm = lineList[11]
                dp1ppm = lineList[12]
                h2ppm = lineList[13]
                dh2ppm = lineList[14]
                p2ppm = lineList[15]
                dp2ppm = lineList[16]
            elif lineList[0] == 'a':
                #get the data:
                peakType = lineList[1]
                fomAria = lineList[2]
                curWeight = lineList[3]
                lowerBound = lineList[4]
                upperBound = lineList[5]
                sumDistanceAve = lineList[6]
                sumDistanceStd = lineList[7]
                backVolumeAria = lineList[8]
                backVolumeAriaStd = lineList[9]
                allAssi = lineList[10]
                nta = lineList[11]
            elif lineList[0] == 'c':
                #get the data for one contribution:
                fomContribution = lineList[1]
                contribution = lineList[2]
                distanceAve = lineList[3]
                distanceStd = lineList[4]
                backVolume = lineList[5]
                backVolumeStd = lineList[6]
                segid1 = lineList[7]
                residue1 = lineList[8]
                aa1 = lineList[9]
                atomnameh1 = (lineList[10],) #atomname is a tuple
                assih1ppm = lineList[11]
                assidh1ppm = lineList[12]
                atomnamep1 = (lineList[13],) #atomname is a tuple
                assip1ppm = lineList[14]
                assidp1ppm = lineList[15]
                segid2 = lineList[16]
                residue2 = lineList[17]
                aa2 = lineList[18]
                atomnameh2 = (lineList[19],) #atomname is a tuple
                assih2ppm = lineList[20]
                assidh2ppm = lineList[21]
                atomnamep2 = (lineList[22],) #atomname is a tuple
                assip2ppm = lineList[23]
                assidp2ppm  = lineList[24]
                #put it into contribution instance:
                CONT = NoeContribution(spectrumName = spectrumName,\
                                       peakNumber = peakNumber,\
                                       testSet = testSet,\
                                       inWeight = inWeight,\
                                       volume = volume,\
                                       volumeError = volumeError,\
                                       intensity = intensity,\
                                       intensityError = intensityError,\
                                       h1ppm = h1ppm,\
                                       dh1ppm = dh1ppm,\
                                       p1ppm = p1ppm,\
                                       dp1ppm = dp1ppm,\
                                       h2ppm = h2ppm,\
                                       dh2ppm = dh2ppm,\
                                       p2ppm = p2ppm,\
                                       dp2ppm = dp2ppm,\
                                       peakType = peakType,\
                                       fomAria = fomAria,\
                                       curWeight = curWeight,\
                                       lowerBound = lowerBound,\
                                       upperBound = upperBound,\
                                       sumDistanceAve = sumDistanceAve,\
                                       sumDistanceStd = sumDistanceStd,\
                                       backVolumeAria = backVolumeAria,\
                                       backVolumeAriaStd = backVolumeAriaStd,\
                                       allAssi = allAssi,\
                                       nta = nta,\
                                       contribution = contribution,\
                                       fomContribution = fomContribution,\
                                       distanceAve = distanceAve,\
                                       distanceStd = distanceStd,\
                                       backVolume = backVolume,\
                                       backVolumeStd = backVolumeStd,\
                                       segid1 = segid1,\
                                       residue1 = residue1,\
                                       aa1 = aa1,\
                                       atomnameh1 = atomnameh1,\
                                       assih1ppm = assih1ppm,\
                                       assidh1ppm = assidh1ppm,\
                                       atomnamep1 = atomnamep1,\
                                       assip1ppm = assip1ppm,\
                                       assidp1ppm = assidp1ppm,\
                                       segid2 = segid2,\
                                       residue2 = residue2,\
                                       aa2 = aa2,\
                                       atomnameh2 = atomnameh2,\
                                       assih2ppm = assih2ppm,\
                                       assidh2ppm = assidh2ppm,\
                                       atomnamep2 = atomnamep2,\
                                       assip2ppm = assip2ppm,\
                                       assidp2ppm = assidp2ppm)
                
                #add contribution to peak:
                PEAK.AddContribution(CONT)
            else:
                print 'WARNING: incorrect format'

        #add the last peak to peakslist:
        if PEAK:
            self.AddPeak(PEAK)
    
    def ReadListPeaksProt(self, listfile, peaksfile, protfile, het1, pro1, het2, pro2):
        """
        reads Xeasy .peaks and .prot file and an ARIA .list file
        this method is used to read the whole data set after a
        complete ARIA iteration

        het1 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro1 = '1', '2', '3', '4' or 'N' (column number or not used)
        het2 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro2 = '1', '2', '3', '4' or 'N' (column number or not used)
        """
        self.fileNames.append(peaksfile)
        self.fileNames.append(protfile)
        self.fileNames.append(listfile)

        #for XEASY I need:
        #from Aria.DataformatConversion import ReadXeasy
        import ReadXeasy
        
        #internal XEASY data stores:
        self._PEAKS = ReadXeasy.XeasyPeaks()   #internal store of XEASY .peaks data
        self._PROT = ReadXeasy.XeasyProt()     #internal store of XEASY .prot data

        #read .peaks and .prot:
        self._PEAKS.ReadPeaks(peaksfile)
        self._PROT.ReadProt(protfile)
        
        #take dimensionality from .peaks:
        self.dimension = self._PEAKS.dimension

        #XEASY inames for the .peaks header:
        self._iname1 = self._PEAKS.iname1
        self._iname2 = self._PEAKS.iname2
        self._iname3 = self._PEAKS.iname3
        self._iname4 = self._PEAKS.iname4

        #important: delete all the existing peaks
        self._cleanUp()
        
        #check, if files exist:
        _DoesFileExist(peaksfile)
        _DoesFileExist(protfile)
        _DoesFileExist(listfile)

        #insert spectrum from the .list file:
        self.ReadList(listfile)
        
        #update the data with the .peaks and .prot data:
        for PEAK in self.peakslist:
            for CON in PEAK.contributions:
##                 #get the peakType:
##                 CON.peakType = self._PEAKS.peaksdic[CON.peakNumber].peakType

##                 #if the spectrumName doesn't fit, don't use it!!!
##                 if CON.spectrumName != inSpectrumName:
##                     continue

                #test whether the peakNumber exists in the .peaks file:
                if not self._PEAKS.peaksdic.has_key(CON.peakNumber):
                    print 'WARNING: peak number', CON.peakNumber, 'does exist in:'
                    print '        ', listfile
                    print '         but NOT in:'
                    print '        ', peaksfile
                    print '         skipping this peak!'
                    continue

                #get the XEASY integration method:
                CON.xeasyinte = self._PEAKS.peaksdic[CON.peakNumber].integration

                #get the XEASY user-defined spectrum-type:
                CON.xeasyudst = self._PEAKS.peaksdic[CON.peakNumber].spectrumType

                #get the missing heteronucleus atomnumbers from the .peaks file:
                if het1 != 'N':
                    CON.atomnumberh1 = getattr(self._PEAKS.peaksdic[CON.peakNumber], 'atomnumber' + het1)
#                    print  CON.atomnumberh1 #test
                    if CON.atomnumberh1 == '0':
                        CON.atomnumberh1 = None
                if het2 != 'N':
                    CON.atomnumberh2 = getattr(self._PEAKS.peaksdic[CON.peakNumber], 'atomnumber' + het2)
                    if CON.atomnumberh2 == '0':
                        CON.atomnumberh2 = None
                
                #get the missing proton atomnumbers, look up the atoms in the
                #.prot file, otherwise try pseudoatoms:
                #pro1:
                if self._PROT.atomdicfa.has_key((CON.residue1, CON.atomnamep1)): #normal case
                    CON.atomnumberp1 = self._PROT.atomdicfa[(CON.residue1, CON.atomnamep1)].atomnumber
                else: #now the pseudoatom case
                    pseudos = PseudoAtom.Atom2Pseudo(CON.atomnamep1[0], CON.aa1)
                    for pseudo in pseudos:
                        if pseudo in self._PROT.atomdicre[CON.residue1]:
                            CON.atomnumberp1 = self._PROT.atomdicfx[(CON.residue1, pseudo)].atomnumber
                            continue #29.1.99 always get the first hit
                #pro2:
                if self._PROT.atomdicfa.has_key((CON.residue2, CON.atomnamep2)): #normal case
                    CON.atomnumberp2 = self._PROT.atomdicfa[(CON.residue2, CON.atomnamep2)].atomnumber
                else: #now the pseudoatom case
                    pseudos = PseudoAtom.Atom2Pseudo(CON.atomnamep2[0], CON.aa2)
                    for pseudo in pseudos:
                        if pseudo in self._PROT.atomdicre[CON.residue2]:
                            CON.atomnumberp2 = self._PROT.atomdicfx[(CON.residue2, pseudo)].atomnumber
                            continue #29.1.99 always get the first hit
        self.RemoveDoubleQuotes()
        
    def ReadOldListPeaksProt(self, listfile, peaksfile, protfile, het1, pro1, het2, pro2):
        """
        insert a spectrum, reading XEASY .peaks, .prot and ARIA's old .list files
        the data from the old .list file are complemented with the data from .peaks
        and .prot
        """
        #some messages:
        print 'reading the files:'
        print '  ' + listfile
        print '  ' + peaksfile
        print '  ' + protfile

        #for XEASY I need:
        #from Aria.DataformatConversion import ReadXeasy
        import ReadXeasy

        #internal XEASY data stores:
        self._ASSIGN = ReadXeasy.XeasyAssign() #internal store of XEASY .assign data
        self._PEAKS = ReadXeasy.XeasyPeaks()   #internal store of XEASY .peaks data
        self._PROT = ReadXeasy.XeasyProt()     #internal store of XEASY .prot data
        
        #check if the files exist:
        _DoesFileExist(listfile)        
        _DoesFileExist(peaksfile)
        _DoesFileExist(protfile)
        
        #for the documentation:
        self.fileNames.append(listfile)
        self.fileNames.append(peaksfile)
        self.fileNames.append(protfile)
        
        #read the old .list file:
        self._ReadOldList(listfile)

        #read .peaks and .prot:
        self._PEAKS.ReadPeaks(peaksfile)
        self._PROT.ReadProt(protfile)

        #get the dimensionality from the .peaks file:
        self.dimension = self._PEAKS.dimension
        
        #update the data with the .peaks and .prot data:
        for PEAK in self.peakslist:
            for CON in PEAK.contributions:
                #get the heteronucleus frequency:
                if het1 != 'N':
                    CON.h1ppm = getattr(self._PEAKS.peaksdic[CON.peakNumber], 'w' +\
                                        het1)
                if het2 != 'N':
                    CON.h2ppm = getattr(self._PEAKS.peaksdic[CON.peakNumber], 'w' +\
                                        het2)
                    
                #get the peakType:
                CON.peakType = self._PEAKS.peaksdic[CON.peakNumber].peakType

                #get the missing heteronucleus atomnumbers from the .peaks file:
                if het1 != 'N':
                    CON.atomnumberh1 = getattr(self._PEAKS.peaksdic[CON.peakNumber], 'atomnumber' + het1)
                    if CON.atomnumberh1 == '0':
                        CON.atomnumberh1 = None
                if het2 != 'N':
                    CON.atomnumberh2 = getattr(self._PEAKS.peaksdic[CON.peakNumber], 'atomnumber' + het2)
                    if CON.atomnumberh2 == '0':
                        CON.atomnumberh2 = None
                
                #get the missing proton atomnumbers, look up the atoms in the
                #.prot file, otherwise try pseudoatoms:
                #pro1:
                if self._PROT.atomdicfa.has_key((CON.residue1, CON.atomnamep1)): #normal case
                    CON.atomnumberp1 = self._PROT.atomdicfa[(CON.residue1, CON.atomnamep1)].atomnumber
                else: #now the pseudoatom case
                    pseudos = PseudoAtom.Atom2Pseudo(CON.atomnamep1[0], CON.aa1)
                    for pseudo in pseudos:
                        if pseudo in self._PROT.atomdicre[CON.residue1]:
                            CON.atomnumberp1 = self._PROT.atomdicfx[(CON.residue1, pseudo)].atomnumber
                            continue #29.1.99 always get the first hit
                #pro2:
                if self._PROT.atomdicfa.has_key((CON.residue2, CON.atomnamep2)): #normal case
                    CON.atomnumberp2 = self._PROT.atomdicfa[(CON.residue2, CON.atomnamep2)].atomnumber
                else: #now the pseudoatom case
                    pseudos = PseudoAtom.Atom2Pseudo(CON.atomnamep2[0], CON.aa2)
                    for pseudo in pseudos:
                        if pseudo in self._PROT.atomdicre[CON.residue2]:
                            CON.atomnumberp2 = self._PROT.atomdicfx[(CON.residue2, pseudo)].atomnumber
                            continue #29.1.99 always get the first hit

                #get the missing heteronucleus atomnames from the .prot file:
                if het1 != 'N' and self._PROT.atomdican.has_key(CON.atomnumberh1):
                    CON.atomnameh1 = (self._PROT.atomdican[CON.atomnumberh1].xeasyatomname,)
                if het2 != 'N' and self._PROT.atomdican.has_key(CON.atomnumberh2):
                    CON.atomnameh2 = (self._PROT.atomdican[CON.atomnumberh2].xeasyatomname,)
        self.RemoveDoubleQuotes()
    
    def ReadPeaksProt(self, peaksfile, protfile, het1, pro1, het2, pro2):
        """
        reads Xeasy .peaks and .prot files
        
        I use the same convention as in the .html file for XEASY:
        het1 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro1 = '1', '2', '3', '4' or 'N' (column number or not used)
        het2 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro2 = '1', '2', '3', '4' or 'N' (column number or not used)
        """
        #some messages:
        print 'reading the XEASY files:'
        print '  ' + peaksfile
        print '  ' + protfile

        #check if the files exist:
        _DoesFileExist(peaksfile)
        _DoesFileExist(protfile)
        
        #for XEASY I need:
        #from Aria.DataformatConversion import ReadXeasy
        import ReadXeasy
        
        #internal XEASY data stores:
        self._ASSIGN = ReadXeasy.XeasyAssign() #internal store of XEASY .assign data
        self._PEAKS = ReadXeasy.XeasyPeaks()   #internal store of XEASY .peaks data
        self._PROT = ReadXeasy.XeasyProt()     #internal store of XEASY .prot data
        
        #for the documentation:
        self.fileNames.append(peaksfile)
        self.fileNames.append(protfile)
        
        #important - delete all the existing peaks:
        self._cleanUp()
        
        #read in .peaks and .prot:
        self._PEAKS.ReadPeaks(peaksfile)
        self._PROT.ReadProt(protfile)

        #get the dimensionality from the .peaks file:
        self.dimension = self._PEAKS.dimension
        
        #put all the information together (like described above):
        #start with the peaks from the .peaks file:
        for EACHP in self._PEAKS.peakslist:
            NEWP = NoeContribution()
            NEWP.peakNumber = EACHP.peakNumber

            #get the right columns for the atomnumbers (0 => None):
            if het1 != 'N':
                NEWP.atomnumberh1 = getattr(EACHP, 'atomnumber' + het1)
                if NEWP.atomnumberh1 == '0':
                    NEWP.atomnumberh1 = None
            if pro1 != 'N':
                NEWP.atomnumberp1 = getattr(EACHP, 'atomnumber' + pro1)
                if NEWP.atomnumberp1 == '0':
                    NEWP.atomnumberp1 = None
            if het2 != 'N':
                NEWP.atomnumberh2 = getattr(EACHP, 'atomnumber' + het2)
                if NEWP.atomnumberh2 == '0':
                    NEWP.atomnumberh2 = None
            if pro2 != 'N':
                NEWP.atomnumberp2 = getattr(EACHP, 'atomnumber' + pro2)
                if NEWP.atomnumberp2 == '0':
                    NEWP.atomnumberp2 = None

            NEWP.peakType = EACHP.peakType
            NEWP.xeasyudst = EACHP.spectrumType
            NEWP.xeasyinte = EACHP.integration

            NEWP.volume = EACHP.volume
            NEWP.volumeError = EACHP.volumeError

            # BARDAIUX 11/04/06
            # Patch for Xeasy spectrum comming from Sparky
            NEWP.intensity = EACHP.volume
            NEWP.intensityError = EACHP.volumeError            
            #
            
            #get the right columns for the frequencies:
            if het1 != 'N':
                NEWP.h1ppm = getattr(EACHP, 'w' + het1)
            if pro1 != 'N':
                NEWP.p1ppm = getattr(EACHP, 'w' + pro1)
            if het2 != 'N':
                NEWP.h2ppm = getattr(EACHP, 'w' + het2)
            if pro2 != 'N':
                NEWP.p2ppm = getattr(EACHP, 'w' + pro2)
            
            #now copy this contribution to the peak:
            WHOLEP = NoePeak()
            WHOLEP.AddContribution(NEWP)
            
            #now copy the whole peak to the self.peakslist attribute:
            self.AddPeak(WHOLEP)
            
        #now get the assigned shifts, shifterrors, atomnames and residue numbers
        #from .prot, loop through the whole NoeList:
        for EPEAK in self.peakslist:
            for ECON in EPEAK.contributions:
                if self._PROT.atomdican.has_key(ECON.atomnumberh1):
                    ECON.ah1ppm = self._PROT.atomdican[ECON.atomnumberh1].shift
                    ECON.dah1ppm = self._PROT.atomdican[ECON.atomnumberh1].shifterror
                    ECON.atomnameh1 = self._PROT.atomdican[ECON.atomnumberh1].ariaatomname
                    ECON.residue1 = self._PROT.atomdican[ECON.atomnumberh1].fragmentnumber
                else:
                    if not ((ECON.atomnumberh1 == None) or (ECON.atomnumberh1 == '0')):
                        print 'WARNING: atomnumber', ECON.atomnumberh1, 'not found in .prot'
                if self._PROT.atomdican.has_key(ECON.atomnumberp1):
                    ECON.ap1ppm = self._PROT.atomdican[ECON.atomnumberp1].shift
                    ECON.dap1ppm = self._PROT.atomdican[ECON.atomnumberp1].shifterror
                    ECON.atomnamep1 = self._PROT.atomdican[ECON.atomnumberp1].ariaatomname
                    ECON.residue1 = self._PROT.atomdican[ECON.atomnumberp1].fragmentnumber
                else:
                    if not ((ECON.atomnumberp1 == None) or (ECON.atomnumberp1 == '0')):
                        print 'WARNING: atomnumber', ECON.atomnumberp1, 'not found in .prot'
                if self._PROT.atomdican.has_key(ECON.atomnumberh2):
                    ECON.ah2ppm = self._PROT.atomdican[ECON.atomnumberh2].shift
                    ECON.dah2ppm = self._PROT.atomdican[ECON.atomnumberh2].shifterror
                    ECON.atomnameh2 = self._PROT.atomdican[ECON.atomnumberh2].ariaatomname
                    ECON.residue2 = self._PROT.atomdican[ECON.atomnumberh2].fragmentnumber
                else:
                    if not ((ECON.atomnumberh2 == None) or (ECON.atomnumberh2 == '0')):
                        print 'WARNING: atomnumber', ECON.atomnumberh2, 'not found in .prot'
                if self._PROT.atomdican.has_key(ECON.atomnumberp2):
                    ECON.ap2ppm = self._PROT.atomdican[ECON.atomnumberp2].shift
                    ECON.dap2ppm = self._PROT.atomdican[ECON.atomnumberp2].shifterror
                    ECON.atomnamep2 = self._PROT.atomdican[ECON.atomnumberp2].ariaatomname
                    ECON.residue2 = self._PROT.atomdican[ECON.atomnumberp2].fragmentnumber
                else:
                    if not ((ECON.atomnumberp2 == None) or (ECON.atomnumberp2 == '0')):
                        print 'WARNING: atomnumber', ECON.atomnumberp2, 'not found in .prot'
        self.RemoveDoubleQuotes()


    def ReadPippPck(self, noeFile, het1, pro1, het2, pro2, assign1, assign2):
        """
        reads an PIPP .PCK file
        
        I use the same convention as in the .html file for PIPP:
        het1 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro1 = '1', '2', '3', '4' or 'N' (column number or not used)
        het2 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro2 = '1', '2', '3', '4' or 'N' (column number or not used)

        assign1 and assign2 are used to determine which 'Assign'
        column belongs to residue 1 or residue 2:
        assign1 = '1' or '2'
        assign2 = '1' or '2'

        all columns are numbered from left to right
        """
        #some messages:
        print 'reading the PIPP .PCK file:'
        print '  ' + noeFile

        #check if the files exist:
        _DoesFileExist(noeFile)
        self.fileNames.append(noeFile)
        
        #important - delete all the existing peaks:
        self._cleanUp()
        
        #read the general NOE file:
        noeHandle = TextFile.TextFile(noeFile)

        #unavailable data set to default:
        peakType = None
        segid1 = None
        segid2 = None
        volumeError = None
        distance = None
        distanceError = None
        lowerBound = None
        upperBound = None
        h1ppm = None
        h2ppm = None
        dh1ppm = None
        dp1ppm = None
        dh2ppm = None
        dp2ppm = None
        residue1 = None
        residue2 = None
        aa1 = None
        aa2 = None
        atomnameh1 = None
        atomnamep1 = None
        atomnameh2 = None
        atomnamep2 = None
            
        pipe = re.compile('|')
        pound = re.compile('#')
        semiColon = re.compile(';')
        startList = 0
        formatCol = []
        varsCol = []
        
        for eachLine in noeHandle:
            if pound.match(eachLine):
                continue
#            print eachLine, #test
            lineList = string.split(eachLine)
            if len(lineList) == 0:
                continue
            elif lineList[0] == 'FORMAT':
                formatCol = lineList[1:]
                continue
            elif lineList[0] == 'VARS':
                startList = 1
                varsCol = lineList[1:]
                continue
            elif lineList[0] == 'DATA' and lineList[1] == 'DIMCOUNT':
                self.dimension = lineList[2]
                continue
            elif startList == 0:
                continue

            #pkID (=peakNumber):
            peakNumber = lineList[varsCol.index('PkID')]
            
            #Sl.Z (=the Z slice the data was found on):
            sliceZ = lineList[varsCol.index('Sl.Z')]
            
            #X, Y, Z, A (=ppm 4D):
            if self.dimension == '4':
            	ppmList = [0, lineList[varsCol.index('X')], lineList[varsCol.index('Y')], \
               	       lineList[varsCol.index('Z')], lineList[varsCol.index('A')]]
            
            #X, Y, Z (=ppm 3D):
            if self.dimension == '3':
            	ppmList = [0, lineList[varsCol.index('X')], lineList[varsCol.index('Y')], \
               	       lineList[varsCol.index('Z')]]
            
            if het1 != 'N':
                exec 'h1ppm = ppmList[' + het1 + ']'
            if het2 != 'N':
                exec 'h2ppm = ppmList[' + het2 + ']'
            exec 'p1ppm = ppmList[' + pro1 + ']'
            exec 'p2ppm = ppmList[' + pro2 + ']'
            
            #Intensity:
            volume = lineList[varsCol.index('Intensity')]

            asterisk = re.compile('\*')
            
            #Assign1, Assign2 (not mandatory!):
            if len(varsCol) > 6:
                #get the whole entry (e.g. 'E68.CB,HB1|HB2'):
                get1 = semiColon.sub('', lineList[varsCol.index('Assign1')])
                try:
                    get2 = lineList[varsCol.index('Assign2')]
                except IndexError:
                    get2 = '****'
                if not asterisk.search(get1):
                    #for the first column (Assign1):
                    #split after the '.' (e.g. E68.CB,HB1|HB2' => ['E68', 'CB,HB1|HB2']):
                    get1List = string.split(get1, '.')
                    #get 3-letter code and residue number:
                    outS = 'aa' + assign1 + ' = AminoAcid.AminoAcid(get1List[0][0])[1];' +\
                           'residue' + assign1 + ' = get1List[0][1:];'
                    #split after the ',' (e.g. 'CB,HB1|HB2' => ['CB', 'HB1|HB2']):
                    atom1List = string.split(get1List[1], ',')
                    #form tuples to get all the possibilities within the same list
                    #e.g. ['CB', 'HB1|HB2']  => [('CB'), ('HB1', 'HB2')]:
                    atom1List[0] = tuple(string.split(atom1List[0], '|'))
                    if len(atom1List) == 2:
                        atom1List[1] = tuple(string.split(atom1List[1], '|'))
                    #now assign the stuff according to the specified columns:
                    outS = outS + 'atomnamep' + assign1 + ' = atom1List[-1];'
                    if len(atom1List) == 2:
                        outS = outS + 'atomnameh' + assign1 + ' = atom1List[-2]'
                    exec outS
                else:
                    outS = 'atomnamep' + assign1 + ' = (None,);'
                    exec outS
                    
                if not asterisk.search(get2):
                    #for the second column (Assign2):
                    #split after the '.' (e.g. E68.CB,HB1|HB2' => ['E68', 'CB,HB1|HB2']):
                    get2List = string.split(get2, '.')
                    #get 3-letter code and residue number:
                    outS = 'aa' + assign2 + ' = AminoAcid.AminoAcid(get2List[0][0])[1];' +\
                           'residue' + assign2 + ' = get2List[0][1:];'
                    #split after the ',' (e.g. 'CB,HB1|HB2' => ['CB', 'HB1|HB2']):
                    atom2List = string.split(get2List[1], ',')
                    #form tuples to get all the possibilities within the same list
                    #e.g. ['CB', 'HB1|HB2']  => [('CB'), ('HB1', 'HB2')]:
                    atom2List[0] = tuple(string.split(atom2List[0], '|'))
                    if len(atom2List) == 2:
                        atom2List[1] = tuple(string.split(atom2List[1], '|'))
                    #now assign the stuff according to the specified columns:
                    outS = outS + 'atomnamep' + assign2 + ' = atom2List[-1];'
                    if len(atom2List) == 2:
                        outS = outS + 'atomnameh' + assign2 + ' = atom2List[-2]'
                    exec outS
                

            #put it in a contribution instance:
            CONT = NoeContribution(peakNumber = peakNumber,\
                                   peakType = peakType,\
                                   residue1 = residue1,\
                                   aa1 = aa1,\
                                   segid1 = segid1,\
                                   atomnameh1 = atomnameh1,\
                                   atomnamep1 = atomnamep1,\
                                   residue2 = residue2,\
                                   aa2 = aa2,\
                                   segid2 = segid2,\
                                   atomnameh2 = atomnameh2,\
                                   atomnamep2 = atomnamep2,\
                                   volume = volume,\
                                   volumeError = volumeError,\
                                   distanceAve = distance,\
                                   distanceStd = distanceError,\
                                   lowerBound = lowerBound,\
                                   upperBound = upperBound,\
                                   h1ppm = h1ppm,\
                                   dh1ppm = dh1ppm,\
                                   p1ppm = p1ppm,\
                                   dp1ppm = dp1ppm,\
                                   h2ppm = h2ppm,\
                                   dh2ppm = dh2ppm,\
                                   p2ppm = p2ppm,\
                                   dp2ppm = dp2ppm)
            
            #create a peak instance and add contribution:
            NOE = NoePeak()
            NOE.AddContribution(CONT)

            #add peak to the NOE list:
            self.AddPeak(NOE)

        
    def ReadRegine(self, fileName):
        """
        reads an Regine-derived 2D, 3D or 4D peaks file

        the number of columns specifies the dimensionality
        
        2D format: whitespace separated fields
                   peakNumber ppmProt1 ppmProt2 Intensity Volume
    
        3D format: whitespace separated fields
                   peakNumber ppmHet1 ppmProt1 ppmProt2 Intensity Volume
    
        4D format: whitespace separated fields
                   peakNumber ppmHet1 ppmProt1 ppmHet2 ppmProt2 Intensity Volume

        2D example (5 columns):
        1      8.535    5.528     13449000    310650000
        2      8.534    1.953     16301000    385720000
        3      9.405    5.329      7526900    160030000

        3D example (6 columns):
        1  117.566    8.535    5.528     13449000    310650000
        2  117.582    8.534    1.953     16301000    385720000
        3  123.610    9.405    5.329      7526900    160030000

        4D example (7 columns):
        1  117.566    8.535    51.435   5.528     13449000    310650000
        2  117.582    8.534    45.235   1.953     16301000    385720000
        3  123.610    9.405    50.435   5.329      7526900    160030000
        

        We always use the volumes for the calibration in ARIA!
        Edit the code if you want to change this behaviour.
        """
        self.fileNames.append(fileName)

        #check, if file exists:
        _DoesFileExist(fileName)
        
        #important: delete all the existing peaks
        self._cleanUp()

        #fileHandle:
        fileHandle = TextFile.TextFile(fileName)

        getDim = 0
        lineCounter = 0
        for eachLine in fileHandle:
            lineCounter = lineCounter + 1
            lineList = string.split(eachLine)
            if len(lineList) < 4:
                continue
            #get the dimensionality from the first data line:
            if not getDim:
                if len(lineList) == 7:
                    getDim = 4
                elif len(lineList) == 6:
                    getDim = 3
                else:
                    getDim = 2
            #print warnings:
            if getDim != (len(lineList) - 3):
                print 'WARNING: number of columns not correct in shifts file!'
                print '         row no.', str(lineCounter), 'incorrect'
            #new peak instance:
            PEAK = NoePeak()
            #get the data:
            if getDim == 2:
                peakNumber = lineList[0]
                h1ppm = None
                p1ppm = lineList[1]
                h2ppm = None
                p2ppm = lineList[2]
                intensity = lineList[3]
                volume = lineList[4]
            elif getDim == 3:
                peakNumber = lineList[0]
                h1ppm = lineList[1]
                p1ppm = lineList[2]
                h2ppm = None
                p2ppm = lineList[3]
                intensity = lineList[4]
                volume = lineList[5]
            elif getDim == 4:
                peakNumber = lineList[0]
                h1ppm = lineList[1]
                p1ppm = lineList[2]
                h2ppm = lineList[3]
                p2ppm = lineList[4]
                intensity = lineList[5]
                volume = lineList[6]
            #put it into contribution instance:
            CONT = NoeContribution(spectrumName = None,\
                                   peakNumber = peakNumber,\
                                   testSet = None,\
                                   inWeight = None,\
                                   volume = volume,\
                                   volumeError = None,\
                                   intensity = intensity,\
                                   intensityError = None,\
                                   h1ppm = h1ppm,\
                                   dh1ppm = None,\
                                   p1ppm = p1ppm,\
                                   dp1ppm = None,\
                                   h2ppm = h2ppm,\
                                   dh2ppm = None,\
                                   p2ppm = p2ppm,\
                                   dp2ppm = None,\
                                   peakType = None,\
                                   fomAria = None,\
                                   curWeight = None,\
                                   lowerBound = None,\
                                   upperBound = None,\
                                   sumDistanceAve = None,\
                                   sumDistanceStd = None,\
                                   backVolumeAria = None,\
                                   backVolumeAriaStd = None,\
                                   allAssi = None,\
                                   nta = None,\
                                   contribution = None,\
                                   fomContribution = None,\
                                   distanceAve = None,\
                                   distanceStd = None,\
                                   backVolume = None,\
                                   backVolumeStd = None,\
                                   segid1 = None,\
                                   residue1 = None,\
                                   aa1 = None,\
                                   atomnameh1 = (None,),\
                                   assih1ppm = None,\
                                   assidh1ppm = None,\
                                   atomnamep1 = (None,),\
                                   assip1ppm = None,\
                                   assidp1ppm = None,\
                                   segid2 = None,\
                                   residue2 = None,\
                                   aa2 = None,\
                                   atomnameh2 = (None,),\
                                   assih2ppm = None,\
                                   assidh2ppm = None,\
                                   atomnamep2 = (None,),\
                                   assip2ppm = None,\
                                   assidp2ppm = None)

            #add contribution to peak:
            PEAK.AddContribution(CONT)
            
            #create a peak instance:
            NOE = NoePeak()
            
            #add contribution to NOE instance:
            NOE.AddContribution(CONT)

            #add peak to the NOE list:
            self.AddPeak(NOE)
    
    def ReadTbl(self, noeFile):
        """
        reads an ARIA .tbl file

        typical lines can look like that:

        ASSI {    2}
           (( segid "    " and resid 3    and name HN  ))
           (( segid "    " and resid 3    and name HA  ))
           2.600 0.800 0.800 peak 2 weight 0.10000E+01 volume 0.96210E+07 ppm1 7.997 ppm2 4.846

        assi ( attr store5 < 2.177 and attr store6 > 2.097 )
             ( attr store5 < 2.407 and attr store6 > 2.367 )
             6.0 0.1 0.1 peak 12 volume 0.000e+00 ppm1 2.137 ppm2 2.387

        assign (resid 1 and name ha# ) (resid 2 and name HN ) 2.55 0.75 0.75
        assign (resid 1 and name ha# ) (resid 2 and name HA ) 3.65 1.85 1.85

        the 'OR' is also understood.
        It's always a better idea to read the corresponding .list files because
        they contain more information than the .tbl files.
        The contents of the .list file are stored in the attribute _tblFileContent
        This is important for the free R factor refinement.
        """
        #some messages:
        print 'reading the ARIA .tbl file:'
        print '  ' + noeFile

        #check if the files exist:
        _DoesFileExist(noeFile)
        self.fileNames.append(noeFile)
        
        #important - delete all the existing peaks:
        self._cleanUp()
        
        #read the general NOE file:
        bigString = DeleteCnsComments.GetString(noeFile)
        assiList = string.split(bigString, 'assign')[1:]
        if len(assiList) == 0:
            assiList = string.split(bigString, 'ASSI')[1:]
        if len(assiList) == 0:
            assiList = string.split(bigString, 'assi')[1:]
        if len(assiList) == 0:
            assiList = string.split(bigString, 'Assi')[1:]

#        print assiList #test

        #IMPORTANT:
        #store the original content of the peak list in the attribute
        #_tblFileContent (e.g. important for free R factor refinement): 
        self._tblFileContent = assiList
        
        #compile some patterns:
        segidP = re.compile('segid\s*"(.{0,4})"\s*')
        residP = re.compile('resid\s*([0-9]+)[\s\)]+')
        nameP = re.compile('name\s*([\w!@#$%&*?|]+)\s*\)')
        store1LessP = re.compile('store5\s*<\s*(\d+)')
        store1GreatP = re.compile('store6\s*>\s*(\d+)')
        
## ARIA 1.0 version used the three following lines:
##         nameP = re.compile('name\s*([\w!@#$%&*?|\'\"]+)\s*\)')
##         store1LessP = re.compile('store1\s*<\s*(\d+)')
##         store1GreatP = re.compile('store1\s*>\s*(\d+)')
        
        distancesP = re.compile('\)\s+([0-9.\-eE+]+)\s+([0-9.\-eE+]+)\s+([0-9.\-eE+]+)\s+')
        peakP = re.compile('\s+peak\s+([0-9.\-eE+]+)\s+')
        weightP = re.compile('\s+weight\s+([0-9.\-eE+]+)\s+')
        volumeP = re.compile('\s+volume\s+([0-9.\-eE+]+)\s+')
        ppm1P = re.compile('\s+ppm1\s+([0-9.\-eE+]+)\s+')
        ppm2P = re.compile('\s+ppm2\s+([0-9.\-eE+]+)\s+')
        # BARDIAUX test
        ppm1H = re.compile('\s+hpm1\s+([0-9.\-eE+]+)\s+')
        ppm2H = re.compile('\s+hpm2\s+([0-9.\-eE+]+)\s+')
        #for creating peaknumbers (if there are no peaknumbers already):
        countPeaks = 1

        for eachAssi in assiList:
##             print eachAssi #test
            #unavailable data set to default:
            peakNumber = None
            peakType = None
            segid1 = None
            segid2 = None
            volume = None
            volumeError = None
            distance = None
            distanceError = None
            lowerBound = None
            upperBound = None
            h1ppm = None
            h2ppm = None
            p1ppm = None
            p2ppm = None
            dh1ppm = None
            dp1ppm = None
            dh2ppm = None
            dp2ppm = None
            residue1 = None
            residue2 = None
            aa1 = None
            aa2 = None
            atomnameh1 = (None,)
            atomnamep1 = (None,)
            atomnameh2 = (None,)
            atomnamep2 = (None,)
            
            #consider the 'OR':
            eachList = string.split(eachAssi, 'OR ')
            if len(eachList) == 1:
                eachList = string.split(eachAssi, 'or ')

            #create a peak instance:
            NOE = NoePeak()

            for eachCont in eachList:
                #mandatory:
                distancesL = distancesP.search(eachCont)
                if distancesL:
                    distance = distancesL.group(1)
                    upperBound = distancesL.group(2)
                    lowerBound = distancesL.group(3)
                peakM = peakP.search(eachCont)
                if peakM:
                    peakNumber = peakM.group(1)
                else:
                    peakNumber = str(countPeaks)
                    countPeaks = countPeaks + 1
                volumeM = volumeP.search(eachCont)
                if volumeM:
                    volume = volumeM.group(1)
                ppm1M = ppm1P.search(eachCont)
                if ppm1M:
                    p1ppm = ppm1M.group(1)
                ppm2M = ppm2P.search(eachCont)
                if ppm2M:
                    p2ppm = ppm2M.group(1)
                # BARDIAUX
                ppm1HM = ppm1H.search(eachCont)
                if ppm1HM:
                    h1ppm = ppm1HM.group(1)
                ppm2HM = ppm2H.search(eachCont)
                if ppm2HM:
                    h2ppm = ppm2HM.group(1)
                    
                #others:
                segidL = segidP.findall(eachCont)
                if len(segidL) > 0:
                    segid1 = segidL[0]
                    segid2 = segidL[1]
                residL = residP.findall(eachCont)
                if len(residL) > 0:
                    residue1 = residL[0]
                    residue2 = residL[1]
                nameL  = nameP.findall(eachCont)
#                print nameL #test
                if len(nameL) > 0:
                    atomnamep1 = (nameL[0],)
                    atomnamep2 = (nameL[1],)
                weightMatch = weightP.search(eachCont)
                if weightMatch:
                    weight = weightMatch.group(1)

                #put it in a contribution instance:
                CONT = NoeContribution(peakNumber = peakNumber,\
                                       peakType = peakType,\
                                       residue1 = residue1,\
                                       aa1 = aa1,\
                                       segid1 = segid1,\
                                       atomnameh1 = atomnameh1,\
                                       atomnamep1 = atomnamep1,\
                                       residue2 = residue2,\
                                       aa2 = aa2,\
                                       segid2 = segid2,\
                                       atomnameh2 = atomnameh2,\
                                       atomnamep2 = atomnamep2,\
                                       volume = volume,\
                                       volumeError = volumeError,\
                                       distanceAve = distance,\
                                       distanceStd = distanceError,\
                                       lowerBound = lowerBound,\
                                       upperBound = upperBound,\
                                       h1ppm = h1ppm,\
                                       dh1ppm = dh1ppm,\
                                       p1ppm = p1ppm,\
                                       dp1ppm = dp1ppm,\
                                       h2ppm = h2ppm,\
                                       dh2ppm = dh2ppm,\
                                       p2ppm = p2ppm,\
                                       dp2ppm = dp2ppm)

                #add contribution to NOE instance:
                NOE.AddContribution(CONT)

            #add peak to the NOE list:
            self.AddPeak(NOE)


    def ReadAureliaUserPeak(self):
        """
        test
        """
        return
    

    def ReadAureliaTenFields(self, noeFile, pro1, pro2):
        """
        reads an AURELIA file with a format like that:
        hypo  peak_id  cluster  center(w1,w2)  match(xy,z)  intensity(abs)  name

     1     1      428     2.732, 10.138   1.00, 1.00      277766   40HB2 41HE1
     2     1       79     7.281, 10.138   1.00, 1.00      119324   
     3     1       71     7.431, 10.138   1.00, 1.00     1462283   41HZ2 41HE1
     4     1      497     2.367, 10.138   1.00, 1.00      324854   40HB3 41HE1.s  

     the .s symmetry identifiers are thrown away
     pro1: proton1 found in column '1' or '2'
     pro2: proton2 found in column '1' or '2'
     """
        #some messages:
        print 'reading the file:'
        print '  ' + noeFile

        #check if the files exist:
        _DoesFileExist(noeFile)
        self.fileNames.append(noeFile)
        
        #important - delete all the existing peaks:
        self._cleanUp()
        
        #read the general NOE file:
        noeHandle = TextFile.TextFile(noeFile)

        #unavailable data set to default:
        peakType = None
        segid1 = None
        segid2 = None
        volumeError = None
        distance = None
        distanceError = None
        lowerBound = None
        upperBound = None
        h1ppm = None
        h2ppm = None
        dh1ppm = None
        dp1ppm = None
        dh2ppm = None
        dp2ppm = None
        residue1 = None
        residue2 = None
        aa1 = None
        aa2 = None
        atomnameh1 = None
        atomnamep1 = None
        atomnameh2 = None
        atomnamep2 = None
            
        startList = 0
        formatCol = []
        varsCol = []

        number = re.compile('(^\d+)(\S+)')
        
        for eachLine in noeHandle:
            if string.strip(eachLine) == '':
                continue
            lineList = string.split(eachLine)
            if len(lineList) < 8 or lineList[0] == 'hypo':
                continue
            
            peakNumber = lineList[0]
            peakId= lineList[1]
            cluster= lineList[2]
            outS = ''
            outS = 'p' + pro1 + 'ppm =  lineList[3][:-1];' #without trailing ','
            outS = outS + 'p' + pro2 + 'ppm = lineList[4]'
            exec outS
            matchxy = lineList[5][:-1]  #without trailing ','
            matchz = lineList[6]
            volume = lineList[7]

            #defaults for the atomnames:
            atomnamep1 = (None,)
            atomnamep2 = (None,)
            if len(lineList) > 9:
                name1 = lineList[8]
                name2 = lineList[9]
                if name1[-2:] == '.s': name1 = name1[:-2]
                if name2[-2:] == '.s': name2 = name2[:-2]
                #get the residue numbers and atomname tuples:
                outS = 'getNumber1 = number.search(name' + pro1 + ');'
                outS = outS + 'getNumber2 = number.search(name' + pro2 + ')'
                exec outS
                residue1 = getNumber1.group(1)
                residue2 = getNumber2.group(1)
                atomnamep1 = PseudoAtom.Pseudo2Tuple(getNumber1.group(2))
                atomnamep2 = PseudoAtom.Pseudo2Tuple(getNumber2.group(2))
                
            #put it in a contribution instance:
            CONT = NoeContribution(peakNumber = peakNumber,\
                                   peakType = peakType,\
                                   residue1 = residue1,\
                                   aa1 = aa1,\
                                   segid1 = segid1,\
                                   atomnameh1 = atomnameh1,\
                                   atomnamep1 = atomnamep1,\
                                   residue2 = residue2,\
                                   aa2 = aa2,\
                                   segid2 = segid2,\
                                   atomnameh2 = atomnameh2,\
                                   atomnamep2 = atomnamep2,\
                                   volume = volume,\
                                   volumeError = volumeError,\
                                   distanceAve = distance,\
                                   distanceStd = distanceError,\
                                   lowerBound = lowerBound,\
                                   upperBound = upperBound,\
                                   h1ppm = h1ppm,\
                                   dh1ppm = dh1ppm,\
                                   p1ppm = p1ppm,\
                                   dp1ppm = dp1ppm,\
                                   h2ppm = h2ppm,\
                                   dh2ppm = dh2ppm,\
                                   p2ppm = p2ppm,\
                                   dp2ppm = dp2ppm)
            
            #create a peak instance and add contribution:
            NOE = NoePeak()
            NOE.AddContribution(CONT)

            #add peak to the NOE list:
            self.AddPeak(NOE)
            


    def ReadAureliaNineFields(self, noeFile, pro1, pro2):
        """
        reads an AURELIA file with a format like that:
!  hypothesis     w1     w2   intensity      volume     di    dv    sym_di sym_dv
          10   9.047  9.735     159232    17308568   3.18  3.0/1
          12   8.918  9.599     163728    17384984   3.17  3.0/1   3.14  3.10/1
          22   8.522  8.867     125160    16040688   3.31  3.0/1   3.07  3.01/1
          23   8.467  8.935     112024    13173216   3.37  3.1/1   3.08  3.06/1

     the .s symmetry identifiers are thrown away
     pro1: proton1 found in column '1' or '2'
     pro2: proton2 found in column '1' or '2'
     """
        #some messages:
        print 'reading the file:'
        print '  ' + noeFile

        self.dimension = 2
        
        #check if the files exist:
        _DoesFileExist(noeFile)
        self.fileNames.append(noeFile)
        
        #important - delete all the existing peaks:
        self._cleanUp()
        
        #read the general NOE file:
        noeHandle = TextFile.TextFile(noeFile)

        #unavailable data set to default:
        peakType = None
        segid1 = None
        segid2 = None
        volumeError = None
        distance = None
        distanceError = None
        lowerBound = None
        upperBound = None
        h1ppm = None
        h2ppm = None
        dh1ppm = None
        dp1ppm = None
        dh2ppm = None
        dp2ppm = None
        residue1 = None
        residue2 = None
        aa1 = None
        aa2 = None
        atomnameh1 = None
        atomnamep1 = None
        atomnameh2 = None
        atomnamep2 = None
            
        startList = 0
        formatCol = []
        varsCol = []

        number = re.compile('(^\d+)(\S+)')
        
        for eachLine in noeHandle:
            if string.strip(eachLine) == '':
                continue
            lineList = string.split(eachLine)
            if len(lineList) < 6 or lineList[0] == 'hypo' or \
               lineList[0] == '!' or lineList[0][1] == '!':
                continue
            
            peakNumber = lineList[0]
            peakId= '1'
            cluster= '0'
            outS = ''
            outS = 'p' + pro1 + 'ppm =  lineList[1];'
            outS = outS + 'p' + pro2 + 'ppm = lineList[2]'
            exec outS
            matchxy = '0' 
            matchz = '0'
            volume = lineList[4]

            #defaults for the atomnames:
            atomnamep1 = (None,)
            atomnamep2 = (None,)
                
            #put it in a contribution instance:
            CONT = NoeContribution(peakNumber = peakNumber,\
                                   peakType = peakType,\
                                   residue1 = residue1,\
                                   aa1 = aa1,\
                                   segid1 = segid1,\
                                   atomnameh1 = atomnameh1,\
                                   atomnamep1 = atomnamep1,\
                                   residue2 = residue2,\
                                   aa2 = aa2,\
                                   segid2 = segid2,\
                                   atomnameh2 = atomnameh2,\
                                   atomnamep2 = atomnamep2,\
                                   volume = volume,\
                                   volumeError = volumeError,\
                                   distanceAve = distance,\
                                   distanceStd = distanceError,\
                                   lowerBound = lowerBound,\
                                   upperBound = upperBound,\
                                   h1ppm = h1ppm,\
                                   dh1ppm = dh1ppm,\
                                   p1ppm = p1ppm,\
                                   dp1ppm = dp1ppm,\
                                   h2ppm = h2ppm,\
                                   dh2ppm = dh2ppm,\
                                   p2ppm = p2ppm,\
                                   dp2ppm = dp2ppm)
            
            #create a peak instance and add contribution:
            NOE = NoePeak()
            NOE.AddContribution(CONT)
            NOE.dimension = 2
            #add peak to the NOE list:
            self.AddPeak(NOE)

    def ReadAureliaStefania(self, noeFile):
        """
        The AURELIA file formats are a major catastrophy
        here comes yet another method to read those files
        
        reads an AURELIA file with a format like that:
!peak        w1             w2             w3         intensity    name
!----------------------------------------------------------------------
# 2367    0.248 (466)  117.053 ( 66)    6.924 (315)      153951  52_HN/51_HD1
#  788    6.416 (181)  127.687 ( 17)    5.745 (424)      169412  112_HD21/112_HD22
#  726    6.848 (161)  115.533 ( 73)    8.958 (127)       97188  49_HN/51_HN
#  727    6.848 (161)  117.704 ( 63)    8.774 (144)      233523  727
#  724    6.848 (161)  110.976 ( 94)    6.697 (336)      238623  724
#  725    6.848 (161)  110.976 ( 94)    7.011 (307)      312605  725
#    2   10.289 (  2)  126.385 ( 23)    8.060 (210)       62500  2

         empty lines and lines starting with a '!' are not parsed
      
         """
        #some messages:
        print 'reading the file:'
        print '  ' + noeFile

        self.dimension = 3
        
        #check if the files exist:
        _DoesFileExist(noeFile)
        self.fileNames.append(noeFile)
        
        #important - delete all the existing peaks:
        self._cleanUp()
        
        #read the general NOE file:
        noeHandle = TextFile.TextFile(noeFile)

        #unavailable data set to default:
        peakType = None
        segid1 = None
        segid2 = None
        volumeError = None
        distance = None
        distanceError = None
        lowerBound = None
        upperBound = None
        h1ppm = None
        h2ppm = None
        dh1ppm = None
        dp1ppm = None
        dh2ppm = None
        dp2ppm = None
        residue1 = None
        residue2 = None
        aa1 = None
        aa2 = None
        atomnameh1 = None
        atomnamep1 = None
        atomnameh2 = None
        atomnamep2 = None

        reLine = re.compile('(#)?\s*([\d\.]+)\s+([\d\.,Ee-]+)\s+\(\s*(\d*)\)\s+([\d\.,Ee-]+)\s+\(\s*(\d*)\)\s+([\d\.,Ee-]+)\s+\(\s*(\d*)\)\s+([\d\.,Ee-]+)\s+(\S*)')
        reAssi = re.compile('(\d)+_(\w#\*%)+/(\d)+_(\w#\*%)+')
        reExcla = re.compile('\s*!')
        for eachLine in noeHandle:
            if string.strip(eachLine) == '':
                continue
            if reExcla.match(eachLine):
                continue
#            print eachLine #test
            lineList = string.split(eachLine)
            if len(lineList) < 6 or lineList[0] == 'peaks':
                continue
            meLine = reLine.search(eachLine)
            if not meLine:
                print 'could not read the following line:'
                print eachLine
                continue

            print meLine.group(2)
            print meLine.group(3)
            print meLine.group(4)
            print meLine.group(5)
            print meLine.group(6)
            print meLine.group(7)
            print meLine.group(8)
            print meLine.group(9)
            

            peakNumber = meLine.group(2)
            volume = meLine.group(9)

            #ppms:
            h1ppm = meLine.group(5)
            h2ppm = None
            p1ppm = meLine.group(7)
            p2ppm = meLine.group(3)
            
            #defaults for the atomnames:
            atomnamep1 = (None,)
            atomnamep2 = (None,)

            #try to understand the assignments:
            assiS = meLine.group(10)
            seAssi = reAssi.search(assiS)
            if seAssi:
                aa1=seAssi.group(1)
                atomnamep1 =(seAssi.group(2),)
                aa2=seAssi.group(3)
                atomnamep2 =(seAssi.group(4),)
            else:
                pass #could not understand assignment

                
            #put it in a contribution instance:
            CONT = NoeContribution(peakNumber = peakNumber,\
                                   peakType = peakType,\
                                   residue1 = residue1,\
                                   aa1 = aa1,\
                                   segid1 = segid1,\
                                   atomnameh1 = atomnameh1,\
                                   atomnamep1 = atomnamep1,\
                                   residue2 = residue2,\
                                   aa2 = aa2,\
                                   segid2 = segid2,\
                                   atomnameh2 = atomnameh2,\
                                   atomnamep2 = atomnamep2,\
                                   volume = volume,\
                                   volumeError = volumeError,\
                                   distanceAve = distance,\
                                   distanceStd = distanceError,\
                                   lowerBound = lowerBound,\
                                   upperBound = upperBound,\
                                   h1ppm = h1ppm,\
                                   dh1ppm = dh1ppm,\
                                   p1ppm = p1ppm,\
                                   dp1ppm = dp1ppm,\
                                   h2ppm = h2ppm,\
                                   dh2ppm = dh2ppm,\
                                   p2ppm = p2ppm,\
                                   dp2ppm = dp2ppm)
            
            #create a peak instance and add contribution:
            NOE = NoePeak()
            NOE.AddContribution(CONT)
            NOE.dimension = 2
            #add peak to the NOE list:
            self.AddPeak(NOE)
        


    def ReadNMRViewPeaksOld(self, fileName):
        """
        reads a NMRView-derived 2D, 3D or 4D peaks file 
        the number of columns specifies the dimensionality

        2D format: whitespace separated fields
                   peakNumber ppmProt1 ppmProt2 Volume  ass1   ass2

        3D format: whitespace separated fields
                   peakNumber ppmHet1 ppmProt1 ppmProt2 Volume  ass1  ass2  ass3

        4D format: whitespace separated fields
                   peakNumber ppmHet1 ppmProt1 ppmHet2 ppmProt2 Volume  ass1  ass2  ass3

        2D example (8 columns):

        3D example (11 columns):

        4D example (14 columns):

        """
        self.fileNames.append(fileName)

        #check, if file exists:
        _DoesFileExist(fileName)

        #important: delete all the existing peaks
        self._cleanUp()

        #fileHandle:
        fileHandle = TextFile.TextFile(fileName)

        getDim = 0
        lineCounter = 0
        for eachLine in fileHandle:
            lineCounter = lineCounter + 1
            lineList = string.split(eachLine)
            if len(lineList) < 4:
                continue
            #get the dimensionality from the first data line:
            if not getDim:
                if len(lineList) == 14:
                    getDim = 4
                    print 'Assuming 4D data'
                elif len(lineList) == 11:
                    getDim = 3
                    print 'Assuming 3D data'
                else:
                    getDim = 2
                    print 'Assuming 2D data'
            #print warnings:
            if (len(lineList) - (getDim*3)) != 2:
                print 'WARNING: number of columns not correct in xpk file!'
                print '         row no.', str(lineCounter), 'incorrect'
            self.dimension = getDim
            #new peak instance:
            PEAK = NoePeak()
            #get the data:
            rassh1 = None
            nassh1 = None
            rassp1 = None
            nassp1 = None
            rassh2 = None
            nassh2 = None
            rassp2 = None
            nassp2 = None
            if getDim == 2:
                peakNumber = lineList[0]
                h1ppm = None
                p1ppm = lineList[1]
                h2ppm = None
                p2ppm = lineList[2]
                volume = lineList[3]
                rass1 = lineList[4]
                nass1 = PseudoAtom.Pseudo2Atom(lineList[5])
                rass2 = lineList[6]
                nass2 = PseudoAtom.Pseudo2Atom(lineList[7])
            elif getDim == 3:
                peakNumber = lineList[0]
                h1ppm = lineList[1]
                p1ppm = lineList[2]
                h2ppm = None
                p2ppm = lineList[3]
                volume = lineList[4]
                rassh1 = lineList[5]
                nassh1 = lineList[6]
                rass1 = lineList[7]
                nass1 = PseudoAtom.Pseudo2Atom(lineList[8])
                rass2 = lineList[9]
                nass2 = PseudoAtom.Pseudo2Atom(lineList[10])
            elif getDim == 4:
                peakNumber = lineList[0]
                h1ppm = lineList[1]
                p1ppm = lineList[2]
                h2ppm = lineList[3]
                p2ppm = lineList[4]
                volume = lineList[5]
                rassh1 = lineList[6]
                nassh1 = lineList[7]
                rass1 = lineList[8]
                nass1 = PseudoAtom.Pseudo2Atom(lineList[9])
                rassh2 = lineList[10]
                nassh2 = lineList[11]
                rass2 = lineList[12]
                nass2 = PseudoAtom.Pseudo2Atom(lineList[13])
            #put it into contribution instance:
        if rassh1 == '999':
            rassh1 = None
            nassh1 = None
        if rass1 == '999':
            rass1 = None
            nass1 = None
        if rassh2 == '999':
            rassh2 = None
            nassh2 = None
        if rass2 == '999':
            rass1 = None
            nass1 = None
            CONT = NoeContribution(spectrumName = None,\
                                   peakNumber = peakNumber,\
                                   testSet = None,\
                                   volume = volume,\
                                   volumeError = None,\
                                   intensity = volume,\
                                   intensityError = None,\
                                   h1ppm = h1ppm,\
                                   dh1ppm = None,\
                                   p1ppm = p1ppm,\
                                   dp1ppm = None,\
                                   h2ppm = h2ppm,\
                                   dh2ppm = None,\
                                   p2ppm = p2ppm,\
                                   dp2ppm = None,\
                                   peakType = None,\
                                   fomAria = None,\
                                   lowerBound = None,\
                                   upperBound = None,\
                                   sumDistanceAve = None,\
                                   sumDistanceStd = None,\
                                   backVolumeAria = None,\
                                   backVolumeAriaStd = None,\
                                   allAssi = None,\
                                   nta = None,\
                                   contribution = None,\
                                   fomContribution = None,\
                                   distanceAve = None,\
                                   distanceStd = None,\
                                   backVolume = None,\
                                   backVolumeStd = None,\
                                   segid1 = None,\
                                   residue1 = rass1,\
                                   aa1 = None,\
                                   atomnameh1 = (nassh1,),\
                                   assih1ppm = None,\
                                   assidh1ppm = None,\
                                   atomnamep1 = (nass1,),\
                                   assip1ppm = None,\
                                   assidp1ppm = None,\
                                   segid2 = None,\
                                   residue2 = rass2,\
                                   aa2 = None,\
                                   atomnameh2 = (nassh2,),\
                                   assih2ppm = None,\
                                   assidh2ppm = None,\
                                   atomnamep2 = (nass2,),\
                                   assip2ppm = None,\
                                   assidp2ppm = None)

            #add contribution to peak:
            PEAK.AddContribution(CONT)

            #create a peak instance:
            NOE = NoePeak()

            #add contribution to NOE instance:
            NOE.AddContribution(CONT)

            #add peak to the NOE list:
            self.AddPeak(NOE)


    def ReadNMRViewPeaks(self, fileName, het1, pro1, het2, pro2,\
                         readOnlyHeader=None):
        """
reads a NMRView 2D, 3D or 4D peaks file 

The headers always have 6 lines.

The second line contains the names for the 4 dimensions, defined by the
user.

The sixth line contains the colums of data that follow.  The X.P columns
will contain ppm data and the int column will contain the intensity.  X is
one of the labels defined in the 2nd line of the header.

Note that the data lines first contain the peaknumber before the data so
it is 'off by one' compared to the definition line.  For example in the
4d data H.P is the second data listed in the header but it is colum 3 in
the data.

Last item, lines 4 and 5 of the header contain the spectral sweepwidths
and frequencies for the dimensions so that one could unfold folded spectra
until one has tested all possible values for the range of a nucleus.  An
additional option would be for the user to set an f1180 flag so that
positive and negative peaks would be treated differently.

        2D format:
label dataset sw sf
H1 H2
noe120_h2oref_add_res.nv
{6982.63 } {6982.63 }
{599.9090 } {599.9090 }
 H1.L  H1.P  H1.W  H1.B  H1.E  H1.J  H1.U  H2.L  H2.P  H2.W  H2.B  H2.E  H2.J  H2.U  vol  int  stat  comment  flag0
0  {42.hn}   9.428   0.037   0.052   ?   0.000   {?}   {42.hb1}   2.353   0.050   0.071   ++   0.000   {?}  73432.78906 10080670.00000 0 {?}
0
1  {42.hn}   9.428   0.035   0.049   ++   0.000   {?}   {42.hb2}   1.456   0.061   0.086   ++   0.000   {?}  268296.25000 8933272.00000 0 {?} 0
2  {42.hn}   9.429   0.031   0.044   ++   0.000   {?}   {42.ha}   5.399   0.049   0.070   ++   0.000   {?}  88934.37500 10483782.00000 0 {?}
0
3  {33.hn}   10.028   0.043   0.061   ++   0.000   {?}   {33.ha1}   5.027   0.062   0.088   ++   0.000   {?}  50151.48438 3166743.00000 0 {?} 0


        3D format:
label dataset sw sf 
H HA C13 
hhc.nv
{8000.00 } {8000.00 } {8000.00 }
{750.1071 } {750.1071 } {188.6273 }
 H.L  H.P  H.W  H.B  H.E  H.J  H.U  HA.L  HA.P  HA.W  HA.B  HA.E  HA.J  HA.U  C13.L  C13.P  C13.W  C13.B  C13.E  C13.J  C13.U  vol  int  stat  comment  flag0 
0  {59.ha}   3.863   0.039   0.039   ++   0.000   {?}   {59.hn}   7.971   0.147   0.080   ++   0.000   {?}   {59.ca}   52.865   0.459   0.247   ++   0.000   {?}  0.28357 0.28357 0 {?} 0
1  {59.HA}   3.864   0.035   0.099   ++   0.000   {?}   {59.HA}   3.867   0.081   0.186   ++   0.000   {?}   {59.CA}   52.888   0.545   1.272   ?   0.000   {?}  0.00000 5.24471 0 {?} 0
2  {59.ha 59.ha}   3.863   0.040   0.045   ++   0.000   {?}   {69.hb# 70.hb}   1.442   0.111   0.124   ++   0.000   {?}   {59.ca 59.ca}   52.908   0.560   0.627   ++   0.000   {?}  0.47090 0.47090 0 {?} 0
3  {59.ha}   3.863   0.048   0.068   ++   0.000   {?}   {59.HB#}   1.253   0.092   0.136   ++   0.000   {?}   {59.ca}   52.882   0.534   1.272   ?   0.000   {?}  0.93336 0.93336 0 {?} 0
4  {59.ha}   3.864   0.045   0.041   ++   0.000   {?}   {58.hg2#}   0.780   0.197   0.166   ++   0.000   {?}   {59.ca}   52.906   0.576   0.507   ++   0.000   {?}  0.37023 0.37023 0 {?} 0
5  {64.hd1#}   0.560   0.027   0.027   ?   0.000   {?}   {64.HA}   4.146   0.184   0.073   ++   0.000   {?}   {64.cd1}   10.030   0.532   0.217   ++   0.000   {?}  0.26050 0.26050 0 {?} 0
6  {64.hd1#}   0.557   0.045   0.045   ++   0.000   {?}   {64.HG12}   1.201   0.000   0.080   ?   0.000   {?}   {64.cd1}   9.900   0.417   1.272   ?   0.000   {?}  0.29701 0.29701 0 {?} 0
7  {64.HD1#}   0.555   0.029   0.044   ++   0.000   {?}   {64.HG2#}   0.801   0.042   0.098   ?   0.000   {?}   {64.CD1}   9.933   0.502   1.272   ?   0.000   {?}  0.81043 0.81043 0 {?} 0
8  {64.HD11}   0.554   0.037   0.054   ++   0.000   {?}   {72.HD1#}   0.722   0.042   0.061   ?   0.000   {?}   {64.CD1}   9.928   0.509   1.272   ?   0.000   {?}  0.83489 0.83489 0 {?} 0
9  {64.HD1#}   0.558   0.030   0.090   ++   0.000   {?}   {64.HD1#}   0.554   0.077   0.193   ++   0.000   {?}   {64.CD1}   9.937   0.491   1.272   ?   0.000   {?}  17.32395 17.32395 0 {?} 0
10  {64.hd1#}   0.558   0.048   0.048   ++   0.000   {?}   {64.HG11}   1.132   0.069   0.069   ?   0.000   {?}   {64.cd1}   9.966   0.450   1.272   ?   0.000   {?}  0.36150 0.36150 0 {?} 0
11  {53.HA}   4.529   0.008   0.008   ?   0.000   {?}   {53.HB#}   1.591   0.020   0.015   ++   0.000   {?}   {53.CA}   9.666   0.444   0.089   ++   0.000   {?}  0.23040 0.23040 0 {?} 0

        4D format:

        """
        self.fileNames.append(fileName)

        #check, if file exists:
        _DoesFileExist(fileName)

        #important: delete all the existing peaks
        self._cleanUp()

        #fileHandle:
        fileHandle = open(fileName)

        #compile some patterns:
        labelDatasetPA = re.compile('label dataset')
        openingBracketPA = re.compile('{')
        closingBracketPA = re.compile('}')
        whiteSpacePA = re.compile('\s')
        getDim = 0
        lineCounter = 0

        def __splitNMRViewLine(lineString):
            iii = 0
            withinBrackets = 0
            outLine = ''
            while 1:
                if iii > len(lineString)-1:
                    break
                character = lineString[iii]
                if character == '{':
                    withinBrackets = 1
                elif character == '}':
                    withinBrackets = 0
                    
                if withinBrackets and character in [' ', ';', ',']:
                    outLine = outLine + '_'
                else:
                    outLine = outLine + character
                iii = iii + 1
            return string.split(outLine)

        for eachLine in fileHandle.readlines():
            if lineCounter == 0:
                if not labelDatasetPA.search(eachLine):
                    print 'WARNING: first line does not start with "label dataset"'
                    print '         Is it really a NMRView file???'
                    print '         aborting.'
                    return 0
            elif lineCounter == 1:
                # get the names for all dimensions:
                dimNameList = string.split(eachLine)
                print 'NMRView file ' + fileName + ' contains ' + str(len(dimNameList)) + ' dimensions: ' + str(dimNameList)
                self.dimension = len(dimNameList)
            elif lineCounter == 2:
                # get the filename:
                nameThirdLine = string.strip(eachLine)
            elif lineCounter == 3:
                # split after '}'
                eachLine = whiteSpacePA.sub('',eachLine)
                fourthLineList = string.split(openingBracketPA.sub('',eachLine), '}')
                fourthLineList.remove('')
            elif lineCounter == 4:
                # split after '}'
                eachLine = whiteSpacePA.sub('',eachLine)
                fifthLineList = string.split(openingBracketPA.sub('',eachLine), '}')
                fifthLineList.remove('')
            elif lineCounter == 5:
                # get the short description of the data:
                descriptionList = string.split(eachLine)
                # get the column positions of the assignments, frequencies, int, vol:
                assignmentPositions = []
                frequencyPositions = []
                intPosition = None
                volPosition = None
                for eachDim in dimNameList:
                    assignmentPositions.append(descriptionList.index(eachDim + '.P'))
                    frequencyPositions.append(descriptionList.index(eachDim + '.W'))
                    intPosition = descriptionList.index('int')
                    volPosition = descriptionList.index('vol')
#                print assignmentPositions, frequencyPositions, intPosition, volPosition #test
            elif lineCounter == 6:
                self._nmrviewDimNameList = dimNameList
                self._nmrviewNameThirdLine = nameThirdLine
                self._nmrviewfourthLineList = fourthLineList
                self._nmrviewfifthLineList = fifthLineList
                self._nmrviewAssignmentPositions = assignmentPositions
                self._nmrviewFrequencyPositions = frequencyPositions
                self._nmrviewIntPosition = intPosition
                self._nmrviewVolPosition = volPosition
                self._nmrviewOriginalFileName = fileName
                self._nmrviewHet1 = het1
                self._nmrviewPro1 = pro1
                self._nmrviewHet2 = het2
                self._nmrviewPro2 = pro2
                #if you are just interested in the first six lines, stop here:
                if readOnlyHeader:
                    return

            # take care of the whitespace in {12.hn 21.hn} etc:
            lineList = __splitNMRViewLine(eachLine)
            if lineCounter < 6 or len(lineList) < 4:
                lineCounter = lineCounter + 1
                continue

##             print eachLine, lineList #test
            
            #new peak instance:
            PEAK = NoePeak()
                
            #get the data:
            assignments = []
            for eachAS in self._nmrviewAssignmentPositions:
                # extract residue numbers & atomnames (take care of "?" assignments):
                rawAssignment = lineList[eachAS]
                rawAssignment = openingBracketPA.sub('', rawAssignment)
                rawAssignment = closingBracketPA.sub('', rawAssignment)
                possibleList = string.split(rawAssignment,'_')
##                 print eachLine, possibleList #test
                outList = []
                for eachP in possibleList:
                    oneAssiList = string.split(eachP,'.')
                    if len(oneAssiList) < 2:
                        # no assignment in this case, e.g for '?' 
                        outList.append([None,None])
#                        print eachLine #test
                    else:
                        outList.append(oneAssiList)
                # check for several copies of one assignment:
                indices = []
                for kkk in range(len(outList)):
                    for lll in range(kkk+1,len(outList)):
                        if (outList[kkk] == outList[lll]):
                            indices.append(lll)
                cleanList = []
                for kkk in range(len(outList)):
                    if not kkk in indices:
                        cleanList.append(outList[kkk])
                assignments.append(cleanList)
            frequencies = []
            for eachAS in self._nmrviewFrequencyPositions:
                frequencies.append(lineList[eachAS])
            intensity = lineList[self._nmrviewIntPosition+1]
            volume = lineList[self._nmrviewVolPosition+1]
            peakNumber = lineList[0]

##             print assignments #test
##             print frequencies #test
##             print intensity #test
##             print volume #test
##             print peakNumber #test

            # for the data of one peak:
            assignmentList = []
            for xxx in range(len(self._nmrviewDimNameList)):
                assignmentList.append([assignments[xxx], frequencies[xxx]])

#            print '# ' + str(lineList) #test
#            print assignmentList #test
            

            #now get the assignments for the hets and pros:
            if het1 == 'N' or not het1:
                h1ppm = None
                atomnameh1 = (None,)
            else:
                h1ppm = assignmentList[string.atoi(het1)-1][1]
            if het2 =='N' or not het2:
                h2ppm = None
                atomnameh2 = (None,)
            else:
                h2ppm = assignmentList[string.atoi(het2)-1][1]

##             print assignmentList #test
##             print het1, pro1, het2, pro2 #test
##             print str(string.atoi(pro1)-1) #test
            
            p1ppm = assignmentList[string.atoi(pro1)-1][1]
            p2ppm = assignmentList[string.atoi(pro2)-1][1]
#            print p1ppm, p2ppm #test
#            print assignmentList #test

            #create a peak instance:
            NOE = NoePeak()

            for eachAssi1 in  assignmentList[string.atoi(pro1)-1][0]:
                for eachAssi2 in  assignmentList[string.atoi(pro2)-1][0]:
#                    print '1: ' + str(eachAssi1) + '  2: ' + str(eachAssi2) #test
                    atomnamep1 = (eachAssi1[1],)
                    atomnamep2 = (eachAssi2[1],)
                    residue1 = eachAssi1[0]
                    residue2 = eachAssi2[0]
##                     print residue1 #test
##                     print residue2 #test

                    if het1 == 'N' or not het1:
                        atomnameh1 = (None,)
                    else:
##                         print assignmentList #test
##                         print string.atoi(het1)-1 #test
##                         print assignmentList[string.atoi(het1)-1] #test
##                         print assignmentList[string.atoi(het1)-1][0] #test
                        atomnameh1 = (assignmentList[string.atoi(het1)-1][0][0][1],)
                    if het2 =='N' or not het2:
                        atomnameh2 = (None,)
                    else:
                        atomnameh2 = (assignmentList[string.atoi(het2)-1][0][0][1],)

                    #put it into contribution instance:
                    CONT = NoeContribution(spectrumName = None,\
                                           peakNumber = peakNumber,\
                                           testSet = None,\
                                           volume = volume,\
                                           volumeError = None,\
                                           intensity = intensity,\
                                           intensityError = None,\
                                           h1ppm = h1ppm,\
                                           dh1ppm = None,\
                                           p1ppm = p1ppm,\
                                           dp1ppm = None,\
                                           h2ppm = h2ppm,\
                                           dh2ppm = None,\
                                           p2ppm = p2ppm,\
                                           dp2ppm = None,\
                                           peakType = None,\
                                           fomAria = None,\
                                           lowerBound = None,\
                                           upperBound = None,\
                                           sumDistanceAve = None,\
                                           sumDistanceStd = None,\
                                           backVolumeAria = None,\
                                           backVolumeAriaStd = None,\
                                           allAssi = None,\
                                           nta = None,\
                                           contribution = None,\
                                           fomContribution = None,\
                                           distanceAve = None,\
                                           distanceStd = None,\
                                           backVolume = None,\
                                           backVolumeStd = None,\
                                           segid1 = None,\
                                           residue1 = residue1,\
                                           aa1 = residue1,\
                                           atomnameh1 = atomnameh1,\
                                           assih1ppm = None,\
                                           assidh1ppm = None,\
                                           atomnamep1 = atomnamep1,\
                                           assip1ppm = None,\
                                           assidp1ppm = None,\
                                           segid2 = None,\
                                           residue2 = residue2,\
                                           aa2 = residue2,\
                                           atomnameh2 = atomnameh2,\
                                           assih2ppm = None,\
                                           assidh2ppm = None,\
                                           atomnamep2 = atomnamep2,\
                                           assip2ppm = None,\
                                           assidp2ppm = None)

                    #add contribution to NOE instance:
                    NOE.AddContribution(CONT)

            #add peak to the NOE list:
            self.AddPeak(NOE)
            
            lineCounter = lineCounter + 1



    def ReadNMRViewPeaksAriaList(self, nmrViewFileName, ariaListFileName,\
                                 het1, pro1, het2, pro2):
        """
        use this function when reading a spectrum from merged.list or 15N.list
        when you also want to have the header information of the .xpk file
        """
        #1. read only the header information first:
        self.ReadNMRViewPeaks(nmrViewFileName,het1, pro1, het2, pro2,1)
        #2. read the ARIA .list file:
        self.ReadList(ariaListFileName)
        
        
    def ReadLolUpl(self, lolFileName, uplFileName):
        """
        reads DYANA .lol and .upl files

        ...NOT YET IMPLEMENTED!

        
        """
        #some messages:
        print 'reading the ANSIG general NOE file:'
        print '  ' + noefile

        #check if the files exist:
        _DoesFileExist(noefile)
        self.fileNames.append(noefile)
        
        #important - delete all the existing peaks:
        self._cleanUp()
        
        #read the .upl file:
        noeHandle = TextFile.TextFile(noefile)

        print 'NOT YET IMPLEMENTED!'

    def ReadSparky(self, noefile,het1,pro1,het2,pro2):
        """
        I use the same convention as in the .html file for ANSIG
        het1= '1' , '2' , '3' , '4', or 'N' (column number or not used)
        pro1= '1' , '2' , '3' , '4', or 'N' (column number or not used)
        het2= '1' , '2' , '3' , '4', or 'N' (column number or not used)
        pro2= '1' , '2' , '3' , '4', or 'N' (column number or not used)

        """
        
        #finds the dimensionality of spectrum
        DM = 0

        if het1 == 'N':
            DM = DM + 1
        if pro1 == 'N':
            DM = DM + 1
        if het2 == 'N':
            DM = DM + 1
        if pro2 == 'N':
            DM = DM + 1

        D = 4 - DM
        print 'the spectrum dimension is: ' + str(D)
        self.dimension = D
        
        #check if rhe files exist
        if  _DoesFileExist(noefile) == 0 :
            return

        self._cleanUp()

        #read the general noefile
        noehandle = open(noefile)
        assignmentPattern = re.compile('([a-zA-Z]\d+){0,1}(\S+)')

        iii = 0
        for line in noehandle.readlines():
            iii = iii + 1
            linelist = string.split(line)
##             print linelist
            #It jumps the first line and the empty one
            if len(linelist) < 4 or linelist[0] =='Assignment':
                continue

            # start analysing the assignments:
            linelist0 = string.split(linelist[0],'-')


            
            # assignmentList contains 3-letter-code,residuenumber,atomname
            assignmentList = [[None,None,None],[None,None,None],\
                             [None,None,None],[None,None,None]]
            
            # get all the assignments in the order how they appear in the file:
            if len(linelist0) > 0 and  linelist0[0] != '?':
                matched = assignmentPattern.match(linelist0[0])
                aminoacidplusnumber = matched.group(1)
                if aminoacidplusnumber:
                    aa = AminoAcid.AminoAcid(aminoacidplusnumber[0])[1]
                    residue = aminoacidplusnumber[1:]
                else:
                    aa = None
                    residue = None
                atomname = matched.group(2)
                assignmentList[0] = [aa,residue,atomname]
            if len(linelist0) > 1 and  linelist0[1] != '?':
                matched = assignmentPattern.match(linelist0[1])
                aminoacidplusnumber = matched.group(1)
                if aminoacidplusnumber:
                    aa = AminoAcid.AminoAcid(aminoacidplusnumber[0])[1]
                    residue = aminoacidplusnumber[1:]
                else:
                     aa =  assignmentList[0][0]
                     residue = assignmentList[0][1]
                atomname = matched.group(2)
                assignmentList[1] = [aa,residue,atomname]
            if len(linelist0) > 2 and  linelist0[2] != '?':
                matched = assignmentPattern.match(linelist0[2])
                aminoacidplusnumber = matched.group(1)
                if aminoacidplusnumber:
                    aa = AminoAcid.AminoAcid(aminoacidplusnumber[0])[1]
                    residue = aminoacidplusnumber[1:]
                else:
                     aa = assignmentList[1][0]
                     residue = assignmentList[1][1]
                atomname = matched.group(2)
                assignmentList[2] = [aa,residue,atomname]
            if len(linelist0) > 3 and  linelist0[3] != '?':
                matched = assignmentPattern.match(linelist0[3])
                aminoacidplusnumber = matched.group(1)
                if aminoacidplusnumber:
                    aa = AminoAcid.AminoAcid(aminoacidplusnumber[0])[1]
                    residue = aminoacidplusnumber[1:]
                else:
                     aa = assignmentList[2][0]
                     residue = assignmentList[2][1]
                atomname = matched.group(2)
                assignmentList[3] = [aa,residue,atomname]
    
            # get the right order with het1, etc....
            aa1 = None
            aa2 = None
            residue1 = None
            residue2 = None
            atomnameh1 = (None,)
            atomnameh2 = (None,)
 
            if het1 != 'N':
                if assignmentList[(string.atoi(het1)-1)][0]:
                    aa1 = assignmentList[(string.atoi(het1)-1)][0]
                if assignmentList[(string.atoi(het1)-1)][1]:
                    residue1 = assignmentList[(string.atoi(het1)-1)][1]
                atomnameh1 = (assignmentList[(string.atoi(het1)-1)][2],)
            if het2 != 'N':
                if assignmentList[(string.atoi(het2)-1)][0]:
                    aa2 = assignmentList[(string.atoi(het2)-1)][0]
                if assignmentList[(string.atoi(het2)-1)][1]:
                    residue2 = assignmentList[(string.atoi(het2)-1)][1]
                atomnameh2 = (assignmentList[(string.atoi(het2)-1)][2],)

            if assignmentList[(string.atoi(pro1)-1)][0]:
                aa1 = assignmentList[(string.atoi(pro1)-1)][0]
            if assignmentList[(string.atoi(pro1)-1)][1]:
                residue1 = assignmentList[(string.atoi(pro1)-1)][1]
            atomnamep1 = (assignmentList[(string.atoi(pro1)-1)][2],)

            if assignmentList[(string.atoi(pro2)-1)][0]:
                aa2 = assignmentList[(string.atoi(pro2)-1)][0]
            if assignmentList[(string.atoi(pro2)-1)][1]:
                residue2 = assignmentList[(string.atoi(pro2)-1)][1]
            atomnamep2 = (assignmentList[(string.atoi(pro2)-1)][2],)
            

            # for the ppms, volumes etc.:
            if het1 and het1 != 'N':
                h1ppm = linelist[string.atoi(het1)]
            else:
                h1ppm = None
            p1ppm = linelist[string.atoi(pro1)]
            if het2 and het2 != 'N':
                h2ppm = linelist[string.atoi(het2)]
            else:
                h2ppm = None
            p2ppm = linelist[string.atoi(pro2)]
            volume = linelist[-3]
            sparkylittlecomment = linelist[-2]
            intensity = linelist[-1]
            
            
            #put it in a contribution instance
            CONT = NoeContribution(peakNumber = str(iii),\
                                   residue1 = residue1,\
                                   aa1 = aa1,\
                                   atomnameh1 = atomnameh1,\
                                   atomnamep1 = atomnamep1,\
                                   residue2 = residue2,\
                                   aa2 = aa2,\
                                   atomnameh2 = atomnameh2,\
                                   atomnamep2 = atomnamep2,\
                                   volume = volume,\
                                   h1ppm = h1ppm,\
                                   p1ppm = p1ppm,\
                                   h2ppm = h2ppm,\
                                   p2ppm = p2ppm,\
                                   sparkylittlecomment = sparkylittlecomment,\
                                   intensity = intensity)
                                   
                                   
            

            #create a peak instance and add contribution
            NOE = NoePeak()
            NOE.AddContribution(CONT)

            #add peak to the NOE list
            self.AddPeak(NOE)

            
    
    def RemovePeak(self, spectrumName, peakNumber):
        """
        removes the whole peakobject for the peak with a given peakNumber
        """
        if self.peaksdic.has_key(spectrumName, peakNumber):
            #get the object from the dictionary:
            whichobject = self.peaksdic[spectrumName, peakNumber]
            #remove it from the list and the dictionary:
            self.peakslist.remove(whichobject)
            del(self.peaksdic[(spectrumName, peakNumber)])
        
    def RemoveDoubleQuotes(self):
        """
        removes double quotation marks from all the atomnames
        this is useful e.g. for XEASY RNA atomnames
        all the double quotation marks are replaced by two single
        quotation marks
        """
        doubleQuote =  re.compile('"')
        for eachP in self.peakslist:
            for eachC in eachP.contributions:
                tmpList = []
                for eachA in eachC.atomnameh1:
                    if eachA:
                        tmpList.append(doubleQuote.sub("''", eachA))
                    else:
                        tmpList.append(eachA)
                eachC.atomnameh1 = tuple(tmpList)
                tmpList = []
                for eachA in eachC.atomnameh2:
                    if eachA:
                        tmpList.append(doubleQuote.sub("''", eachA))
                    else:
                        tmpList.append(eachA)
                eachC.atomnameh2 = tuple(tmpList)
                tmpList = []
                for eachA in eachC.atomnamep1:
                    if eachA:
                        tmpList.append(doubleQuote.sub("''", eachA))
                    else:
                        tmpList.append(eachA)
                eachC.atomnamep1 = tuple(tmpList)
                tmpList = []
                for eachA in eachC.atomnamep2:
                    if eachA:
                        tmpList.append(doubleQuote.sub("''", eachA))
                    else:
                        tmpList.append(eachA)
                eachC.atomnamep2 = tuple(tmpList)

    def RestraintCombination(self):
        #author: Michele Fossi
        #this methode works on LongList:it must be applied
        #after the methode IntraSeqMediumLong.
        #it takes in a random way
        #its elements two by two, adds their volumes together,
        #creates a new element and deletes the old two ones.
        
        
        N=len(self.longList)
        
#        print  N
        RandomList=[]
        i=0
        import whrandom
        generator =whrandom.whrandom()
        while  i < N:
            
            Randomnumber=generator.randint(0,N-1)
           
            if Randomnumber in RandomList:
                pass
            else:
                RandomList.append(Randomnumber)
                i=i+1
#        print 'lenght of randomlist', len(RandomList)
        A=[]
        if len(RandomList)%2 == 0:
            #if the list is composed by an even number
            #of elements
            k=0
            
        else:
            #if the list is composed by an odd number
            #of elements
            k=1
            o=RandomList[0]
            A.append(self.longList[o])
           
        while k < len(RandomList):
            #m,n are a couple of adjacentnumbers in the random
            #list
            m=RandomList[k]
            n=RandomList[k+1]
            if self.longList[m].contributions[0].volume:
                if self.longList[n].contributions[0].volume:
                    vm = self.longList[m].contributions[0].volume
#                    print 'volume1',vm
                    vn = self.longList[n].contributions[0].volume 
#                    print 'volume 2',vn
                    vt = string.atof(vm)+string.atof(vn)
            
                    for contribution in self.longList[n].contributions:
                        #I append to the first peak all the contribution
                        #of the second one too
                        self.longList[m].contributions.append(contribution)
                        #the volume has to be changed for every contribution of
                        #the new peak
                        for CONT in self.longList[m].contributions:
                            vt = str(vt)
                            CONT.volume=vt
#                            print 'the new volume',vt
                         
                                               

                       
                        #I append the new peak in A, the list of new peaks                      
                        A.append(self.longList[m])    
                        k=k+2
                        
                        if self.longList[m].contributions[0].intensity:
                            if self.longList[n].contributions[0].intensity:
                                intm = self.longList[m].contributions[0].intensity
#                                print 'int 1',intm
                                intn = self.longList[n].contributions[0].intensity
#                                print 'int 2',intn
                                it = string.atof(intm)+string.atof(intn)
                                for CONT in self.longList[m].contributions:
                                    it = str(it)
                                    CONT.intensity=it
#                                    print 'the new intensity', it
                else:
                    k=k+2
                    continue
            else:
#                print 'the volume in longList[',m,'],contributions[0] doesn t exist!'
                k=k+2
                continue

        #the new peakslist contains only half of the old elements,
        #each of them with a new volume which is the sum of the volumes
        #of two old peaks
#        print ' len(A)',len(A)
        
        self.longList=A
        


        self.peakslist=self.intraList +\
                       self.seqList+\
                       self.mediumList+\
                       self.longList+\
                       self.unassignedList+\
                       self.ambiguousList
#        print len(self.peakslist)
        return self.peakslist                                                                                                        

    
    

    def SetSpectrumName(self, spectrumName):
        """sets the spectrum name for all contributions"""
        for eachP in self.peakslist:
            for eachC in eachP.contributions:
                eachC.spectrumName = spectrumName
                
    def Sort(self):
        """
        sorts all the data:   1. after spectrumName
                              2. after peakNumber
        """
        allLists = ['peakslist',
                    'intraList',
                    'seqList',
                    'mediumList',
                    'longList',
                    'unassignedList',
                    'ambiguousList']
        for eachList in allLists:
            getattr(self, eachList).sort(Comparisons.CmpComposite(Comparisons.CmpAttr('spectrumName', 0),\
                                                                  Comparisons.CmpAttr('peakNumber', 1)))

        
    def Stdout(self):
        """writes the attributes to stdout in the new .list format"""
        print 'name:', self.name, '\ndimension:', self.dimension,\
              '\ncomment:', self.comment
        print 'p spec pnum test WIn volu volE int intE h1ppm dh1ppm p1ppm dp1ppm h2ppm dh2ppm p2ppm dp2ppm'
        print 'a type fom Wcur low up SDA SDS bVol bVolS all nta'
        print 'c con fom dAve dStd bVol bVolS seg1 res1 aa1 Nh1 ah1ppm adh1ppm Np1 ap1ppm adp1ppm seg2 res2 aa2 Nh2 ah2ppm adh2ppm Np2 ap2ppm adp2ppm'
        for PEAK in self.peakslist:
            forTheFirstCon = 1
            for CON in PEAK.contributions:
                if forTheFirstCon:
                    #the 'p' line:
                    print 'p',
                    for each in (CON.spectrumName, CON.peakNumber,\
                                 CON.testSet, CON.inWeight,\
                                 CON.volume, CON.volumeError,\
                                 CON.intensity, CON.intensityError,\
                                 CON.h1ppm,  CON.dh1ppm,\
                                 CON.p1ppm, CON.dp1ppm,\
                                 CON.h2ppm, CON.dh2ppm,\
                                 CON.p2ppm, CON.dp2ppm):
                        if each == None:
                            print '-',
                        elif each == (None,):
                            print ('-',),
                        else:
                            print each,
                    print ''
                    #the 'a' line:
                    print 'a',
                    for each in (CON.peakType, CON.fomAria,\
                                 CON.curWeight, CON.lowerBound,\
                                 CON.upperBound, CON.sumDistanceAve,\
                                 CON.sumDistanceStd, CON.backVolume,\
                                 CON.backVolumeStd, CON.allAssi,\
                                 CON.nta):
                        if each == None:
                            print '-',
                        elif each == (None,):
                            print ('-',),
                        else:
                            print each,
                    print ''
                #the 'c' line:
                forTheFirstCon = 0
                print 'c',
                for each in (CON.contribution, CON.fomContribution,\
                             CON.distanceAve, CON.distanceStd,\
                             CON.backVolume, CON.backVolumeStd,\
                             CON.segid1, CON.residue1,\
                             CON.aa1, CON.atomnameh1,\
                             CON.assih1ppm, CON.assidh1ppm,\
                             CON.atomnamep1, CON.assip1ppm,\
                             CON.assidp1ppm, CON.segid2,\
                             CON.residue2, CON.aa2,\
                             CON.atomnameh2, CON.assih2ppm,\
                             CON.assidh2ppm, CON.atomnamep2,\
                             CON.assip2ppm, CON.assidp2ppm):
##                              CON.atomnumberh1, CON.atomnumberp1,\
##                              CON.atomnumberh2, CON.atomnumberp2):
                    if each == None:
                        print '-',
                    elif each == (None,):
                        print ('-',),
                    else:
                        print each,
                print ''


    def WriteAqua(self, fileName):
        """
        a method for writing AQUA restraint file
        
The restraint files are line-oriented, giving one restraint per line.
Every restraint line starts with a keyword identifier.
Lines beginning with # are ignored.
In the Aqua routines empty or blank lines force end of reading, and
lines starting with an unknown keyword generate an error. 

1.The first line is a title header (not interpreted). 
2.The second line is: 
  count NUM type FORMAT MR
  where 
  NUM is integer value equal to the number of valid restraint lines in file 
  FORMAT is string indicating restraint format (DISGEO, X-PLOR, BIOSYM etc.;
  is used to specify AtomLIB, see Names)
  MR is string with value mr (MR type file), or empty string 
3.The third and following lines contain the restraint information in a format
  that is described below. 

A restraint line contains: 
keyword  residue_name  residue_number  atom_name
         ...residue_name  residue_number  atom_name  bound_1  [ bound_2 ]

The keyword can be: 

   keyword          data                        BOUND_1        BOUND_2
   ======================================================================
    NOEUPP      NOE restraints                 upper bound        n.a.
    NOELOW      NOE restraints                 lower bound        n.a.
    NOEUPLO     NOE restraints                 upper bound    lower bound
    HBUPP       H-bond restraints              upper bound        n.a.
    HBLOW       H-bond restraints              lower bound        n.a.
    HBUPLO      H-bond restraints              upper bound    lower bound
    SSUPP       Disulphide restraints          upper bound        n.a.
    SSLOW       Disulphide restraints          lower bound        n.a.
    SSUPLO      Disulphide restraints          upper bound    lower bound
    DISUPP      Generic distance restraints    upper bound        n.a.
    DISLOW      Generic distance restraints    lower bound        n.a.
    DISUPLO     Generic distance restraints    upper bound    lower bound

The C-format used to write the file in the conversion scripts is: 

       %-7s %-4s %3i %-5s %-4s %3i %-5s %10.3f %10.3f\n

The Fortran equivalent of this should be (didn't test it): 

       (A7,1X,A4,1X,I3,1X,A5,1X,A4,1X,I3,1X,A5,1X,F10.3,1X,F10.3) 

Note: the format should be treated as field- rather than column-oriented,
since chain identifiers and residue insertion codes lead to extra fields
being inserted.
The extended format is: 
keyword  [CHAIN id]  residue_name  residue_number  [INSERT code]  atom_name
         ...[CHAIN id]  residue_name  residue_number  [INSERT code]  atom_name
         ...bound_1  [ bound_2 ]
        """
        noeHandle = open(fileName, 'w')
        print 'writing  an AQUA general NOE file:\n  ', fileName
        if len(self.peakslist) == 0:
            print 'peaklist empty. WriteAqua method aborted.'
            return
        noeHandle.write('ARIA generated AQUA NOE restraint file\n')
        secondLineString = "count %i type %s \n" %(len(self.peakslist), 'DISGEO')
        noeHandle.write(secondLineString)
        for PEAK in self.peakslist:
            CON = PEAK.contributions[0]
            #initializing:
            outAa1 = ''
            outResidue1 = ''
            outAtomnamep1 = ''
            outAa2 = ''
            outResidue2 = ''
            outAtomnamep2 = ''
            outUpperBound = ''
            outLowerBound = ''

            if CON.aa1 and len(CON.aa1) == 3:
                outAa1 = CON.aa1
            if CON.residue1:
                outResidue1 = CON.residue1
            if CON.atomnamep1 and len(CON.atomnamep1) > 0:
                if CON.atomnamep1[0][-1] == '%':
                    outAtomnamep1 = PseudoAtom.Atom2Pseudo(CON.atomnamep1[0],CON.aa1)[0]
                else:
                    outAtomnamep1 = CON.atomnamep1[0]
            if CON.aa2 and len(CON.aa2) == 3:
                outAa2 = CON.aa2
            if CON.residue2:
                outResidue2 = CON.residue2
            if CON.atomnamep2 and len(CON.atomnamep2) > 0:
                if CON.atomnamep2[0][-1] == '%':
                    outAtomnamep2 = PseudoAtom.Atom2Pseudo(CON.atomnamep2[0],CON.aa2)[0]
                else:
                    outAtomnamep2 = CON.atomnamep2[0]
            if CON.upperBound:
                outUpperBound = CON.upperBound
            if CON.lowerBound:
                outLowerBound = CON.lowerBound

            lineString = """%-7s %-4s %3s %-5s %-4s %3s %-5s %10.3f %10.3f\n""" %\
                         ('NOEUPLO', outAa1, outResidue1, outAtomnamep1,\
                          outAa2, outResidue2, outAtomnamep2,\
                          string.atof(outUpperBound), string.atof(outLowerBound))
            noeHandle.write(lineString)
                
    
    def WriteAnsigNoe(self, fileName, het1, pro1, het2, pro2):
        """
        writes an ANSIG general NOE file
        
        I use the same convention as in the .html file for ANSIG:
        het1 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro1 = '1', '2', '3', '4' or 'N' (column number or not used)
        het2 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro2 = '1', '2', '3', '4' or 'N' (column number or not used)):
        """
        noeHandle = TextFile.TextFile(fileName, 'w')
        print 'writing  an ANSIG general NOE file:\n  ', fileName
        
        if len(self.peakslist) == 0:
            print 'peaklist empty. WriteAnsigNoe method aborted.'
            return
        
        #the first line contains only comments:
        noeHandle.write('! ANSIG v3.3 general format NOE file\n')

        #the output format as explained in the ANSIG manual:
        format2 = FortranFormat.FortranFormat('I1,A1,3A4,F8.3,3A4,F8.3,I6,A12,F6.1,2E11.4')
        format3 = FortranFormat.FortranFormat('I1,A1,3A4,F8.3,3A4,F8.3,3A4,F8.3,I6,A12,F6.1,2E11.4')
        format4 = FortranFormat.FortranFormat('I1,A1,3A4,F8.3,3A4,F8.3,3A4,F8.3,3A4,F8.3,I6,A12,F6.1,2E11.4')

        #order the dimensions with a list of tuples called outlist:
        outlist = range(5)  #zero-element is not used
        
        for PEAK in self.peakslist:
            for CON in PEAK.contributions:
                #convert 3-lettercode from XXX to Xxx and create outlist:
                if het1 != 'N':
                    outlist[string.atoi(het1)] = (CON.residue1, CON.aa1[:1] + \
                                                  string.lowerBound(CON.aa1[1:3]), \
                                                  CON.atomnameh1, CON.h1ppm)
                outlist[string.atoi(pro1)] = (CON.residue1, CON.aa1[:1] + \
                                              string.lowerBound(CON.aa1[1:3]), \
                                              CON.atomnamep1, CON.p1ppm)
                if het2 != 'N':
                    outlist[string.atoi(het2)] = (CON.residue2, CON.aa2[:1] + \
                                                  string.lowerBound(CON.aa2[1:3]), \
                                                  CON.atomnameh2, CON.h2ppm)
                outlist[string.atoi(pro2)] = (CON.residue2, CON.aa2[:1] + \
                                              string.lowerBound(CON.aa2[1:3]), \
                                              CON.atomnamep2, CON.p2ppm)

                #write the stuff with FortranFormat:
                #set some default values if field is None:
                if CON.peakNumber == None:
                    outpeakNumber = 0
                else:
                    outpeakNumber = string.atoi(CON.peakNumber)
                if CON.mixtime == None:
                    outmixtime = 0.0
                else:
                    outmixtime = string.atof(CON.mixtime)
                if CON.volume == None:
                    outvolume = 0.0
                else:
                    outvolume = string.atof(CON.volume)
                if CON.relativeint == None:
                    outrelativeint = 0.0
                else:
                    outrelativeint = string.atof(CON.relativeint)
                
                if self.dimension == 2:
                    #use attributes spectrumName, mixtime, relative,
                    #take the first atomname from the tuple:
                    line = FortranFormat.FortranLine([self.dimension, ' ',\
                                                      outlist[1][0], outlist[1][1],\
                                                      outlist[1][2][0], string.atof(outlist[1][3]),\
                                                      outlist[2][0], outlist[2][1],\
                                                      outlist[2][2][0], string.atof(outlist[2][3]),\
                                                      outpeakNumber, CON.spectrumName,\
                                                      outmixtime, outvolume,\
                                                      outrelativeint], format2)
                elif self.dimension == 3:
                    #use attributes spectrumName, mixtime, relative,
                    #take the first atomname from the tuple:
                    line = FortranFormat.FortranLine([self.dimension, ' ',\
                                                      outlist[1][0], outlist[1][1],\
                                                      outlist[1][2][0], string.atof(outlist[1][3]),\
                                                      outlist[2][0], outlist[2][1],\
                                                      outlist[2][2][0], string.atof(outlist[2][3]),\
                                                      outlist[3][0], outlist[3][1],\
                                                      outlist[3][2][0], string.atof(outlist[3][3]),\
                                                      outpeakNumber, CON.spectrumName,\
                                                      outmixtime, outvolume,\
                                                      outrelativeint], format3)
                elif self.dimension == 4:
                    #use attributes spectrumName, mixtime, relative,
                    #take the first atomname from the tuple:
                    line = FortranFormat.FortranLine([self.dimension, ' ',\
                                                      outlist[1][0], outlist[1][1],\
                                                      outlist[1][2][0], string.atof(outlist[1][3]),\
                                                      outlist[2][0], outlist[2][1],\
                                                      outlist[2][2][0], string.atof(outlist[2][3]),\
                                                      outlist[3][0], outlist[3][1],\
                                                      outlist[3][2][0], string.atof(outlist[3][3]),\
                                                      outlist[4][0], outlist[4][1],\
                                                      outlist[4][2][0], string.atof(outlist[4][3]),\
                                                      outpeakNumber, CON.spectrumName,\
                                                      outmixtime, outvolume,\
                                                      outrelativeint], format4)
                #write the stuff to disk:
                noeHandle.write(str(line) +'\n')
        noeHandle.close()


    def WriteAmber(self, pdbFN, outputFN):
        """
        writes a list in AMBER distance restraint format
        """
        file = TextFile.TextFile(outputFN, 'w')
        print 'writing the data in AMBER format to', outputFN
        if len(self.peakslist) == 0:
            print 'peaklist empty. Method aborted.'
            return
        # import the MMTK/Scientific PDB stuff:
        from Scientific.IO import PDB
        conf = PDB.Structure(pdbFN)
#        for residue in conf.residues:
#            for atom in residue:
#                print atom.name, residue.name
#        return atom,residue
        bigString = ''

###############################################################################
        atomnameDic = {"H5'": "H5'1",\
                       "H5''": "H5'2",\
                       "O2'": "HO'2",\
                       "H2'": "H2'1"}
##         residueDic = {"GUA": "RG",\
##                       "ADE": "RA",\
##                       "CYT": "RC",\
##                       "URI": "RU"}
##         # 5' end of RNA is called RG5, RA5, RC5, RU5
##         # 3' end of RNA is called RG3, RA3, etc.

###############################################################################        
        
        for PEAK in self.peakslist:
            if len(PEAK.contributions) < 1:
                continue
            if len(PEAK.contributions) > 1:
                print 'CAUTION: ' + PEAK.contributions[0].peakNumber + ' has ' +\
                      len(PEAK.contributions) + ' contributions -> not used for AMBER restraints!'
                continue
            
            # get the peak numbers from the PDB file:
            lookupAtomname1 = PEAK.contributions[0].atomnamep1[0]
            lookupAtomname2 = PEAK.contributions[0].atomnamep2[0]
            if atomnameDic.has_key(lookupAtomname1):
                lookupAtomname1 = atomnameDic[lookupAtomname1]
            if atomnameDic.has_key(lookupAtomname2):
                lookupAtomname2 = atomnameDic[lookupAtomname2]
            iat1 = str(conf.residues[string.atoi(PEAK.contributions[0].residue1) - 1].atoms[lookupAtomname1].properties['serial_number'])
#            except:
#                print 'WARNING: could not find PDB atom number for: ' + PEAK.contributions[0].residue1 + ' ' + PEAK.contributions[0].atomnamep1[0]
#                continue
            iat2 = str(conf.residues[string.atoi(PEAK.contributions[0].residue2) - 1].atoms[lookupAtomname2].properties['serial_number'])
#            except:
#                print 'WARNING: could not find PDB atom number for: ' + PEAK.contributions[0].residue2 + ' ' + PEAK.contributions[0].atomnamep2[0]
#                continue
                
            # create the AMBER potential:
            lowerLimit = string.atof(PEAK.contributions[0].lowerBound) - 0.5
            upperLimit = string.atof(PEAK.contributions[0].upperBound) + 0.5
            if lowerLimit > upperLimit:
                average = (lowerLimit + upperLimit)/2.0
                lowerLimit = average
                upperLimit = average
            r1 = str(lowerLimit)
            r2 = PEAK.contributions[0].lowerBound
            r3 = PEAK.contributions[0].upperBound
            r4 = str(upperLimit)

            #some default values for AMBER restraints:
            rk2 = '1'
            rk3 = '1'
            ir6 = '1'
            
            restraintString = """ &rst iat= %s, %s, r3= %s, r4= %s,
       r1= %s, r2= %s, rk2= %s, rk3= %s, ir6= %s, &end
""" % (iat1, iat2, r3, r4, r1, r2, rk2, rk3, ir6)
            file.write(restraintString)
        file.close()
        

    def WriteAureliaUserPeak(self, fileName):
        """
        writes a list in AURELIA user peak list format
        """
        file = TextFile.TextFile(fileName, 'w')
        print 'writing the data in AURELIA user peak list format to', fileName
        print 'We always use the following format for the user peak lists:'
        print '  # 110.23 4.56 15.8 2.05 N HN CA HC 3 76 ABCD EFGH'
        print '  # h1ppm p1ppm h2ppm p2ppm h1 p1 h2 p2 res1 res2 segid1 segid2'
        print '  segid should contain 4 letters or 4 spaces'
        print '  other formats can not be read in by this method!'
        if len(self.peakslist) == 0:
            print 'peaklist empty. WriteAureliaUserPeak method aborted.'
            return
        for PEAK in self.peakslist:
            for CON in PEAK.contributions:
                file.write('# ' + CON.h1ppm + ' ' + CON.p1ppm + ' ' +\
                           CON.h2ppm + ' ' + CON.p2ppm + ' ' +\
                           CON.atomnameh1[0] + ' ' + CON.atomnamep1[0] + ' ' +\
                           CON.atomnameh2[0] + ' ' + CON.atomnamep2[0] + ' ' +\
                           CON.residue1 + ' ' + CON.residue2 +\
                           CON.segid1 + ' ' + CON.segid2 + '\n')
        file.close()
        

    def WriteList(self, fileName):
        """
        writes a list in ARIA's .list format
        """
        file = TextFile.TextFile(fileName, 'w')
        print 'writing the data in .list format to', fileName
        if len(self.peakslist) == 0:
            print 'peaklist empty. WriteList method aborted.'
            return
        for PEAK in self.peakslist:
            if len(PEAK.contributions) == 0:
                continue
            #peak line:
            if PEAK.spectrumName:
                outspectrumName = PEAK.contributions[0].spectrumName
            else:
                outspectrumName = '-'
            if PEAK.contributions[0].peakNumber:
                outpeakNumber = PEAK.contributions[0].peakNumber
            else:
                outpeakNumber = '-'
            if PEAK.contributions[0].testSet:
                outtestSet = PEAK.contributions[0].testSet
            else:
                outtestSet = '-'
            if PEAK.contributions[0].inWeight:
                outinWeight = PEAK.contributions[0].inWeight
            else:
                outinWeight = '-'
            if PEAK.contributions[0].volume:
                outvolume = PEAK.contributions[0].volume
            else:
                outvolume = '-'
            if PEAK.contributions[0].volumeError:
                outvolumeError = PEAK.contributions[0].volumeError
            else:
                outvolumeError = '-'
            if PEAK.contributions[0].intensity:
                outintensity = PEAK.contributions[0].intensity
            else:
                outintensity = '-'
            if PEAK.contributions[0].intensityError:
                outintensityError = PEAK.contributions[0].intensityError
            else:
                outintensityError = '-'
            if PEAK.contributions[0].h1ppm:
                outh1ppm = PEAK.contributions[0].h1ppm
            else:
                outh1ppm = '-'
            if PEAK.contributions[0].dh1ppm:
                outdh1ppm = PEAK.contributions[0].dh1ppm
            else:
                outdh1ppm = '-'
            if PEAK.contributions[0].p1ppm:
                outp1ppm = PEAK.contributions[0].p1ppm
            else:
                outp1ppm = '-'
            if PEAK.contributions[0].dp1ppm:
                outdp1ppm = PEAK.contributions[0].dp1ppm
            else:
                outdp1ppm = '-'
            if PEAK.contributions[0].h2ppm:
                outh2ppm = PEAK.contributions[0].h2ppm
            else:
                outh2ppm = '-'
            if PEAK.contributions[0].dh2ppm:
                outdh2ppm = PEAK.contributions[0].dh2ppm
            else:
                outdh2ppm = '-'
            if PEAK.contributions[0].p2ppm:
                outp2ppm = PEAK.contributions[0].p2ppm
            else:
                outp2ppm = '-'
            if PEAK.contributions[0].dp2ppm:
                outdp2ppm = PEAK.contributions[0].dp2ppm
            else:
                outdp2ppm = '-'
            
            outString = 'p ' +\
            outspectrumName + ' ' +\
            outpeakNumber + ' ' +\
            outtestSet + ' ' +\
            outinWeight + ' ' +\
            outvolume + ' ' +\
            outvolumeError + ' ' +\
            outintensity + ' ' +\
            outintensityError + ' ' +\
            outh1ppm + ' ' +\
            outdh1ppm + ' ' +\
            outp1ppm + ' ' +\
            outdp1ppm + ' ' +\
            outh2ppm + ' ' +\
            outdh2ppm + ' ' +\
            outp2ppm + ' ' +\
            outdp2ppm + ' \n'
            
            file.write(outString)

            iii = 0
            for CON in PEAK.contributions:
                if iii == 0:
                    if CON.peakType:
                        outpeakType = CON.peakType
                    else:
                        outpeakType = '-'
                    if CON.fomAria:
                        outfomAria = CON.fomAria
                    else:
                        outfomAria = '-'
                    if CON.curWeight:
                        outcurWeight = CON.curWeight
                    else:
                        outcurWeight = '-'
                    if CON.lowerBound:
                        outlowerBound = CON.lowerBound
                    else:
                        outlowerBound = '-'
                    if CON.upperBound:
                        outupperBound = CON.upperBound
                    else:
                        outupperBound = '-'
                    if CON.sumDistanceAve:
                        outsumDistanceAve = CON.sumDistanceAve
                    else:
                        outsumDistanceAve = '-'
                    if CON.sumDistanceStd:
                        outsumDistanceStd = CON.sumDistanceStd
                    else:
                        outsumDistanceStd = '-'
                    if CON.backVolumeAria:
                        outbackVolumeAria = CON.backVolumeAria
                    else:
                        outbackVolumeAria = '-'
                    if CON.backVolumeAriaStd:
                        outbackVolumeAriaStd = CON.backVolumeAriaStd
                    else:
                        outbackVolumeAriaStd = '-'
                    if CON.allAssi:
                        outallAssi = CON.allAssi
                    else:
                        outallAssi = '-'
                    if CON.nta:
                        outnta = CON.nta
                    else:
                        outnta = '-'

                    outString = 'a ' +\
                    outpeakType + ' ' +\
                    outfomAria + ' ' +\
                    outcurWeight + ' ' +\
                    outlowerBound + ' ' +\
                    outupperBound + ' ' +\
                    outsumDistanceAve + ' ' +\
                    outsumDistanceStd + ' ' +\
                    outbackVolumeAria + ' ' +\
                    outbackVolumeAriaStd + ' ' +\
                    outallAssi + ' ' +\
                    outnta + ' \n'
                    file.write(outString)
                iii = iii + 1


                if CON.contribution:
                    outcontribution = CON.contribution
                else:
                    outcontribution = '-'
                if CON.fomContribution:
                    outfomContribution = CON.fomContribution
                else:
                    outfomContribution = '-'
                if CON.distanceAve:
                    outdistanceAve = CON.distanceAve
                else:
                    outdistanceAve = '-'
                if CON.distanceStd:
                    outdistanceStd = CON.distanceStd
                else:
                    outdistanceStd = '-'
                if CON.backVolume:
                    outbackVolume = CON.backVolume
                else:
                    outbackVolume = '-'
                if CON.backVolumeStd:
                    outbackVolumeStd = CON.backVolumeStd
                else:
                    outbackVolumeStd = '-'
                if CON.segid1:
                    outsegid1 = CON.segid1
                else:
                    outsegid1 = '-'
                if CON.residue1:
                    outresidue1 = CON.residue1
                else:
                    outresidue1 = '-'
                if CON.aa1:
                    outaa1 = CON.aa1
                else:
                    outaa1 = '-'
                if len(CON.atomnameh1) > 0:
                    if CON.atomnameh1[0]:
                        outatomnameh1 = CON.atomnameh1[0]
                    else:
                        outatomnameh1 = '-'
                else:
                    outatomnameh1 = '-'
                if CON.assih1ppm:
                    outassih1ppm = CON.assih1ppm
                else:
                    outassih1ppm = '-'
                if CON.assidh1ppm:
                    outassidh1ppm = CON.assidh1ppm
                else:
                    outassidh1ppm = '-'
                if len(CON.atomnamep1) > 0:
                    if CON.atomnamep1[0]:
                        outatomnamep1 = CON.atomnamep1[0]
                    else:
                        outatomnamep1 = '-'
                else:
                    outatomnamep1 = '-'
                if CON.assip1ppm:
                    outassip1ppm = CON.assip1ppm
                else:
                    outassip1ppm = '-'
                if CON.assidp1ppm:
                    outassidp1ppm = CON.assidp1ppm
                else:
                    outassidp1ppm = '-'
                if CON.segid2:
                    outsegid2 = CON.segid2
                else:
                    outsegid2 = '-'
                if CON.residue2:
                    outresidue2 = CON.residue2
                else:
                    outresidue2 = '-'
                if CON.aa2:
                    outaa2 = CON.aa2
                else:
                    outaa2 = '-'
                if len(CON.atomnameh2) > 0:
                    if CON.atomnameh2[0]:
                        outatomnameh2 = CON.atomnameh2[0]
                    else:
                        outatomnameh2 = '-'
                else:
                    outatomnameh2 = '-'
                if CON.assih2ppm:
                    outassih2ppm = CON.assih2ppm
                else:
                    outassih2ppm = '-'
                if CON.assidh2ppm:
                    outassidh2ppm = CON.assidh2ppm
                else:
                    outassidh2ppm = '-'
                if len(CON.atomnamep2) > 0:
                    if CON.atomnamep2[0]:
                        outatomnamep2 = CON.atomnamep2[0]
                    else:
                        outatomnamep2 = '-'
                else:
                    outatomnamep2 = '-'
                if CON.assip2ppm:
                    outassip2ppm = CON.assip2ppm
                else:
                    outassip2ppm = '-'
                if CON.assidp2ppm:
                    outassidp2ppm = CON.assidp2ppm
                else:
                    outassidp2ppm = '-'


                outString = 'c ' +\
                outfomContribution + ' ' +\
                outcontribution + ' ' +\
                outdistanceAve + ' ' +\
                outdistanceStd + ' ' +\
                outbackVolume + ' ' +\
                outbackVolumeStd + ' ' +\
                outsegid1 + ' ' +\
                outresidue1 + ' ' +\
                outaa1 + ' ' +\
                outatomnameh1 + ' ' +\
                outassih1ppm + ' ' +\
                outassidh1ppm + ' ' +\
                outatomnamep1 + ' ' +\
                outassip1ppm + ' ' +\
                outassidp1ppm + ' ' +\
                outsegid2 + ' ' +\
                outresidue2 + ' ' +\
                outaa2 + ' ' +\
                outatomnameh2 + ' ' +\
                outassih2ppm + ' ' +\
                outassidh2ppm + ' ' +\
                outatomnamep2 + ' ' +\
                outassip2ppm + ' ' +\
                outassidp2ppm + ' \n'
                file.write(outString)

        file.close()


    def WriteTbl(self, fileName, ppmdhet1='0.5', ppmdpro1='0.03',\
                 ppmdhet2='0.5', ppmdpro2='0.03', whichpeaks='all'):
        """
        writes a list in ARIA's .tbl format

        you have to specify the deltas for the frequency windows:
        ppmdhet1
        ppmdpro1
        ppmdhet2
        ppmdpro2
        The given values are just default values. It's advisable
        to find the optimal values for each dataset.
        
        
        whichpeaks offers the choices:
        all         all peaks
        allw1       all peaks, set weight > 1
        allassiw1   all peaks, set weight > 1 for fully assigned peaks 
        onlyassi    only fully assigned peaks
        onlyassiw1  only fully assigned peaks, set weight > 1

        The 'hydrogen' statement is pretty useful for some NMR experiments involving 15N or 13C.

        this method looks if the chemical shifts are existent.
        then the assi statement is written out
        """
        print 'writing the .tbl file', fileName
        tblhandle = TextFile.TextFile(fileName, 'w')
        tblhandle.write('! derived from the ' + str(self.dimension) +\
                        'D file(s):\n')
        for eachfile in self.fileNames:
            tblhandle.write('! ' + eachfile + '\n')
        for EACHP in self.peakslist:
            contCounter = 0
            for EACHC in EACHP.contributions:
                #set lowerBound and upper bounds:
                if EACHC.h1ppm != None:
                    lowerBoundh1 = str(float(EACHC.h1ppm) - float(ppmdhet1))
                    higherh1 = str(float(EACHC.h1ppm) + float(ppmdhet1))
                if EACHC.p1ppm != None:
                    lowerBoundp1 = str(float(EACHC.p1ppm) - float(ppmdpro1))
                    higherp1 = str(float(EACHC.p1ppm) + float(ppmdpro1))
                if EACHC.h2ppm != None:   
                     lowerBoundh2 = str(float(EACHC.h2ppm) - float(ppmdhet2))
                     higherh2 = str(float(EACHC.h2ppm) + float(ppmdhet2))
                if EACHC.p2ppm != None:   
                    lowerBoundp2 = str(float(EACHC.p2ppm) - float(ppmdpro2))
                    higherp2 = str(float(EACHC.p2ppm) + float(ppmdpro2))

                #for the volumes:
                if EACHC.volume == None:
                    outvolume = '0.000'
                else:
                    outvolume = EACHC.volume

                #for the spectrumName:
                if EACHC.spectrumName:
                    if len(EACHC.spectrumName) > 4:
#                        print 'WARNING: spectrum name contains more than 4 characters!'
#                        print '         only use first four characters...'
                        useSN = EACHC.spectrumName[:4]
                    else:
                        useSN = EACHC.spectrumName
                    outSpectrumName = ' spectrum ' + useSN
                else:
                    outSpectrumName = ''
                
                #distance is always set to 6 by default:
                dist = '6.0'
                if EACHC.volume == 0.0:
                    errp = '0.0'
                    errm = '6.0'
                else:
                    errp = '0.1'
                    errm = '0.1'
                #set the segids (with keyword, braces and 'and') if they exist:
                segid1 = ''
                segid2 = ''
                outsegid1 = ''
                outsegid2 = ''
                if EACHC.segid1 != None:
                    segid1 = ' and segid "' + EACHC.segid1 + '"'
                if EACHC.segid2 != None:
                    segid2 = ' and segid "' + EACHC.segid2 + '"'

                #either 'ASSI' or 'OR':
                assiOr = 'ASSI'
                if contCounter > 0:
                    assiOr = '  OR'

                #for empty residue numbers:
                if EACHC.residue1 == None:
                    outresidue1 = '-'
                else:
                    outresidue1 = EACHC.residue1
                if EACHC.residue2 == None:
                    outresidue2 = '-'
                else:
                    outresidue2 = EACHC.residue2
                
                #check all the cases for 2D, 3D and 4D:
                #4D case:
                if (EACHC.h1ppm != None) and\
                   (EACHC.p1ppm != None) and\
                   (EACHC.h2ppm != None) and\
                   (EACHC.p2ppm != None):
                    #look if it's fully assigned:
                    if EACHC.atomnameh1 != (None,) and\
                       EACHC.atomnamep1 != (None,) and\
                       EACHC.atomnameh2 != (None,) and\
                       EACHC.atomnamep2 != (None,):
                        #for the fully assigned peaks:
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassi':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassiw1':
                            weight = ' weight 1.1'                            
                        #write all the stuff to disk (assigned format):
                        tblhandle.write(assiOr + ' ( resid ' + outresidue1 + \
                                        ' and name ' + EACHC.atomnamep1[0] + \
                                        segid1 + ' )\n')
                        tblhandle.write('     ( resid ' + outresidue2 + \
                                        ' and name ' + EACHC.atomnamep2[0] + \
                                        segid2 + ' )\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' + \
                                            errp + ' peak ' + \
                                            EACHC.peakNumber + \
                                            outSpectrumName + weight + ' volume ' + \
                                            outvolume + \
                                            ' hpm1 ' + EACHC.h1ppm + \
                                            ' ppm1 ' + EACHC.p1ppm + \
                                            ' hpm2 ' + EACHC.h2ppm + \
                                            ' ppm2 ' + EACHC.p2ppm + '\n')
                    else:
                        #for the other peaks (not fully assigned):
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ''
                        elif whichpeaks == 'onlyassi':
                            continue  #don't use this contribution!
                        elif whichpeaks == 'onlyassiw1':
                            continue  #don't use this contribution!
                        #write all the stuff to disk (frequency window format):
                        tblhandle.write(assiOr + ' (( attr store1 < ' + higherp1 + \
                                        ' and attr store1 > ' + lowerBoundp1 + ' )\n')
                        tblhandle.write('     and bondedto ( attr store1 < ' + \
                                        higherh1 + ' and attr store1 > ' + \
                                        lowerBoundh1 + ' ) ' + segid1 + ' )\n')
                        tblhandle.write('     (( attr store1 < ' + higherp2 + \
                                        ' and attr store1 > ' + lowerBoundp2 + ' )\n')
                        tblhandle.write('     and bondedto ( attr store1 < ' + \
                                        higherh2 + ' and attr store1 > ' + \
                                        lowerBoundh2 + ' ) ' + segid2 + ' )\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume +\
                                            ' hpm1 ' + EACHC.h1ppm +\
                                            ' ppm1 ' + EACHC.p1ppm +\
                                            ' hpm2 ' + EACHC.h2ppm +\
                                            ' ppm2 ' + EACHC.p2ppm + '\n')
                
                #3D case, heteronucleus in het1:
                elif (EACHC.h1ppm != None) and\
                     (EACHC.p1ppm != None) and\
                     (EACHC.h2ppm == None) and\
                     (EACHC.p2ppm != None):
                    #look if it's fully assigned:
                    if EACHC.atomnameh1 != (None,) and\
                       EACHC.atomnamep1 != (None,) and\
                       EACHC.atomnamep2 != (None,):
                        #for the fully assigned peaks:
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassi':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassiw1':
                            weight = ' weight 1.1'
                        #write all the stuff to disk (assigned format):
                        tblhandle.write(assiOr + ' ( resid ' + outresidue1 + \
                                        ' and name ' + EACHC.atomnamep1[0] + \
                                        segid1 + ' )\n')
                        tblhandle.write('     ( resid ' + outresidue2 + \
                                        ' and name ' + EACHC.atomnamep2[0] + \
                                        segid2 + ' )\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume +\
                                            ' hpm1 ' + EACHC.h1ppm +\
                                            ' ppm1 ' + EACHC.p1ppm +\
                                            ' ppm2 ' + EACHC.p2ppm + '\n')
                    else:
                        #for the other peaks (not fully assigned):
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ''
                        elif whichpeaks == 'onlyassi':
                            continue  #don't use this contribution!
                        elif whichpeaks == 'onlyassiw1':
                            continue  #don't use this contribution!
                        #write all the stuff to disk (frequency window format):
                        tblhandle.write(assiOr + ' (( attr store1 < ' + higherp1 + \
                                        ' and attr store1 > ' + lowerBoundp1 + ' )\n')
                        tblhandle.write('     and bondedto ( attr store1 < ' + \
                                        higherh1 + ' and attr store1 > ' + \
                                        lowerBoundh1 + ' ) ' + segid1 + ' )\n')
                        tblhandle.write('     ( attr store1 < ' + higherp2 + \
                                        ' and attr store1 > ' + lowerBoundp2 + ' )' + segid2 + '\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume +\
                                            ' hpm1 ' + EACHC.h1ppm +\
                                            ' ppm1 ' + EACHC.p1ppm +\
                                            ' ppm2 ' + EACHC.p2ppm + '\n')
                #3D case, heteronucleus in het2:
                elif (EACHC.h1ppm == None) and\
                     (EACHC.p1ppm != None) and\
                     (EACHC.h2ppm != None) and\
                     (EACHC.p2ppm != None):
                    #look if it's fully assigned:
                    if EACHC.atomnamep1 != (None,) and\
                       EACHC.atomnameh2 != (None,) and\
                       EACHC.atomnamep2 != (None,):
                        #for the fully assigned peaks:
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassi':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassiw1':
                            weight = ' weight 1.1'
                        #write all the stuff to disk (assigned format):
                        tblhandle.write(assiOr + ' ( resid ' + outresidue1 + \
                                        ' and name ' + EACHC.atomnamep1[0] + \
                                        segid1 + ' )\n')
                        tblhandle.write('     ( resid ' + outresidue2 + \
                                        ' and name ' + EACHC.atomnamep2[0] + \
                                        segid2 + ' )\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume + ' ppm1 ' + EACHC.p1ppm +\
                                            ' hpm2 ' + EACHC.h2ppm +\
                                            ' ppm2 ' + EACHC.p2ppm + '\n')
                    else:
                        #for the other peaks (not fully assigned):
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ''
                        elif whichpeaks == 'onlyassi':
                            continue  #don't use this contribution!
                        elif whichpeaks == 'onlyassiw1':
                            continue  #don't use this contribution!
                        #write all the stuff to disk (frequency window format):
                        tblhandle.write(assiOr + '     ( attr store1 < ' + higherp1 + \
                                            ' and attr store1 > ' + lowerBoundp1 + ' )' + segid1 + '\n')
                        tblhandle.write('     (( attr store1 < ' + higherp2 + \
                                        ' and attr store1 > ' + lowerBoundp2 + ' )\n')
                        tblhandle.write('     and bondedto ( attr store1 < ' + \
                                        higherh2 + ' and attr store1 > ' + \
                                        lowerBoundh2 + ' ) ' + segid2 + ' )\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume + ' ppm1 ' + EACHC.p1ppm +\
                                            ' hpm2 ' + EACHC.h2ppm +\
                                            ' ppm2 ' + EACHC.p2ppm + '\n')
                #3D case, het1 & het2 & pro1:
                elif (EACHC.h1ppm != None) and\
                     (EACHC.p1ppm != None) and\
                     (EACHC.h2ppm != None) and\
                     (EACHC.p2ppm == None):
                    #look if it's fully assigned:
                    if EACHC.atomnamep1 != (None,) and\
                       EACHC.atomnameh1 != (None,) and\
                       EACHC.atomnameh2 != (None,):
                        #for the fully assigned peaks:
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassi':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassiw1':
                            weight = ' weight 1.1'
                        #write all the stuff to disk (assigned format):
                        tblhandle.write(assiOr + ' ( resid ' + outresidue1 + \
                                        ' and name ' + EACHC.atomnamep1[0] + \
                                        segid1 + ' )\n')
                        tblhandle.write('     ( (hydrogen) and bondedto (resid ' + outresidue2 + \
                                        ' and name ' + EACHC.atomnameh2[0] + \
                                        segid2 + ' ))\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume + ' ppm1 ' + EACHC.p1ppm +\
                                            ' hpm1 ' + EACHC.h1ppm +\
                                            ' hpm2 ' + EACHC.h2ppm + '\n')
                    else:
                        #for the other peaks (not fully assigned):
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ''
                        elif whichpeaks == 'onlyassi':
                            continue  #don't use this contribution!
                        elif whichpeaks == 'onlyassiw1':
                            continue  #don't use this contribution!
                        #write all the stuff to disk (frequency window format):
                        tblhandle.write(assiOr + ' ( attr store1 < ' + higherp1 + \
                                            ' and attr store1 > ' + lowerBoundp1 + ' )' + segid1 + '\n')
                        tblhandle.write('     ((hydrogen)\n')
                        tblhandle.write('     and bondedto ( attr store1 < ' + \
                                        higherh2 + ' and attr store1 > ' + \
                                        lowerBoundh2 + ' ) ' + segid2 + ' )\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume + ' ppm1 ' + EACHC.p1ppm +\
                                            ' hpm1 ' + EACHC.h1ppm +\
                                            ' hpm2 ' + EACHC.h2ppm + '\n')
                #3D case, het1 & het2 & pro2:
                elif (EACHC.h1ppm != None) and\
                     (EACHC.p1ppm == None) and\
                     (EACHC.h2ppm != None) and\
                     (EACHC.p2ppm != None):
                    #look if it's fully assigned:
                    if EACHC.atomnamep2 != (None,) and\
                       EACHC.atomnameh1 != (None,) and\
                       EACHC.atomnameh2 != (None,):
                        #for the fully assigned peaks:
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassi':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassiw1':
                            weight = ' weight 1.1'
                        #write all the stuff to disk (assigned format):
                        tblhandle.write(assiOr + ' ( (hydrogen) and bondedto (resid ' + outresidue1 + \
                                        ' and name ' + EACHC.atomnameh1[0] + \
                                        segid1 + ' ))\n')
                        tblhandle.write('     ( resid ' + outresidue2 + \
                                        ' and name ' + EACHC.atomnamep2[0] + \
                                        segid2 + ' )\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume + ' ppm2 ' + EACHC.p2ppm +\
                                            ' hpm1 ' + EACHC.h1ppm +\
                                            ' hpm2 ' + EACHC.h2ppm + '\n')
                    else:
                        #for the other peaks (not fully assigned):
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ''
                        elif whichpeaks == 'onlyassi':
                            continue  #don't use this contribution!
                        elif whichpeaks == 'onlyassiw1':
                            continue  #don't use this contribution!
                        #write all the stuff to disk (frequency window format):
                        tblhandle.write(assiOr + ' ((hydrogen)' +\
                                        '     and bondedto ( attr store1 < ' + \
                                        higherh1 + ' and attr store1 > ' + \
                                        lowerBoundh1 + ' ) ' + segid1 + ' )\n')
                        tblhandle.write(' ( attr store1 < ' + higherp2 + \
                                        ' and attr store1 > ' + lowerBoundp2 + ' )' + segid2 + '\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume + ' ppm2 ' + EACHC.p2ppm +\
                                            ' hpm1 ' + EACHC.h1ppm +\
                                            ' hpm2 ' + EACHC.h2ppm + '\n')

                #2D case:
                elif (EACHC.h1ppm == None) and\
                     (EACHC.p1ppm != None) and\
                     (EACHC.h2ppm == None) and\
                     (EACHC.p2ppm != None):
                    #look if it's fully assigned:
                    if EACHC.atomnamep1 != (None,) and\
                       EACHC.atomnamep2 != (None,):
                        #for the fully assigned peaks:
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassi':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'onlyassiw1':
                            weight = ' weight 1.1'
                        #write all the stuff to disk (assigned format):
##                        print outresidue1,EACHC.atomnamep1[0], segid1 #test
                        tblhandle.write(assiOr + ' ( resid ' + outresidue1 + \
                                        ' and name ' + EACHC.atomnamep1[0] + \
                                        segid1 + ' )\n')
                        tblhandle.write('     ( resid ' + outresidue2 + \
                                        ' and name ' + EACHC.atomnamep2[0] + \
                                        segid2 + ' )\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume + ' ppm1 ' + EACHC.p1ppm +\
                                            ' ppm2 ' + EACHC.p2ppm + '\n')
                    else:
                        #for the other peaks (not fully assigned):
                        if whichpeaks == 'all':
                            weight = ''
                        elif whichpeaks == 'allw1':
                            weight = ' weight 1.1'
                        elif whichpeaks == 'allassiw1':
                            weight = ''
                        elif whichpeaks == 'onlyassi':
                            continue  #don't use this contribution!
                        elif whichpeaks == 'onlyassiw1':
                            continue  #don't use this contribution!
                        #write all the stuff to disk (frequency window format):
                        tblhandle.write(assiOr + ' ( attr store1 < ' + higherp1 + \
                                        ' and attr store1 > ' + lowerBoundp1 + ' )' + segid1 + '\n')
                        tblhandle.write('     ( attr store1 < ' + higherp2 + \
                                        ' and attr store1 > ' + lowerBoundp2 + ' )' + segid2 + '\n')
                        if contCounter == 0:
                            tblhandle.write('     ' + dist + ' ' + errm + ' ' +\
                                            errp + ' peak ' + \
                                            EACHC.peakNumber +\
                                            outSpectrumName + weight + ' volume ' +\
                                            outvolume + ' ppm1 ' + EACHC.p1ppm +\
                                            ' ppm2 ' + EACHC.p2ppm + '\n')
                contCounter = contCounter + 1
        #closing file handles and bye:
        tblhandle.write('\n')
        tblhandle.close()

    def WriteNMRViewPeaksKnownHeader(self,fileName):
        """
        use this function when you have already read the header of the original
        NMRView file
        """
        if self._nmrviewHet1 and self._nmrviewHet1 != 'N':
            het1width = self._nmrviewfourthLineList[string.atoi(self._nmrviewHet1)-1]
            het1frequency = self._nmrviewfifthLineList[string.atoi(self._nmrviewHet1)-1]
            het1name = self._nmrviewDimNameList[string.atoi(self._nmrviewHet1)-1]
        else:
            het1width = None
            het1frequency = None
            het1name = None
            
        if self._nmrviewHet2 and self._nmrviewHet2 != 'N':
            het2width = self._nmrviewfourthLineList[string.atoi(self._nmrviewHet2)-1]
            het2frequency = self._nmrviewfifthLineList[string.atoi(self._nmrviewHet2)-1]
            het2name = self._nmrviewDimNameList[string.atoi(self._nmrviewHet2)-1]
        else:
            het2width = None
            het2frequency = None
            het2name = None

        pro2width = self._nmrviewfourthLineList[string.atoi(self._nmrviewPro2)-1]
        pro2frequency = self._nmrviewfifthLineList[string.atoi(self._nmrviewPro2)-1]
        pro1width = self._nmrviewfourthLineList[string.atoi(self._nmrviewPro1)-1]
        pro1frequency = self._nmrviewfifthLineList[string.atoi(self._nmrviewPro1)-1]
        pro1name = self._nmrviewDimNameList[string.atoi(self._nmrviewPro1)-1]
        pro2name = self._nmrviewDimNameList[string.atoi(self._nmrviewPro2)-1]

#        print pro1width, pro2frequency #test
#        print self._nmrviewfifthLineList #test
        self.WriteNMRViewPeaks(fileName=fileName,\
                               het1=self._nmrviewHet1,\
                               pro1=self._nmrviewPro1,\
                               het2=self._nmrviewHet2,\
                               pro2=self._nmrviewPro2,\
                               het1name=het1name,\
                               pro1name=pro1name,\
                               het2name=het2name,\
                               pro2name=pro2name,\
                               het1width=het1width,\
                               pro1width=pro1width,\
                               het2width=het2width,\
                               pro2width=pro2width,\
                               het1frequency=het1frequency,\
                               pro1frequency=pro1frequency,\
                               het2frequency=het2frequency,\
                               pro2frequency=pro2frequency,\
                               spectrumFileName=self._nmrviewNameThirdLine)
        
    def WriteNMRViewPeaks(self, fileName,\
                          het1,pro1,het2,pro2,\
                          het1name,pro1name,het2name,pro2name,\
                          het1width,pro1width,het2width,pro2width,\
                          het1frequency,pro1frequency,\
                          het2frequency, pro2frequency,\
                          spectrumFileName):
        """
        writes a list in nmrView xpk format

        .L assignment
        .P chemical shift
        .W 0.000
        .B 0.000
        .E ++
        .J 0.000
        .U {?}
        """

        #default values: TODO: parse the defaultValues as well!!!
        defaultW = '0.000'
        defaultB = '0.000'
        defaultE = '++'
        defaultJ = '0.000'
        defaultU = '{?}'
        defaultstat = '0'
        defaultflag0 = '0'
        defaultcomment = '{?}'
        print 'writing the NMRView xpk file', fileName
        xpkHandle = open(fileName, 'w')

        # the header - first line:
        firstLine = "label dataset sw sf \n"
        xpkHandle.write(firstLine)

        # get dimension names, widths and frequencies in correct order:
        namesOrdered = [None,None,None,None]
        widthsOrdered = [None,None,None,None]
        frequenciesOrdered = [None,None,None,None]
        if het1 and het1 != 'N':
            namesOrdered[string.atoi(het1)-1] = het1name
            widthsOrdered[string.atoi(het1)-1] = het1width
            frequenciesOrdered[string.atoi(het1)-1] = het1frequency
        if het2 and het2 != 'N':
            namesOrdered[string.atoi(het2)-1] = het2name
            widthsOrdered[string.atoi(het2)-1] = het2width
            frequenciesOrdered[string.atoi(het2)-1] = het2frequency
        namesOrdered[string.atoi(pro1)-1] = pro1name
        namesOrdered[string.atoi(pro2)-1] = pro2name
        widthsOrdered[string.atoi(pro1)-1] = pro1width
        widthsOrdered[string.atoi(pro2)-1] = pro2width
        frequenciesOrdered[string.atoi(pro1)-1] = pro1frequency
        frequenciesOrdered[string.atoi(pro2)-1] = pro2frequency
        
        #second line contains the names of all dimensions:
        secondLine = ''
        for eachN in namesOrdered:
            if eachN:
                secondLine = secondLine + eachN + ' '
        secondLine =  secondLine + '\n'
        xpkHandle.write(secondLine)

        #third line contains the spectrum filename:
        xpkHandle.write(spectrumFileName + '\n')

        #fourth line contains the sweep widths:
        fourthLine = ''
        for eachN in widthsOrdered:
            if eachN:
                fourthLine = fourthLine + '{' + eachN + ' } '
        fourthLine =  fourthLine[:-1] + '\n'
        xpkHandle.write(fourthLine)


        #fifth line contains the frequencies:
        fifthLine = ''
        for eachN in frequenciesOrdered:
            if eachN:
                fifthLine = fifthLine + '{' + eachN + ' } '
        fifthLine =  fifthLine[:-1] + '\n'
        xpkHandle.write(fifthLine)

        #sixth line contains the data description:
        tags = ['.L', '.P', '.W', '.B', '.E', '.J', '.U']
        theRest = ['vol', 'int', 'stat', 'comment', 'flag0']

        sixthLine = ''
        for eachN in namesOrdered:
            if eachN:
                for eachT in tags:
                    sixthLine = sixthLine + ' ' + eachN + eachT + ' '
        for eachTR in theRest:
            sixthLine = sixthLine + ' ' + eachTR + ' '
        xpkHandle.write(sixthLine + '\n')

        bracketsPA = re.compile('{\s*}')
        for EACHP in self.peakslist:
            contCounter = 0
            lineString = EACHP.peakNumber + '  '
            #get all the atomnames:
            het1atomnames = []
            pro1atomnames = []
            het2atomnames = []
            pro2atomnames = []
            aa1s = []
            aa2s = []
            for EACHC in EACHP.contributions:
                #always take the first of the tuple:
                het1atomnames.append(EACHC.atomnameh1[0])
                pro1atomnames.append(EACHC.atomnamep1[0])
                het2atomnames.append(EACHC.atomnameh2[0])
                pro2atomnames.append(EACHC.atomnamep2[0])
                aa1s.append(EACHC.residue1)
                aa2s.append(EACHC.residue2)
            #get the frequencies, volumes, intensities:
            if len(EACHP.contributions) > 0:
                h1ppm = EACHP.contributions[0].h1ppm
                h2ppm = EACHP.contributions[0].h2ppm
                p1ppm = EACHP.contributions[0].p1ppm
                p2ppm = EACHP.contributions[0].p2ppm
                #for the volumes:
                outVolume = EACHP.contributions[0].volume
                outIntensity= EACHP.contributions[0].intensity

            if outVolume == None:
                outVolume = '0.000'
            if outIntensity == None:
                outIntensity = '0.000'


            #ordering the atomnames, residue numbers, ppms:
            atomnamesOrdered = [(None,),(None,),(None,),(None,)]
            ppmsOrdered = [None,None,None,None]
            residuesOrdered = [None,None,None,None]
            if het1 and het1 != 'N':
                atomnamesOrdered[string.atoi(het1)-1] = het1atomnames
                ppmsOrdered[string.atoi(het1)-1] = h1ppm
                residuesOrdered[string.atoi(het1)-1] = aa1s
            if het2 and het2 != 'N':
                atomnamesOrdered[string.atoi(het2)-1] = het2atomnames
                ppmsOrdered[string.atoi(het2)-1] = h2ppm
                residuesOrdered[string.atoi(het2)-1] = aa2s
            atomnamesOrdered[string.atoi(pro1)-1] = pro1atomnames
            atomnamesOrdered[string.atoi(pro2)-1] = pro2atomnames
            ppmsOrdered[string.atoi(pro1)-1] = p1ppm
            ppmsOrdered[string.atoi(pro2)-1] = p2ppm
            residuesOrdered[string.atoi(pro1)-1] = aa1s
            residuesOrdered[string.atoi(pro2)-1] = aa2s

            # output:
            for iii in range(0,4):
                if not namesOrdered[iii]: continue
                #1. assignment (.L):
                outAssi = '{'

                www = 0
                for eachAssi in atomnamesOrdered[iii]:
                    if residuesOrdered[iii][www]:
                        outRes = residuesOrdered[iii][www]
                    else:
                        outRes = '99999' #default
                    if not eachAssi:
                        continue
                    outAssi = outAssi + outRes + '.' + eachAssi + ' '
                    www = www +1
                if not outAssi[-1] == '{':
                    outAssi=outAssi[:-1]
                outAssi = outAssi + '}  '
                outAssi = bracketsPA.sub('{?}',outAssi)
                lineString = lineString + outAssi
                
                #2. ppm (.P):
                lineString = lineString + ' ' +str(ppmsOrdered[iii]) + '   '
                #3. .W:
                lineString = lineString + defaultW + '   '
                #4. .B:
                lineString = lineString + defaultW + '   '
                #5. .E:
                lineString = lineString + defaultE + '   '
                #6. .J:
                lineString = lineString + defaultJ + '   '
                #7. .U: 
                lineString = lineString + defaultU + '   '
            #8. volume:
            lineString = lineString[:-1] + outVolume + ' '
            #9. intensity:
            lineString = lineString + outIntensity + ' '
            #10. stat
            lineString = lineString + defaultstat + ' '
            #11. comment:
            lineString = lineString + defaultcomment + ' '
            #12. flag:
            lineString = lineString + defaultflag0 
                
            xpkHandle.write(lineString + '\n')
        
    def WriteOldList(self, fileName):
        """
        writes a list in ARIA's old .list format
        """
        file = TextFile.TextFile(fileName, 'w')
        print 'writing the data in .list format to', fileName
        if len(self.peakslist) == 0:
            print 'peaklist empty. WriteList method aborted.'
            return
        for PEAK in self.peakslist:
            for CON in PEAK.contributions:
                if (CON.peakType == '2') or (CON.peakType == '4'):
                    sign = '-'
                else:
                    sign = '+'
                file.write(sign + ' ' + CON.peakNumber + ' ' +\
                           CON.residue1 + ' ' + CON.aa1 + ' ' +\
                           CON.atomnamep1[0] + ' ' +\
                           CON.residue2 + ' ' + CON.aa2 + ' ' +\
                           CON.atomnamep2[0] + ' ' + CON.distanceAve + ' ' +\
                           CON.volume + ' ' + CON.lowerBound + ' ' + CON.upperBound +\
                           CON.p1ppm + ' ' + CON.p2ppm + ' ' +\
                           CON.atomnumberp1 + ' ' + CON.atomnumberp2 + ' ' +\
                           CON.contribution + ' ' + CON.weight + ' ' +\
                           CON.nta + '\n')
        file.close()

    def WriteSparky(self, fileName):
    
        """
        writes a list in Sparky .list format

           Assignment         w1         w2         w3        Volume    Fit Height 

                 ?-?-?      3.977     45.152      3.982   9.89e+10 lo    211314157 
         I41HD#-CD-HD#      0.818     12.622      0.820   4.10e+10 lo    138032890 
         A66HB#-CB-HB#     -0.053     15.971     -0.050   2.60e+10 lo     98784762
       """

        File= TextFile.TextFile(fileName, 'w')
        print 'writing the data in .list format to', fileName
        if len(self.peakslist) == 0:

             print 'peaklist empty. WriteList method aborted.'


        
        iii = 0
        for PEAK in self.peakslist:
            for CON  in PEAK.contributions:
                iii = iii + 1
                h1ppm=CON.h1ppm
                p1ppm=CON.p1ppm
                h2ppm=CON.h2ppm
                p2ppm=CON.p2ppm

                
                #calculates spectrum dimension
                DM = 0
              

                if h1ppm == None:
                    DM = DM + 1
                if p1ppm == None:
                    DM = DM + 1
                if h2ppm == None:
                    DM = DM + 1
                if p2ppm == None:
                    DM = DM + 1

                Dimension = 4 - DM
    
                #Writes the first line

                #4Dcase
                if Dimension==4:
                    File.write('       %s        %s         %s         %s         %s       %s   %s' %('Assignment','w1','w2','w3','w4','volume','Data Height'))
                #3Dcase
                if Dimension==3:
                    File.write('       %s        %s         %s         %s       %s   %s' %('Assignment','w1','w2','w3','volume','Fit Height'))
                #2Dcase
                if Dimension==2:
                     File.write('       %s        %s         %s       %s   %s' %('Assignment','w1','w2','volume','Data Height'))
            
                File.write('\n')
                File.write('\n')

            
            break
                
        
        
        for PEAK in self.peakslist:
            for CON  in PEAK.contributions:
                
                if CON.aa1 != None:
                    aa1 = AminoAcid.AminoAcid(CON.aa1)[0]
                else:
                    aa1=CON.aa1
                residue1=CON.residue1
                atomnameh1=CON.atomnameh1
                atomnamep1=CON.atomnamep1
                residue2=CON.residue2
                if CON.aa2 != None:
                    aa2=AminoAcid.AminoAcid(CON.aa2)[0]
                else:
                    aa2=CON.aa2
                #Maybe I should directly transform here tuples in strings    
                atomnameh2=CON.atomnameh2
                atomnamep2=CON.atomnamep2
                
                volume=string.atof(CON.volume)
                
                h1ppm=CON.h1ppm
                p1ppm=CON.p1ppm
                h2ppm=CON.h2ppm
                p2ppm=CON.p2ppm
                intensity = CON.intensity

                sparkylittlecomment=CON.sparkylittlecomment

                #IMPAGINATION
                
                p1ppmlength=len(str(p1ppm))
                numWhiteSpaces1 =11-p1ppmlength             
                WhiteSpaces1=numWhiteSpaces1* ' '
                p1ppm             =WhiteSpaces1+p1ppm
                
                if h1ppm:
                    h1ppmlength=len(str(h1ppm))
                    numWhiteSpaces2 =11-h1ppmlength             
                    WhiteSpaces2=numWhiteSpaces2* ' '
                    h1ppm             =WhiteSpaces2+h1ppm
                else:
                    h1ppm = ''
                    
                if h2ppm:
                    h2ppmlength=len(str(h2ppm))
                    numWhiteSpaces3 =11-h2ppmlength             
                    WhiteSpaces3=numWhiteSpaces3* ' '
                    h2ppm             =WhiteSpaces3+h2ppm
                else:
                    h2ppm = ''
                    
                
                p2ppmlength=len(str(p2ppm))
                numWhiteSpaces4 =11-p2ppmlength             
                WhiteSpaces4=numWhiteSpaces4* ' '
                p2ppm             =WhiteSpaces4+p2ppm

               

                #in case sparkylittlecomment doesn't exist
                #it has to write 3 whitespaces.
                if sparkylittlecomment:
                    sparkylittlecomment = ' '+sparkylittlecomment
                else:
                    sparkylittlecomment = '   '
                    
                intensitylength = len(str(intensity))
                numWhiteSpacesInt =13-intensitylength          
                WhiteSpacesInt=numWhiteSpacesInt* ' '
                intensity          =WhiteSpacesInt+intensity

                #peakdescription:2D,3D and 4D are different.

                #
                if aa1 or  aa2 != None:
                    #interresidualcase
                    if aa1 != aa2:
                        if Dimension==4:
                             peakdescription=aa1+str(residue1)+atomnamep1[0]+'-'+atomnameh1[0]+'-'+aa2+str(residue2)+atomnameh2[0]+'-'+atomnamep2[0]
                        if Dimension==3:
                            if atomnameh1[0]:
                                 peakdescription=aa2+str(residue2)+atomnamep2[0]+'-'+aa1+str(residue1)+atomnameh1[0]+'-'+atomnamep1[0]
                                 h2ppm=''
                            if atomnameh2[0]:
                                 peakdescription=aa1+str(residue1)+atomnamep1[0]+'-'+aa2+str(residue2)+atomnameh2[0]+'-'+atomnamep2[0]
                                 h1ppm=''
                        if Dimension==2:
                            peakdescription=aa1+str(residue1)+atomnamep1[0]+'-'+aa2+str(residue2)+'-'+atomnamep2[0]
                            h1ppm=''
                            h2ppm=''
                    #intraresidual case
                    if aa1 == aa2:
                        if Dimension==4:
                            peakdescription=aa1+str(residue1)+atomnamep1[0]+'-'+atomnameh1[0]+'-'+atomnameh2[0]+'-'+atomnamep2[0]
                        if Dimension==3:
                            if atomnameh1[0]:
                                peakdescription=aa1+str(residue1)+atomnamep2[0]+'-'+atomnameh1[0]+'-'+atomnamep1[0]
                                h2ppm=''
                            if atomnameh2[0]:
                                peakdescription=aa1+str(residue1)+atomnamep1[0]+'-'+atomnameh2[0]+'-'+atomnamep2[0]
                                h1ppm=''
                        if Dimension==2:
                            peakdescription=aa1+str(residue1)+atomnamep1[0]+'-'+atomnamep2[0]
                            h1ppm=''
                            h2ppm=''
                            
                    peakdescriptionlength=len(peakdescription)
                    numWhiteSpaces0=18-peakdescriptionlength
                    WhiteSpaces0=numWhiteSpaces0* ' '
                    peakdescription=WhiteSpaces0+peakdescription
                else:
                    if Dimension == 4:
                        peakdescription='           ?-?-?-?'
                    if Dimension == 3:
                        peakdescription='             ?-?-?'
                        if atomnameh1 == None:
                            h1ppm = ''
                        if atomnameh2 == None:
                            h2ppm = ''
                    if Dimension == 2:
                        peakdescription='               ?-?'
                        h1ppm = ''
                        h2ppm = ''

                     
                
                    
             
                print ( peakdescription+p1ppm+h1ppm+h2ppm+p2ppm+'   %1.2e' % (volume)+ sparkylittlecomment+intensity+'\n') # %(volume)
                File.write( peakdescription+p1ppm+h1ppm+h2ppm+p2ppm+'   %1.2e'  %(volume)+ sparkylittlecomment+intensity+'\n')
                File.write('\n')
                
        File.close()                          
                
            
    def WriteUplLol(self, fileroot):
        """
        writes out .upl and .lol files for MolMol
        pseudoatoms are always used!!! all atomnames are converted to pseudoatoms
        important:
        I use the contribution to the peak as relative limit in these files!
        If more than one atomname is possible, all the possibilities are
        written out!
        
        accepted peaks:
          _amb_acc.upl
          _unamb_acc.upl
          _amb_acc.lol
          _unamb_acc.lol
        rejected peaks:
          _amb_rej.upl
          _unamb_rej.upl
          _amb_rej.lol
          _unamb_rej.lol

        e.g. WriteUplLol('/home/linge/vasp')
             will create vasp.upl and vasp.lol in /home/linge

        """
        uplAmbAccFile = fileroot + '_amb_acc.upl'
        lolAmbAccFile = fileroot + '_amb_acc.lol'
        uplUnambAccFile = fileroot + '_unamb_acc.upl'
        lolUnambAccFile = fileroot + '_unamb_acc.lol'
        uplAmbRejFile = fileroot + '_amb_rej.upl'
        lolAmbRejFile = fileroot + '_amb_rej.lol'
        uplUnambRejFile = fileroot + '_unamb_rej.upl'
        lolUnambRejFile = fileroot + '_unamb_rej.lol'
        print '  writing all data to the following MolMol files:\n   ',\
              uplAmbAccFile, '\n   ', lolAmbAccFile, '\n   ',\
              uplUnambAccFile, '\n   ', lolUnambAccFile, '\n   ',\
              uplAmbRejFile, '\n   ', lolAmbRejFile, '\n   ',\
              uplUnambRejFile, '\n   ', lolUnambRejFile
        uplAmbAccHandle = TextFile.TextFile(uplAmbAccFile, 'w')
        lolAmbAccHandle = TextFile.TextFile(lolAmbAccFile, 'w')
        uplUnambAccHandle = TextFile.TextFile(uplUnambAccFile, 'w')
        lolUnambAccHandle = TextFile.TextFile(lolUnambAccFile, 'w')
        uplAmbRejHandle = TextFile.TextFile(uplAmbRejFile, 'w')
        lolAmbRejHandle = TextFile.TextFile(lolAmbRejFile, 'w')
        uplUnambRejHandle = TextFile.TextFile(uplUnambRejFile, 'w')
        lolUnambRejHandle = TextFile.TextFile(lolUnambRejFile, 'w')
        
        for PEAK in self.peakslist:
##             if PEAK.contributions == []:  #test
##                 continue                  #test
            for CON in PEAK.contributions:
                #I have to loop over the proton atomname tuples:
                for each1 in CON.atomnamep1:
                    for each2 in CON.atomnamep2:
                        if CON.upperBound == None:
                            outupperBound = '6'  #default: 6 Angstrom
                        else:
                            outupperBound = CON.upperBound
                        if CON.lowerBound == None:
                            outlowerBound = '0'  #default: 0 Angstrom
                        else:
                            outlowerBound = CON.lowerBound
                        if CON.contribution == None:
                            outcontribution = '1'
                        else:
                            outcontribution = CON.contribution
                        if CON.residue1 == None:
                            outResidue1 = 0
                        else:
                            outResidue1 = string.atoi(CON.residue1)
                        if CON.residue2 == None:
                            outResidue2 = 0
                        else:
                            outResidue2 = string.atoi(CON.residue2)

                        #convert the atomnames back to pseudoatoms if possible:
                        if CON.residue1:
                            converted1 = PseudoAtom.Atom2Pseudo(each1, CON.aa1)
                        else:
                            converted1 = None
#                        print converted1, (each1, CON.aa1) #test
#                        if converted1 and each1[-1] == '%':
                        if converted1: #only those ending with a '%'
                            outAtomname1 = converted1[0]
                        else:
                            outAtomname1 = each1
                        if CON.residue2:
                            converted2 = PseudoAtom.Atom2Pseudo(each2, CON.aa2)
                        else:
                            converted2 = None
#                        print converted2, (each2, CON.aa2) #test
#                        if converted2 and each2[-1] == '%': #only those ending with a '%'
                        if converted2:
                            outAtomname2 = converted2[0]
                        else:
                            outAtomname2 = each2
                                                   
                        if len(PEAK.contributions) == 1 and CON.peakType == '1':
                            uplUnambAccHandle.write('%3i %-4s %-5s %3i %-4s %-7s %-5.2f %-5.2f \n' %\
                                                    (outResidue1,\
                                                     CON.aa1,\
                                                     outAtomname1,\
                                                     outResidue2,\
                                                     CON.aa2,\
                                                     outAtomname2,\
                                                     string.atof(outupperBound),\
                                                     string.atof(outcontribution)))
                            lolUnambAccHandle.write('%3i %-4s %-5s %3i %-4s %-7s %-5.2f %-5.2f \n' %\
                                                    (outResidue1,\
                                                     CON.aa1,\
                                                     outAtomname1,\
                                                     outResidue2,\
                                                     CON.aa2,\
                                                     outAtomname2,\
                                                     string.atof(outlowerBound),\
                                                     string.atof(outcontribution)))
                        if len(PEAK.contributions) > 1 and CON.peakType == '3':
                            uplAmbAccHandle.write('%3i %-4s %-5s %3i %-4s %-7s %-5.2f %-5.2f \n' %\
                                                  (outResidue1,\
                                                   CON.aa1,\
                                                   outAtomname1,\
                                                   outResidue2,\
                                                   CON.aa2,\
                                                   outAtomname2,\
                                                   string.atof(outupperBound),\
                                                   string.atof(outcontribution)))
                            lolAmbAccHandle.write('%3i %-4s %-5s %3i %-4s %-7s %-5.2f %-5.2f \n' %\
                                                  (outResidue1,\
                                                   CON.aa1,\
                                                   outAtomname1,\
                                                   outResidue2,\
                                                   CON.aa2,\
                                                   outAtomname2,\
                                                   string.atof(outlowerBound),\
                                                   string.atof(outcontribution)))
                        if len(PEAK.contributions) == 1 and CON.peakType == '2':
                            uplUnambRejHandle.write('%3i %-4s %-5s %3i %-4s %-7s %-5.2f %-5.2f \n' %\
                                                    (outResidue1,\
                                                     CON.aa1,\
                                                     outAtomname1,\
                                                     outResidue2,\
                                                     CON.aa2,\
                                                     outAtomname2,\
                                                     string.atof(outupperBound),\
                                                     string.atof(outcontribution)))
                            lolUnambRejHandle.write('%3i %-4s %-5s %3i %-4s %-7s %-5.2f %-5.2f \n' %\
                                                    (outResidue1,\
                                                     CON.aa1,\
                                                     outAtomname1,\
                                                     outResidue2,\
                                                     CON.aa2,\
                                                     outAtomname2,\
                                                     string.atof(outlowerBound),\
                                                     string.atof(outcontribution)))
                        if len(PEAK.contributions) > 1 and CON.peakType == '4':
                            uplAmbRejHandle.write('%3i %-4s %-5s %3i %-4s %-7s %-5.2f %-5.2f \n' %\
                                                  (outResidue1,\
                                                   CON.aa1,\
                                                   outAtomname1,\
                                                   outResidue2,\
                                                   CON.aa2,\
                                                   outAtomname2,\
                                                   string.atof(outupperBound),\
                                                   string.atof(outcontribution)))
                            lolAmbRejHandle.write('%3i %-4s %-5s %3i %-4s %-7s %-5.2f %-5.2f \n' %\
                                                  (outResidue1,\
                                                   CON.aa1,\
                                                   outAtomname1,\
                                                   outResidue2,\
                                                   CON.aa2,\
                                                   outAtomname2,\
                                                   string.atof(outlowerBound),\
                                                   string.atof(outcontribution)))

        #close all the filehandles:
        uplAmbAccHandle.close()
        lolAmbAccHandle.close()
        uplUnambAccHandle.close()
        lolUnambAccHandle.close()
        uplAmbRejHandle.close()
        lolAmbRejHandle.close()
        uplUnambRejHandle.close()
        lolUnambRejHandle.close()


    def WriteXeasyAssignPeaks(self, fileroot, het1, pro1, het2, pro2):
        """
        writes out the XEASY .assign and .peaks files
        fileroot should be a name without any extensions
        (no .prot file is created, just copy the .prot file of the last iteration)
        
        I use the same convention as in the .html file for XEASY:
        het1 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro1 = '1', '2', '3', '4' or 'N' (column number or not used)
        het2 = '1', '2', '3', '4' or 'N' (column number or not used)
        pro2 = '1', '2', '3', '4' or 'N' (column number or not used)

        I use pro1 and pro2 for the matrix -> binary format conversion

        IMPORTANT: it is only possible to write out these files if the
                   atomnumbers have been read in from other XEASY files!
                   also the spectrum type, integration method must be
                   available!

        NOTE: the matrix is set up with the protons!
              more than 10 ambiguous assignments are not possible in XEASY
        """
##         #have a look if there are atomnumbers defined:
##         print 'checking the atomnumbers:',
##         print 'o.k.' #test - I have to work on that!
        
        #creating the fileNames from fileroot:
        peaksfile = fileroot + '.peaks'
        assignfile = fileroot + '.assign'

        #only for XEASY stuff I need:
        try:
            import numpy
        except:
            print 'WARNING: COULD NOT IMPORT numpyc Python module'
            print '         NoeList.WriteXeasyAssignPeaks() aborted!'
            print '-'*60
            traceback.print_exc(file=sys.stdout)
            print '-'*60
            sys.exit()
        #from Aria.DataformatConversion import ReadXeasy
        import ReadXeasy
        
        #message to stdout:
        print 'writing all data to the following XEASY files:\n ',\
              peaksfile, '\n ', assignfile
        
        #opening filehandles:
        peakshandle = TextFile.TextFile(peaksfile, 'w')
        assignhandle = TextFile.TextFile(assignfile, 'w')
        
        #write some comments (dimension and INAME) for .peaks:
        peakshandle.write('# Number of dimensions ' + str(self.dimension) + '\n')
        if self.dimension >= 1:
            peakshandle.write('#INAME 1 ' + self._iname1 + '\n')
        if self.dimension >= 2:
            peakshandle.write('#INAME 2 ' + self._iname2 + '\n')
        if self.dimension >= 3:
            peakshandle.write('#INAME 3 ' + self._iname3 + '\n')
        if self.dimension >= 4:
            peakshandle.write('#INAME 4 ' + self._iname4 + '\n')
            
        #write first line with the dimension for .assign:
        assignhandle.write(str(self.dimension) + '\n')
        
        #first loop over all the peaks, only one contribution possible in .peaks:
        for PEAK in self.peakslist:
            #for empty fields:
            if PEAK.contributions == []:
                print 'NOTE: peak without contributions are rejected.'
                continue
                
            #get the assigned atom numbers from the first contribution:
            if het1 != 'N':
                outnumberh1 = PEAK.contributions[0].atomnumberh1
            outnumberp1 = PEAK.contributions[0].atomnumberp1
            if het2 != 'N':
                outnumberh2 = PEAK.contributions[0].atomnumberh2
            outnumberp2 = PEAK.contributions[0].atomnumberp2

            #looking whether the assigned atomnumbers are the same for every contribution,
            #because only one contribution is allowed in .peaks:
            for iii in range(len(PEAK.contributions) - 1):
                if het1 != 'N':
                    if PEAK.contributions[iii].atomnumberh1 == PEAK.contributions[iii+1].atomnumberh1:
                        outnumberh1 = PEAK.contributions[iii].atomnumberh1
                    else:
                        outnumberh1 = '0'   #default
                if PEAK.contributions[iii].atomnumberp1 == PEAK.contributions[iii+1].atomnumberp1:
                    outnumberp1 = PEAK.contributions[iii].atomnumberp1
                else:
                    outnumberp1 = '0'   #default
                if het2 != 'N':
                    if PEAK.contributions[iii].atomnumberh2 == PEAK.contributions[iii+1].atomnumberh2:
                        outnumberh2 = PEAK.contributions[iii].atomnumberh2
                    else:
                        outnumberh2 = '0'   #default
                if PEAK.contributions[iii].atomnumberp2 == PEAK.contributions[iii+1].atomnumberp2:
                    outnumberp2 = PEAK.contributions[iii].atomnumberp2
                else:
                    outnumberp2 = '0'   #default    

            #arrange the outnumbers in a list in the correct order (for the columns):
            outnumberlist = [0, None, None, None, None]  #zero-element is not used
            if het1 != 'N':
                outnumberlist[string.atoi(het1)] = outnumberh1
            outnumberlist[string.atoi(pro1)] = outnumberp1
            if het2 != 'N':
                outnumberlist[string.atoi(het2)] = outnumberh2
            outnumberlist[string.atoi(pro2)] = outnumberp2

            #arrange the ppms in a list in the correct order (for the columns):
            outppmlist = [0, None, None, None, None]  #zero-element is not used
            if het1 != 'N':
                outppmlist[string.atoi(het1)] = PEAK.contributions[0].h1ppm
            outppmlist[string.atoi(pro1)] = PEAK.contributions[0].p1ppm
##            print string.atoi(pro1), PEAK.contributions[0].p1ppm #test
            if het2 != 'N':
                outppmlist[string.atoi(het2)] = PEAK.contributions[0].h2ppm
            outppmlist[string.atoi(pro2)] = PEAK.contributions[0].p2ppm

            #replace None by 0 for every ppm in outppmlist, for atof():
            for www in range(0, 4):
                if None in outppmlist:
                    index = outppmlist.index(None)
                    outppmlist.remove(None)
                    outppmlist.insert(index, '0')
                else:
                    break
            
##            print outppmlist #test
            
            #replace None by 0 for volume and volumeError, for atof():
            outvolume = PEAK.contributions[0].volume
            outvolumeError = PEAK.contributions[0].volumeError
            if outvolume == None:
                outvolume = '0'
            if outvolumeError == None:
                outvolumeError = '0'
            if PEAK.contributions[0].peakNumber == None:
                outPeakNumber = 0
            else:
                outPeakNumber = string.atoi(PEAK.contributions[0].peakNumber)
            if PEAK.contributions[0].peakType == None:
                outPeakType = 6
            else:
                outPeakType = string.atoi(PEAK.contributions[0].peakType)   

            if outnumberlist[1] == None:
                outNumber1 = 0
            else:
                outNumber1 = string.atoi(outnumberlist[1])
            if outnumberlist[2] == None:
                outNumber2 = 0
            else:
                outNumber2 = string.atoi(outnumberlist[2])
            if outnumberlist[3] == None:
                outNumber3 = 0
            else:
                outNumber3 = string.atoi(outnumberlist[3])
            if outnumberlist[4] == None:
                outNumber4 = 0
            else:
                outNumber4 = string.atoi(outnumberlist[4])

            #replace integration method by 'm' (if None):
            integrationM = PEAK.contributions[0].xeasyinte
            if not integrationM:
                integrationM = 'm'

            #write out with the correct order:
            if self.dimension == 4:
                   peakshandle.write('%4i %7.3f %7.3f %7.3f %7.3f %1i %-9s %10.3e %9.2e %1s %3i %4i %4i %4i %4i %1i\n' %\
                                  (outPeakNumber,\
                                   string.atof(outppmlist[1]),\
                                   string.atof(outppmlist[2]),\
                                   string.atof(outppmlist[3]),\
                                   string.atof(outppmlist[4]),\
                                   outPeakType,\
                                   PEAK.contributions[0].xeasyudst,\
                                   string.atof(outvolume),\
                                   string.atof(outvolumeError),\
                                   integrationM,\
                                   0,          #unused1, always 0
                                   outNumber1,\
                                   outNumber2,\
                                   outNumber3,\
                                   outNumber4,\
                                   0))         #unused2, always 0
            elif self.dimension == 3:
                peakshandle.write('%4i %7.3f %7.3f %7.3f %1i %-9s %10.3e %9.2e %1s %3i %4i %4i %4i %1i\n' %\
                                  (outPeakNumber,\
                                   string.atof(outppmlist[1]),\
                                   string.atof(outppmlist[2]),\
                                   string.atof(outppmlist[3]),\
                                   outPeakType,\
                                   PEAK.contributions[0].xeasyudst,\
                                   string.atof(outvolume),\
                                   string.atof(outvolumeError),\
                                   integrationM,\
                                   0,          #unused1, always 0
                                   outNumber1,\
                                   outNumber2,\
                                   outNumber3,\
                                   0))         #unused2, always 0
            elif self.dimension == 2:
               peakshandle.write('%4i %7.3f %7.3f %1i %-9s %10.3e %9.2e %1s %3i %4i %4i %1i\n' %\
                                  (outPeakNumber,\
                                   string.atof(outppmlist[1]),\
                                   string.atof(outppmlist[2]),\
                                   outPeakType,\
                                   PEAK.contributions[0].xeasyudst,\
                                   string.atof(outvolume),\
                                   string.atof(outvolumeError),\
                                   integrationM,\
                                   0,          #unused1, always 0
                                   outNumber1,\
                                   outNumber2,\
                                   0))         #unused2, always 0
               
            #initialize the lists of possible assignments in all dimensions:
            possibleh1 = []
            possiblep1 = []
            possibleh2 = []
            possiblep2 = []
                
            #loop through all the contributions, set up the matrix:
            if len(PEAK.contributions) > 1:
                for CON in PEAK.contributions:
                    if PEAK.contributions == []:
                        continue
                    if not (CON.atomnumberh1 in possibleh1):
                        possibleh1.append(CON.atomnumberh1)
                    if not (CON.atomnumberp1 in possiblep1):
                        possiblep1.append(CON.atomnumberp1)
                    if not (CON.atomnumberh2 in possibleh2):
                        possibleh2.append(CON.atomnumberh2)
                    if not (CON.atomnumberp2 in possiblep2):
                        possiblep2.append(CON.atomnumberp2)
                if possibleh1.count(None) != 0:
                    possibleh1.remove(None)
                if possiblep1.count(None) != 0:
                    possiblep1.remove(None)
                if possibleh2.count(None) != 0:
                    possibleh2.remove(None)
                if possiblep2.count(None) != 0:
                    possiblep2.remove(None)
##                 print possibleh1, possiblep1, possibleh2, possiblep2  #test
                if (len(possiblep1)) > 0 and (len(possiblep2) > 0):
                    matrix = numpy.zeros((len(possiblep1), len(possiblep2)))
                for CON in PEAK.contributions:
                    if PEAK.contributions == []:  #test
                        continue  #test
##                        print matrix, CON.atomnumberp1, CON.atomnumberp2  #test
                    matrix[possiblep1.index(CON.atomnumberp1),\
                           possiblep2.index(CON.atomnumberp2)] = 1

                #calculate the 4 binary numbers from the matrix:
                XAP = ReadXeasy.XeasyAssignPeak()
                #print warnings, if a peak has more than 10 contributions:
                if len(possibleh1) > 11:
                    print 'NOTE: peak', PEAK.peakNumber, 'has', len(possibleh1), 'contributions,',\
                          'only 10 can be used in XEASY'
                if len(possiblep1) > 11:
                    print 'NOTE: peak', PEAK.peakNumber, 'has', len(possiblep1), 'contributions,',\
                          'only 10 can be used in XEASY'
                if len(possibleh2) > 11:
                    print 'NOTE: peak', PEAK.peakNumber, 'has', len(possiblep2), 'contributions,',\
                          'only 10 can be used in XEASY'
                if len(possiblep2) > 11:
                    print 'NOTE: peak', PEAK.peakNumber, 'has', len(possiblep2), 'contributions,',\
                          'only 10 can be used in XEASY'
                #just use the first 10 contributions,
                #use the protons (set: 2 and 4) to set up the matrix:
                XAP.AddMatrix(PEAK.peakNumber, matrix, self.dimension,\
                              '2', '4',\
                              possibleh1[0:10], possiblep1[0:10],\
                              possibleh2[0:10], possiblep2[0:10])
                
                #write output to .assign file:
                assignhandle.write('# ' + CON.peakNumber + '\n')
                for qqq in range(1, 5):
                    qqq = str(qqq)
                    if het1 == qqq:
                        assignhandle.write(str(len(possibleh1[0:10])))
                        for eachatom in possibleh1[0:10]:
                            assignhandle.write(' ' + eachatom)
                        assignhandle.write('\n')
                    if pro1 == qqq:
                        assignhandle.write(str(len(possiblep1[0:10])))
                        for eachatom in possiblep1[0:10]:
                            assignhandle.write(' ' + eachatom)
                        assignhandle.write('\n')
                    if het2 == qqq:
                        assignhandle.write(str(len(possibleh2[0:10])))                        
                        for eachatom in possibleh2[0:10]:
                            assignhandle.write(' ' + eachatom)
                        assignhandle.write('\n')
                    if pro2 == qqq:
                        assignhandle.write(str(len(possiblep2[0:10])))
                        for eachatom in possiblep2[0:10]:
                            assignhandle.write(' ' + eachatom)
                        assignhandle.write('\n')
                for eachbin in XAP.binlist:
                    eachbin = str(eachbin)
                    eachbin = eachbin[:-1]   #get rid of the 'L' of the long integer
                    assignhandle.write(str(eachbin) + ' ')
                assignhandle.write('\n')
            else:
                #now for the case that you only have one contribution:
                assignhandle.write('# ' + PEAK.peakNumber + '\n')
                for qqq in range(1, 5):
                    qqq= str(qqq)
                    #for the case that there are no atomnumbers, use 0 instead,
                    #I don't know how XEASY copes with it:
                    if PEAK.contributions[0].atomnumberh1:
                       outatomnumberh1 = PEAK.contributions[0].atomnumberh1
                    else:
                        outatomnumberh1 = '0'
                    if PEAK.contributions[0].atomnumberp1:
                       outatomnumberp1 = PEAK.contributions[0].atomnumberp1
                    else:
                        outatomnumberp1 = '0'
                    if PEAK.contributions[0].atomnumberh2:
                       outatomnumberh2 = PEAK.contributions[0].atomnumberh2
                    else:
                        outatomnumberh2 = '0'
                    if PEAK.contributions[0].atomnumberp2:
                       outatomnumberp2 = PEAK.contributions[0].atomnumberp2
                    else:
                        outatomnumberp2 = '0'

                    
                    if het1 == qqq:
                        if PEAK.contributions[0].atomnumberh1 == None:
                            assignhandle.write('0\n')
                        else:
                            assignhandle.write('1 ' + outatomnumberh1 + '\n')
                    if pro1 == qqq:
                        if PEAK.contributions[0].atomnumberp1 == None:
                            assignhandle.write('0\n')
                        else:
                            assignhandle.write('1 ' + outatomnumberp1 + '\n')
                    if het2 == qqq:
                        if PEAK.contributions[0].atomnumberh2 == None:
                            assignhandle.write('0\n')
                        else:
                            assignhandle.write('1 ' + outatomnumberh2 + '\n')
                    if pro2 == qqq:
                        if PEAK.contributions[0].atomnumberp2 == None:
                            assignhandle.write('0\n')
                        else:
                            assignhandle.write('1 ' + outatomnumberp2 + '\n')
                            
                #now write the binary numbers for fully assigned peaks:
                assignhandle.write('-1 -1 -1 -1\n')

        #close all the filehandles:
        peakshandle.close()
        assignhandle.close()

        
     
    def WriteXML2String(self, xmlVersion="1.0", encoding="UTF-8",\
                        dtdFileName="noe1.3.dtd", dtdVersionName="1.3"):
        """
        returns a string containing the NOE data in XML format

        xmlVersion       appears in the XML declaration
        encoding         appears in the XML declaration
        dtdFileName      points to the corresponding DTD
        dtdVersionName   version name/number of the DTD
        """
        def __None2Empty(inObject):
            """converts None to ''"""
            if inObject == None:
                return ''
            return inObject

        N2E = __None2Empty #shortcut
        
        #beginning and end:
        getSpectrumName = N2E(self.peakslist[0].contributions[0].spectrumName)
        beginSection = """<?xml version="%s" encoding="%s"?>
<!DOCTYPE aria_noe:list SYSTEM "%s">
<aria_noe:list>
  <aria_noe:version>%s</aria_noe:version>
  <aria_noe:spectrum>
    <aria_noe:spectrum_name>%s</aria_noe:spectrum_name>
""" % (N2E(xmlVersion), N2E(encoding),\
       N2E(dtdFileName), N2E(dtdVersionName), getSpectrumName)
        outString = beginSection
        endPeak = """    </aria_noe:peak>
"""
        endList = """  </aria_noe:spectrum>
</aria_noe:list>
"""
        
        for eachP in self.peakslist:
            beginString = """    <aria_noe:peak>
      <aria_noe:peak_number>%s</aria_noe:peak_number>
      <aria_noe:hetero1_ppm>%s</aria_noe:hetero1_ppm>
      <aria_noe:hetero1_ppm_error>%s</aria_noe:hetero1_ppm_error>
      <aria_noe:proton1_ppm>%s</aria_noe:proton1_ppm>
      <aria_noe:proton1_ppm_error>%s</aria_noe:proton1_ppm_error>
      <aria_noe:hetero2_ppm>%s</aria_noe:hetero2_ppm>
      <aria_noe:hetero2_ppm_error>%s</aria_noe:hetero2_ppm_error>
      <aria_noe:proton2_ppm>%s</aria_noe:proton2_ppm>
      <aria_noe:proton2_ppm_error>%s</aria_noe:proton2_ppm_error>
      <aria_noe:volume>%s</aria_noe:volume>
      <aria_noe:volume_error>%s</aria_noe:volume_error>
      <aria_noe:intensity>%s</aria_noe:intensity>
      <aria_noe:intensity_error>%s</aria_noe:intensity_error>
      <aria_noe:aria_peak_parameters>
        <aria_noe:peak_type>%s</aria_noe:peak_type>
        <aria_noe:initial_weight>%s</aria_noe:initial_weight>
        <aria_noe:current_weight>%s</aria_noe:current_weight>
        <aria_noe:lower_bound>%s</aria_noe:lower_bound>
        <aria_noe:upper_bound>%s</aria_noe:upper_bound>
        <aria_noe:calculated_volume>%s</aria_noe:calculated_volume>
        <aria_noe:calculated_volume_std>%s</aria_noe:calculated_volume_std>
        <aria_noe:test_set_flag>%s</aria_noe:test_set_flag>
        <aria_noe:average_figure_of_merit>%s</aria_noe:average_figure_of_merit>
        <aria_noe:average_summed_distance>%s</aria_noe:average_summed_distance>
        <aria_noe:average_summed_distance_std>%s</aria_noe:average_summed_distance_std>
      </aria_noe:aria_peak_parameters>
""" % (N2E(eachP.contributions[0].peakNumber), N2E(eachP.contributions[0].h1ppm),\
       N2E(eachP.contributions[0].dh1ppm), N2E(eachP.contributions[0].p1ppm), N2E(eachP.contributions[0].dp1ppm),\
       N2E(eachP.contributions[0].h2ppm), N2E(eachP.contributions[0].dh2ppm), N2E(eachP.contributions[0].p2ppm),\
       N2E(eachP.contributions[0].dp2ppm), N2E(eachP.contributions[0].volume), N2E(eachP.contributions[0].volumeError),\
       N2E(eachP.contributions[0].intensity), N2E(eachP.contributions[0].intensityError),\
       N2E(eachP.contributions[0].peakType), N2E(eachP.contributions[0].inWeight), N2E(eachP.contributions[0].curWeight),\
       N2E(eachP.contributions[0].lowerBound), N2E(eachP.contributions[0].upperBound),\
       N2E(eachP.contributions[0].backVolumeAria), N2E(eachP.contributions[0].backVolumeAriaStd),\
       N2E(eachP.contributions[0].testSet), N2E(eachP.contributions[0].fomAria), N2E(eachP.contributions[0].sumDistanceAve),\
       N2E(eachP.contributions[0].sumDistanceStd))
            outString = outString + beginString
            for eachC in eachP.contributions:
                contributionString = "      <aria_noe:contribution>\n"
                proton1 = eachC.atomnamep1[0]
                contributionString = contributionString + """        <aria_noe:atom_list>
          <aria_noe:atom>
            <aria_noe:segment_name>%s</aria_noe:segment_name>
            <aria_noe:residue_number>%s</aria_noe:residue_number>
            <aria_noe:hetero_atom_name>%s</aria_noe:hetero_atom_name>
            <aria_noe:proton_atom_name>%s</aria_noe:proton_atom_name>
          </aria_noe:atom>
        </aria_noe:atom_list>""" % (N2E(eachC.segid1), N2E(eachC.residue1),\
                                  N2E(eachC.atomnameh1[0]), N2E(proton1))
                proton2 = eachC.atomnamep2[0]
                contributionString = contributionString + """
        <aria_noe:atom_list>
          <aria_noe:atom>
            <aria_noe:segment_name>%s</aria_noe:segment_name>
            <aria_noe:residue_number>%s</aria_noe:residue_number>
            <aria_noe:hetero_atom_name>%s</aria_noe:hetero_atom_name>
            <aria_noe:proton_atom_name>%s</aria_noe:proton_atom_name>
          </aria_noe:atom>
        </aria_noe:atom_list>""" % (N2E(eachC.segid2), N2E(eachC.residue2),\
                                  N2E(eachC.atomnameh2[0]), N2E(proton2))
                contributionString = contributionString + """
        <aria_noe:aria_contribution_parameters>
          <aria_noe:figure_of_merit>%s</aria_noe:figure_of_merit>
          <aria_noe:contribution_percentage>%s</aria_noe:contribution_percentage>
          <aria_noe:distance>%s</aria_noe:distance>
          <aria_noe:distance_std>%s</aria_noe:distance_std>
          <aria_noe:spindiffusion_corrected_distance>%s</aria_noe:spindiffusion_corrected_distance>
        </aria_noe:aria_contribution_parameters>  
      </aria_noe:contribution>
""" %  (N2E(eachC.fomContribution), N2E(eachC.contribution), N2E(eachC.distanceAve),\
        N2E(eachC.distanceStd), N2E(eachC.backVolume))
                outString = outString + contributionString
            outString = outString + endPeak
        
        return outString + endList

    
    def WriteXML2File(self, fileName):
        """writes the NOE data to an XML file"""
        outString = self.WriteXML2String()
        try:
            outhandle = TextFile.TextFile(fileName, 'w')
        except IOError:
            print 'could not open the file', fileName
            print 'Abort WriteXML2File method.'
            print '-'*60
            traceback.print_exc(file=sys.stdout)
            print '-'*60
            return
        print 'writing to the file:', fileName
        outhandle.write(outString)
        outhandle.close()

        
    
    def _ReadOldList(self, listfile):
        """
        reads in an old ARIA .list file and writes a list of peak
        objects to the attribute .peakslist
        """
        #get a clean list:
        self.peakslist = []
        NP = NoePeak()
        CON = NoeContribution()
        oldpeakNumber = 0
        withsegid = 0      #if 0, segid is empty
        #set some variables:
        segid1 = None
        segid2 = None
        atomnameh1 = (None,)
        atomnameh2 = (None,)
        atomnamep1 = []
        atomnamep2 = []
        volumeError = None
        distanceError = None
        for line in TextFile.TextFile(listfile):
            #for the negative volumes, add one space character:
            #add spaces for heteronucleus - this is a ugly hack!!!
            line = line[:54] + ' ' + line[54:]
            line = line[:77] + ' ' + line[77:]
            line = line[:87] + ' ' + line[87:]

            linelist = string.split(line)

            #for non-data lines:
            if len(linelist) < 5:
                continue

            #for non-data lines:
            sign = linelist[0]
            if sign != '+' and sign != None:
                continue   #throws out all non-data lines

            #for the segid field:
            try:
                testatoi = string.atoi(linelist[2])
            except ValueError:
                withsegid = 1
                
            if withsegid == 0:
                peakNumber = linelist[1]
                residue1 = linelist[2]
                aa1 = linelist[3]
                atomnamep1 = (linelist[4], )
                residue2 = linelist[5]
                aa2 = linelist[6]
                atomnamep2 = (linelist[7], )
                distance = linelist[8]
                volume = linelist[9]
                lowerBound = linelist[10]
                upperBound = linelist[11]
                hppm1 = linelist[12]
                hppm2 = linelist[13]
                pppm1 = linelist[14]
                pppm2 = linelist[15]
                contribution = linelist[16]
                weight = linelist[17]
                nontrivial = '1' #linelist[18]
            else:
                peakNumber = linelist[1]
                segid1 = linelist[2]
                residue1 = linelist[3]
                aa1 = linelist[4]
                atomnamep1 = (linelist[5], )
                segid2 = linelist[6]
                residue2 = linelist[7]
                aa2 = linelist[8]
                atomnamep2 = (linelist[9], )
                distance = linelist[10]
                volume = linelist[11]
                lowerBound = linelist[12]
                upperBound = linelist[13]
                hppm1 = linelist[14]
                hppm2 = linelist[15]
                pppm1 = linelist[16]
                pppm2 = linelist[17]
                contribution = linelist[18]
                weight = linelist[19]
                nontrivial = linelist[20]

            oldpeakNumber = CON.peakNumber
            
            #convert +/- to peakTypes:
            if string.atoi(nontrivial) > 1:
                if sign == '+':
                    peakType = '3'
                else:
                    peakType = '4'
            else:
                if sign == '+':
                    peakType = '1'
                else:
                    peakType = '2'
                    
            CON = NoeContribution(peakNumber = peakNumber,\
                                  peakType = peakType,\
                                  residue1 = residue1,\
                                  aa1 = aa1,\
                                  segid1 = segid1,\
                                  atomnameh1 = atomnameh1,\
                                  atomnamep1 = atomnamep1,\
                                  residue2 = residue2,\
                                  aa2 = aa2,\
                                  segid2 = segid2,\
                                  atomnameh2 = atomnameh2,\
                                  atomnamep2 = atomnamep2,\
                                  volume = volume,\
                                  volumeError = volumeError,\
                                  distanceAve = distance,\
                                  distanceStd = distanceError,\
                                  lowerBound = lowerBound,\
                                  upperBound = upperBound,\
                                  pppm1 = pppm1,\
                                  pppm2 = pppm2,\
                                  inWeight = weight,\
                                  contribution = contribution,\
                                  nta = nontrivial)

            #add contribution to peak:
            if oldpeakNumber == peakNumber:
                NP.AddContribution(CON)
            else:
                #the output:
                self.AddPeak(NP)
                #create a new instance of a NOE peak:
                NP = NoePeak()
                NP.AddContribution(CON)
        #for the last peak:
        self.AddPeak(NP)
        
###############################################################################
class NoePeak:
    """
    Represents one NOE peak with all possible contributions (like several
    lines for one peak in an ARIA .list file)
    
    attributes:
        contributions    a list of NoeContribution objects
        peakNumber       peakNumber (e.g. from XEASY)
        spectrumName     spectrum name
    methods:
        AddContribution  adds a contribution object to the peak
    """
    def __init__(self):
        self.contributions = []
        self.peakNumber = None
        self.spectrumName = None
    
    def AddContribution(self, newContribution):
        self.contributions.append(newContribution)
        if self.peakNumber == None:
            self.peakNumber = newContribution.peakNumber
        if self.spectrumName == None:
            self.spectrumName = newContribution.spectrumName

    
###############################################################################
class NoeContribution:
    """
    Represents one contribution to an NOE peak (like one line in an ARIA
    .list file).
    This object just contains all the data as attributes, there are no
    methods available.    
    """
    def __init__(self,\
                 spectrumName = None,\
                 peakNumber = None,\
                 testSet = None,\
                 inWeight = None,\
                 volume = None,\
                 volumeError = None,\
                 intensity = None,\
                 intensityError = None,\
                 h1ppm = None,\
                 dh1ppm = None,\
                 p1ppm = None,\
                 dp1ppm = None,\
                 h2ppm = None,\
                 dh2ppm = None,\
                 p2ppm = None,\
                 dp2ppm = None,\
                 peakType = None,\
                 fomAria = None,\
                 curWeight = None,\
                 lowerBound = None,\
                 upperBound = None,\
                 sumDistanceAve = None,\
                 sumDistanceStd = None,\
                 backVolumeAria = None,\
                 backVolumeAriaStd = None,\
                 allAssi = None,\
                 nta = None,\
                 contribution = None,\
                 fomContribution = None,\
                 distanceAve = None,\
                 distanceStd = None,\
                 backVolume = None,\
                 backVolumeStd = None,\
                 segid1 = None,\
                 residue1 = None,\
                 aa1 = None,\
                 atomnameh1 = (None, ),\
                 assih1ppm = None,\
                 assidh1ppm = None,\
                 atomnamep1 = (None, ),\
                 assip1ppm = None,\
                 assidp1ppm = None,\
                 segid2 = None,\
                 residue2 = None,\
                 aa2 = None,\
                 atomnameh2 = (None, ),\
                 assih2ppm = None,\
                 assidh2ppm = None,\
                 atomnamep2 = (None, ),\
                 assip2ppm = None,\
                 assidp2ppm = None,\
                 xeasyudst = 'T',\
                 xeasyinte = None,\
                 atomnumberh1 = None,\
                 atomnumberp1 = None,\
                 atomnumberh2 = None,\
                 atomnumberp2 = None,\
                 sparkylittlecomment= None,\
                 comment = None):

## the old constructor contained:                 
##                 peakNumber = None, peakType = None,\
##                 residue1 = None, aa1 = None, segid1 = None,\
##                 atomnameh1 = (None, ), atomnamep1 = (None, ),\
##                 residue2 = None, aa2 = None, segid2 = None,\
##                 atomnameh2 = (None, ), atomnamep2 = (None, ),\
##                 volume = None, volumeError = None,\
##                 distance = None, distanceError = None,\
##                 lowerBound = None, upperBound = None,\
##                 h1ppm = None, dh1ppm = None, p1ppm = None, dp1ppm = None,\
##                 h2ppm = None, dh2ppm = None, p2ppm = None, dp2ppm = None,\
##                 weight = None, contribution = None, nta = None, comment = None,\
##                 ah1ppm = None, dah1ppm = None, ap1ppm = None, dap1ppm = None,\
##                 ah2ppm = None, dah2ppm = None, ap2ppm = None, dap2ppm = None,\
##                 xeasyudst = 'T', xeasyinte = None):

        #peak line:
        self.spectrumName = spectrumName #spectrum name (not longer than 4 characters)
        self.peakNumber = peakNumber     #peak number
        self.testSet = testSet           #either 1 or 0 (1=used for the calculation, 0=not used)
        self.inWeight = inWeight         #initial weight of the peak
        self.volume = volume             #NOE volume
        self.volumeError = volumeError   #NOE volume error
        self.intensity = intensity       #NOE intensity
        self.intensityError = intensityError       #NOE intensity error
        self.h1ppm = h1ppm               #heteronucleus1 chemical shift
        self.dh1ppm = dh1ppm             #delta of heteronucleus1 chemical shift
        self.p1ppm = p1ppm               #proton1 chemical shift
        self.dp1ppm = dp1ppm             #delta of proton1 chemical shift
        self.h2ppm = h2ppm               #heteronucleus2 chemical shift
        self.dh2ppm = dh2ppm             #delta of heteronucleus2 chemical shift
        self.p2ppm = p2ppm               #proton2 chemical shift
        self.dp2ppm = dp2ppm             #delta of proton2 chemical shift
        
        #aria line:
        self.peakType = peakType         #colour coding in the range [1, 6]
        self.fomAria = fomAria           #figure of merit
        self.curWeight = curWeight       #current weight
        self.lowerBound = lowerBound               #lower distance bound
        self.upperBound = upperBound               #upper distance bound
        self.sumDistanceAve = sumDistanceAve       #sum of all distances
        self.sumDistanceStd = sumDistanceStd       #std of the sum of all distances
        self.backVolumeAria = backVolumeAria       #backcalculated volume
        self.backVolumeAriaStd = backVolumeAriaStd #std of the backcalculated volume
        self.allAssi = allAssi           #number of all assignments
        self.nta = nta                   #number of non-trivial assignments
        
        #contribution line:
        self.fomContribution = fomContribution     #figure of merit
        self.contribution = contribution #contribution to the whole peak volume
        self.distanceAve = distanceAve   #average distance
        self.distanceStd = distanceStd   #standard deviation of the sum of the average distance
        self.backVolume = backVolume     #back calculated volume
        self.backVolumeStd = backVolumeStd #standard deviation of the backcalculated volume
        self.segid1 = segid1             #segid 1
        self.residue1 = residue1         #residue number 1
        self.aa1 = aa1                   #three-letter code 1
        self.atomnameh1 = atomnameh1     #heteronucleus atomname 1
        self.assih1ppm = assih1ppm       #ppm of heteronucleus 1
        self.assidh1ppm = assidh1ppm     #delta of ppm of heteronucleus 1
        self.atomnamep1 = atomnamep1     #atomname of proton 1
        self.assip1ppm = assip1ppm       #ppm of proton 1
        self.assidp1ppm = assidp1ppm     #delta of ppm of proton 1
        self.segid2 = segid2             #segid 2
        self.residue2 = residue2         #residue number 2
        self.aa2 = aa2                   #three-letter code 2
        self.atomnameh2 = atomnameh2     #heteronucleus atomname 2
        self.assih2ppm = assih2ppm       #ppm of heteronucleus 2
        self.assidh2ppm = assidh2ppm     #delta of ppm of heteronucleus 2
        self.atomnamep2 = atomnamep2     #atomname of proton 2
        self.assip2ppm = assip2ppm       #ppm of proton 2
        self.assidp2ppm = assidp2ppm     #delta of ppm of proton 2
        
        #xeasy specific stuff:
        self.xeasyudst = xeasyudst       #user-defined spectrum type
        self.xeasyinte = xeasyinte       #integration method
        self.atomnumberh1 = atomnumberh1 #XEASY atomnumber hetero 1
        self.atomnumberp1 = atomnumberp1 #XEASY atomnumber proton 1
        self.atomnumberh2 = atomnumberh2 #XEASY atomnumber hetero 2
        self.atomnumberp2 = atomnumberp2 #XEASY atomnumber proton 2

        #comments:
        self.comment = comment           #each contribution may be commented!
        self.sparkylittlecomment = sparkylittlecomment
###############################################################################
# Below there are some private classes and functions                          #
###############################################################################

def _DoesFileExist(fileName):
    if os.path.exists(fileName) == 0:
        print 'WARNING:', fileName, 'does not exist.'
        return 0
    return 1


###############################################################################
#test code:
if __name__ == "__main__":
    print 'testing module...'
    NC1 = NoeContribution('test spectrum', '1234', '1', '0.9', '123456789',\
                          '111', '1.2345', '0.00123', '1.111', '0.001', '2.222',\
                          '0.002', '3.333', '0.003', '4.444', '0.004', '6',\
                          '1.98475', '0.5134', '1.0', '5.0', '3.0', '0.1',\
                          '1134543', '123', '1', '1', '1', '0.2456', '3.62',\
                          '0.341', '156346', '12', 'QWER', '62', 'ARG', ('CA',),\
                          '1.123', '0.001', ('HA',), '2.234', '0.002',\
                          'BKG', '63', 'LEU', ('CA',),\
                          '1.414', '0.005', ('HA',), '3.34455', '0.05',\
                          ('CB',), ('HB2',))
    NC2 = NoeContribution('test spectrum', '1234', '1', '0.9', '123456789',\
                          '111', '1.2345', '0.00123', '1.111', '0.001', '2.222',\
                          '0.002', '3.333', '0.003', '4.444', '0.004', '6',\
                          '1.98475', '0.5134', '1.0', '5.0', '3.0', '0.1',\
                          '1134543', '123', '1', '1', '1', '0.2456', '3.62',\
                          '0.341', '156346', '12', 'QWER', '62', 'ARG', ('CB',),\
                          '1.123', '0.001', ('HB2','HB3'), '2.234', '0.002',\
                          'BKG', '63', 'LEU', ('CA',),\
                          '1.414', '0.005', ('HA',), '3.34455', '0.05',\
                          ('CB',), ('HB2',))
    NP = NoePeak()
    NP.AddContribution(NC1)
    NP.AddContribution(NC2)
    NL = NoeList()
    NL.dimension = 4
    NL.name = 'testing - 1 peak with 2 identical contributions'
    NL.AddPeak(NP)
    
    print '\nThe Stdout() method prints:\n'
    NL.Stdout()
    print '\n\nThe XML output is:\n'
    print NL.WriteXML2String()
    print '\ntest done. bye.'


