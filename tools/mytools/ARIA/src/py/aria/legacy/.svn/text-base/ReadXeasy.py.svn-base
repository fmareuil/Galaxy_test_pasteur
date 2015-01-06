"""
The ReadXeasy module contains three classes to read in XEASY .assign,
.peaks and .prot files and supply the data in objects which can easily
be accessed by other modules, e.g. this module is used by NoeList.
"""
__author__   = "$Author: bardiaux $"
__revision__ = "$Revision: 1.1.1.1 $"
__date__     = "$Date: 2010/03/23 15:27:24 $"

import os, re, string
import numpy

import PseudoAtom
import DictWithDefault, TextFile

###############################################################################
class XeasyAssign:
    """
    reads in all the assignments of a XEASY .assign file

    this class basically provides a list of peaks with
    possible assignments
    
    attributes:
        assignments    list of AssignPeak objects
        assipeaksdic   (key: peakNumber, value: AssignPeak object) 
        comment
        filename       name of the .assignfile
        first          column number of the first dimension of the matrix
        second         column number of the second dimension of the matrix
        name
        
    methods:
        AddAssign      adds an assignment
        ReadAssign     reads an XEASY .assign file
        WriteAssign    writes an XEASY .assign file
        Stdout         writes the data to stdout
        Test           test for matrix -> binary format -> matrix

    all the data are strings, except 'dimension' which is a integer!
    """
    def __init__(self):
        self.assignments = []
        self.assipeaksdic = {}

    def AddAssign(self, object):
        self.assignments.append(object)
        self.assipeaksdic[object.peakNumber] = object
        
    
    def ReadAssign(self, assignfile, first, second):
        """
        reads a XEASY .assign file
        
        use proton1 for first
        use proton2 for second
        """
        self.filename = assignfile

        self.first = first
        self.second = second

        # look for the .assign file and open it:
        _DoesFileExist(assignfile)
        assignhandle = TextFile.TextFile(assignfile, 'r')

        # get the dimensionality:
        try:
            dimension = string.atoi(string.strip(assignhandle.readline()))
        except:
            AriaError = 'first line in .assign has to contain the dimension!'
            raise AriaError
        if dimension != 2 and dimension != 3 and dimension != 4:
            AriaError = 'DimensionError:  No dimension specified'
            raise AriaError
        print 'dimension:', dimension, ' (read from', assignfile + ')'

        # loop over the rest of the file:
        while 1:
            #read a line:
            line = string.strip(assignhandle.readline())
            #end of loop condition:
            if string.strip(line) == '':
                break
            #get a list of fields:
            linelist = string.split(line)
            #initializing the lists of ambiguous assignments (atomnumbers):
            possible1 = []
            possible2 = []
            possible3 = []
            possible4 = []
            #get the lines with # to read the peakNumber:
            if linelist[0] == '#':
                peakNumber = linelist[1]
                #get the 2, 3 or 4 lines for the ambiguous assignments:
                for xxx in range(dimension):
                    linelist = string.split(assignhandle.readline())
                    for yyy in range(string.atoi(linelist[0])):
                        exec 'possible' + str(xxx + 1) +\
                             '.append(linelist[' + str(yyy + 1) + '])'
##                 print possible1, #test
##                 print possible2, #test
##                 print possible3, #test
##                 print possible4, #test
                #get the four binary numbers:
                binlist = []
                linelist = string.split(assignhandle.readline())
                for bbb in range(4):
                    binlist.append(long(linelist[bbb]))
##                 print 'binlist: ', binlist #test
                #create an instance of XeasyAssignPeak for the binary numbers:
                XAP = XeasyAssignPeak()
                XAP.AddBin(peakNumber, binlist, dimension, first, second,\
                           possible1, possible2, possible3, possible4)
                #add this objects to the list in the assignments attribute:
                self.AddAssign(XAP)
        

    def Stdout(self):
        """
        writes out all the data from the .assign file to stdout
        """
        print 'write all the data from the XEASY.assign class to stdout:'
        print 'peakNumber dimension binlist first second 1 2 3 4 matrix'
        for AP in self.assignments:
            print AP.peakNumber, AP.dimension, AP.binlist, AP.first,\
                  AP.second, AP.possible1, AP.possible2, AP.possible3,\
                  AP.possible4, '\n', AP.assimatrix
        

    def Test(self):
        """just a test: matrix -> binary format -> matrix"""
        matrix = numpy.zeros((4, 5))
        matrix[1,3] = 1
        dimension = 3
        first = '1'
        second = '2'
        possible1 = [1, 2, 3, 4]
        possible2 = [11, 12, 13, 14, 15]
        possible3 = [21, 22, 23, 24, 25, 26, 27, 28]
        possible4 = [31, 32, 33, 34, 35, 36, 37, 38, 39]
        print 'Test: matrix -> binary format'
        AP = XeasyAssignPeak()
        AP.AddMatrix(123456, matrix, dimension,\
                     first, second,\
                     possible1, possible2, possible3, possible4)
        matrixandbintest = AP.binlist
        print AP.assimatrix
        print AP.binlist
        print 'Test: binary format -> matrix'
        AP = XeasyAssignPeak()
        AP.AddBin(123456, matrixandbintest, dimension,\
                  first, second,\
                  possible1, possible2, possible3, possible4)
        print AP.binlist
        print AP.assimatrix


    def WriteAssign(self, fileName):
        outHandle = open(fileName, 'w')
        if len(self.assignments) > 0:
            outHandle.write(str(self.assignments[0].dimension)+'\n')
        for eachAP in self.assignments:
            outHandle.write('# ' + eachAP.peakNumber+'\n')
            outHandle.write(str(len(eachAP.possible1)) + ' ' + string.join(eachAP.possible1)+'\n')
            outHandle.write(str(len(eachAP.possible2)) + ' ' + string.join(eachAP.possible2)+'\n')
            if eachAP.dimension == 3:
                outHandle.write(str(len(eachAP.possible3)) + ' ' + string.join(eachAP.possible3)+'\n')
            if eachAP.dimension == 4:
                outHandle.write(str(len(eachAP.possible4)) + ' ' + string.join(eachAP.possible4)+'\n')
            outHandle.write(str(eachAP.binlist[0])[:-1] + ' ' + str(eachAP.binlist[1])[:-1] + ' ' +\
                            str(eachAP.binlist[2])[:-1] + ' ' + str(eachAP.binlist[3])[:-1] + ' \n')
                           

###############################################################################    
class XeasyAssignPeak:
    """
    Represents a peak with all possible assignments (as given by an XEASY
    .assign file


    IMPORTANT: XEASY uses 32bit integers for the four assignment numbers
               in the .assign file
               

    attributes:
        peakNumber        XEASY peakNumber
        dimension         dimensionality of the spectrum
        assimatrix
        binlist           LIST of the 4 XEASY integers
        first
        second
        possible1         LIST of atomnumbers = possible assignments in w1
        possible2         LIST of atomnumbers = possible assignments in w2
        possible3         LIST of atomnumbers = possible assignments in w3
        possible4         LIST of atomnumbers = possible assignments in w4
    methods:
        AddBin            add possible assignments using binary format
        AddMatrix         add possible assignments using matrix format

        """
    def __init__(self):
        self.assimatrix = numpy.ones((10, 10))
        
    
    def AddBin(self, peakNumber, binlist, dimension, first, second,\
               possible1, possible2, possible3, possible4):
        """
        uses the rows with the numbers specified by
        'first' and 'second'
        sets up the matrix with the specified lists 'possiblex'
        
        NOTE: binlist has to be a list of long integers
        
        """
        self.peakNumber = peakNumber
        self.binlist = binlist
        self.dimension = dimension
        
        self.possible1 = possible1
        self.possible2 = possible2
        self.possible3 = possible3
        self.possible4 = possible4
        self.plusone = [0L, 0L, 0L, 0L]
        #for the correct setup of the matrix:
        if string.atoi(first) > string.atoi(second):
            self.first = second
            self.second = first
        else:
            self.first = first
            self.second = second
        #add one to each long integer in binlist:
        for counter in range(0, 4):
            self.plusone[counter] = self.binlist[counter] + 1
        #get the right 'possible' lists:
        exec 'possiblea = possible' + str(self.first)
        exec 'possibleb = possible' + str(self.second)
        #set up the empty matrix:
        self.assimatrix = numpy.ones((len(possiblea), len(possibleb)))
##         print 'matrix a b', len(possiblea), len(possibleb)  #test
        #now the dirty part for the calculation of the matrix:
        thirtytwo = range(32)
        thirtytwo.reverse()
        for kkk in range(0, 4):
            for expo in thirtytwo:
                if self.plusone[kkk] == 0:
                    continue
                #if-clause for the integer overflow, create negative values:
                if self.plusone[kkk] > 0:
                    self.plusone[kkk] = - (2147483648L - self.plusone[kkk])\
                                        - 2147483648L
##                    print 'simulate the overflow:', self.plusone[kkk] #test
                if self.plusone[kkk] <= -(2L ** expo):
                    self.plusone[kkk] = self.plusone[kkk] + 2L ** expo
                    totalexpo = 32 * kkk + expo
                    aaa = totalexpo % 10
                    bbb = totalexpo / 10
##                     print aaa,bbb #test
                    #0 means not assigned (button not pressed):
                    try:
                        self.assimatrix[aaa, bbb] = 0
#                        print '#########'
#                        print self.assimatrix, aaa, bbb
                    except:
                        print 'WARNING: problem with assignment matrix!'
                        print self.assimatrix, aaa, bbb

##         print self.assimatrix  #test
        
    
    def AddMatrix(self, peakNumber, assimatrix, dimension,\
                  first, second,\
                  possible1, possible2, possible3, possible4):
        """
        adds a matrix and calculates the binary numbers
        """
        self.peakNumber = peakNumber
        self.assimatrix = assimatrix
        self.dimension = dimension
        self.possible1 = possible1
        self.possible2 = possible2
        self.possible3 = possible3
        self.possible4 = possible4
        #for the correct setup of the matrix:
        if string.atoi(first) > string.atoi(second):
            self.first = second
            self.second = first
        else:
            self.first = first
            self.second = second
        #use long integers, be careful with 32 bit:
        self.binlist = [-1L, -1L, -1L, -1L]
##         print self.assimatrix   #test
        
        #get the right 'possible' lists:
        exec 'possiblea = possible' + str(self.first)
        exec 'possibleb = possible' + str(self.second)
        
        #now comes the dirty part for calculating the binary numbers:
        for aaa in range(0, len(possiblea)):
            for bbb in range(0, len(possibleb)):
                if self.assimatrix[aaa, bbb] == 0:  #1 means assigned
                    which = 10 * bbb + aaa
                    if which < 32:
                        kkk = 0
                    elif which < 64:
                        kkk = 1
                    elif which < 96:
                        kkk = 2
                    elif which < 128:
                        kkk = 3
                    exponent = which - 32 * kkk
                    # now come 4 tricky lines to overcome the integer overflow
                    # of XEASY 32bit integers, just simulate the overflow:
                    if self.binlist[kkk] - 2L ** exponent > -2147483469L:
                        self.binlist[kkk] = self.binlist[kkk] - 2L ** exponent
                    else:
                        self.binlist[kkk] = self.binlist[kkk] + 2L ** exponent
##         print self.binlist  #test
        

###############################################################################
class XeasyPeaks:
    """
    contains the whole data of a Xeasy .peaks file
    
    attributes:
        dimension    dimension (read from the first line of the .peaks file)
        iname1       INAME1 (from the header)
        iname2       INAME2 (from the header)
        iname3       INAME3 (from the header)
        iname4       INAME4 (from the header)
        peakslist    list of XeasyPeak objects
        peaksdic     dictionary of XeasyPeak objects (key: peakNumber,
                     value: XeasyPeak objects)
        
    public methods:
        AddPeak      adds one peak
        ReadPeaks    reads a Xeasy .peaks file
        WritePeaks   write a Xeasy .peaks file
        Stdout       writes the whole list to stdout        
    """
    def __init__(self):
        self.dimension = 0
        self.iname1 = '?'  #default
        self.iname2 = '?'  #default
        self.iname3 = '?'  #default
        self.iname4 = '?'  #default
        self.peakslist = []
        self.peaksdic = {}
        
     
    def AddPeak(self, pobj):
        self.peakslist.append(pobj)
        self.peaksdic[pobj.peakNumber] = pobj
        
    
    def ReadPeaks(self, peaksfile):
        self.peaksfile = peaksfile
        peakshandle = TextFile.TextFile(self.peaksfile)
        scratch = re.compile('#')
        numberofdim = re.compile('Number of dimension')
        iname1 = re.compile('INAME 1')
        iname2 = re.compile('INAME 2')
        iname3 = re.compile('INAME 3')
        iname4 = re.compile('INAME 4')
        for line in peakshandle:
            peakslist= string.split(line)
            if scratch.match(line):
                if numberofdim.search(line):
                    self.dimension = string.atoi(peakslist[4])
                if iname1.search(line):
                    self.iname1 = peakslist[-1]
                if iname2.search(line):  
                    self.iname2 = peakslist[-1]
                if iname3.search(line):
                    self.iname3 = peakslist[-1]
                if iname4.search(line):
                    self.iname4 = peakslist[-1]
                continue
            XP = XeasyPeaksPeak()
            XP.peakNumber = peakslist[0]
            XP.w1 = peakslist[1]
            XP.w2 = peakslist[2]
            if self.dimension == 0:
                print 'need "Number of dimensions" in .peaks file,',\
                      'ReadPeaks method aborted.'
                return
            if self.dimension == 4:
                XP.w3 = peakslist[3]
                XP.w4 = peakslist[4]
                XP.peakType = peakslist[5]
                XP.spectrumType = peakslist[6]
                XP.volume = peakslist[7]
                XP.volumeError = peakslist[8]
                XP.integration = peakslist[9]
                XP.unused1 = peakslist[10]
                XP.atomnumber1 = peakslist[11]
                XP.atomnumber2 = peakslist[12]
                XP.atomnumber3 = peakslist[13]
                XP.atomnumber4 = peakslist[14]
                XP.unused2 = peakslist[15]
            elif self.dimension == 3:
                XP.w3 = peakslist[3]
                XP.peakType = peakslist[4]
                XP.spectrumType = peakslist[5]
                XP.volume = peakslist[6]
                XP.volumeError = peakslist[7]
                XP.integration = peakslist[8]
                XP.unused1 = peakslist[9]
                XP.atomnumber1 = peakslist[10]
                XP.atomnumber2 = peakslist[11]
                XP.atomnumber3 = peakslist[12]
                XP.unused2 = peakslist[13]
            elif self.dimension == 2:
                XP.peakType = peakslist[3]
                XP.spectrumType = peakslist[4]
                XP.volume = peakslist[5]
                XP.volumeError = peakslist[6]
                XP.integration = peakslist[7]
                XP.unused1 = peakslist[8]
                XP.atomnumber1 = peakslist[9]
                XP.atomnumber2 = peakslist[10]
                XP.unused2 = peakslist[11]
            else:
                NodimensionError = 'No dimension specified'
                raise NodimensionError
            self.AddPeak(XP)
        peakshandle.close()

    def WritePeaks2String(self):
        """returns a string which contains the whole .peaks file"""
        #write some comments (dimension and INAME) for .peaks:
        outString = '# Number of dimensions ' + str(self.dimension) + '\n'
        if self.dimension >= 1:
            outString = outString + '#INAME 1 ' + self.iname1 + '\n'
        if self.dimension >= 2:
            outString = outString + '#INAME 2 ' + self.iname2 + '\n'
        if self.dimension >= 3:
            outString = outString + '#INAME 3 ' + self.iname3 + '\n'
        if self.dimension >= 4:
            outString = outString + '#INAME 4 ' + self.iname4 + '\n'

        #loop over all the peaks:
        for EP in self.peakslist:
            #write out with the correct order:
            if self.dimension == 4:
                lineString = '%4i %7.3f %7.3f %7.3f %7.3f %1i %-9s %10.3e %9.2e %1s %3i %4i %4i %4i %4i %1i\n' %\
                             (string.atoi(EP.peakNumber),\
                              string.atof(EP.w1),\
                              string.atof(EP.w2),\
                              string.atof(EP.w3),\
                              string.atof(EP.w4),\
                              string.atoi(EP.peakType),\
                              EP.spectrumType,\
                              string.atof(EP.volume),\
                              string.atof(EP.volumeError),\
                              EP.integration,\
                              string.atoi(EP.unused1),\
                              string.atoi(EP.atomnumber1),\
                              string.atoi(EP.atomnumber2),\
                              string.atoi(EP.atomnumber3),\
                              string.atoi(EP.atomnumber4),\
                              string.atoi(EP.unused2))
                outString = outString + lineString
            elif self.dimension == 3:
                lineString = '%4i %7.3f %7.3f %7.3f %1i %-9s %10.3e %9.2e %1s %3i %4i %4i %4i %1i\n' %\
                             (string.atoi(EP.peakNumber),\
                              string.atof(EP.w1),\
                              string.atof(EP.w2),\
                              string.atof(EP.w3),\
                              string.atoi(EP.peakType),\
                              EP.spectrumType,\
                              string.atof(EP.volume),\
                              string.atof(EP.volumeError),\
                              EP.integration,
                              string.atoi(EP.unused1),\
                              string.atoi(EP.atomnumber1),\
                              string.atoi(EP.atomnumber2),\
                              string.atoi(EP.atomnumber3),\
                              string.atoi(EP.unused2))
                outString = outString + lineString
            elif self.dimension == 2:
                lineString = '%4i %7.3f %7.3f %1i %-9s %10.3e %9.2e %1s %3i %4i %4i %1i\n' %\
                             (string.atoi(EP.peakNumber),\
                              string.atof(EP.w1),\
                              string.atof(EP.w2),\
                              string.atoi(EP.peakType),\
                              EP.spectrumType,\
                              string.atof(EP.volume),\
                              string.atof(EP.volumeError),\
                              EP.integration,\
                              string.atoi(EP.unused1),\
                              string.atoi(EP.atomnumber1),\
                              string.atoi(EP.atomnumber2),\
                              string.atoi(EP.unused2))
                outString = outString + lineString
        return outString
            

    def WritePeaks(self, fileName):
        """writes a Xeasy .peaks file to a file"""
        try:
            outHandle = open(fileName, 'w')
        except:
            print 'WARNING: could not open', fileName
            return
        outString = self.WritePeaks2String()
        outHandle.write(outString)

    def Stdout(self):
        """writes all data to stdout"""
        if self.peakslist == '':
            print 'no peaks read in'
        else:
            print 'pn w1 w2 w3 w4 pt st vo ve in u1 a1 a2 a3 a4 u2'
            for PEAK in self.peakslist:
                print PEAK.peakNumber, PEAK.w1, PEAK.w2, PEAK.w3, PEAK.w4,\
                      PEAK.peakType, PEAK.spectrumType, PEAK.volume,\
                      PEAK.volumeError, PEAK.integration, PEAK.unused1,\
                      PEAK.atomnumber1, PEAK.atomnumber2, PEAK.atomnumber3,\
                      PEAK.atomnumber4, PEAK.unused2
        
    
###############################################################################
class XeasyPeaksPeak:
    """
    represents one line in a Xeasy .peaks file
    
    This object just contains all the data as attributes, there are no
    methods available.

    
    attributes:
        peakNumber       Xeasy peakNumber
        w1               frequency 1
        w2               frequency 2
        w3               frequency 3
        w4               frequency 4
        peakType         peakType (colour code from 1-6)
        spectrumType     user defined type of the spectrum
        volume           NOE volume
        volumeError      NOE volume error
        integration      integration method
        unused1          unused field 1
        atomnumber1      Xeasy atom number 1
        atomnumber2      Xeasy atom number 2
        atomnumber3      Xeasy atom number 3
        atomnumber4      Xeasy atom number 4
        unused2          unused field 2
    """
    def __init__(self):
        self.peakNumber = '-'
        self.w1 = '-'
        self.w2 = '-'
        self.w3 = '-'
        self.w4 = '-'
        self.peakType = '-'
        self.spectrumType = 'N'   #use N as default!
        self.volume = '-'
        self.volumeError = '-'
        self.integration = '-'
        self.unused1 = 0          #unused number = 0
        self.atomnumber1 = '-'
        self.atomnumber2 = '-'
        self.atomnumber3 = '-'
        self.atomnumber4 = '-'
        self.unused2 = 0          #unused number = 0
        

###############################################################################
class XeasyProt:
    """
    contains the whole data of a Xeasy .prot file
    
    attributes:
        atomlist      list of XeasyProtAtom objects
        atomdican     dictionary (key: atomnumber,
                      value:XeasyProtAtom object)
        atomdicfa     dictionary (key: (fragmentnumber, ariaatomname),
                      value:XeasyProtAtom object)
        atomdicfx     dictionary (key: (fragmentnumber, xeasyatomname),
                      value:XeasyProtAtom object
        atomdicre     dictionary (key: residuenumber, value: list of
                      xeasyatomnames)
                      uses DictWithDefault instead of a simple dictionary
    public methods:
        AddProt       adds one line in a .prot file
        ReadProt      read a Xeasy .prot file
        Stdout        write all data to stdout
    """
    def __init__(self):
        self.atomlist = []
        self.atomdican = {}
        self.atomdicfa = {}
        self.atomdicfx = {}
        self.atomdicre = DictWithDefault.DictWithDefault([])
        
         
    def AddProt(self, pobj):
        self.atomlist.append(pobj)
        self.atomdican[pobj.atomnumber] = pobj
        self.atomdicfa[(pobj.fragmentnumber, pobj.ariaatomname)] = pobj
        self.atomdicfx[(pobj.fragmentnumber, pobj.xeasyatomname)] = pobj
        self.atomdicre[pobj.fragmentnumber].append(pobj.xeasyatomname)
        
    
    def ReadProt(self, protfile, verbose = 0):
        """reads in a Xeasy .prot file"""
        self.protfile = protfile
        prothandle = TextFile.TextFile(self.protfile)
        scratch = re.compile('#')
        for line in prothandle:
            if scratch.match(line):
                continue      #for the comments
            XPA = XeasyProtAtom()
            protlist = string.split(line)
            if len(protlist) < 4:
                continue
            XPA.atomnumber = protlist[0]
            XPA.shift = protlist[1]
            XPA.shifterror = protlist[2]
            XPA.xeasyatomname = string.upper(protlist[3])
            XPA.fragmentnumber = protlist[4]
            XPA.ariaatomname= PseudoAtom.Pseudo2Tuple(XPA.xeasyatomname)
            self.AddProt(XPA)

            if verbose:
                print XPA.fragmentnumber, XPA.xeasyatomname,\
                      XPA.ariaatomname, XPA.atomnumber
            
        prothandle.close()
        

    def Stdout(self):
        """writes all the data of the .prot list to stdout"""
        if self.atomlist == []:
            print 'no atoms read in yet'
        else:
            print 'atomnumber shift shifterror xeasyatomname fragmentnumber',\
                  'ariaatomname'
            for ATOM in self.atomlist:
                print ATOM.atomnumber, ATOM.shift, ATOM.shifterror,\
                      ATOM.xeasyatomname, ATOM.fragmentnumber,\
                      ATOM.ariaatomname


###############################################################################                
class XeasyProtAtom:
    """
    represents the chemical shift assignment of one Xeasy atom (like one
    line in a Xeasy .prot file)
    
    This object just contains all the data as attributes, there are no
    methods available.
    
    attributes:
       atomnumber            Xeasy atom number
       shift                 chemical shift in ppm
       shifterror            deviation from the mean value
       xeasyatomname         Xeasy atom name
       fragmentnumber        Xeasy fragment number (ideally the residue number)
       ariaatomname          TUPLE of the corresponding Aria atom name(s)
    """
    def __init__(self):
        self.atomnumber = '-' 
        self.shift = '-'
        self.shifterror = '-'
        self.xeasyatomname = '-'
        self.fragmentnumber = '-'
        self.ariaatomname = '-'
        

###############################################################################
def _DoesFileExist(filename):
    if os.path.exists(filename) == 0:
        print 'WARNING:', filename, 'does not exist.'
        return 0
    return 1

###############################################################################

