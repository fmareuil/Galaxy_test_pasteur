"""
A module to deal with the sequence information from XPLOR or CNS files.
can also read in files in the DIANA .seq format or .pdb files
uses the AminoAcid and DeleteComments modules
"""
__author__   = "$Author: bardiaux $"
__revision__ = "$Revision: 1.1.1.1 $"
__date__     = "$Date: 2010/03/23 15:27:24 $"

import os, re, string
import AminoAcid
import DeleteComments
import TextFile

###############################################################################
class SequenceList:
    """
    contains a list of aminoacids in 3-letter code notation

    attributes:
        aalist     list of aminoacids in 3-letter code notation
                   (from N- to C-terminus)
        fileName   name of the input sequence file

    methods:
        ReadAnsig  reads a sequence file from ANSIG
        ReadFasta  reads a sequence in FASTA format
        ReadOne    reads a sequence in 1-letter format
        ReadPdb    simple method that reads in a pdb file to get the sequence
        ReadPipp   reads a PIPP .par file (SEQUENCE tag)
        ReadSeq    reads a sequence file (CNS, DIANA or XEASY)
        Stdout     writes the list of aminoacids to stdout (3-letter code)
        StdoutOne  writes the list of aminoacids to stdout (1-letter code,
                   all the stuff  in one line
        WriteFasta writes a FASTA sequence file
        WriteSeq   writes a 'clean' file without comments in 3-letter code
        WriteSeq1  writes a 'clean' file without comments in 1-letter code
        WriteXML2String  returns a string containing the sequence in XML format
        WriteXML2File    writes the sequence to an XML file
    """
    def __init__(self, fileName = ''):
        self.aalist = []
        self.fileName = fileName
        

    def __repr__(self):
        return self.fileName
    
    def ReadAnsig(self, fileName):
        """
        reads in an ANSIG sequence file
        uses DeleteComments
       """
        if os.path.exists(fileName):
            print 'reading sequence', fileName
            self.fileName = fileName
        else:
            print 'the file', fileName
            print 'does not exist. Abort ReadANSIG method.'
            return
        #important: delete the aalist aatribute:
        self.aalist = []
        #get the file without the comments:
        bigString = DeleteComments.GetString(fileName)
        #delete the word 'residue' and the line (non-greedy) with 'sequence'
        #and the line with 'end_sequence':
        resword = re.compile('residue')
        bigString = resword.sub('', bigString)
        bigString = re.sub('end_sequence', '', bigString)
        bigString = re.sub('sequence.*?\012', '', bigString)
#        print bigString #test
        #split the string in lines:
        lines = string.split(bigString, '\n')
        for line in lines:
            if string.strip(line) == '':
                continue
            linelist = string.split(line)
            residuenumber = linelist[0]
            getamino = string.upper(linelist[1])
            #append in to the list attribute:
            self.aalist.append(getamino)

###############################################################################
    def ReadFasta(self, fileName):
        """
        reads FASTA files
        The lines which begin with an '>' are comments,
        the rest is in 1-letter code, e.g.:
>BRC2_HUMAN
MPIGSKERPTFFEIFKTRCNKADLGPISLNWFEELSSEAPPYNSEPAEESEHKNNNYEPN
LFKTPQRKPSYNQLASTPIIFKEQGLTLPLYQSPVKELDKFKLDLGRNVPNSRHKSLRTV
KTKMDQADDVSCPLLNSCLSESPVVLQCTHVTPQRDKSVVCGSLFHTPKFVKGRQTPKHI
SESLGAEVDPDMSWSSSLATPPTLSSTVLIVRNEEASETVFPHDTTANVKSYFSNHDESL
        """
        if os.path.exists(fileName):
            print 'reading FASTA file', fileName
            self.fileName = fileName
        else:
            print 'the file', fileName
            print 'does not exist. Abort ReadFasta method.'
            return
        #important: delete the aalist aatribute!
        self.aalist = []
        self.fileName = fileName

        greaterThan = re.compile('>')
        for line in TextFile.TextFile(fileName):
            #FASTA comments start with an '>':
            if greaterThan.match(line):
                continue
            for iii in range(len(line)):
                oneLetter = string.strip(line[iii])
                if oneLetter:
                    self.aalist.append(AminoAcid.AminoAcid(oneLetter)[1])
    
###############################################################################    
    def ReadOne(self, fileName):
        """
        reads a file in 1-letter code. e.g.
        MPIGSKERPTFFEIFKTRCNKADLGPISLNWFEELSSEAPPYNSEPAEESEHKNNNYEPN

        uses DeleteComments to get rid of standard comments like
        '# bla' or '! bla'
        
        """
        if os.path.exists(fileName):
            print 'reading sequence file (one-letter format)', fileName
            self.fileName = fileName
        else:
            print 'the file', fileName
            print 'does not exist. Abort ReadOne method.'
            return
        #important: delete the aalist aatribute!
        self.aalist = []
        self.fileName = fileName

        #get the file without the comments:
        bigString = DeleteComments.GetString(fileName)
        
        for iii in range(len(bigString)):
            oneLetter = string.strip(bigString[iii])
            if oneLetter:
                self.aalist.append(AminoAcid.AminoAcid(oneLetter)[1])
                  
    
###############################################################################
    def ReadPdb(self, pdbfile):
        """
        simple method that reads in a pdb file to get the sequence
        it just reads all the lines which begin with 'ATOM'
        all the other lines are irrelevant!
        """
        if os.path.exists(pdbfile):
            print 'reading pdb file', pdbfile
            self.fileName = pdbfile
        else:
            print 'the file', pdbfile
            print 'does not exist. Abort ReadPdb method.'
            return
        #important: delete the aalist aatribute!
        self.aalist = []
        self.fileName = pdbfile
        #initialize:
        seqlist = {}
        aanumber = ''
        aatype = ''
        #look for all the lines which begin with 'ATOM':
        atom = re.compile('ATOM')
        for line in TextFile.TextFile(pdbfile):
            #ENDMDL means: the sequence is finished:
            if re.search('ENDMDL', line):
                break
            #get the lines with ATOM:
            if atom.match(line):
                allfields = string.split(line)
                aanumber = allfields[4];       #get the amino acid number   
                aatype = allfields[3];         #get the amino acid type
                if seqlist.has_key(aanumber):  #uses hash table (see below)
                    pass
                else:
                    self.aalist.append(aatype)
                    seqlist[str(aanumber)] = aatype;  #defines hash table
        
    
###############################################################################
    def ReadPipp(self, fileName):
        """
        reads an PIPP stapp.par file
        reads in the 1-letter code after the SEQUENCE tag:
        SEQUENCE        GGGGGGGGGG AAAAAAAAAA \
                        VVVVVVVVVV LLLLLLLLLL \
                        IIIIIIIIII
        stops after the first linebreak withouf a '\' in front of '\n'
       """
        if _DoesFileExist(fileName) == 0:
            return
        print 'reading a PIPP sequence file', fileName
        #important: delete the aalist aatribute:
        self.aalist = []
        fileHandle = open(fileName)
        startNow = 0
        sequence = ''
        backSlash = re.compile(r'\\\s*\n')
        whiteSpace = re.compile('\s')
        for line in fileHandle.readlines():
            if line[:8] == 'SEQUENCE':
                line = line[8:]
                startNow = 1
            if startNow and backSlash.search(line):
                sequence = sequence + backSlash.sub('', line)
            elif startNow and not backSlash.search(line):
                sequence = sequence + line
                break
        sequence = whiteSpace.sub('', sequence)
        for eachChar in sequence:
            self.aalist.append(AminoAcid.AminoAcid(eachChar)[1])
        
    
    
###############################################################################
    def ReadSeq(self, fileName):
        """
        -reads in a sequence file which contains:
          -sequence in 3-letter or 1-letter notation
          -spaces or tabs or lineends
          -uses DeleteComments, the following comments are recognized:
            { bla }
            { bla { bla } bla }
            ! bla until end of line
            # bla until end of line
        -if an aminoacid has more than 4 characters (like CYSS or ASP- in
         the DIANA format), only the first three characters will be used
        -numbers and aminoacids in non-standard notation are not used!
        -any non-alphabetic characters are removed automatically
        """
        if os.path.exists(fileName):
            print 'reading sequence', fileName
            self.fileName = fileName
        else:
            print 'the file', fileName
            print 'does not exist. Abort ReadSeq method.'
            return
        #important: delete the aalist aatribute:
        self.aalist = []
        #get the file without the comments:
        bigString = DeleteComments.GetString(fileName)
        #delete all non-alphabetic characters:
        notABC = re.compile('[^a-zA-Z\n\t]')
        bigString = notABC.sub(' ', bigString)
        #split the string in lines and loop over the lines:
        lines = string.split(bigString, '\n')
        for line in lines:
            linelist = string.split(line)
            for element in linelist:
                #get the 3-letter code from the first 3 characters:
                aaoutlist = AminoAcid.AminoAcid(element[0:3])
                #catch all the numbers and non-standard amino acids:
                if string.strip(aaoutlist[2]) == '':
                    print 'could not understand:', element[0:3]
                    print '=> this is not included in the sequence!'
                    continue
                #append in to the list attribute:
                self.aalist.append(aaoutlist[1])
        
    
###############################################################################
    def Stdout(self):
        """
        writes the list of aminoacids to stdout
        """
        print 'the file', self.fileName, 'contains the sequence:'
        for eachaa in self.aalist:
            print eachaa
        
###############################################################################    
    def StdoutOne(self):
        """
        writes the list of aminoacids to stdout (1-letter code,
        all the stuiff in one line)
        """
        print 'the file', self.fileName, 'contains the sequence:'
        outS = ''
        for eachaa in self.aalist:
            outS = outS + AminoAcid.AminoAcid(eachaa)[0]
        print outS
        print 'the sequence includes ' + str(len(self.aalist)) + ' residues.'


###############################################################################
    def WriteFasta(self, outfile):
        """
        writes the list of aminoacids to the specified fileName
        in FASTA format
        """
        try:
            outhandle = TextFile.TextFile(outfile, 'w')
        except IOError:
            print 'could not open the file', outfile
            print 'Abort WriteFasta method.'
            return
        print 'writing to the file:', outfile

        #for the comment use the fileName
        #(if it's empty use 'SEQUENCELIST_OUTPUT'):
        if string.strip(self.fileName):
            outhandle.write('>' + string.strip(self.fileName)+ '\n')
        else:
            outhandle.write('>SEQUENCELIST_OUTPUT\n')

        iii = 0
        for eachaa in self.aalist:
            outhandle.write(AminoAcid.AminoAcid(eachaa)[0])
            iii = iii + 1
            if iii == 60:
                outhandle.write('\n')
                iii=0
        
        outhandle.write('\n')
        outhandle.close()



###############################################################################
    def WriteSeq(self, outfile):
        """
        writes the list of aminoacids to the specified file
        only amino acids in 3-letter notation will appear there
        """
        try:
            outhandle = TextFile.TextFile(outfile, 'w')
        except IOError:
            print 'could not open the file', outfile
            print 'Abort WriteSeq method.'
            return
        print 'writing to the file:', outfile
        for eachaa in self.aalist:
            outhandle.write(eachaa + '\n')
        outhandle.close()
        

###############################################################################
    def WriteSeqSmall(self, outfile):
        """
        writes the list of aminoacids to the specified file
        only amino acids in 3-letter notation will appear there
        """
        try:
            outhandle = TextFile.TextFile(outfile, 'w')
        except IOError:
            print 'could not open the file', outfile
            print 'Abort WriteSeq method.'
            return
        print 'writing to the file:', outfile
        for eachaa in self.aalist:
            outhandle.write(string.lower(eachaa) + '\n')
        outhandle.close()
        
###############################################################################
    def WriteSeq1(self, outfile):
        """
        writes the list of aminoacids to the specified file
        1-letter code
        """
        try:
            outhandle = TextFile.TextFile(outfile, 'w')
        except IOError:
            print 'could not open the file', outfile
            print 'Abort WriteSeq1 method.'
            return
        print 'writing to the file:', outfile
        for eachaa in self.aalist:
            outhandle.write(AminoAcid.AminoAcid(eachaa)[0])
        outhandle.write('\n')
        outhandle.close()


    def WriteXML2String(self, startNumber=1, xmlVersion="1.0", encoding="UTF-8",\
                        dtdFileName="sequence1.1.dtd", dtdVersionName="1.1",\
                        segmentName="", segmentType=""):
        """
        returns a string containing the sequence in XML format

        startNumber      sequence number the whole sequence start with
        xmlVersion       appears in the XML declaration
        encoding         appears in the XML declaration
        dtdFileName      points to the corresponding DTD
        dtdVersionName   version name/number of the DTD
        segmentName      segment name of your protein, RNA, DNA, ...
        """
        #beginning and end:
        beginSection = """<?xml version="%s" encoding="%s"?>
<!DOCTYPE aria_seq:sequence SYSTEM "%s">
<aria_seq:sequence>
  <aria_seq:version>%s</aria_seq:version>
  <aria_seq:segment>""" %\
        (xmlVersion, encoding, dtdFileName, dtdVersionName)

        endSection = """  </aria_seq:segment>
</aria_seq:sequence>
"""
        #loop over all the residues:
        outString = ''
        for eachS in self.aalist:
            eachEntry = """
    <aria_seq:segment_name>%s</aria_seq:segment_name>
    <aria_seq:segment_type>%s</aria_seq:segment_type>
    <aria_seq:residue>
      <aria_seq:residue_number>%s</aria_seq:residue_number>
      <aria_seq:residue_type>%s</aria_seq:residue_type>
    </aria_seq:residue>""" % (segmentName, segmentType, str(startNumber), eachS)
            outString = outString + eachEntry
            startNumber = startNumber + 1
        outString = beginSection + outString + endSection
        return outString
    
    def WriteXML2File(self, fileName):
        """writes the sequence to an XML file"""
        outString = self.WriteXML2String()
        try:
            outhandle = TextFile.TextFile(fileName, 'w')
        except IOError:
            print 'could not open the file', fileName
            print 'Abort WriteXML2File method.'
            return
        print 'writing to the file:', fileName
        outhandle.write(outString)
        outhandle.close()
        
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
    SL=SequenceList()
    SL.aalist = ['ARG', 'GLY', 'HIS', 'PRO', 'GLU', 'TYR', 'ASP', 'GLN']
    print '1. The test sequence is:'
    print SL.aalist
    print '\n2. The Stdout() method prints:'
    SL.Stdout()
    print '\n3. The XML output is:'
    print SL.WriteXML2String()
    print '\ntest done. bye.'
