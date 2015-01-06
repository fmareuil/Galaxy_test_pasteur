#  This file is a part of COV2HTML software
#  COV2HTML <http://sourceforge.net/projects/cov2html/> 
#  A visualization and analysis tool of Bacterial NGS data for biologists.
# 
#  ---------------------------------------------------------------------- 
#  Copyright (C) 2012 Marc Monot
#  
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version. 
#   <http://www.gnu.org/licenses/gpl-3.0.html>
# 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#  ---------------------------------------------------------------------- 
#  
#  @author Marc Monot <marc.monot@pasteur.fr>
#  @author Mickael Orgeur


### PYTHON Module
# Tkinter
from Tkinter import *
from tkFileDialog import *
from tkMessageBox import *

# Project
from functions import * 

# Add from Python Package Index
from ordereddict import *


class LaunchInterface(Tk):
    """Create the Tk interface"""

    def __init__(self, interface):
        Tk.__init__(self, interface)
        self.interface = interface
        self.initialize()

    def initialize(self):
        self.grid()
                
        ### VARIABLES initialization
        # Stockage
        self.gbk = StringVar()
        self.embl = StringVar()
        self.fasta = StringVar()
        self.cds = StringVar()
        self.sam = StringVar()
        self.eland = StringVar()
        self.wig = StringVar()
        self.datatype = StringVar()
        self.secret = StringVar()
        self.length = StringVar()
        self.strand = StringVar()
        self.coverage = StringVar()        
        self.outdir = StringVar()
        self.help = '(i) Browse file(s) associated to the reference;\
\n(ii) Browse your alignment file;\n(iii) Select the corresponding data type;\
\n(iv) Provide a directory for the output file;\n(v) Set eventually optional parameters.'
        
        
        ### INTERFACE
        # Reference information
        self.frame = Frame(self, width=300, height=50, bg='white', bd=5, relief='raised')
        self.frame.grid(column=1, row=0, columnspan=2)
        self.label = Label(self, 
                            anchor='center', bg='white', fg='black', font=('Sans', 13, 'bold'),
                            text='Reference Information')
        self.label.grid(column=1, row=0, columnspan=2)

        # GenBank annotation file
        self.frame = Frame(self, width=300, height=80, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=1, row=1, columnspan=2, rowspan=2, padx=10, pady=5)
        self.label = Label(self, bg='white', fg='blue', text='GenBank file')
        self.label.grid(column=1, row=1, columnspan=2)
        self.button = Button(self, text='Browse', command=self.browseGenbank)
        self.button.grid(column=1, row=2)
        self.entry = Entry(self, textvariable=self.gbk, state='readonly')
        self.entry.grid(column=2, row=2)

        # EMBL annotation file
        self.frame = Frame(self, width=300, height=80, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=1, row=3, columnspan=2, rowspan=2)
        self.label = Label(self, bg='white', fg='blue', text='EMBL file')
        self.label.grid(column=1, row=3, columnspan=2)
        self.button = Button(self, text='Browse', command=self.browseEmbl)
        self.button.grid(column=1, row=4)
        self.entry = Entry(self, textvariable=self.embl, state='readonly')
        self.entry.grid(column=2, row=4)

        # Fasta sequence and CDS region files
        self.frame = Frame(self, width=300, height=180, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=1, row=5, columnspan=2, rowspan=4, padx=10, pady=5)
        # Reference sequence
        self.label = Label(self, bg='white', fg='blue', text='Reference Sequence (Fasta)')
        self.label.grid(column=1, row=5, columnspan=2)
        self.button = Button(self, text='Browse', command=self.browseFasta)
        self.button.grid(column=1, row=6)
        self.entry = Entry(self, textvariable=self.fasta, state='readonly')
        self.entry.grid(column=2, row=6)
        # CDS regions
        self.label = Label(self, bg='white', fg='blue', text='CDS regions (Fasta or GFF3)')
        self.label.grid(column=1, row=7, columnspan=2)
        self.button = Button(self, text='Browse', command=self.browseCds)
        self.button.grid(column=1, row=8)
        self.entry = Entry(self, textvariable=self.cds, state='readonly')
        self.entry.grid(column=2, row=8)

        # Alignment file
        self.frame = Frame(self, width=300, height=50, bg='white', bd=5, relief='raised')
        self.frame.grid(column=3, row=0, columnspan=2)
        self.label = Label(self, 
                            anchor='center', bg='white', fg='black', font=('Sans', 13, 'bold'),
                            text='Alignment File')
        self.label.grid(column=3, row=0, columnspan=2)

        # SAM/BAM file
        self.frame = Frame(self, width=300, height=80, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=3, row=1, columnspan=2, rowspan=2, padx=10, pady=5)
        self.label = Label(self, bg='white', fg='blue', text='SAM/BAM file')
        self.label.grid(column=3, row=1, columnspan=2)
        self.button = Button(self, text='Browse', command=self.browseSam)
        self.button.grid(column=3, row=2)
        self.entry = Entry(self, textvariable=self.sam, state='readonly')
        self.entry.grid(column=4, row=2)

        # ELAND file
        self.frame = Frame(self, width=300, height=80, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=3, row=3, columnspan=2, rowspan=2)
        self.label = Label(self, bg='white', fg='blue', text='ELAND file')
        self.label.grid(column=3, row=3, columnspan=2)
        self.button = Button(self, text='Browse', command=self.browseEland)
        self.button.grid(column=3, row=4)
        self.entry = Entry(self, textvariable=self.eland, state='readonly')
        self.entry.grid(column=4, row=4)

        # WIG file
        self.frame = Frame(self, width=300, height=80, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=3, row=5, columnspan=2, rowspan=2, padx=10, pady=5, sticky='n')
        self.label = Label(self, bg='white', fg='blue', text='WIG file')
        self.label.grid(column=3, row=5, columnspan=2)
        self.button = Button(self, text='Browse', command=self.browseWig)
        self.button.grid(column=3, row=6)
        self.entry = Entry(self, textvariable=self.wig, state='readonly')
        self.entry.grid(column=4, row=6)

        # Data type
        self.frame = Frame(self, width=300, height=80, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=3, row=7, columnspan=2, rowspan=2, padx=10, pady=5, sticky='s')
        self.label = Label(self, bg='white', fg='blue', text='Data type')
        self.label.grid(column=3, row=7, columnspan=2, sticky='s')
        # DNA data
        self.radiobutton = Radiobutton(self, text='DNA', variable=self.datatype, value='dna', 
                                        bg='white', highlightbackground='white')
        self.radiobutton.select()
        self.radiobutton.grid(column=3, row=8, columnspan=2, sticky='w', padx=20)
        # RNA data
        self.radiobutton = Radiobutton(self, text='RNA', variable=self.datatype, value='rna',
                                        bg='white', highlightbackground='white')
        self.radiobutton.grid(column=3, row=8, columnspan=2, sticky='w', padx=95)
        # TSS data
        self.radiobutton = Radiobutton(self, text='TSS', variable=self.datatype, value='tss',
                                        bg='white', highlightbackground='white')
        self.radiobutton.grid(column=3, row=8, columnspan=2, sticky='e', padx=95)
        # CHip data
        self.radiobutton = Radiobutton(self, text='ChIP', variable=self.datatype, value='chip',
                                        bg='white', highlightbackground='white')
        self.radiobutton.grid(column=3, row=8, columnspan=2, sticky='e', padx=20)
      
        # Output directory
        self.frame = Frame(self, width=300, height=80, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=1, row=9, columnspan=4, rowspan=2, padx=10, pady=5)
        self.label = Label(self, bg='white', fg='blue', text='Output directory')
        self.label.grid(column=1, row=9, columnspan=4)
        self.button = Button(self, text='Browse', command=self.browseOutputDir)
        self.button.grid(column=1, row=10, columnspan=3, sticky='w', padx=180)
        self.entry = Entry(self, textvariable=self.outdir, state='readonly')
        self.entry.grid(column=1, row=10, columnspan=3, sticky='e')
        
        # Optional parameters
        self.frame = Frame(self, width=300, height=50, bg='white', bd=5, relief='raised')
        self.frame.grid(column=5, row=0, columnspan=2)
        self.label = Label(self, 
                            anchor='center', bg='white', fg='black', font=('Sans', 13, 'bold'),
                            text='Optional Parameters')
        self.label.grid(column=5, row=0, columnspan=2)

        # Confidential data
        self.frame = Frame(self, width=300, height=80, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=5, row=1, columnspan=2, rowspan=2, padx=10, pady=5)
        self.label = Label(self, bg='white', fg='blue', text='Anonymize data')
        self.label.grid(column=5, row=1, columnspan=2)
        # Do not anonymize data
        self.radiobutton = Radiobutton(self, text='No', variable=self.secret, value='no',
                                        bg='white', highlightbackground='white')
        self.radiobutton.select()
        self.radiobutton.grid(column=5, row=2)
        # Anonymize data
        self.radiobutton = Radiobutton(self, text='Yes', variable=self.secret, value='yes', 
                                        bg='white', highlightbackground='white')
        self.radiobutton.grid(column=6, row=2)
        
        # Read threshold
        self.frame = Frame(self, width=300, height=80, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=5, row=3, columnspan=2, rowspan=2, padx=10, pady=5)
        self.label = Label(self, bg='white', fg='blue', text='Read-length threshold (0 to disable)')
        self.label.grid(column=5, row=3, columnspan=2)
        self.length.set('20')
        self.entry = Entry(self, textvariable=self.length, width=10)
        self.entry.grid(column=5, row=4, columnspan=2)
        
        # Strand-specific data
        self.frame = Frame(self, width=300, height=80, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=5, row=5, columnspan=2, rowspan=2, padx=10, pady=5)
        self.label = Label(self, bg='white', fg='blue', text='Strand-specific data')
        self.label.grid(column=5, row=5, columnspan=2)
        # Stranded data
        self.radiobutton = Radiobutton(self, text='No', variable=self.strand, value='no',
                                        bg='white', highlightbackground='white')
        self.radiobutton.select()
        self.radiobutton.grid(column=5, row=6)
        # Not stranded data
        self.radiobutton = Radiobutton(self, text='Yes', variable=self.strand, value='yes', 
                                        bg='white', highlightbackground='white')
        self.radiobutton.grid(column=6, row=6)

        # Genetic Elements coverage
        self.frame = Frame(self, width=300, height=90, bg='white', bd=1.5, relief='solid')
        self.frame.grid(column=5, row=7, columnspan=2, rowspan=2, padx=10, pady=5)
        self.label = Label(self, bg='white', fg='blue', text='Genetic Elements Coverage Style')
        self.label.grid(column=5, row=7, columnspan=2)
        # Mean Coverages
        self.radiobutton = Radiobutton(self, text='Raw Counts', variable=self.coverage, value='counts',
                                        bg='white', highlightbackground='white')
        self.radiobutton.select()
        self.radiobutton.grid(column=5, row=8)
        # Raw Counts
        self.radiobutton = Radiobutton(self, text='Mean Coverage', variable=self.coverage, value='mean',
                                        bg='white', highlightbackground='white')
        self.radiobutton.grid(column=6, row=8)

        # Analysis
        self.button = Button(self, text='Launch Analysis', command=self.launchAnalysis,
                                bg='white', fg='red', font=('Sans', 13, 'bold'),
                                height=2, width=15, bd=3, relief='raised')
        self.button.grid(column=5, row=9, rowspan=2, columnspan=2)
        
        # Quit
        self.button = Button(self, text='Quit', command=self.destroy,
                                bg='white', fg='black', font=('Sans', 13),
                                height=1, width=5, bd=3, relief='raised')
        self.button.grid(column=4, row=9, rowspan=2, padx=10, sticky='e')
        
        # Fix the windows size
        self.grid_columnconfigure(0, weight=1)
        self.resizable(False, False)
    
    
    ### BROWSE functions
    def browseGenbank(self):
        """Get GenBank input file"""
        self.gbk.set(askopenfilename())
        if self.embl.get():
            self.embl.set('')
        if self.fasta.get():
            self.fasta.set('')
        if self.cds.get():
            self.cds.set('')

    def browseEmbl(self):
        """Get EMBL input file"""
        self.embl.set(askopenfilename())
        if self.gbk.get():
            self.gbk.set('')
        if self.fasta.get():
            self.fasta.set('')
        if self.cds.get():
            self.cds.set('')

    def browseFasta(self):
        """Get Fasta input file"""
        self.fasta.set(askopenfilename())
        if self.gbk.get():
            self.gbk.set('')
        if self.embl.get():
            self.embl.set('')
        
    def browseCds(self):
        """Get CDS-regions input-file"""
        self.cds.set(askopenfilename())
        if self.gbk.get():
            self.gbk.set('')
        if self.embl.get():
            self.embl.set('')
        
    def browseSam(self):
        """Get SAM/BAM input file"""
        self.sam.set(askopenfilename())
        if self.eland.get():
            self.eland.set('')
        if self.wig.get():
            self.wig.set('')
        
    def browseEland(self):
        """Get ELAND input file"""
        self.eland.set(askopenfilename())
        if self.sam.get():
            self.sam.set('')
        if self.wig.get():
            self.wig.set('')

    def browseWig(self):
        """Get WIG input file"""
        self.wig.set(askopenfilename())
        if self.sam.get():
            self.sam.set('')
        if self.eland.get():
            self.eland.set('')

    def browseOutputDir(self):
        """Get output directory"""
        self.outdir.set(askdirectory())
    
    
    # VARIABLES storage
    def storeInputs(self, x, input):
        inputs[x] = input
    
    
    ### FEEDBACK
    def popup(self):
        """Create popup processing window"""
        win = Toplevel()
        win.title('MAP2COV - Analysis Follow-Up')
        win.resizable(False, False)
        return win

    def display(self, win, message, widget):
        """Display processing messages"""
        if widget == 'label':
            label = Label(win, bg='white', fg='black', width=50,
                            text=message, anchor='w').pack()
        elif widget == 'text':
            text = Text(win, bg='white', fg='black', height=1, width=57, wrap='none')
            text.insert(END, message)
            text.pack()
        win.update()
        win.after(100) # Display bugs happen for some fast parsing steps
    
    
    ### CHECKPOINT
    def error(self, checkpoint, message):
        """Stop the mainloop in case of error and re-initialize it"""
        if checkpoint:
            showerror('ERROR', message)
            return 1
    
    def warning(self, checkpoint, message):
        """Information message which does not stop the processing"""
        if checkpoint:
            showwarning('WARNING', message)
            return 1
            
    ### LAUNCH function
    def launchAnalysis(self):
        """Execute input files analysis"""
        
        ### INPUTS
        """[0] GenBank, [1] EMBL, [2] Fasta, [3] CDS, [4] SAM/BAM, [5] ELAND, [6] WIG, \
        [7] datatype, [8] outname, [9] secret, [10] length, [11] strand, [12] coverage"""
        global inputs
        inputs = ['', '', '', '', '', '', '', '', '', '', '', '', '']
        
        ### INPUTS test
        # Reference sequence test
        checkpoint = 0
        if self.gbk.get() != '':
            self.storeInputs(0, self.gbk.get())
            pass
        elif self.embl.get() != '':
            self.storeInputs(1, self.embl.get())
            pass
        elif self.fasta.get() != '':
            self.storeInputs(2, self.fasta.get())
            pass
        else:
            checkpoint = 1
        if self.error(checkpoint, 'ERROR: No reference sequence provided'):
            return
        
        # Fasta reference and CDS regions test
        checkpoint = 1
        if self.fasta.get() != '' and self.cds.get() == '':
            self.error(checkpoint, 'ERROR: No CDS Region file provided with Fasta reference')
            return
        elif self.fasta.get() == '' and self.cds.get() != '':
            self.error(checkpoint, 'ERROR: No Reference Sequence provided with CDS file')
            return
        elif self.fasta.get() != '' and self.cds.get() != '':
            self.storeInputs(3, self.cds.get())
        
        # Alignment file test
        checkpoint = 0
        if self.sam.get() != '':
            self.storeInputs(4, self.sam.get())
            pass
        elif self.eland.get() != '':
            self.storeInputs(5, self.eland.get())
            pass
        elif self.wig.get() != '':
            self.storeInputs(6, self.wig.get())
            pass
        else:
            checkpoint = 1
        if self.error(checkpoint, 'ERROR: No alignment file provided'):
            return
        
        # Datatype
        self.storeInputs(7, self.datatype.get())

        # Confidential data
        self.storeInputs(9, self.secret.get())
        
        # Read-length threshold test
        try:
            length = self.length.get()
            length = int(length)
            if length < 0:
                checkpoint = 1
        except ValueError:
            checkpoint = 1
        if checkpoint == 1:
            self.error(checkpoint, 'ERROR: Read-length threshold must be a positive integer')
            return
        self.storeInputs(10, length)
        
        # Strand-specific data
        self.storeInputs(11, self.strand.get())
        
        # Genetic Elements coverage
        self.storeInputs(12, self.coverage.get())
        
        # Tests
        checkpoint = 1
        if inputs[6] != '' and inputs[12] == 'counts':
            self.error(checkpoint, 'The Coverage style Raw counts cannot be used with the WIG format, no sequence information')
            return
        if inputs[6] != '' and inputs[10] != 0 and inputs[10] != 20:
            self.warning(checkpoint, 'The read-length thresholf is ignored with WIG format, no sequence information')
            self.storeInputs(10, 0)
        if inputs[6] != '' and inputs[11] == 'yes':
            self.warning(checkpoint, 'The strand-specific processing is ignored with WIG format, no sequence information')
            self.storeInputs(11, 'no')
        
        # Output directory test
        if self.outdir.get() != '':
            dir = self.outdir.get() + '/' 
        else:
            dir = './'
        if inputs[11] == 'no':
            outname = dir + 'export_map2cov.txt.bz2'
        else:
            outname = dir + 'export_map2cov_stranded.txt.bz2'
            
        # Output file test
        if os.path.isfile(outname):
            i = 1
            basename = outname[:-8]
            outname = '%s(%i).txt.bz2' % (basename, i)
            while os.path.isfile(outname):
                i += 1
                outname = '%s(%i).txt.bz2' % (basename, i)
        self.storeInputs(8, outname)
        
        ### ANALYSIS
        self.compute(inputs)

    ### COMPUTE
    def compute(self, inputs):
        """Data processing"""
        
        ### INPUTS
        gbk, embl, fasta, cds, sam, eland, wig, datatype, outname, secret, length, strand, coverage = inputs[:]
        genome_forward = []
        genome_reverse = []
        genome_coverage = []
        dict_cds = {}
        checkpoint = 1
        reads = positions = 0
        
        ### GENOMIC FILES: GBK, EMBL or FASTA / GFF3
        # GBK FORMAT
        if gbk != '':
            try:
                filin = open(gbk, 'rU')
            except IOError:
                self.error(checkpoint, 'ERROR: No such file %s' % gbk)
                return
            
            # Check GenBank format
            if checkGenbankFile(filin):
                self.error(checkpoint, 'ERROR: Annotation file is not in GenBank format')
                return

            # Processing popup
            win = self.popup()
            
            # CDS infos: name, CDS product, CDS first, CDS end, CDS way, CDS coverage (0)
            self.display(win, '> START Importing CDS Information', 'label')
            dict_cds, exon = cdsGenbankRecovery(filin, secret)
            self.display(win, '%i CDS regions extracted\n' % len(dict_cds), 'label')
            if exon > 0:
                self.display(win, 'Number of excluded CDS: %i\n' % exon, 'label')
            
            # Genome sequence forward
            self.display(win, '> START Importing Reference Sequence File', 'label')
            genome_forward = genomeGenbankRecovery(filin)
            if genome_forward == []:
                self.display(win, 'ERROR: Processing interrupted\n', 'label')
                self.error(checkpoint, 'ERROR: No genome sequence in GenBank file')
                return
            self.display(win, 'Length of the sequence: %i bp\n' % len(genome_forward), 'label')
            
            filin.close()
        
        # EMBL FORMAT
        elif embl != '':
            try:
                filin = open(embl, 'rU')
            except IOError:
                self.error(checkpoint, 'ERROR: No such file %s' % embl)
                return
            
            # Check EMBL format
            if checkEmblFile(filin):
                self.error(checkpoint, 'ERROR: Annotation file is not in EMBL format')
                return
            
            # Processing popup
            win = self.popup()

            # CDS infos: name, CDS product, CDS first, CDS end, CDS way, CDS coverage (0)
            self.display(win, '> START Importing CDS Information', 'label')
            dict_cds, exon = cdsEmblRecovery(filin, secret)
            self.display(win, '%i CDS regions extracted\n' % len(dict_cds), 'label')
            if exon > 0:
                self.display(win, 'Number of excluded CDS: %i\n' % exon, 'label')
            
            # Genome sequence forward
            self.display(win, '> START Importing Reference Sequence', 'label')
            genome_forward = genomeEmblRecovery(filin)
            if genome_forward == []:
                self.display(win, 'ERROR: Processing interrupted\n', 'label')
                self.error(checkpoint, 'ERROR: No genome sequence in EMBL file')
                return
            self.display(win, 'Length of the sequence: %i bp\n' % len(genome_forward), 'label')
            
            filin.close()
        
        # FASTA / GFF3 FORMAT
        elif fasta != '' and cds != '':
            # Check genome and CDS files format (fasta / GFF3)
            if checkFastaFile(fasta):
                self.error(checkpoint, 'ERROR: Reference file is not in Fasta format')
                return
                               
            filin = open(cds, 'rU')     
            fileformat = checkFastaOrGFF(filin)
            if fileformat == 0:
                self.error(checkpoint, 'ERROR: CDS file is not in Fasta or GFF3 format')
                return
                
            filin.close()
            
            
            # Processing popup
            win = self.popup()

            # Genome sequence forward and reverse
            self.display(win, '> START Importing Genome File', 'label')
            genome_forward, genome_reverse = genomeFastaRecovery(fasta)
            self.display(win, 'Length of the sequence: %i bp\n' % len(genome_forward), 'label')
            
            # CDS infos: name, CDS product, CDS first, CDS end, CDS way, CDS coverage (0)
            self.display(win, '> START Importing CDS File', 'label')
            if fileformat == 'F':
                dict_cds = cdsFastaRecovery(cds, genome_forward, genome_reverse, secret)
            if fileformat == 'G':
                dict_cds = cdsGFFRecovery(cds, genome_forward, genome_reverse, secret)
            self.display(win, '%i CDS regions extracted\n' % len(dict_cds), 'label')
        
        
        ### READING FILES : SAM/BAM, ELAND or WIG
        # SAM/BAM FORMAT
        if sam:
            # Check SAM/BAM alignment format
            format = checkSamFile(sam)
            if format == 1:
                self.display(win, 'ERROR: Processing interrupted\n', 'label')
                self.error(checkpoint, 'ERROR: Alignment file is not in SAM/BAM format')
                return
            
            self.display(win, '> START Reading SAM/BAM File', 'label')
            self.display(win, '* * * Data Processing * * *', 'label')
            genome_coverage, reads, tss_coverage = readSamFile(sam, format, length, genome_forward, datatype)
            if genome_coverage == 2:
                self.display(win, 'ERROR: Processing interrupted\n', 'label')
                self.error(checkpoint, 'ERROR: Length of reference sequence provided does not correspond to the mapping reference')
                return
            if genome_coverage == 3:
                self.display(win, 'ERROR: Processing interrupted\n', 'label')
                self.error(checkpoint, 'ERROR: Mapping position out of the reference sequence')
                return
            self.display(win, 'Number of mapped reads: %i\n' % reads, 'label')
        
        # ELAND FORMAT
        elif eland:
            # Check ELAND alignment format
            format = checkElandFile(eland)
            if format == 1:
                self.display(win, 'ERROR: Processing interrupted\n', 'label')
                self.error(checkpoint, 'ERROR: Alignment file is not in ELAND format')
                return
            
            self.display(win, '> START Reading ELAND File', 'label')
            self.display(win, '* * * Data Processing * * *', 'label')
            genome_coverage, reads, tss_coverage = readElandFile(eland, format, length, genome_forward, datatype)
            if genome_coverage == 1:
                self.display(win, 'ERROR: Processing interrupted\n', 'label')
                self.error(checkpoint, 'ERROR: Mapping position out of the reference sequence')
                return
            self.display(win, 'Number of mapped reads: %i\n' % reads, 'label')
        
        # WIG FORMAT
        elif wig:
            # Check WIG coverage format
            if checkWigFile(wig):
                self.display(win, 'ERROR: Processing interrupted\n', 'label')
                self.error(checkpoint, 'ERROR: Coverage file is not in WIG format')
                return
            
            self.display(win, '> START Reading WIG File', 'label')
            self.display(win, '* * * Data Processing * * *', 'label')
            genome_coverage, positions, tss_coverage = readWigFile(wig, genome_forward)
            if genome_coverage == 1:
                self.display(win, 'ERROR: Processing interrupted\n', 'label')
                self.error(checkpoint, 'ERROR: Coverage position out of the reference sequence')
                return
            self.display(win, 'Number of coverage positions: %i\n' % positions, 'label')
            
            
        ### RESULTS FILE
        # CDS coverage
        self.display(win, '> START Creating CDS Coverage', 'label')
        if coverage == 'mean':
            dict_cds = cdsCoverage(genome_coverage, dict_cds, datatype, coverage)
        if coverage == 'counts':
            dict_cds = cdsCoverage(tss_coverage, dict_cds, datatype, coverage)
        self.display(win, '%i CDS coverages calculated\n' % len(dict_cds), 'label')

        # Intergenic region coverage
        self.display(win, '> START Creating Intergenic-Region Coverage', 'label')
        if coverage == 'mean':
            dict_ig = intergenicCoverage(genome_coverage, dict_cds, coverage)
        if coverage == 'counts':
            dict_ig = intergenicCoverage(tss_coverage, dict_cds, coverage)
        self.display(win, '%i intergenic-region coverages calculated\n' % len(dict_ig), 'label')
        
        
        # Output editing
        self.display(win, '> START Writing Result File', 'label')
        output = writeResults(genome_forward, genome_coverage, tss_coverage, dict_cds, dict_ig, datatype, strand)
        if len(output) == 1:
            compressResults(output[0], outname)
            self.display(win, 'Path to the output file:', 'label')
            self.display(win, outname, 'text')
            self.display(win, '\nData analysis: DONE', 'label')
        else:
            compressResults(output[0] + "\n" + output[1], outname)
            self.display(win, 'Path to the stranded output file: ', 'label')
            self.display(win, outname, 'text')
            self.display(win, '\nData analysis: DONE', 'label')
