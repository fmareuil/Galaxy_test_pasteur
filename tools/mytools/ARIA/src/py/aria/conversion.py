"""

          ARIA -- Ambiguous Restraints for Iterative Assignment

                 A software for automated NOE assignment

                               Version 2.3


Copyright (C) Benjamin Bardiaux, Michael Habeck, Therese Malliavin,
              Wolfgang Rieping, and Michael Nilges

All rights reserved.


NO WARRANTY. This software package is provided 'as is' without warranty of
any kind, expressed or implied, including, but not limited to the implied
warranties of merchantability and fitness for a particular purpose or
a warranty of non-infringement.

Distribution of substantively modified versions of this module is
prohibited without the explicit permission of the copyright holders.

$Author: bardiaux $
$Revision: 1.1.1.1 $
$Date: 2010/03/23 15:27:24 $
"""



## TODO: is the residue number checked in PPM list conversion?

import os, re

from aria.ariabase import AriaBaseClass
from aria.TypeChecking import *
from aria.Settings import Settings
from aria.xmlutils import XMLBasePickler as _XMLBasePickler

from aria.Topology import load_topology, TYPE_AMINO_ACID, TYPE_RNA_BASE, \
     TYPE_DNA_BASE, TYPE_NONBASE
from aria.Chain import TYPE_PROTEIN, TYPE_DNA, TYPE_RNA, TYPE_NONPOLYMER
from aria.ConversionTable import ConversionTable, NAMING_CONVENTIONS, \
     IUPAC_CONVENTION, CNS_CONVENTION, DYANA_CONVENTION
from time import clock

## FOR TIMING ...

## AriaBaseClass.display_debug=1
## AriaBaseClass.display_warnings=0

SEQUENCE_FORMATS = ('seq', 'pdb')
## CHEMICAL_SHIFT_FORMATS = ('xeasy', 'ansig', 'nmrview', 'sparky')
## CROSSPEAK_FORMATS = ('xeasy', 'ansig', 'nmrview', 'sparky')
CHEMICAL_SHIFT_FORMATS = ('xeasy', 'ansig', 'nmrview', 'sparky', 'chem')
CROSSPEAK_FORMATS = ('xeasy', 'ansig', 'nmrview', 'sparky', 'tbl')

CONVENTIONS = {'ansig': CNS_CONVENTION,
               'xeasy': DYANA_CONVENTION,
               'nmrview': CNS_CONVENTION,
               'sparky': CNS_CONVENTION,
               'tbl': CNS_CONVENTION,
               'chem': CNS_CONVENTION}

## BARDIAUX 2.2
SPECTRUM_AMBIGUITY_INTRA = 'intra'
SPECTRUM_AMBIGUITY_INTER = 'inter'
SPECTRUM_AMBIGUITY_ALL = 'all'

SPECTRUM_AMBIGUITIES = [SPECTRUM_AMBIGUITY_INTRA,
                        SPECTRUM_AMBIGUITY_INTER,
                        SPECTRUM_AMBIGUITY_ALL]

CCPN_XML_TEMPLATE = \
'''<!DOCTYPE conversion SYSTEM "conversion1.0.dtd">

<!--##################### Conversion ##################################

  Template xml file for data conversion. Fill out the necessary fields
  between the quotation marks "". Leave optional fields that you do
  not want to set as they are.

-->

<conversion>

<!--#################### Project definition ###########################

  Defines the CCP project name and file.
  
  If you are using the conversion the first time and do not have a
  CCP .xml project file, then one will be created based on the data
  you enter lower down. User interaction is necessary in this case
  to resolve ambiguities in atom names, ...
  
  If the CCP .xml project file already exists you can just give its
  location: molecule data that is already present will not be changed,
  peak lists and chemical shift lists will be updated.

  The namingSystem can be used for setting the naming system used for
  the atom/resonance names for the assignments. Recommended setting is
  "auto" (the best match will be automatically determined this way).

  All fields are mandatory.
  
-->  

  <project
     name = "%(project_name)s"
     filename = "%(project_filename)s"
	 naming_system = "auto">
  </project>       

<!--################### Molecular system definition ####################

  Defines the molecular system (e.g. a monomer, a dimer, a ligand-
  protein complex, ...), and the molecules that compose it.
  
  Add new molecules by copying the <molecule> block within the
  <molecularSystem> block. Have to do this to create a dimer (use
  a different segID in this case).
  
  #### Molecule #####

  Type of the molecule can be either "protein", "DNA", "RNA".
  Currently only polymers are handled.
  
  Change firstResidueNumber if the numbering in the chemical shift file
  is different from the sequence (e.g. first residue in sequence is
  numbered as 25 in chemical shift file).
  
  sequenceFilename gives the location of the input sequence file.
  
  sequenceFormat has to be one of (please copy exact name):
  ["Ansig","NmrStar","NmrView","PDB","Sparky","XEasy"]

  Optional fields:

  - segID (default set to " ")

-->

  <molecular_system
  
    name = "%(mol_system_name)s">
	
    <molecule
	
      name = "%(mol_name)s"
      type = "%(mol_type)s"
      segID = "%(mol_segid)s"
      first_residue_number = "%(mol_first_residue_number)d"

      sequence_filename = "%(mol_input)s"
      sequence_format = "%(mol_format)s"
	 
	 />

  </molecular_system>

%(spectra)s  
  
</conversion>
'''

CCPN_SPECTRUM_BLOCK = \
'''<!--################### Spectrum definition ###########################

  Defines the spectra and the chemical shift list and peak list
  associated with it.
  
  ### spectrum ###
  
  Add/remove spectra by copying/deleting the <spectrum> block.

  type has to be: 2D homonuclear noesy: "noesy.hh"
                  3D noesy (carbon):    "noesy_hsqc_HCH.hhc"
                  3D noesy (nitrogen):  "noesy_hsqc_HNH.hhn"
                  4D noesy (carbon):    "noesy_hsqc_HCCH.hhcc"
                  4D noesy (nitrogen):  "noesy_hsqc_HNNH.hhnn"
                  4D noesy (carb/nitr): "noesy_hsqc_HCNH.hhcn"
                  4D noesy (nitr/carb): "noesy_hsqc_HNCH.hhnc"
  
  dataFilesFormat has to be one of:
  ["Ansig","NmrView","Pronto","Sparky","XEasy"]
  
  ### chemicalShifts ###
  
  Note that if filename = "" the values will be calculated from the
  assignments in the peak list (notably for Ansig format this is
  necessary). This will only work for fully assigned peak lists.
  
  ### peakList ###
  
  The entry for the 1st and the 2nd proton should always be present
  and should be a number (1, 2, 3, or 4) that corresponds to the relevant
  dimension in the original peak list file. For multi-dimensional
  experiments the corresponding hetero dimensions have to be set in the
  same way (note that hetero1 corresponds to the hetero dimension linked
  to proton1).
  
  ### peakListAssignments ###
  
  Can only be used with .assign files from XEasy. Set filename to "" for
  all other cases.

-->

  <spectrum
  
    name = "%(name)s"
    type = "%(type)s"
    data_files_format = "%(format)s">
	
    <chemical_shifts
	  filename = "%(shifts_filename)s"
    />
	
    <cross_peaks
      
      filename = "%(peaks_filename)s"
	  
      proton1 = "%(proton1)s"
      hetero1 = "%(hetero1)s"
      proton2 = "%(proton2)s"
      hetero2 = "%(hetero2)s"
	  
	/>
	
	<cross_peaks_assignments
	  filename = "%(assignments)s"
	/>

  </spectrum>
'''

TEMPLATE_XML = \
'''<!DOCTYPE conversion SYSTEM "conversion1.0.dtd">

<!--===================== Conversion ==================================

  Template XML file for converting NMR data from various formats into
  ARIA XML format. Fill in the fields between the quotation marks "".
  Leave optional fields unchanged.

  ARIAs conversion routines supports the following formats:

  Ansig, NMRView, XEasy

  If data conversion is performed with the CCPNMR software suite (in this case
  use the ARIA command line option "convert_ccpn" instead of "convert" to
  launch the conversion) the following formats are supported:

  Ansig, NmrView, Pronto, Sparky, XEasy

-->

<conversion>

<!--==================== Project definition ===========================

  If you want ARIA to automatically create a project XML file, specify the
  desired filename and the name of your project in the fields "filename" and
  "name", respectively. During the conversion process, your data files
  are converted into ARIA XML format (the original files still exist after
  conversion) and referenced by the project file. Use ARIAs graphical user
  interface (GUI), an XML editor, or a text editor to edit the project file.

  Alternatively, you may also leave the fields empty and create an empty
  project XML file by yourself by invoking ARIA from the command line. 

  Optional fields:

  - filename

-->  

  <project name="my_struct">
     <output
       filename="my_project.xml"/>
  </project>       

<!--=================== Molecule definition ===========================

  Supported input formats are "seq" (a sequence of three letter codes) or
  "pdb". If your sequence is in PDB format, a naming convention for the
  residues and atoms needs to be specified; ARIA supports "iupac", "dyana",
  and "cns" format. Furthermore, the "molecule_segid" entry will be ignored
  and instead read from the PDB file. The type of the molecule can
  be set to "PROTEIN", "DNA" or "RNA". The additional attribute
  "first_residue_number" specifies where the residue numbering starts
  (in case of SEQ format). Fill in the attribute "molecule_name" only if you
  want to use CCPNMR conversion.

  Optional fields:

  - molecule_type        (in case of PDB input files)
  - molecule_name        
  - molecule_segid       (in case of PDB input files)
  - first_residue_number (in case of PDB input files)
  
  Conversion via CCPNMR:

  Supported formats are "Ansig", "Fasta", "NmrStar", "NmrView", 
                        "Pdb", "Sparky", "XEasy"

  When working on symmetric multimeric protein, you can specify the segids of each
  monomer by seprating the segid by an "/". Example molecule_segid="A/B"

-->

  <molecule
  
    molecule_type="PROTEIN"
    molecule_name="my_molecule"
    molecule_segid=""
    first_residue_number="1">
    
    <input
      filename=""
      format="seq"
      naming_convention=""/>
      
    <output
      filename=""/>
  </molecule>

<!--=================== Spectrum definition ===========================

  Each spectrum is defined by a list of chemical shifts,
  and a list of NOESY cross peaks. For converting the cross peaks, you need
  to specify the dimensions that correspond to the resonances. The entry for
  the 1st and the 2nd proton dimension must be a number (1, 2, 3, or 4).
  Specify the dimension of the linked heavy nuclei in "hetero1" and "hetero2",
  respectively. If you want to specify the segment(s) for which you have
  chemical shifts or cross-peaks assignments, use the "segid"-field in the
  header of the spectrum-block. In case of multiple segments the corresponding
  segids are separated by a slash (e.g., segids="A" or segids="A/B").

  Optional fields:

  - spectrum_name
  - segids

  Data conversion via CCPNMR

  If you want to use the CCPNMR software suite and the format converter for
  data conversion, you need to specify the type of your NOESY experiment.
  The following experiments are supported:

  spectrum_type:  2D homonuclear noesy: noesy.hh
                  3D noesy (carbon):    noesy_hsqc_HCH.hhc
                  3D noesy (nitrogen):  noesy_hsqc_HNH.hhn
                  4D noesy (carbon):    noesy_hsqc_HCCH.hhcc
                  4D noesy (nitrogen):  noesy_hsqc_HNNH.hhnn
                  4D noesy (carb/nitr): noesy_hsqc_HCNH.hhcn
                  4D noesy (nitr/carb): noesy_hsqc_HNCH.hhnc
                  
  Otherwise, leave the field empty. Also, the format of chemical shifts and
  cross-peaks must be identical.

  Symmetric dimer:
  
  When working on symmetric dimers, you can specify the ambiguity level of
  your NOESY experiment. For example, if your NOESY spectrum comes from an
  asymmetric labeling experiment, you must specify spectrum_ambiguity as
  "inter". ARIA will then consider all the noes from this spectrum as intre-
  molecular ones.
  This field is not compatible with the conversion via CCPNMR.

  If the spectrum_ambiguity is set ton "inter" or "all" you must sepecify
  the segids of the involved monomers.(e.g., segids="A/B").
  
  If you are working on a monomeric protein, you can just leave this field empty
  (spectrum_ambiguity set to "intra" by default)
  
  spectrum_ambiguity:  intra    (noes involving atoms from one monomer only)
                       inter    (noes involving atoms atoms from different monomers)
                       all      (no known information, all noes are ambigous in terms of monomer)

                       


  
  Supported formats: "Ansig", "NmrView", "Pronto", "Sparky", "XEasy"

  For XEasy you may specify your cross peak assignments in a separate
  ".assign" file. In that case specify the filename (attribute "filename" in
  element "assignments"). Set filename to "" for all other cases.

-->

<!--      If you want to convert more than one spectrum, copy the block
          <spectrum> ... </spectrum>, i.e. for n spectra
          would you need n blocks. -->

  <spectrum
  
    spectrum_name=""
    spectrum_type=""
    spectrum_ambiguity=""
    segids="">
    
    <chemical_shifts>
      <input
        filename=""
        format=""/>
      <output
        filename=""/>
    </chemical_shifts>
    
    <cross_peaks>
      <input
        filename=""
        format=""
        proton1=""
        hetero1=""
        proton2=""
        hetero2=""/>
        
      <output
        filename=""/>
        
      <assignments
        filename=""/>

      <peaks_ambiguity
        filename=""/>        
    </cross_peaks>
  </spectrum>
  
<!-- ============================================================== -->

</conversion>'''

## auxiliary functions

def convert_value(x):
    if type(x) == type(''):
        return eval(x)
    else:
        return x

def is_pseudo(atom_name):

    if re.search('[%#\*Q]', atom_name):
        return 1
    else:
        return 0

def convert_atomname(residue_type, atom_name, naming_convention):

    from legacy.PseudoAtom import Pseudo2IupacTuple

    atom_name = atom_name.upper()

    if is_pseudo(atom_name):
        atom_name = Pseudo2IupacTuple(residue_type, atom_name)
        naming_convention = IUPAC_CONVENTION

    else:
        atom_name = (atom_name,)
        
    ## legacy patch: Pseudo2IupacTuple creates an additional
    ## HZ atom in TYR 
        
    if residue_type == 'TYR':
        atom_name = list(atom_name)
        if 'HZ' in atom_name: atom_name.remove('HZ')
        atom_name = tuple(atom_name)

    return atom_name, naming_convention

def abspath(filename):
    if not len(filename): return filename
    from os.path import expanduser, abspath
    return abspath(expanduser(filename))

## wrapper classes for conversion

class Base(AriaBaseClass):

    topology = None
    table = None
    types = [TYPE_AMINO_ACID, TYPE_RNA_BASE, TYPE_DNA_BASE, TYPE_NONBASE]

    def __init__(self, type, number):

        self.check_type(type)

        AriaBaseClass.__init__(self)

        self.type = type
        self.number = number
        
        self.name = None
        self.atom_names = []
        self.terminus = None
        self.residue = None
        self.segid = None
        
    def check_type(self, type):
        if not type in self.types:
            self.error(TypeError, 'Base type "%s" not supported.' % str(type))

    def set_name(self, base_name, format = IUPAC_CONVENTION):

        if not self.type is TYPE_NONBASE:
            args = (base_name, format, IUPAC_CONVENTION, self.type, 0)
            name = self.table.convert_residue(*args)
            if name is not None: base_name = name

        self.name = base_name

    def get_name(self, format = IUPAC_CONVENTION):

        if not self.type is TYPE_NONBASE:
            args = (self.name, IUPAC_CONVENTION, format, self.type, 0)
            name = self.table.convert_residue(*args)
            if name is not None:
                return name

        return self.name
        
    def get_definition(self):

        if self.type == TYPE_AMINO_ACID:
            find = lambda x, f = self.topology.find: f(amino_acid = x)

        elif self.type == TYPE_RNA_BASE:
            find = lambda x, f = self.topology.find: f(rna_base = x)

        elif self.type == TYPE_DNA_BASE:
            find = lambda x, f = self.topology.find: f(dna_base = x)

        try:
            return find(self.name)
        except:
            return None

    def get_atom_definitions(self):

        from aria.Topology import TERMINUS_N_STANDARD, TERMINUS_C_STANDARD, \
             TERMINUS_N_AMINYL, TERMINUS_C_CARBOXYL, \
             TERMINUS_C5_PRIME_PHOSPHATE, TERMINUS_C5_PRIME_HYDROXYL, \
             TERMINUS_C3_PRIME_HYDROXYL

        residue_def = self.get_definition()

        if residue_def is None: return []

        atom_defs = []
        atom_defs += residue_def.getBackboneAtoms()
        atom_defs += residue_def.getSidechainAtoms()
            
        if self.type == TYPE_AMINO_ACID:

            n_terminus = TERMINUS_N_STANDARD
            c_terminus = TERMINUS_C_STANDARD
            
            if self.terminus == 'N':
                n_terminus = TERMINUS_N_AMINYL

            elif self.terminus == 'C':
                c_terminus = TERMINUS_C_CARBOXYL

        elif self.type in [TYPE_DNA_BASE, TYPE_RNA_BASE]:

            n_terminus = TERMINUS_C5_PRIME_PHOSPHATE
            c_terminus = None
            
            if self.terminus == 'N':
                n_terminus = TERMINUS_C5_PRIME_HYDROXYL

            elif self.terminus == 'C':
                c_terminus = TERMINUS_C3_PRIME_HYDROXYL

        elif self.type == TYPE_NONBASE:
            return []
        
        ## get terminus definitions

        if n_terminus and residue_def.getTermini().has_key(n_terminus):

            terminus = residue_def.getTermini()[n_terminus]
            atom_defs += terminus.getAtoms()

        if c_terminus and residue_def.getTermini().has_key(c_terminus):

            terminus = residue_def.getTermini()[c_terminus]
            atom_defs += terminus.getAtoms()

        return atom_defs

    def __str__(self):

        return '%s(%s, %d, %s)' % (self.__class__.__name__, self.type,
                                   self.number, self.name)
        
    def create_atom(self, atom_name, factory, naming_convention):
        """
        Creates the atom instances.
        """
        names, convention = convert_atomname(self.get_name(CNS_CONVENTION),
                                             atom_name, naming_convention)

        if not self.type is TYPE_NONBASE:
            args = (self.get_name(convention), names, convention,
                    IUPAC_CONVENTION, self.type, 0)
            new_atoms = self.table.convert_atoms(*args)
            if new_atoms is not None: names = new_atoms

        atoms = []

        for name in names:

            try:
                atoms.append(factory.createAtom(self.segid, self.number, name))
            except:
                if self.terminus == 'N' and self.residue is not None:
                    from aria.Topology import EQUIV_NTERMINUS
                    g = [g for g in self.residue.getEquivalentGroups()
                         if g.getType() == EQUIV_NTERMINUS]
                    if len(g):
                        atoms += g[0].getAtoms()

        return atoms

    def create_residue(self, factory):
        """
        Converts the base into a Residue-instance.
        """

        if self.residue is not None: return self.residue
        
        from aria.Residue import Residue
        from aria.Topology import ATOM_TYPE_UNKNOWN
        
        ## get atom definitions from topology

        atom_defs = self.get_atom_definitions()

        ## create residue

        residue = Residue(self.number, self.get_name(IUPAC_CONVENTION))

        ## create atoms

        atom_names = self.atom_names

        ## if no atoms specified in the amino acid use the
        ## definition from the library

        if not atom_names and atom_defs:
            atom_names = [a.getName() for a in atom_defs]

        atoms = []
        
        for atom_name in atom_names:

            atom = self.create_atom(atom_name, factory, IUPAC_CONVENTION)

            if not len(atom):
                m = 'No entry for atom(s) with specifier "%s" could be ' + \
                    'created in residue "%s%d"; entry will be skipped.'
                self.warning(m % (atom_name, self.name, self.number))

            else:
                atoms += atom
                
        [residue.addAtom(a) for a in atoms]

        ## complete atoms if possible
            
        for a in residue.getAtoms():

            d = [b for b in atom_defs if a.getName() == b.getName()]

            if len(d) == 1:
                a.getSettings().update(d[0].getSettings())

            else:
                type_entity = a.getSettings().getEntity('type')
                if type_entity.is_valid(a.getName()[0]):
                    m = 'No information on the atom-type exists for atom ' + \
                        '"%s" in residue "%s" of segment "%s"; assuming ' + \
                        'type "%s".'
                    self.warning(m % (a.getName(), residue.getName(),
                                      str(a.getSegid()), a.getName()[0]))
                    a.setType(a.getName()[0])
                else:
                    a.setType(ATOM_TYPE_UNKNOWN)
                
        ## complete residue if possible

        residue_def = self.get_definition()

        if residue_def is not None:

            for group in residue_def.getEquivalentGroups(): 

                if not 0 in [residue.hasAtom(a) for a in group.getAtomNames()]:

                    residue.addEquivalentGroup(group)

        residue.link()

        self.residue = residue

        return residue

class SequenceList(AriaBaseClass):

    def __init__(self, chain_types = None, offset = 1):

        AriaBaseClass.__init__(self)

        check_type(chain_types, NONE, DICT)
        self.chain_types = chain_types
        self.bases = []
        self.offset = offset

    def parse(self, file, format, naming_convention):

        if format == 'seq':
            self._parse_seq(file)
        elif format == 'pdb':
            self._parse_pdb(file, naming_convention)
            
    def _parse_seq(self, file):

        from legacy.SequenceList import SequenceList as Reader
        from aria.PDBReader import BASE_TYPES

        import sys, StringIO
        
        s = Reader()

        stdout = sys.stdout
        sys.stdout = StringIO.StringIO()
        s.ReadSeq(file)
        sys.stdout.seek(0)
        m = sys.stdout.read()
        sys.stdout = stdout
        self.message(m.strip())

        if self.chain_types is None:
            m = 'The molecule type has not been specified, i.e., either ' + \
                '%s, %s, %s or %s.' % (TYPE_PROTEIN, TYPE_DNA, TYPE_RNA,
                                       TYPE_NONPOLYMER)
            self.error(TypeError, m)

        ## only molecules consisting of one chain can be stored in
        ## seq format

        ## BARDIAUX 2.2 multimers

        big_bases = []
        
        for segid, chain_type in self.chain_types.items():
            
            bases = []  
            
            base_type = BASE_TYPES[chain_type]

            for i in range(len(s.aalist)):

                base = Base(base_type, i + self.offset)
                base.segid = segid
                base.set_name(s.aalist[i])

                ## modif BARDIAUX
                #base.set_struc(s.aalist_struc[i])

                bases.append(base)

            
            bases[0].terminus = 'N'
            bases[-1].terminus = 'C'
            
            big_bases += bases

            
        self.bases = big_bases
        
##         segid, chain_type = self.chain_types.items()[0]
##         base_type = BASE_TYPES[chain_type]
##         bases = []

##         for i in range(len(s.aalist)):

##             base = Base(base_type, i + self.offset)
##             base.segid = segid
##             base.set_name(s.aalist[i])
##             bases.append(base)

##         bases[0].terminus = 'N'
##         bases[-1].terminus = 'C'

##         self.bases = bases
    
    def _parse_pdb(self, file, naming_convention):

        from aria.PDBReader import BASE_TYPES, PDBReader

        p = PDBReader()
        chains = p.read(file, self.chain_types, naming_convention)
        
        bases = []

        chain_types = {}

        for segid in chains.keys():

            if not len(segid) == 4: continue

            chain = chains[segid]
            chain_type = chain['chain_type']
            if self.chain_types is None:
                chain_types[segid] = chain_type
            base_type = BASE_TYPES[chain_type]
            numbers = chain.keys()
            numbers.remove('chain_type')
            numbers.sort()

            for n in numbers:

                residue = chain[n]

                base = Base(base_type, n)
                base.set_name(residue['residue_type'], IUPAC_CONVENTION)
                base.atom_names = [a for a in residue.keys()
                                   if a != 'residue_type']
                base.segid = segid
                bases.append(base)

        if self.chain_types is None: self.chain_types = chain_types

        bases[0].terminus = 'N'
        bases[-1].terminus = 'C'

        self.bases = bases

    def create_chains(self, factory):

        from aria.Chain import Chain, ChainSettings

        chains = {}

        for base in self.bases:

            segid = base.segid

            if not chains.has_key(segid):
                s = ChainSettings()
                s['type'] = self.chain_types[segid]
                chains[segid] = Chain(s, segid)

            chains[segid].addResidue(base.create_residue(factory))

        return chains

class PpmList(AriaBaseClass):

    def __init__(self):

        AriaBaseClass.__init__(self)

        self.assignments = []

    def parse(self, file, format):

        from legacy.PpmList import PpmList as Reader

        import sys, StringIO
        
        r = Reader()

        stdout = sys.stdout
        sys.stdout = StringIO.StringIO()

        if format == 'nmrview':
            r.ReadNmrView(file)
        elif format == 'ansig':
            r.ReadAnsig(file)
        elif format == 'xeasy':
            r.ReadXeasyProt(file)
        elif format == 'sparky':
            r.ReadSparky(file)
        # test
        elif format == 'chem':
            r.ReadChem(file)
            
        sys.stdout.seek(0)
        m = sys.stdout.read()
        sys.stdout = stdout
        self.message(m.strip())
        
        self.assignments = r.atomlist

    def _convert_ppm(self, x):

        if x is not None: return float(x)

    def _sort(self, assignments, bases, factory, naming_convention, segids):
        """
        Extract all chemical shifts and group them according to the
        residue number.
        """

        ## store information in dictionary

        chain = {}
        for base in bases:

            if segids is not None and base.segid not in segids:
                continue
##             if chain.has_key(base.number):
##                 m = 'Chemical shift assignments are not uniquely ' + \
##                     'resolvable: residues with same residue number ' + \
##                     '"%d" occur (at least) in the chains with segids ' + \
##                     '"%s" and "%s".'
##                 self.error(ValueError, m % (base.number,
##                                             chain[base.number].segid,
##                                             base.segid))
##             chain[base.number] = base
            
            ## BARDIAUX 2.2 multimers
            chain.setdefault(base.segid,{})
            chain[base.segid][base.number] = base                

        ## create all residues that occur in the chemical shift list

        residues = {}

        ## BARDIAUX 2.2 multimers
        for ss in segids:
            
            for assignment in assignments:

                n = int(assignment.residuenumber)
                assignment.segid = ss
                seg = assignment.segid 


                if not n in residues:

                    if not chain[ss].has_key(n):
                        m = 'Chemical shift assignment entry for atom "%s" ' + \
                            'of residue "%d" could not be converted; residue ' +\
                            'with this number not found%s.'

                        if segids is not None:
                            s = ' in segments ' + '/'.join(segids)
                        else:
                            s = ' in molecule'
                        self.warning(m % (str(assignment.atomname), n, s))
                        continue

                    base = chain[seg][n]
                    residue = base.residue

                    residues.setdefault(seg, {})
                    residues[seg][n] = {'residue': base, 'atoms': {}, 'groups': {}}

                    ## create fake assignments for all atoms and groups in residue

                    d = residues[seg][n]['atoms']

                    for atom in residue.getAtoms():
                        d[atom.getName()] = (None, None)

                    d = residues[seg][n]['groups']

                    for equivalent_group in residue.getEquivalentGroups():
                        d[equivalent_group] = (None, None)

            ## collect assignments for each residue

            for i in range(len(assignments)):

                assignment = assignments[i]
                seg = assignment.segid
                n = int(assignment.residuenumber)

                if not residues[seg].has_key(n): continue

                base = residues[seg][n]['residue']
                atom_assignments = residues[seg][n]['atoms']
                group_assignments = residues[seg][n]['groups']

                ## get chemical shift

                chemical_shift = (self._convert_ppm(assignment.shift),
                                  self._convert_ppm(assignment.shifterror))

                ## resolve atom names and create atoms

                atoms = []

                for atom_name in assignment.atomname:
                    atom = base.create_atom(atom_name, factory, naming_convention)
                    if not len(atom):
                        m = 'Chemical shift assignment entry for atom with ' + \
                            'specifier "%s" in residue "%s%d" of segment ' + \
                            '"%s" could not be processed: atom specifier ' + \
                            'is not resolvable. Entry will be skipped.'
                        self.warning(m % (atom_name, base.name, base.number,
                                          str(base.segid)))
                    else:
                        atoms += atom

                ## store assignment, check wether chemical shift was assigned
                ## to a group of atoms

                if not len(atoms):

                    ## do nothing, no atoms could be created
                    continue

                elif len(atoms) == 1:

                    ## stereo-specific assignment                        
                    atom_assignments[atoms[0].getName()] = chemical_shift

                elif len(atoms) > 1:

                    ## assignment for a group of atoms (averaging or
                    ## floating chirality assignment)
                    atom_names = [a.getName() for a in atoms]
                    atom_names.sort()

                    ## legacy patch

                    if base.residue.getType() == 'TYR' and \
                           atom_names == ['HD1', 'HD2', 'HE1', 'HE2', 'HZ']:

                        atom_names.remove('HZ')

                    is_equivalent_group = 0

                    for equivalent_group in base.residue.getEquivalentGroups():

                        equivalent_atoms = equivalent_group.getAtomNames()
                        equivalent_atoms = list(equivalent_atoms)
                        equivalent_atoms.sort()

                        if atom_names == equivalent_atoms:                        
                            is_equivalent_group = 1
                            break

                    if not is_equivalent_group:

                        message = 'Inconsistency in conversion script: ' + \
                                  'atom group "%s" (residue %d) in ' + \
                                  'assignment %d is not complete.' 
                        message %= (str(assignment.atomname), n, i)

    ##                    self.error(ValueError, message)
                        self.warning(message)
                    group_assignments[equivalent_group] = chemical_shift

        return residues

    def _methylene_assignment(self, residue, group, group_shift,
                              atom_assignments):

        from aria.ShiftAssignment import ShiftAssignment, SpinSystem, \
             AVERAGING_METHOD_FAST, ASSIGNMENT_METHOD_STEREO_SPECIFIC, \
             ASSIGNMENT_METHOD_EQUIVALENT, ASSIGNMENT_METHOD_FLOATING, \
             AVERAGING_METHOD_NONE
        from aria.Topology import EQUIV_METHYLENE, EQUIV_AROMATIC, \
             EQUIV_METHYL, EQUIV_NTERMINUS, EQUIV_ISOPROPYL
        from aria.Datum import ChemicalShift

        atom_shifts = [atom_assignments[a] for a in group.getAtomNames()] 
        atom_ppms = [s[0] for s in atom_shifts]
        atom_ppms = [s for s in atom_ppms if s is not None]

        group_ppm = group_shift[0]

        assignment_method = None

        if group_ppm is not None:

            if not atom_ppms:
                assignment_method = ASSIGNMENT_METHOD_EQUIVALENT

            elif not group_ppm in atom_ppms:
                message = ('Residue "%s": inconsistent assignment '+\
                           'for group "%s". Group chemical shift: %s, ' +\
                           'atom chemical-shifts: %s') \
                          % (residue.getName(), str(group), str(group_ppm),
                             str([round(p, 3) for p in atom_ppms]))
                self.warning(message)
                return None
            
            else:
                assignment_method = ASSIGNMENT_METHOD_FLOATING

        elif atom_ppms:

            if len(atom_ppms) == 2 and atom_ppms[0] == atom_ppms[1]:
                assignment_method = ASSIGNMENT_METHOD_EQUIVALENT
                ppm_value = atom_ppms[0]
                ppm_error = max([atom_shifts[0][1], atom_shifts[1][1]])
                group_shift = (ppm_value, ppm_error)

            else:
                assignment_method = ASSIGNMENT_METHOD_FLOATING

        else: return None

        ## create shift assignment

        if assignment_method == ASSIGNMENT_METHOD_EQUIVALENT:

            atoms = [residue.getAtom(n) for n in group.getAtomNames()]

            shift = ChemicalShift(*group_shift)

            spin_system = SpinSystem(AVERAGING_METHOD_FAST)
            spin_system.setAtoms(tuple(atoms))
            spin_system.setChemicalShifts((shift,))

            shift_assignment = ShiftAssignment(assignment_method)
            shift_assignment.setSpinSystems((spin_system,))

        elif assignment_method == ASSIGNMENT_METHOD_FLOATING:

            shifts = [ChemicalShift(*atom_assignments[a]) for
                      a in group.getAtomNames()]
            shifts = tuple(shifts)

            atoms = [residue.getAtom(n) for n in group.getAtomNames()]

            spin_systems = []

            for atom in atoms:

                spin_system = SpinSystem(AVERAGING_METHOD_NONE)
                spin_system.setAtoms((atom,))
                spin_system.setChemicalShifts(shifts)

                spin_systems.append(spin_system)

            shift_assignment = ShiftAssignment(ASSIGNMENT_METHOD_FLOATING)
            shift_assignment.setSpinSystems(tuple(spin_systems))

        ## delete atoms

        if assignment_method is not None:

            for atom_name in group.getAtomNames():

                del atom_assignments[atom_name]

        else:
            shift_assignment = None

        return shift_assignment

    def _methyl_assignment(self, residue, group, group_shift,
                           atom_assignments):

        from aria.ShiftAssignment import ShiftAssignment, SpinSystem, \
             AVERAGING_METHOD_FAST, ASSIGNMENT_METHOD_STEREO_SPECIFIC, \
             ASSIGNMENT_METHOD_EQUIVALENT, ASSIGNMENT_METHOD_FLOATING, \
             AVERAGING_METHOD_NONE
        from aria.Topology import EQUIV_METHYLENE, EQUIV_AROMATIC, \
             EQUIV_METHYL, EQUIV_NTERMINUS, EQUIV_ISOPROPYL
        from aria.Datum import ChemicalShift

        group_ppm = group_shift[0]

        atom_shifts = [atom_assignments[a] for a in group.getAtomNames()]
        atom_shifts = [s for s in atom_shifts if s[0] is not None]
        atom_ppms = [s[0] for s in atom_shifts]
        
        if group_ppm is not None:
            atom_shifts.append(group_shift)

        ## inconsistent shift assignment

        if [s for s in atom_shifts if s[0] <> atom_shifts[0][0]]:

            message = ('Residue "%s": inconsistent assignment '+\
                       'for group "%s". Group chemical shift: %s, ' +\
                       'atom chemical-shifts: %s') \
                      % (residue.getName(), str(group), str(group_ppm),
                         str([round(p, 3) for p in atom_ppms]))
            self.warning(message)
            return None

        ## assign shift to the whole methyl group

        if not atom_shifts:
            ppm_value = None
        else:
            ppm_value = atom_shifts[0][0]

        ppm_errors = [s[1] for s in atom_shifts if s[0] is not None]

        if ppm_errors:
            ppm_error = max(ppm_errors)
        else:
            ppm_error = None

        chemical_shift = ChemicalShift(ppm_value, ppm_error)

        atoms = [residue.getAtom(n) for n in group.getAtomNames()]

        spin_system = SpinSystem(AVERAGING_METHOD_FAST)
        spin_system.setAtoms(tuple(atoms))
        spin_system.setChemicalShifts((chemical_shift, ))

        shift_assignment = ShiftAssignment(ASSIGNMENT_METHOD_EQUIVALENT)
        shift_assignment.setSpinSystems((spin_system,))

        ## delete atoms

        for atom_name in group.getAtomNames():
            del atom_assignments[atom_name]

        return shift_assignment

    def _isopropyl_assignment(self, residue, group_assignments,
                              atom_assignments):

        from aria.ShiftAssignment import ShiftAssignment, SpinSystem, \
             AVERAGING_METHOD_FAST, ASSIGNMENT_METHOD_STEREO_SPECIFIC, \
             ASSIGNMENT_METHOD_EQUIVALENT, ASSIGNMENT_METHOD_FLOATING, \
             AVERAGING_METHOD_NONE
        from aria.Topology import EQUIV_METHYLENE, EQUIV_AROMATIC, \
             EQUIV_METHYL, EQUIV_NTERMINUS, EQUIV_ISOPROPYL
        from aria.Datum import ChemicalShift

        methyl_groups = [g for g in residue.getEquivalentGroups()
                         if g.getType() == EQUIV_METHYL]

        if not len(methyl_groups) == 2: return ()

        ## assign first methyl group

        group = methyl_groups[0]
        group_shift = group_assignments[group]
        del group_assignments[group]

        args = (residue, group, group_shift, atom_assignments)

        shift_assignment1 = self._methyl_assignment(*args)

        ## assign second methyl group

        group = methyl_groups[1]
        group_shift = group_assignments[group]
        del group_assignments[group]

        args = (residue, group, group_shift, atom_assignments)

        shift_assignment2 = self._methyl_assignment(*args)

        ## check wether methyls could be assigned

        if shift_assignment1 is None and shift_assignment2 is not None:
            return (shift_assignment2,)

        if shift_assignment2 is None and shift_assignment1 is not None:
            return (shift_assignment1,)

        if shift_assignment1 is None and shift_assignment2 is None:

            if residue.getType() == 'ILE': return ()

            isopropyl_group = [g for g in residue.getEquivalentGroups()
                               if g.getType() == EQUIV_ISOPROPYL]

            if not len(isopropyl_group) == 1: return ()

            isopropyl_group = isopropyl_group[0]

            args = (residue, isopropyl_group,
                    group_assignments[isopropyl_group],
                    atom_assignments)

            del group_assignments[isopropyl_group]

            shift_assignment = self._methyl_assignment(*args)

            if shift_assignment is not None: return (shift_assignment,)

            else: return ()

        ## get isopropyl group if existing

        isopropyl_group = None

        ## TODO: VAL, LEU hard-coded

        if residue.getType() in ['VAL', 'LEU']:

            iso_groups = [g for g in residue.getEquivalentGroups()
                          if g.getType() == EQUIV_ISOPROPYL]

            if len(iso_groups) == 1: isopropyl_group = iso_groups[0]

        ## spin systems and shifts we have to deal with

        spin_system1 = shift_assignment1.getSpinSystems()[0]
        spin_system2 = shift_assignment2.getSpinSystems()[0]

        shift1 = spin_system1.getChemicalShifts()[0]
        shift2 = spin_system2.getChemicalShifts()[0]

        ## floating assignment

        if shift1[0] <> shift2[0]:

            ## check consistency with isopropyl group

            if isopropyl_group:

                isopropyl_shift = group_assignments[isopropyl_group]
                if isopropyl_shift[0] is not None and not \
                   isopropyl_shift[0] in [shift1[0], shift2[0]]:
                    message = ('Residue "%s": inconsistent assignment '+\
                               'for group "%s". Group chemical shift: %s, ' +\
                               'atom chemical-shifts: %s') \
                              % (residue.getName(), str(isopropyl_group),
                                 str(isopropyl_shift[0]),
                                 str([round(shift1[0], 3),
                                      round(shift2[0], 3)]))
                    self.warning(message)
                    return ()

            ## create floating assignment

            shift_assignment = ShiftAssignment(ASSIGNMENT_METHOD_FLOATING)

            spin_system1.setChemicalShifts((shift1, shift2))
            spin_system2.setChemicalShifts((shift1, shift2))

            shift_assignment.setSpinSystems((spin_system1, spin_system2))

            ## delete isopropyl assignment

            if isopropyl_group: del group_assignments[isopropyl_group]

            return (shift_assignment,)

        elif residue.getType() == 'ILE':

            return (shift_assignment1, shift_assignment2)

        else:

            ## check consistency with isopropyl group

            if isopropyl_group:

                isopropyl_shift = group_assignments[isopropyl_group]
                if isopropyl_shift[0] is not None:
                    if shift1[0] is not None and not \
                       isopropyl_shift[0] == shift1[0]:
                        message = ('Residue "%s": inconsistent assignment '+\
                                   'for group "%s". Group chemical shift: '+\
                                   '%s, atom chemical-shifts: %s') \
                                   % (residue.getName(), str(isopropyl_group),
                                      str(isopropyl_shift[0]),
                                      str([round(shift1[0], 3),
                                           round(shift2[0], 3)]))
                        self.warning(message)
                        return ()

            else: isopropyl_shift = [None, None]

            atoms = list(spin_system1.getAtoms()) + \
                    list(spin_system2.getAtoms())

            ppm_values = [isopropyl_shift[0], shift1[0], shift2[0]]
            ppm_values = [p for p in ppm_values if p is not None]

            if ppm_values:                
                ppm_value = ppm_values[0]

                ## TODO: check is probably superfluous here
                if ppm_values.count(ppm_values) == len(ppm_values):
                    message = ('Residue "%s": inconsistent assignment '+\
                               'for group "%s". Group chemical shift: '+\
                               '%s, atom chemical-shifts: %s') \
                               % (residue.getName(), str(isopropyl_group),
                                  str(isopropyl_shift[0]),
                                  str([round(shift1[0], 3),
                                       round(shift2[0], 3)]))
                    self.warning(message)
                    return ()
            else:
                ppm_value = None

            ppm_error = max((isopropyl_shift[1], shift1[1], shift2[1]))

            shift = ChemicalShift(ppm_value, ppm_error)

            spin_system = SpinSystem(AVERAGING_METHOD_FAST)
            spin_system.setAtoms(tuple(atoms))
            spin_system.setChemicalShifts((shift,))

            shift_assignment = ShiftAssignment(ASSIGNMENT_METHOD_EQUIVALENT)
            shift_assignment.setSpinSystems((spin_system,))

            return (shift_assignment,)

    def _aromatic_ring_assignment(self, residue, group_assignments,
                                  atom_assignments):

        from aria.ShiftAssignment import ShiftAssignment, SpinSystem, \
             AVERAGING_METHOD_FAST, ASSIGNMENT_METHOD_STEREO_SPECIFIC, \
             ASSIGNMENT_METHOD_EQUIVALENT, ASSIGNMENT_METHOD_FLOATING, \
             AVERAGING_METHOD_NONE
        from aria.Topology import EQUIV_METHYLENE, EQUIV_AROMATIC, \
             EQUIV_METHYL, EQUIV_NTERMINUS, EQUIV_ISOPROPYL
        from aria.Datum import ChemicalShift

        aromatic_groups = [g for g in residue.getEquivalentGroups()
                           if g.getType() == EQUIV_AROMATIC]

        small_groups = [g for g in aromatic_groups \
                        if len(g.getAtomNames()) == 2]

        if not len(small_groups) == 2: return ()

        ## create shift assignments for methylene-type groups

        ## first pair of aromatic protons

        group = small_groups[0]
        group_shift = group_assignments[group]
        args = (residue, group, group_shift, atom_assignments)
        del group_assignments[group]

        shift_assignment1 = self._methylene_assignment(*args)

        ## second pair of aromatic protons

        group = small_groups[1]
        group_shift = group_assignments[group]
        args = (residue, group, group_shift, atom_assignments)
        del group_assignments[group]

        shift_assignment2 = self._methylene_assignment(*args)

        ## aromatic ring

        ring = [g for g in aromatic_groups if len(g.getAtomNames()) >= 4]

        if not len(ring) == 1:
            return ()

        ring = ring[0]
        ring_shift = group_assignments[ring]
        del group_assignments[ring]

        ## check wether both assignments are valid

        if shift_assignment2 is None and shift_assignment1 is not None: 

            return (shift_assignment1,)

        elif shift_assignment1 is None and shift_assignment2 is not None:

            return (shift_assignment2, )

        elif shift_assignment1 is None and shift_assignment2 is None:

            args = (residue, ring, ring_shift, atom_assignments)
            
            shift_assignment = self._methyl_assignment(*args)

            if shift_assignment: return (shift_assignment, )

            else: return ()

        else:

            spin_systems1 = shift_assignment1.getSpinSystems()
            spin_systems2 = shift_assignment2.getSpinSystems()

            shifts1 = spin_systems1[0].getChemicalShifts()
            shifts2 = spin_systems2[0].getChemicalShifts()

            shifts1 = [s for s in shifts1 if s[0] is not None]
            shifts2 = [s for s in shifts2 if s[0] is not None]            

            if len(shifts1) == 1 and len(shifts2) == 1 and \
               shifts1[0][0] == shifts2[0][0]:

                ppm_value = shifts1[0][0]
                ppm_error = max([shifts1[0][1], shifts2[0][1]])

                if ring_shift[0] is not None and ring_shift[0] <> ppm_value:

                    message = ('Residue "%s": inconsistent assignment '+\
                               'for group "%s". Group chemical shift: %s, ' +\
                               'atom chemical-shift: %s') \
                              % (residue.getName(), str(ring),
                                 str(ring_shift[0]), str(ppm_value))
                    self.warning(message)

                    return ()

                ppm_error = max([ppm_error, ring_shift[1]])

                chemical_shift = ChemicalShift((ppm_value, ppm_error))

                atoms = [residue.getAtom(a) for a in ring.getAtomNames()]

                spin_system = SpinSystem(AVERAGING_METHOD_FAST)
                spin_system.setAtoms(tuple(atoms))
                spin_system.setChemicalShifts((chemical_shift, ))

                shift_assignment = \
                                 ShiftAssignment(ASSIGNMENT_METHOD_EQUIVALENT)
                shift_assignment.setSpinSystems((spin_system,))

                for atom in atoms:

                    if atom_assignments.has_key(atom.getName()):

                        del atom_assignments[atom.getName()]

                return (shift_assignment,)

            else:

                if ring_shift[0] is not None and not \
                   ring_shift[0] in [s[0] for s in shifts1] + \
                   [s[0] for s in shifts2]:

                    message = ('Residue "%s": inconsistent assignment '+\
                               'for group "%s". Group chemical shift: %s' +\
                               'atom chemical-shift: %s') \
                              % (residue.getName(), str(ring),
                                 str(ring_shift[0]),
                                 str([round(p, 3) for p in shifts1 + shifts2]))
                    self.warning(message)
                    return ()

                return (shift_assignment1, shift_assignment2)

    def create_shiftlist(self, chain, factory, naming_convention, segids):

        ## TODO: take SEGID from PpmFile?
        
        ## Information on SEGID's that might occur in the chemical shift
        ## file (although that's not very likely) is not retained. 
        
        from aria.ChemicalShiftList import ChemicalShiftList
        from aria.ShiftAssignment import ShiftAssignment, SpinSystem, \
             AVERAGING_METHOD_FAST, ASSIGNMENT_METHOD_STEREO_SPECIFIC, \
             ASSIGNMENT_METHOD_EQUIVALENT, ASSIGNMENT_METHOD_FLOATING, \
             AVERAGING_METHOD_NONE
        from aria.Topology import EQUIV_METHYLENE, EQUIV_AROMATIC, \
             EQUIV_METHYL, EQUIV_NTERMINUS, EQUIV_ISOPROPYL
        from aria.Datum import ChemicalShift

        ## get a dictionary of shift assignments for each residue

        residues = self._sort(self.assignments, chain, factory,
                              naming_convention, segids)

        ## create shift assignment instances for each residue
        ## by extracting groups of atoms that have been assigned
        ## to the same chemical shift value and could represent
        ## either a group of atoms in fast exchange or with missing
        ## stereo-specific assignment

        shift_list = ChemicalShiftList()

        ## BARDIAUX 2.2 multimers
        for s in segids:
            
            for n in residues[s].keys():

                residue = residues[s][n]['residue'].residue
                atom_assignments = residues[s][n]['atoms']
                group_assignments = residues[s][n]['groups']

                ## treat special cases first

                if residue.getType() in ['VAL', 'ILE', 'LEU']:
                    args = (residue, group_assignments, atom_assignments)
                    shift_assignments = self._isopropyl_assignment(*args)
                    [shift_list.addShiftAssignment(a) for a in shift_assignments]

                elif residue.getType() in ['PHE', 'TYR']:
                    args = (residue, group_assignments, atom_assignments)
                    shift_assignments = self._aromatic_ring_assignment(*args)
                    [shift_list.addShiftAssignment(a) for a in shift_assignments]

                ## handle groups of potentially equivalent or floating
                ## atoms now

                for group, group_shift in group_assignments.items():

                    ## handle only complete groups

                    if 0 in [atom_assignments.has_key(a) for a in
                             group.getAtomNames()]:
                        continue

                    ## methylene group

                    if group.getType() == EQUIV_METHYLENE:
                        args = (residue, group, group_shift, atom_assignments)
                        shift_assignment = self._methylene_assignment(*args)

                        if shift_assignment is not None:
                            shift_list.addShiftAssignment(shift_assignment)

                    ## methyl group

                    elif group.getType() in (EQUIV_METHYL, EQUIV_NTERMINUS):
                        args = (residue, group, group_shift, atom_assignments)
                        shift_assignment = self._methyl_assignment(*args)

                        if shift_assignment is not None:
                            shift_list.addShiftAssignment(shift_assignment)

                    elif len(group.getAtomNames()) > 3 :
                        message = ('Residue "%s": assignment for group "%s" ' + \
                                  'should have already been handled.') \
                                  % (residue.getName(), str(group))
                        self.warning(message)

                ## finally add the stereospecific assignments of the remaining
                ## atoms

                for atom_name, atom_shift in atom_assignments.items():

                    spin_system = SpinSystem(AVERAGING_METHOD_NONE)
                    spin_system.setAtoms((residue.getAtom(atom_name), ))
                    spin_system.setChemicalShifts((ChemicalShift(*atom_shift), ))

                    shift_assignment = ShiftAssignment(\
                        ASSIGNMENT_METHOD_STEREO_SPECIFIC)

                    shift_assignment.setSpinSystems((spin_system, ))

                    shift_list.addShiftAssignment(shift_assignment)

        return shift_list

class NoeList(AriaBaseClass):

    def __init__(self):

        AriaBaseClass.__init__(self)

        self.peaks = []
        self.peaks_ambiguity = {}

    def _convert_dimension(self, dim):

        if dim is not None:
            return str(dim)
        else:
            return 'N'

    def parse_peaks_ambiguity(self, file):

        d = {}

        f = open(file)
        for l in f:
            s = l.split()
            d[int(s[0])] = str(s[1])
            
        f.close()

        self.peaks_ambiguity = d
            
        
    def parse(self, file, format, proton1, hetero1, proton2, hetero2,
              ppm_file = None):

        from legacy.NoeList import NoeList as Reader

        import sys, StringIO
        
        r = Reader()
        a = [file]

        if ppm_file is not None and format == 'xeasy':
            a.append(ppm_file)
            
        kw = {}

        kw['pro1'] = self._convert_dimension(proton1)
        kw['het1'] = self._convert_dimension(hetero1)        
        kw['pro2'] = self._convert_dimension(proton2)
        kw['het2'] = self._convert_dimension(hetero2)        

        stdout = sys.stdout
        sys.stdout = StringIO.StringIO()

        if format == 'nmrview':
            r.ReadNMRViewPeaks(*a, **kw)
        elif format == 'ansig':
            r.ReadAnsigCpe(*a, **kw)
        elif format == 'tbl':
            r.ReadTbl(*a)#, **kw)
        elif format == 'xeasy':
            r.ReadPeaksProt(*a, **kw)
        elif format == 'sparky':
            r.ReadSparky(*a, **kw)

        sys.stdout.seek(0)
        m = sys.stdout.read()
        sys.stdout = stdout
        self.message(m.strip())
        
        self.peaks = r.peakslist

    def create_spectrum(self, bases, factory, naming_convention,
                        segids = None, ambiguity = SPECTRUM_AMBIGUITY_INTRA):

        from aria.NOESYSpectrum import NOESYSpectrum
        from aria.CrossPeak import CrossPeak
        from aria.Datum import ChemicalShift, Datum
        import aria.Assignment as Assi

        spectra = {}

        ## store information in dictionary

        chain = {}
        for base in bases:

            if segids is not None and base.segid not in segids:
                continue
##             if chain.has_key(base.number):
##                 m = 'Cross peak assignments won\'t be uniquely ' + \
##                     'resolvable: residues with same residue number ' + \
##                     '"%d" occur (at least) in the chains with segids ' + \
##                     '"%s" and "%s".'
##                 self.error(ValueError, m % (base.number,
##                                             chain[base.number].segid,
##                                             base.segid))
                
##             chain[base.number] = base

            ## BARDIAUX 2.2 multimers
            chain.setdefault(base.segid, {})
            chain[base.segid][base.number] = base        

        shift_attrs = {'Proton1': 'p1', 'Proton2': 'p2',
                       'Hetero1': 'h1', 'Hetero2': 'h2'}

        for peak in self.peaks:

            shift_attrs = {'Proton1': 'p1', 'Proton2': 'p2',
                           'Hetero1': 'h1', 'Hetero2': 'h2'}

            contributions = peak.contributions
            peak = contributions[0]

            number = int(peak.peakNumber)
            peak_ambiguity = self.peaks_ambiguity.get(number)

            intensity = Datum(convert_value(peak.intensity),
                              convert_value(peak.intensityError))
            volume = Datum(convert_value(peak.volume),
                           convert_value(peak.volumeError))

            
            cross_peak = CrossPeak(number, volume, intensity)
            cross_peak.setAmbiguity(peak_ambiguity)
            
            for method, attr in shift_attrs.items():

                value = convert_value(getattr(peak, '%sppm' % attr))
                error = convert_value(getattr(peak, 'd%sppm' % attr))
                shift = ChemicalShift(value, error)
                getattr(cross_peak, 'set%sChemicalShift' % method)(shift)
                
            for c in contributions:

                ## BARDIAUX 2.2
                ## if spectrum_ambiguity intra = OK
                ## if spectrum_ambiguity inter => p1 = segA, p2 = segA+segB
                ## if spectrum_ambiguity unknown => p1 = segA, p2 = segA+segB
                for s in segids:

                    
                    base1 = None
                    base2 = None

                    if c.residue1 is not None: base1 = chain[s].get(int(c.residue1))
                    if c.residue2 is not None: base2 = chain[s].get(int(c.residue2))

                    if base1 is None and base2 is None: continue
                    elif base1 is None and base2 is not None: base1 = base2
                    elif base2 is None and base1 is not None: base2 = base1
                    
                    base = {'p1': base1, 'h1': base1, 'p2': base2, 'h2': base2}

                    ## BARDIAUX 2.2
                    if ambiguity == SPECTRUM_AMBIGUITY_INTER or peak_ambiguity == SPECTRUM_AMBIGUITY_INTER:
                        if cross_peak.getProton1Assignments():
                            shift_attrs = {'Proton2': 'p2', 'Hetero2': 'h2'}
                        else:
                            shift_attrs = {'Proton1': 'p1', 'Hetero1': 'h1'}
                        
                    elif ambiguity == SPECTRUM_AMBIGUITY_ALL:
                        if cross_peak.getProton1Assignments():
                            if peak_ambiguity == SPECTRUM_AMBIGUITY_INTRA:
                                shift_attrs = {}
                            else:
                                shift_attrs = {'Proton2': 'p2', 'Hetero2': 'h2'}

                    elif ambiguity == SPECTRUM_AMBIGUITY_INTRA or peak_ambiguity == SPECTRUM_AMBIGUITY_INTRA:
                        shift_attrs = {'Proton1': 'p1', 'Proton2': 'p2',
                                       'Hetero1': 'h1', 'Hetero2': 'h2'}


                        
                    for method, attr in shift_attrs.items():

                        name = getattr(c, 'atomname%s' % attr)[0]
                        if not name: continue
                        atoms = base[attr].create_atom(name, factory,
                                                       naming_convention)


                        if len(atoms) == 2:

                            names = [a.getName() for a in atoms]

                            names.sort()

                            print 'Methylene group found. Using %s instead of %s' % (names[0], str(names))

                            atoms = [a for a in atoms if a.getName() == names[0]]

                        if len(atoms):
                            a_type = Assi.ASSIGNMENT_TYPE_MANUAL
                            a = Assi.Assignment(atoms, a_type)
                            if method == 'Proton1':
                                p1assi = [xx.getAtoms()[0] for xx in cross_peak.getProton1Assignments()]
                                if a.getAtoms()[0] not in  p1assi:
                                    cross_peak.addProton1Assignment(a)
                            elif method == 'Proton2':
                                p1assi = [xx.getAtoms()[0] for xx in cross_peak.getProton2Assignments()]
                                if a.getAtoms()[0] not in  p1assi:
                                    cross_peak.addProton2Assignment(a)
                            else:      
                                getattr(cross_peak, 'add%sAssignment' % method)(a)
                        else:
                            m = 'Cross-peak assignment for atom with ' + \
                                'specifier "%s" in residue "%s%d" of ' + \
                                'segment "%s" could not be processed ' + \
                                'and will be skipped.'
                            self.warning(m % (name, base[attr].name,
                                              base[attr].number,
                                              str(base[attr].segid)))

            name = peak.spectrumName
            if not spectra.has_key(name): spectra[name] = []
            spectra[name].append(cross_peak)

        for name in spectra.keys():
            spectra[name] = NOESYSpectrum(name, spectra[name])

        if len(spectra) > 1:
            self.warning('Multiple spectra found.')

        return spectra.values()[0]

## main class for data conversion

class MoleculeSettings(Settings):

    def create(self):

        from aria.Settings import String, FourLetterString, ChoiceEntity, Path, \
             Integer
        from aria.Chain import TYPE_PROTEIN, TYPE_DNA, TYPE_RNA

        c = (TYPE_PROTEIN, TYPE_DNA, TYPE_RNA)
        d = 'Polymer type: ' + ' \ '.join(c)
        type_choice = ChoiceEntity(c, d)

        c = tuple(SEQUENCE_FORMATS)
        d = 'Sequence file has to be in either of the following formats: ' + \
            ' / '.join(c)
        format_choice = ChoiceEntity(c, d)

        d = 'Atom naming convention in case that the molecule is ' + \
            'defined by a PDB file; supported naming conventions ' + \
            'are: ' + ' \ '.join(NAMING_CONVENTIONS)            
        convention = ChoiceEntity(list(NAMING_CONVENTIONS)+[''], d)

        kw = {'name': String(),
              #'segid': FourLetterString(),
              'segid' : String(),
              'type': type_choice,
              'input': Path(exists = 1),
              'format': format_choice,
              'output': Path(exists = 0),
              'naming_convention': convention,
              'first_residue_number': Integer()}

        return kw

class MoleculeSettingsXMLPickler(_XMLBasePickler):

    def load_from_element(self, node):

        from aria.tools import string_to_segid

        d = MoleculeSettings()

        d['name'] = str(node.molecule_name)
        #d['segid'] = string_to_segid(str(node.molecule_segid))
        d['segid'] = str(node.molecule_segid)
        d['type'] = str(node.molecule_type).upper()
        d['input'] = abspath(str(node.input.filename))
        d['format'] = str(node.input.format).lower()
        d['naming_convention'] = str(node.input.naming_convention)
        d['output'] = abspath(str(node.output.filename))
        d['first_residue_number'] = int(node.first_residue_number)

        if d['format'] == 'pdb' and d['naming_convention'] == '':

            m = 'If sequence file is in PDB format the atom naming ' + \
                'convention has to be specified.'
            self.error(ValueError, m)
        
        return d

class ChemicalShiftsSettings(Settings):

    def create(self):

        from aria.Settings import Path, ChoiceEntity

        c = tuple(CHEMICAL_SHIFT_FORMATS)
        m = 'Supported formats for chemical shift file: ' + ' / '.join(c)
        format_choice = ChoiceEntity(c, m)

        kw = {'input': Path(exists = 1),
              'format': format_choice,
              'output': Path(exists = 0)}

        return kw

class ChemicalShiftsSettingsXMLPickler(_XMLBasePickler):

    def load_from_element(self, e):

        if str(e.input.filename).strip() == '':
            return None

        d = ChemicalShiftsSettings()

        d['input'] = abspath(str(e.input.filename))
        d['format'] = str(e.input.format).lower()
        d['output'] = abspath(str(e.output.filename))

        return d

class CrossPeaksSettings(Settings):

    def create(self):

        from aria.Settings import Path, ChoiceEntity

        c = tuple(CROSSPEAK_FORMATS)
        m = 'Supported formats for the cross peak file: ' + ' / '.join(c)
        format_choice = ChoiceEntity(c, m)

        c = [1, 2, 3, 4, None]
        m = 'If the %s dimension has been measured it has to be set to ' + \
            'either ' + ' / '.join([str(x) for x in c[:-1]])

        kw = {'input': Path(exists = 1),
              'format': format_choice,
              'output': Path(exists = 0),
              'assignments': Path(exists = 0),
              'proton1': ChoiceEntity(c[:-1], m % 'proton 1'),
              'proton2': ChoiceEntity(c[:-1], m % 'proton 2'),
              'hetero1': ChoiceEntity(c, m % 'hetero 1'),
              'hetero2': ChoiceEntity(c, m % 'hetero 2'),
              'peaks_ambiguity' : Path(exists = 0),
              }

        return kw

class CrossPeaksSettingsXMLPickler(_XMLBasePickler):

    def convert_dimension(self, d):

        if d: return int(d)
        else: return None

    def load_from_element(self, e):

        if str(e.input.filename).strip() == '':
            return None

        d = CrossPeaksSettings()

        d['input'] = abspath(str(e.input.filename))
        d['format'] = str(e.input.format).lower()
        d['output'] = abspath(str(e.output.filename))

        
        if hasattr(e, 'assignments'):
            d['assignments'] = abspath(str(e.assignments.filename))
        else:
            d['assignments'] = ''
            
        if hasattr(e, 'peaks_ambiguity'):
            d['peaks_ambiguity'] = abspath(str(e.peaks_ambiguity.filename))
        else:
            d['peaks_ambiguity'] = ''            
            
        d['proton1'] = self.convert_dimension(str(e.input.proton1))
        d['proton2'] = self.convert_dimension(str(e.input.proton2))
        d['hetero1'] = self.convert_dimension(str(e.input.hetero1))
        d['hetero2'] = self.convert_dimension(str(e.input.hetero2))

        return d

class SpectrumSettings(Settings):

    def create(self):

        from aria.Settings import TypeEntity, String, NonNegativeFloat, MultiTypeEntity

        kw = {}

        kw['chemical_shifts'] = TypeEntity('ChemicalShiftsSettings')
        kw['cross_peaks'] = TypeEntity('CrossPeaksSettings')
        kw['spectrum_name'] = String()
        kw['spectrum_type'] = String()

        ## BARDIAUX 2.2
        kw['spectrum_ambiguity'] = String()
        kw['segids'] = String()

        
        return kw

class SpectrumSettingsXMLPickler(_XMLBasePickler):

    def load_from_element(self, e):

        if e.chemical_shifts is None or e.cross_peaks is None:
            return None

        d = SpectrumSettings()

        d['chemical_shifts'] = e.chemical_shifts
        d['cross_peaks'] = e.cross_peaks
        d['spectrum_name'] = str(e.spectrum_name)
        d['segids'] = str(e.segids)
        
        ## compatibility with aria2.0

        if hasattr(e, 'spectrum_type'):
            d['spectrum_type'] = str(e.spectrum_type)
        else:
            d['spectrum_type'] = ''

        if hasattr(e, 'spectrum_ambiguity'):
            d['spectrum_ambiguity'] = str(e.spectrum_ambiguity)
        else:
            d['spectrum_ambiguity'] = SPECTRUM_AMBIGUITY_INTRA

            
        
        return d

class ConverterSettings(Settings):

    def create(self):

        from aria.Settings import TypeEntity, String

        kw = {'project_output_filename': String(),
              'project_name': String(),
              'molecule': TypeEntity('MoleculeSettings'),
              'spectrum': TypeEntity(LIST, TUPLE)}

        return kw

    def create_default_values(self):

        return {'spectrum': []}

class ConverterSettingsXMLPickler(_XMLBasePickler):

    def __init__(self):

        _XMLBasePickler.__init__(self)

        self.conversion = self
        self.molecule = MoleculeSettingsXMLPickler()
        self.spectrum = SpectrumSettingsXMLPickler()
        self.chemical_shifts = ChemicalShiftsSettingsXMLPickler()
        self.cross_peaks = CrossPeaksSettingsXMLPickler()
        
    def load_from_element(self, e):

        from aria.tools import as_tuple

        x = ConverterSettings()
        x.reset()

        project = e.project

        x['project_output_filename'] = str(project.output.filename)

        x['molecule'] = e.molecule
        [x['spectrum'].append(d) for d in as_tuple(e.spectrum) \
         if d is not  None]
            
        if hasattr(project, 'name'):

            project_name = str(project.name)

            if project_name == '':
                self.error('Project name (element <project>) missing.')
            
            x['project_name'] = project_name

        else:

            ## Compatibility with older conversion files
            ## Fake project name
        
            if x['molecule']['name'].strip() <> '':
                x['project_name'] = x['molecule']['name']
            else:
                x['project_name'] = 'my_struct'
                
        return x
              
class Converter(AriaBaseClass):

    def __init__(self):

        from aria.ConversionTable import ConversionTable
        from aria.Topology import load_topology
        from aria.AriaXML import AriaXMLPickler

        self.__pickler = AriaXMLPickler()
        self.__sequence = None

        Base.table = ConversionTable()
        Base.topology = load_topology()

        self.__export_ok = 1

    def create_factory(self):

        from aria.Singleton import AtomFactory

        return AtomFactory()

    def _convert_sequence(self):

        from aria.Molecule import Molecule
        from aria.tools import string_to_segid
        
        settings = self.getSettings()['molecule']
        
        # BARDIAUX 2.2 multimers
        segids = settings['segid']
        segids = segids.split('/')

        segids = [string_to_segid(segid) for segid in segids]
            
        if settings['format'] == 'seq':
            # chain_types = {settings['segid']: settings['type']}
            # BARDIAUX 2.2 multimers
            chain_types = {}
            for s in segids:                
                chain_types[s] = settings['type']
                
        elif settings['format'] == 'pdb':
            chain_types = None
            
        sequence = SequenceList(chain_types, settings['first_residue_number'])

        sequence.parse(settings['input'], settings['format'],
                       settings['naming_convention'])
        
        factory = self.create_factory()

        factory.reset()
        factory.unfreeze()

        chains = sequence.create_chains(factory)

        molecule = Molecule(settings['name'])
        # BARDIAUX 2.2
        # sort segids
        #[molecule.add_chain(chain) for chain in chains.values()]
        [molecule.add_chain(chains[seg]) for seg in segids]

        self.__sequence = sequence

        ## If no filename has been specified, skip
        ## pickling of XML file

        if settings['output'] == '':
            self.message('No output filename specified for Molecule ' + \
                         'definition. XML export Skipped.')
            self.__export_ok = 0
            return

        try:

            ## TIMING

            t0=clock()
            
            self.__pickler.dump(molecule, settings['output'])

        except Exception, msg:
            s = 'Could not write molecule XML file %s. Is the filename ' + \
                'correct? Do you have write permissions? Error message was %s.'
            
            self.error(Exception, s % (settings['output'], msg))

        m = 'Sequence file "%s" read, converted to XML and written to "%s".'
        
        self.message(m % (settings['input'], settings['output']))

    def _convert_chemicalshifts(self, number):

        from aria.tools import string_to_segid

        settings = self.getSettings()['spectrum'][number]['chemical_shifts']

        shifts = PpmList()
        shifts.parse(settings['input'], settings['format'])

        factory = self.create_factory()
        factory.freeze()

        naming_convention = CONVENTIONS[settings['format']]

        segids = self.getSettings()['spectrum'][number]['segids']
        segids = segids.split('/')

##         if not segids[0]:
##             segids = None
##         else:
        segids = [string_to_segid(segid) for segid in segids]
            
        shift_list = shifts.create_shiftlist(self.__sequence.bases, factory,
                                             naming_convention, segids)

        ## If no output filename has been specified,
        ## do not write XML file

        if settings['output'] == '':
            self.message('No output filename specified for chemical-' + \
                         'shift list. XML export skipped.')

            self.__export_ok = 0
            return

        try:

            ## TIMING
            
            t0=clock()
            
            self.__pickler.dump(shift_list, settings['output'])
            
            if self.display_debug:
                self.debug('shifts %s: %f' % (settings['input'],clock()-t0))
            
        except Exception, msg:
            s = 'Could not write chemical-shifts XML file %s. Is the ' + \
                'filename correct? Do you have write permissions?'
            self.error(Exception, s % settings['output'])

        m = 'Chemical-shift file "%s" read, converted to XML and written ' + \
            'to "%s".'
        self.message(m % (settings['input'], settings['output']))

    def _convert_crosspeaks(self, number):

        from aria.tools import string_to_segid
        
        settings = self.getSettings()['spectrum'][number]

        input = settings['cross_peaks']['input']
        output = settings['cross_peaks']['output']
        format = settings['cross_peaks']['format']

        ambig = settings['cross_peaks']['peaks_ambiguity']

        proton1 = settings['cross_peaks']['proton1']
        proton2 = settings['cross_peaks']['proton2']
        hetero1 = settings['cross_peaks']['hetero1']
        hetero2 = settings['cross_peaks']['hetero2']

        if format == settings['chemical_shifts']['format']:
            ppm_file = settings['chemical_shifts']['input']
        else:
            ppm_file = None

        noelist = NoeList()
        if os.path.exists(ambig):
            noelist.parse_peaks_ambiguity(ambig)

        noelist.parse(input, format, proton1, hetero1, proton2, hetero2,
                      ppm_file)

        factory = self.create_factory()
        factory.freeze()
        
        naming_convention = CONVENTIONS[format]

        segids = self.getSettings()['spectrum'][number]['segids']
        segids = segids.split('/')

        ## BARDIAUX 2.2
        ambiguity = settings['spectrum_ambiguity']
            
        if settings['spectrum_ambiguity'] not in SPECTRUM_AMBIGUITIES:
            self.message('Spectrum_ambiguity set to %s' % (SPECTRUM_AMBIGUITY_INTRA))
            settings['spectrum_ambiguity'] = SPECTRUM_AMBIGUITY_INTRA

##         if not segids[0]:
##             segids = None
##         else:
        segids = [string_to_segid(segid) for segid in segids]

        if settings['spectrum_ambiguity'] == SPECTRUM_AMBIGUITY_INTRA and len(segids) > 1:
            self.warning('Spectrum ambiguity %s and segids defintion %s are not compatible. ' % (ambiguity, '/'.join(segids))+ \
                         'Segid set to %s.' % segids[0])
            segids = [segids[0]]


        spectrum = noelist.create_spectrum(self.__sequence.bases, factory,
                                           naming_convention, segids, settings['spectrum_ambiguity']) ## BARDIAUX 2.2
        spectrum._setName(settings['spectrum_name'])

        if output == '':
            self.message('No output filename specified for spectrum XML.' + \
                         'XML export skipped.')
            self.__export_ok = 0
            return

        try:

            ## TIMING
            
            t0=clock()
            
            self.__pickler.dump(spectrum, output)

            if self.display_debug:
                self.debug('xpk %s: %f' % (settings['cross_peaks']['input'],clock()-t0))
            
        except Exception, msg:
            s = 'Could not write spectrum XML file %s. Is the filename ' + \
                ' correct? Do you have write permissions?'
            self.error(Exception, s % output)

        m = 'Peak list "%s" read, converted to XML and written to "%s".'
        self.message(m % (input, output))

    def _convert_spectrum(self, number):

        self._convert_chemicalshifts(number)
        self._convert_crosspeaks(number)

    def write_project(self):
        
        from aria.ariabase import PROJECT_TEMPLATE
        import os
        import aria.AriaXML as AriaXML
        import aria.DataContainer as DC
        
        project_filename = self.getSettings()['project_output_filename']
        if project_filename == '':
            return

        if not self.__export_ok:
            self.error('To generate a project XML file, all three ' + \
                       'data-types (molecule, chemical-shift-list and ' + \
                       'spectra) must be converted.')

        ## Read project template.

        template = os.path.join(AriaBaseClass.data_path, PROJECT_TEMPLATE)

        p = AriaXML.AriaXMLPickler()
        project = p.load_relaxed(template)

        ## Fill in filenames of all data-sources.

        settings = self.getSettings()
        
        project.getSettings()['file_root'] = settings['project_name']
        project.getSettings()['name'] = settings['project_name']

        ## Molecule definition

        s_mol = settings['molecule']

        ## Default sequence-data

        seq = project.getData(DC.DATA_SEQUENCE)[0]

        ## Remove default (empty) molecule-definition from project
        project.delData(seq)
        
        seq['filename'] = s_mol['output']
        seq['format'] = 'xml'

        project.addData(seq)

        ## BARDIAUX 2.2 multimers
        segids = self.getSettings()['molecule']['segid']
        n_seg = len(segids.split('/'))
        if n_seg > 1:

            sym_data = project.getData(DC.DATA_SYMMETRY)[0]
            project.delData(sym_data)
            sym_data['enabled'] = 'yes'
            symmetry_type = "C" + str(n_seg)
            if n_seg == "4":
                symmetry_type = "D2"
            sym_data['symmetry_type'] = symmetry_type
            sym_data['n_monomers'] = n_seg
            sym_data['packing_enabled'] = 'yes'
            sym_data['ncs_enabled'] = 'yes'
            project.addData(sym_data)
            
        ## Spectra / Chemical shift list

        for s_spectrum in settings['spectrum']:
            
            s_csl = s_spectrum['chemical_shifts']
            s_peaks = s_spectrum['cross_peaks']

            shifts = DC.ShiftData()
            shifts.reset()
            shifts['filename'] = s_csl['output']
            shifts['format'] = 'xml'
            
            peaks = DC.PeakData()
            peaks.reset()
            peaks['filename'] = s_peaks['output']
            peaks['format'] = 'xml'

            # BARDIAUX 2.2
            exp = DC.ExperimentData()
            exp.reset()
            exp['ambiguity_type'] = s_spectrum['spectrum_ambiguity']

            ## Create spectrum-data with default settings

            SD = DC.SpectrumData()
            SD.reset()
            SD['shifts'] = shifts
            SD['peaks'] = peaks

            # BARDIAUX 2.2
            SD['experiment_data'] = exp

            project.addData(SD)
            
        ## Save project

        if os.path.exists(project_filename):
            
            new_name = raw_input('Project with filename %s exists. Please enter new name (press Ctrl-C-C to cancel): ' % project_filename)

            project_filename = os.path.abspath(new_name)

        try:
            p.dump(project, project_filename)
            s = 'Project XML file %s written.'

        except Exception, msg:
            s = 'Could not write project XML file %s (error message was "%s"). Is the filename ' + \
                'correct? Do you have write permissions? ' % msg
            
        self.message(s % project_filename)

    def check_input_data(self):

        settings = self.getSettings()

        spec_names = {}
        peaks = {}
        shifts = {}
        
        for s in settings['spectrum']:

            name = s['spectrum_name']

            if not name in spec_names:
                spec_names[name] = 1

            else:
                spec_names[name] += 1

            output = s['cross_peaks']['output']

            if not output in peaks:
                peaks[output] = 1

            else:
                peaks[output] += 1

            output = s['chemical_shifts']['output']

            if not output in shifts:
                shifts[output] = 1

            else:
                shifts[output] += 1
        
        err = 0

        for name, n in spec_names.items():

            if n == 1:
                continue

            err = 1

            print 'Error: Multiple use of spectrum name "%s". Names must be unique.' % name

        for name, n in peaks.items():

            if n == 1:
                continue

            err = 1

            print 'Error: Multiple use of "%s" as name of spectrum XML file. Filenames must be unique.' % name

        for name, n in shifts.items():

            if n == 1:
                continue

            err = 1

            print 'Error: Multiple use of "%s" as name of chemical shift list XML file. Filenames must be unique.' % name

        return err
    
    def convert(self, xmlfile):

        settings = self.__pickler.load(xmlfile)
        self.setSettings(settings)

        ## check uniqueness of internal names and output (xml) name
        ## of all spectra

        err = self.check_input_data()

        if err:

            print

            import sys
            
            sys.exit(0)

        ## convert molecule

        self._convert_sequence()

        ## convert all spectra

        [self._convert_spectrum(i) for i in range(len(settings['spectrum']))]

        self.write_project()

    def create_ccpn_template(self):

        import os

        s = self.getSettings()
        
        abs_proj_filename = os.path.abspath(s['project_output_filename'])
        
        path, proj_filename = os.path.split(abs_proj_filename)
        proj_filename, ext = os.path.splitext(proj_filename)

        ccpn_filename = os.path.join(path, proj_filename + '_ccpn' + ext)

        proj_dict = {'project_name': s['project_name'],
                     'project_filename': ccpn_filename}

        translator = {'seq': 'XEasy',
                      'pdb': 'PDB',
                      'xeasy': 'XEasy',
                      'ansig': 'Ansig',
                      'nmrview': 'NmrView',
                      'pronto': 'Pronto',
                      'sparky': 'Sparky'}
        
        s_mol = s['molecule']

        mol_dict = {'mol_system_name': 'my_molecule',
                    'mol_name': s_mol['name'],
                    'mol_type': s_mol['type'],
                    'mol_segid': s_mol['segid'],
                    'mol_first_residue_number': s_mol['first_residue_number'],
                    'mol_input': s_mol['input'],
                    'mol_format': translator[s_mol['format']]}

        values = {'spectra': ''}
        values.update(proj_dict)
        values.update(mol_dict)

        str_spectra = ''

        for spec in s['spectrum']:

            spec_dict = {'name': spec['spectrum_name'],
                         'type': spec['spectrum_type'],
                         'segids': spec['segids'],
                         'format': translator[spec['cross_peaks']['format']],
                         'shifts_filename': spec['chemical_shifts']['input'],
                         'peaks_filename': spec['cross_peaks']['input'],
                         'assignments': spec['cross_peaks']['assignments']}

            for name in ('proton1', 'hetero1', 'proton2', 'hetero2'):

                val = spec['cross_peaks'][name]

                if val is None:
                    val = ''

                spec_dict[name] = str(val)

            str_spectra += CCPN_SPECTRUM_BLOCK % spec_dict

        values['spectra'] = str_spectra

        template = CCPN_XML_TEMPLATE % values

        return template

    def write_xml(self, molecule, peak_shift_lists):
        
        from aria.AriaXML import AriaXMLPickler

        pickler = AriaXMLPickler()

        settings = self.getSettings()

        try:
            pickler.dump(molecule, settings['molecule']['output'])
        except Exception, msg:
            s = 'Could not write molecule XML file %s. Is the filename ' + \
                'correct? Do you have write permissions?'
            self.error(Exception, s % settings['molecule']['output'])

        m = 'Sequence converted to XML and written to "%s".'
        self.message(m % settings['molecule']['output'])

        for i in range(len(settings['spectrum'])):

            peak_list, shift_list = peak_shift_lists[i]

            s = settings['spectrum'][i]['chemical_shifts']

            ## If no output filename has been specified,
            ## do not write XML file

            if s['output'] == '':
                self.message('No output filename specified for chemical-' + \
                             'shift list. XML export skipped.')

                self.__export_ok = 0
                return

            try:
                pickler.dump(shift_list, s['output'])
            except Exception, msg:
                s2 = 'Could not write chemical-shifts XML file "%s". Is ' + \
                    'the filename correct? Do you have write permissions?'
                self.error(Exception, s2 % s['output'])

            m = 'Chemical-shifts converted to XML and written ' + \
                'to "%s".'
            self.message(m % s['output'])

            ## write spectrum

            s = settings['spectrum'][i]['cross_peaks']

            peak_list._setName(settings['spectrum'][i]['spectrum_name'])

            if s['output'] == '':
                self.message('No output filename specified for spectrum ' + \
                             'XML. XML export skipped.')
                self.__export_ok = 0
                return

            try:
                pickler.dump(peak_list, s['output'])
            except Exception, msg:
                s2 = 'Could not write spectrum XML file "%s". Is ' + \
                    'the filename correct? Do you have write permissions?'
                self.error(Exception, s2 % s['output'])

            m = 'Peak list converted to XML and written to "%s".'
            self.message(m % s['output'])

    def convert_ccpn_objects(self, ccpn_project):
        
        from aria.Experiment import Experiment
        import aria.DataContainer as DC
        import aria.importFromCcpn as ccpnImport

        S = self.getSettings()

        molsystem_id = 'my_molecule', ' '

        exp_names = []

        project_name = S['project_name']

        for i in range(len(S['spectrum'])):

            peak_list_keys = project_name, i+1, 1, 1
            shift_list_keys = project_name, i+1

            exp_names.append((peak_list_keys, shift_list_keys))

        if exp_names or molsystem_id:
            
            ccpChain         = ccpnImport.getCcpnChains(ccpn_project, molsystem_id)[0]

            peak_shift_lists = ccpnImport.getCcpnPeakAndShiftLists(ccpn_project, ccpChain.molSystem,
                                                                   exp_names)

            if molsystem_id:
                molecule = ccpnImport.makeAriaMolecule(ccpChain.molSystem)

            aria_peak_shift_lists = []

            n = 0
            for peakList, shiftList in peak_shift_lists:

                aria_shift_list = ccpnImport.makeAriaShiftList(shiftList, ccpChain.molSystem, molecule)
                aria_spectrum   = ccpnImport.makeAriaSpectrum(peakList, molecule)
                # BARDIAUX 2.2
                spec_type = S['spectrum'][n]['ambiguity']
                aria_spectrum.setType(spec_type)
                aria_peak_shift_lists.append((aria_spectrum, aria_shift_list))
                n += 1

        self.write_xml(molecule, aria_peak_shift_lists)
                
    def convert_ccpn(self, xmlfile):

        import os
        from ccpnmr.format.converters import AriaXmlFormat        

        settings = self.__pickler.load(xmlfile)
        self.setSettings(settings)

        ## create XML read by format converter

        template = self.create_ccpn_template()

        template_filename = os.path.abspath('./ccpn_temp.xml')

        try:
            f = open(template_filename, 'w')
            f.write(template)
            f.close()

        except IOError:
            self.error('Unable to write temporary CCPNMR XML file. Do you have write permissions on the current directory?')

        aria_xml_format = AriaXmlFormat.AriaXmlFormat(template_filename)

        try:
            os.unlink(template_filename)
        except:
            pass

        ccpn_project = aria_xml_format.ccpNmrConv.ccpnProject

        self.convert_ccpn_objects(ccpn_project)

        self.write_project()

        self.message('Saving CCPN project ...')

        ccpn_project.saveModified()
        
    def write_template(self, filename):

        filename = os.path.expanduser(filename)

        if os.path.exists(filename):
            self.warning('Template file "%s" conversion exists ' % filename + \
                         'and will be overwritten.')

        try:
            file = open(filename, 'w')
        except:
            self.error('Could not create %s.' % filename)
            
        file.write(TEMPLATE_XML)
        file.close()

        self.message('Template xml file for data conversion has been ' + \
                     'written to "%s".' % filename)
        

def test_sparky():

    import os

    ## molecule settings

    m = MoleculeSettings()
    m['format'] = 'seq'
    m['input'] = os.path.expanduser('~/projects/aria2.0/src/py/legacy/' + \
                                    'example/sparky/sparky.seq')
    m['output'] = '/tmp/sparky.xml'
    m['type'] = 'PROTEIN'
    m['segid'] = '    '
    m['first_residue_number'] = 1
    m['naming_convention'] = 'iupac'
    m['name'] = 'sparky'

    s = ChemicalShiftsSettings()
    s['input'] = os.path.expanduser('~/projects/aria2.0/src/py/legacy/' + \
                                    'example/sparky/sparky_ppmlist.list')
    s['output'] = '/tmp/ppm.xml'
    s['format'] = 'sparky'

    n = CrossPeaksSettings()
    n['input'] = os.path.expanduser('~/projects/aria2.0/src/py/legacy/' + \
                                    'example/sparky/3D_13C_NOESY.list')#2D_NOESY.list')
    n['format'] = 'sparky'
    n['output'] = '/tmp/noesy.xml'
    n['proton1'] = 1
    n['proton2'] = 3
    n['hetero1'] = None
    n['hetero2'] = 2

    c = ConverterSettings()
    c['molecule'] = m
    c['spectrum'] = [SpectrumSettings()]
    c['spectrum'][0]['chemical_shifts'] = s
    c['spectrum'][0]['cross_peaks'] = n
    c['spectrum'][0]['segids'] = '    '    
    c['spectrum'][0]['spectrum_name'] = '3D_NOESY'

    converter = Converter()
    converter.setSettings(c)
    converter._convert_sequence()
    converter._convert_spectrum(0)

    return converter
