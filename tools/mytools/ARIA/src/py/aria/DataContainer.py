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



from aria.ariabase import *
from aria.Settings import Settings
from aria.xmlutils import XMLElement, XMLBasePickler

DATA_DEFAULT = 'default'
DATA_SPECTRUM = 'spectrum'
DATA_SEQUENCE = 'sequence'
DATA_HBONDS = 'hbonds'
DATA_DIHEDRALS = 'dihedrals'
DATA_KARPLUS = 'karplus'
DATA_RDCS = 'rdcs'
DATA_SSBONDS = 'ssbonds'
DATA_SSBRIDGE = 'ssbridge'
DATA_HISPATCH = 'hispatch'
DATA_AMBIGUOUS = 'ambiguous'
DATA_UNAMBIGUOUS = 'unambiguous'
DATA_ITERATION = 'iteration'
DATA_ANNEALING = 'annealing'
DATA_DYNAMICS = 'dynamics'
DATA_SHIFTS = 'shifts'
DATA_PEAKS = 'peaks'
DATA_TEMPLATE_STRUCTURE = 'template_structure'
DATA_INITIAL_STRUCTURE = 'initial_structure'
# BARDIAUX 2.2
DATA_SYMMETRY = 'symmetry'
DATA_EXPERIMENT = 'experiment_data'
DATA_CISPROPATCH = 'cispropatch'
DATA_CYSPATCH = 'cyspatch'
DATA_ZNPATCH = 'znpatch'
# BARDIAUX 2.3
DATA_OTHER = 'other_data'

DATA_CCPN = 'ccpn'

DATA_ANNEALING_AMBIG = 'ambiguous_annealing'
DATA_ANNEALING_UNAMBIG = 'unambiguous_annealing'
DATA_ANNEALING_DIHEDRAL = 'dihedral_annealing'
DATA_ANNEALING_KARPLUS = 'karplus_annealing'
DATA_ANNEALING_RDC = 'rdc_annealing'
DATA_ANNEALING_HBOND = 'hbond_annealing'
DATA_ANNEALING_FBHW = 'fbhw_annealing'
# BARDIAUX
DATA_ANNEALING_SYM = 'symmetry_annealing'
# BERNARD Aymeric
DATA_ANNEALING_LOGHARMONIC = 'logharmonic_annealing'

## experimental data types

DATA_TYPES = (DATA_SPECTRUM, DATA_SEQUENCE, DATA_HBONDS, DATA_DIHEDRALS,
              DATA_KARPLUS, DATA_RDCS, DATA_SSBONDS, DATA_TEMPLATE_STRUCTURE,
              DATA_INITIAL_STRUCTURE,
              DATA_SSBRIDGE, DATA_HISPATCH, DATA_AMBIGUOUS, DATA_UNAMBIGUOUS,
              DATA_SYMMETRY, DATA_EXPERIMENT, DATA_CISPROPATCH, DATA_CYSPATCH, DATA_ZNPATCH, # BARDIAUX 2.2
              DATA_OTHER) # BARDIAUX 2.3

TYPES = DATA_TYPES + (DATA_DEFAULT, DATA_ITERATION, DATA_ANNEALING,
                      DATA_DYNAMICS,
                      DATA_SHIFTS, DATA_PEAKS, DATA_ANNEALING_AMBIG,
                      DATA_ANNEALING_UNAMBIG, DATA_ANNEALING_DIHEDRAL,
                      DATA_ANNEALING_KARPLUS, DATA_ANNEALING_RDC,
                      DATA_ANNEALING_HBOND, DATA_ANNEALING_FBHW,
                      DATA_ANNEALING_SYM, # BARDIAUX 2.2
                      DATA_ANNEALING_LOGHARMONIC )# BERNARD Aymeric

class DataContainer(Settings):

    def __init__(self, type = DATA_DEFAULT, keywords = None,
                 default_settings = None):

        check_string(type)

        if type not in TYPES:
            self.error(TypeError, 'Type "%s" not supported.' % type)

        Settings.__init__(self, keywords, default_settings)

        self.__type = type

    def create(self):
        return {}

    def getType(self):
        return self.__type

class SimpleDataContainer(DataContainer):

    known_formats = ()

    def create(self):

        from aria.Settings import Path, ChoiceEntity, NonEmptyString

        keywords = DataContainer.create(self)
        keywords['filename'] = Path(exists=-1)

        msg = 'Format "%s" not known. Supported formats: ' + \
              ' / '.join(self.known_formats)
        descr = 'Specifies the format in which your data is represented. This can either be ARIA XML format ("xml") or a CCPN data model ("ccpn").'
        keywords['format'] = ChoiceEntity(self.known_formats,
                                          description = descr,
                                          error_message = msg)

        descr = 'If set to "%s", the data will be used.' % str(YES)
        keywords['enabled'] = ChoiceEntity([YES, NO], description = descr)

        keywords['ccpn_id'] = NonEmptyString(description = 'Internal name used to retrieve data from a CCPN data model. In order to read data from a CCPN data model, select "CCPN" as format and specifiy the file that contains the CCPN data model (Node "CCPN").')

        return keywords

    def create_default_values(self):
        return {'ccpn_id': ''}
        
    def __setitem__(self, key, value):
        if key == 'format' and is_type(value, STRING):
            value = value.lower()

        DataContainer.__setitem__(self, key, value)

    def getLocation(self):
        return self['filename'], self['format']

class SequenceData(SimpleDataContainer):

    known_formats = ('xml','ccpn')

    ## Standard DNA/RNA/Protein linkage-
    ## and topology files. Keys: names used in xml-file,
    ## values: labels / descriptions

    USER_DEFINED = 'user_defined'
    AUTOMATIC = 'automatic'
    
    linkage_names = {'dna-rna.link': 'DNA-RNA',
                     'dna-rna-pho.link': 'DNA-RNA-PHO',
                     'topallhdg.pep': 'TOPALLHDG',
                     'topallhdg5.3.pep': 'TOPALLHDG5.3',
                     AUTOMATIC: 'Automatic',
                     USER_DEFINED: 'User defined'}
    
    topology_names = {'dna-rna-allatom.top': 'DNA-RNA',
                      'topallhdg5.0.pro': 'TOPALLHDG5.0',
                      'topallhdg5.1.pro': 'TOPALLHDG5.1',
                      'topallhdg5.2.pro': 'TOPALLHDG5.2',
                      'topallhdg5.3.pro': 'TOPALLHDG5.3',
		      'topallhdg5.3soft.pro': 'TOPALLHDG5.3soft',
                      AUTOMATIC: 'Automatic',
                      USER_DEFINED: 'User defined'}

    parameter_names = {'dna-rna-allatom.param': 'DNA-RNA',
                       'parallhdg5.0.pro': 'PARALLHDG5.0',
                       'parallhdg5.1.pro': 'PARALLHDG5.1',
                       'parallhdg5.2.pro': 'PARALLHDG5.2',
                       'parallhdg5.3.pro': 'PARALLHDG5.3',
		       'parallhdg5.3softbig.pro': 'PARALLHDG5.3soft',
                       AUTOMATIC: 'Automatic',
                       USER_DEFINED: 'User defined'}

    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_SEQUENCE)
        del self['enabled']

    def create(self):
        from aria.Settings import ChoiceEntity, Path
        import aria.Chain as Chain

        d = SimpleDataContainer.create(self)

        ## Linkage Definition
        
        descr = \
'''ARIA supports a number of pre-defined Linkage definitions.

<%(auto)s>: Every chain (specified in the molecule-description XML-file) owns an attribute, "chain_type", which can be set to "PROTEIN", "DNA" or "RNA". The linkage defs are determined automatically from the type:
        
%(protein)s:\ttopallhdg5.3.pep,
%(rna)s, %(dna)s:\tdna-rna.link

<%(user_defined)s>: the user has to specify a proper linkage file. During project setup, the file is then copied to the local CNS topology/parameters directory PROJECT_PATH/RUNxxx/cns/toppar.'''

        names = {'auto': self.linkage_names[SequenceData.AUTOMATIC],
                 'protein': Chain.TYPE_PROTEIN,
                 'rna': Chain.TYPE_RNA,
                 'dna': Chain.TYPE_DNA,
                 'user_defined': self.linkage_names[SequenceData.USER_DEFINED]}
               
        e = ChoiceEntity(elements = self.linkage_names.keys(),
                         description = descr % names)

        d['linkage_name'] = e

        descr = 'If the Linkage definition name has been set to "%(user_defined)s", ARIA copies the user-defined linkage file to the local "toppar" directory of the project.'
        d['linkage_filename'] = Path(description = descr % names,
                                     exists = 0)

        ## Topology Definition

        descr = \
'''ARIA supports a number of pre-defined Topology definitions.

<%(auto)s>: Every chain (specified in the molecule-description XML-file) owns an attribute, "chain_type", which can be set to "PROTEIN", "DNA" or "RNA". The linkage defs are determined automatically from the type:
        
%(protein)s:\ttopallhdg5.3.pro,
%(rna)s, %(dna)s:\tdna-rna-allatom.link

<%(user_defined)s>: the user has to specify a proper topology file. During project setup, the file is then copied to the local CNS topology/parameters directory PROJECT_PATH/RUNxxx/cns/toppar.'''

        names['auto'] = self.topology_names[SequenceData.AUTOMATIC]
        names['user_defined'] = self.topology_names[SequenceData.USER_DEFINED]
               
        e = ChoiceEntity(elements = self.topology_names.keys(),
                         description = descr % names)
        d['topology_name'] = e

        descr = 'If the Topology definition name has been set to "%(user_defined)s", ARIA copies the user-defined topology file to the local "toppar" directory of the project.'
        d['topology_filename'] = Path(description = descr % names,
                                      exists = 0)

        ## Parameter Definition

        descr = \
'''ARIA supports a number of pre-defined Parameter definitions.
<%(auto)s>: Every chain (specified in the molecule-description XML-file) owns an attribute, "chain_type", which can be set to "PROTEIN", "DNA" or "RNA". The linkage defs are determined automatically from the type:
        
%(protein)s:\ttopallhdg5.3.pro,
%(rna)s, %(dna)s:\tdna-rna-allatom.param

<%(user_defined)s>: the user has to specify a proper parameter file. During project setup, the file is then copied to the local topology/parameters directory PROJECT_PATH/RUNxxx/cns/toppar.'''

        names['auto'] = self.parameter_names[SequenceData.AUTOMATIC]
        names['user_defined'] = self.parameter_names[SequenceData.USER_DEFINED]
               
        e = ChoiceEntity(elements = self.parameter_names.keys(),
                         description = descr % names)
        
        d['parameter_name'] = e

        descr = 'If the Parameter definition name has been set to "%(user_defined)s", ARIA copies the user-defined parameter file to the local "toppar" directory of the project.'
        d['parameter_filename'] = Path(description = descr % names,
                                        exists = 0)
        
        return d

    def create_default_values(self):
        d = {}
        
        d['format'] = 'xml'

        d['linkage_name'] = SequenceData.AUTOMATIC
        d['topology_name'] = SequenceData.AUTOMATIC
        d['parameter_name'] = SequenceData.AUTOMATIC
        d['linkage_filename'] = ''
        d['topology_filename'] = ''
        d['parameter_filename'] = ''

        return d

class HBondData(SimpleDataContainer):

    known_formats = ('tbl', 'ccpn',)
    
    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_HBONDS)

    def create(self):

        from aria.Settings import ChoiceEntity

        d = SimpleDataContainer.create(self)

        choices = ['standard', 'csi']
        msg = 'Possible types of H-bonds: ' + ', '.join(choices)
        d['type'] = ChoiceEntity(choices, error_message = msg)

        return d
        
    def create_default_values(self):

        d = SimpleDataContainer.create_default_values(self)

        d.update({'format': 'tbl', 'type': 'standard', 'enabled': YES})

        return d
    
# BARDIAUX 2.2
class AmbiguousDistanceData(SimpleDataContainer):

    known_formats = ('tbl', 'ccpn',)

    def __init__(self):

        SimpleDataContainer.__init__(self, DATA_AMBIGUOUS)

    # BARDIAUX 2.2
    def create(self):
        
        from aria.Settings import ChoiceEntity
        
        d = SimpleDataContainer.create(self)
        
        warn = "(This feature requires restraints from a CCPN project.)"
        
        descr = 'If set to "%s", the restraints will be added to the network of assignments for the Network-Anchoring Analysis. The restraints will stay unchanged. %s' % (str(YES), warn)
        d['add_to_network'] = ChoiceEntity([YES, NO], description = descr)
        
        descr = 'Determine how the restraints should be calibrate by ARIA.\n"all_iterations" : the restraints will be calibrated at each iteration.\n"all_iterations_except_first" : the restraints will  be calibrated at all iteration except the first one.\n"no": the restraints will never be calibrated by ARIA. %s' % warn
        d['calibrate'] = ChoiceEntity(['all_iterations', 'all_iterations_except_first', NO], description = descr)

        descr = 'The restraints will enter the ARIA network-anchoring analysis. %s' % warn
        d['run_network_anchoring'] = ChoiceEntity([YES, NO], description = descr)
        
        descr = 'The restraints will enter the ARIA violation analysis machinery to filter restraints according to the violation analysis and the contributions in order to reduce ambiguity. %s' % warn
        d['filter_contributions'] = ChoiceEntity([YES, NO], description = descr)        
     

        return d
    
    def create_default_values(self):

        d = SimpleDataContainer.create_default_values(self)

        d.update({'format': 'tbl'})

        # BARDIAUX 2.2
        d.update({'enabled': YES}) # Thanks to A. Wilter
        d.update({'add_to_network' : NO})
        d.update({'calibrate' : NO})
        d.update({'run_network_anchoring' : NO})
        d.update({'filter_contributions' : NO})
        
        return d
    
# BARDIAUX 2.2    
class UnambiguousDistanceData(SimpleDataContainer):

    known_formats = ('tbl', 'ccpn',)

    def __init__(self):

        SimpleDataContainer.__init__(self, DATA_UNAMBIGUOUS)

    # BARDIAUX 2.2
    def create(self):
        
        from aria.Settings import ChoiceEntity
        
        d = SimpleDataContainer.create(self)

        warn = "(This feature requires restraints from a CCPN project.)"
        
        descr = 'If set to "%s", the restraints will be added to the network of assignment for the Network-Anchoring Analysis. The restraints will stay unchanged. %s' % (str(YES), warn)
        d['add_to_network'] = ChoiceEntity([YES, NO], description = descr)
        
        descr = 'Determine how the restraints should be calibrate by ARIA.\n"all_iterations" : the restraints will be calibrated at each iteration.\n"all_iterations_except_first" : the restraints will not be calibrated at all iteration except the first one.\n"no": the restraints will never be calibrated by ARIA. %s' % warn
        d['calibrate'] = ChoiceEntity(['all_iterations', 'all_iterations_except_first', NO], description = descr)

        descr = 'The restraints will enter the ARIA network-anchoring analysis. %s' % warn
        d['run_network_anchoring'] = ChoiceEntity([YES, NO], description = descr)
        
        descr = 'The restraints will enter the ARIA violation analysis machinery to filter restraints according to the violation analysis and the contributions in order to reduce ambiguity. %s' % warn
        d['filter_contributions'] = ChoiceEntity([YES, NO], description = descr)  

        
        return d
    
    def create_default_values(self):

        d = SimpleDataContainer.create_default_values(self)

        d.update({'format': 'tbl'})

        # BARDIAUX 2.2
        d.update({'enabled': YES}) # Thanks to A. Wilter
        d.update({'add_to_network' : NO})
        d.update({'calibrate' : NO})
        d.update({'run_network_anchoring' : NO})
        d.update({'filter_contributions' : NO})
        
        return d

class DihedralData(SimpleDataContainer):

    known_formats = ('tbl', 'ccpn',)

    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_DIHEDRALS)

    def create(self):

        from aria.Settings import ChoiceEntity

        d = SimpleDataContainer.create(self)

        choices = ['standard', 'talos', 'csi']
        msg = 'Choices for type of dihedral angle restraints: ' + \
              ', '.join(choices)
        d['type'] = ChoiceEntity(choices, error_message = msg)

        return d

    def create_default_values(self):

        d = SimpleDataContainer.create_default_values(self)

        d.update({'format': 'tbl', 'type': 'standard', 'enabled': YES})

        return d

class KarplusData(SimpleDataContainer):

    known_formats = ('tbl', 'ccpn',)
    
    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_KARPLUS)

    def create(self):

        from aria.Settings import ChoiceEntity

        d = SimpleDataContainer.create(self)

        choices = [1, 2, 3, 4, 5]
        msg = 'Choices for class of Karplus restraint: ' + \
              ', '.join(map(str, choices))
        d['class'] = ChoiceEntity(choices, error_message = msg)

        return d

    def create_default_values(self):

        d = SimpleDataContainer.create_default_values(self)

        d.update({'format': 'tbl', 'class': 1, 'enabled': YES})

        return d

class RDCData(SimpleDataContainer):

    known_formats = ('tbl', 'ccpn', )

    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_RDCS)

    def create(self):

        from aria.Settings import ChoiceEntity

        d = SimpleDataContainer.create(self)

        choices = (1, 2, 3, 4, 5)
        msg = 'Choices for class of RDC restraint: ' + \
              ', '.join(map(str, choices))
        
        descr = """Aria groups residual dipolar coupling (RDC) measurements into 5 classes. Every class has its own parameter set which the user has to specify. It contains the specification of the alignment tensor and settings for the simulated annealing protocol."""
        
        d['class'] = ChoiceEntity(choices, description = descr,
                                  error_message = msg)

        return d

    def create_default_values(self):
        
        d = SimpleDataContainer.create_default_values(self)

        d.update({'format': 'tbl', 'class': 1, 'enabled': YES})

        return d

# BARDIAUX 2.3
class OtherData(SimpleDataContainer):

    known_formats = ('tbl', ) #'ccpn',)

    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_OTHER)

    def create(self):

        from aria.Settings import ChoiceEntity

        d = SimpleDataContainer.create(self)

        choices = ['planarity',]
        msg = 'Choices for type of restraints: ' + \
              ', '.join(choices)
        d['type'] = ChoiceEntity(choices, error_message = msg)

        return d

    def create_default_values(self):

        d = SimpleDataContainer.create_default_values(self)

        d.update({'format': 'tbl',
                  'type': 'planarity',
                  'enabled': YES})

        return d

# BARDIAUX 2.2: add DISN for SSBond and ss_ambigous
class CysPatch(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_CYSPATCH)

    def create(self):

        from aria.Settings import Integer, FourLetterString

        segid_err_msg = 'Segid must be string of length 4'

        d = DataContainer.create(self)
        d['residue'] = Integer()
        d['segid'] = FourLetterString(error_message = segid_err_msg)
        return d

    def create_default_values(self):
        return {'segid': ''}
        
class SSBondData(SimpleDataContainer):

    known_formats = ('tbl', 'ccpn', )

    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_SSBONDS)
        self.reset()

    # BARDIAUX 2.2: add DISU for SSBond and ss_ambigous
    def create(self):

        from aria.Settings import TypeEntity
        
        d = SimpleDataContainer.create(self)
        
        entity = TypeEntity(LIST)
        entity.set([])

        d['cyspatch'] = entity

        return d

    def create_default_values(self):

       d = SimpleDataContainer.create_default_values(self)

       d.update({'format': 'tbl', 'enabled': YES})
       d.update({'cyspatch': []})

       return d

class SSBridge(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_SSBRIDGE)

    def create(self):

        from aria.Settings import Integer, FourLetterString

        segid_err_msg = 'Segid must be string of length 4'

        d = DataContainer.create(self)
        d['residue1'] = Integer()
        d['segid1'] = FourLetterString(error_message = segid_err_msg)
        d['residue2'] = Integer()
        d['segid2'] = FourLetterString(error_message = segid_err_msg)

        return d

    def create_default_values(self):
        return {'segid1': '', 'segid2': ''}

class HisPatch(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_HISPATCH)

    def create(self):

        from aria.Settings import Integer, FourLetterString, ChoiceEntity

        segid_err_msg = 'Segid must be string of length 4'

        d = DataContainer.create(self)
        d['residue'] = Integer()
        d['segid'] = FourLetterString(error_message = segid_err_msg)
        d['proton'] = ChoiceEntity(('HISD', 'HISE'))
        
        return d

    def create_default_values(self):
        return {'segid': '', 'proton': 'HISD'}

## BARDIAUX 2.2
class CisProPatch(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_CISPROPATCH)

    def create(self):

        from aria.Settings import Integer, FourLetterString

        segid_err_msg = 'Segid must be string of length 4'

        d = DataContainer.create(self)
        d['residue'] = Integer()
        d['segid'] = FourLetterString(error_message = segid_err_msg)
        return d

    def create_default_values(self):
        return {'segid': ''}

# BARDIAUX 2.2: ZN patch
class ZnPatch(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_ZNPATCH)

    def create(self):

        from aria.Settings import Integer, FourLetterString, ChoiceEntity

        segid_err_msg = 'Segid must be string of length 4'
        
        choices = ("SSSS", "SSSE", "SSSD")
        msg = 'Choices for the Zinc coordination: ' + \
              ', '.join(map(str, choices))
        descr = "Cys4: 4 Cysteine residues\n" + \
                "Cys3Hisd : 3 Cysteine residues and 1 Histidine HISD\n" + \
                "Cys3Hise : 3 Cysteine residues and 1 Histidine HISE"
        
        d = DataContainer.create(self)
        d['type'] = ChoiceEntity(choices, description = descr, error_message = msg)        
        d['residue_zn'] = Integer()       
        d['segid_zn'] = FourLetterString(error_message = segid_err_msg)        
        d['residue1'] = Integer()       
        d['segid1'] = FourLetterString(error_message = segid_err_msg)
        d['residue2'] = Integer()
        d['segid2'] = FourLetterString(error_message = segid_err_msg)
        d['residue3'] = Integer()
        d['segid3'] = FourLetterString(error_message = segid_err_msg)
        d['residue4'] = Integer()
        d['segid4'] = FourLetterString(error_message = segid_err_msg)

        return d

    def create_default_values(self):
        return {'type' : 'SSSS',
                'segid_zn' :'',
                'segid1': '', 
                'segid2': '', 
                'segid3': '', 
                'segid4': ''}
    
## BARDIAUX  : Symmetry settings for multimers 2.2
class Symmetry(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_SYMMETRY)

    def create(self):

        from aria.Settings import Integer, YesNoChoice, ChoiceEntity


        d = DataContainer.create(self)

        descr = """If enabled, aria will treat this molecule as a symmetric multimer."""
        d['enabled'] =   YesNoChoice(description = descr)

        # BARDIAUX 2.3
        choices = ("standard", )
        msg = 'Choices for the symmery method: ' + \
              ', '.join(map(str, choices))
        descr = """ The method that will be used to ensure the symmetry\n.""" +\
                """ "Standard" : symmetry-ADR (Nilges, Proteins, 1993) """  
        d['method'] =   ChoiceEntity(choices, description = descr,
                                     error_message = msg)
        
        descr = """You can define the number of monomers present in the molecule."""
        d['n_monomers'] = Integer(description = descr)
        
        choices = ("None", "C2", "C3", "C5", "D2")
        msg = 'Choices for the symmetry type: ' + \
              ', '.join(map(str, choices))
        descr = """ Type of Symmetry."""        
        d['symmetry_type'] = ChoiceEntity(choices, description = descr,
                                         error_message = msg)

        descr = """ NCS restraints are applied to minimize the RMSD between monomers during the refinment."""
        d['ncs_enabled'] =   YesNoChoice(description = descr)
        
        descr = """ Packing restraints are applied to keep the monomers close to each-other. If think you have enough inter-monomer restraints, just disable."""
        d['packing_enabled'] = YesNoChoice(description = descr)
        
        return d

    def create_default_values(self):
        d = {}


        d['enabled'] =   NO
        d['method'] =   'standard' 
        d['n_monomers'] = 1
        d['symmetry_type'] = "None"       
        d['ncs_enabled'] =   NO   
        d['packing_enabled'] = NO
        
        return d
    
class ShiftData(SimpleDataContainer):

    known_formats = ('xml', 'ccpn')

    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_SHIFTS)

        del self['enabled']

    def create(self):

        from aria.Settings import NonNegativeFloat

        d = SimpleDataContainer.create(self)

        descr = 'Default error [ppm] of chemical-shift measurements' + \
                ' for which an error is not available.'
        
        d['default_shift_error'] = NonNegativeFloat(description = descr)

        return d
        
    def create_default_values(self):
        
        d = SimpleDataContainer.create_default_values(self)
        
        d['format'] = 'xml'
        d['default_shift_error'] = 0.00

        return d

class BoundCorrection(Settings):

    def create(self):
        from aria.Settings import NonNegativeFloat, YesNoChoice
        
        d = {}
        
        d['value'] = NonNegativeFloat()

        descr = 'If set to "%s", distance bounds will possibly be ' + \
                'changed in the course of structure calculation.' 
        
        d['enabled'] = YesNoChoice(description = descr % str(YES))

        return d

    def create_default_values(self):

        d = {}
        d['enabled'] = NO

        return d

class LowerBoundCorrection(BoundCorrection):
    def create(self):
        d = BoundCorrection.create(self)

        descr = \
"""
When performing a Violation Analysis on an ensemble of structures, ARIAs default approach is to deactivate a restraint for the structure calculation in the subsequent iteration, if the restraint has been violated in more than a certain fraction of all structures.
However, if 'bound-correction' is turned on, ARIA runs a 2-pass Violation Analysis: the 1st pass determines the set of violated restraints which is then analysed further: if the restraint is still violated with respect to the new user-defined LOWER bound, nothing changes. If not, the restraint will be used (with modified lower bound) in the subsequent iteration. 
"""
        
        d['value'].setDescription(descr)

        return d
        
    def create_default_values(self):
        d = BoundCorrection.create_default_values(self)
        d['value'] = 0.

        return d

class UpperBoundCorrection(BoundCorrection):
    def create(self):
        d = BoundCorrection.create(self)

        descr = \
"""
When performing a Violation Analysis on an ensemble of structures, ARIAs default approach is to deactivate a restraint for the structure calculation in the subsequent iteration, if the restraint has been violated in more than a certain fraction of all structures.
However, if 'bound-correction' is turned on, ARIA runs a 2-pass Violation Analysis: the 1st pass determines the set of violated restraints which is then analysed further: if the restraint is still violated with respect to the new user-defined UPPER bound, nothing changes. If not, the restraint will be used (with modified upper bound) in the subsequent iteration. 
"""
        
        d['value'].setDescription(descr)

        return d
        
    def create_default_values(self):
        d = BoundCorrection.create_default_values(self)
        d['value'] = 6.

        return d
    
class PeakData(SimpleDataContainer):

    known_formats = ('xml', 'ccpn')

    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_PEAKS)

    def create(self):
        from aria.Settings import NonNegativeFloat, PeakType, TypeEntity
        
        d = SimpleDataContainer.create(self)

        s = '%s-shift error must be a non-negative float.'

        q = {'proton1_shift_err': 
             NonNegativeFloat(error_message = s % 'Proton1'),
             'proton2_shift_err':
             NonNegativeFloat(error_message = s % 'Proton2'),
             'hetero1_shift_err':
             NonNegativeFloat(error_message = s % 'Hetero1'),
             'hetero2_shift_err':
             NonNegativeFloat(error_message = s % 'Hetero2'),
             'volume_or_intensity': PeakType(),
             'lower_bound_correction':
             TypeEntity('LowerBoundCorrection'),
             'upper_bound_correction':
             TypeEntity('UpperBoundCorrection')}

        d.update(q)
        del d['enabled']

        return d

    def create_default_values(self):

        d = SimpleDataContainer.create_default_values(self)
        
        d['format'] = 'xml'
        d['proton1_shift_err'] = 0.04
        d['hetero1_shift_err'] = 0.5        
        d['proton2_shift_err'] = 0.02
        d['hetero2_shift_err'] = 0.5
        d['volume_or_intensity'] = 'volume'

        l_c = LowerBoundCorrection()
        u_c = UpperBoundCorrection()
        
        l_c.reset()
        u_c.reset()
        
        d['lower_bound_correction'] = l_c
        d['upper_bound_correction'] = u_c

        return d

## BARDIAUX rMat
class ExperimentData(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_EXPERIMENT)
        self.reset()
        
    def create(self):
        
        from aria.Settings import NonNegativeFloat
        from aria.TypeChecking import FLOAT, NONE
        
        d = DataContainer.create(self)

        m = """The rotation correlation time of the molecule, in nano-seconds [ns]. Required for the Spin Diffusion Correction.\n\nMust be zero (0.0) if the user doesn't want to use the  Spin Diffusion Correction."""
        d['molecule_correlation_time'] = NonNegativeFloat(description = m)

        m = """The mixing time for transfer step, in milli-seconds [ms]. Required for the Spin Diffusion Correction.\n\nMust be zero (0.0) if the user doesn't want to use the  Spin Diffusion Correction."""
        
        d['spectrum_mixing_time'] = NonNegativeFloat(description = m)

        m = """The spectrometer frequency, in Mega-Hertz [MHz] where the spectrum was recorded. Required for the Spin Diffusion Correction.\n\nMust be zero (0.0) if the user doesn't want to use the  Spin Diffusion Correction."""
        d['spectrometer_frequency'] = NonNegativeFloat(description = m)

        # inter/intra ambig
        from aria.Settings import ChoiceEntity
        m = """ Ambiguity level of the spectra in terms of inter/intra molecular assignments. Possible values are\n
        - intra    (noes involving atoms from one monomer only)
        - inter    (noes involving atoms from different monomers)
        - all      (no known information, all noes are ambigous in terms of monomer)
        """
        
        c = ('intra', 'inter', 'all')
        d['ambiguity_type'] = ChoiceEntity(c, description = m)

 
        return d
    
    def create_default_values(self):

        d = {}

        d['molecule_correlation_time'] = .0
        d['spectrum_mixing_time'] = .0
        d['spectrometer_frequency'] = .0       

        d['ambiguity_type'] = 'intra'
        
        return d
    
class SpectrumData(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_SPECTRUM)
        self.reset()

    def create(self):
        
        from aria.Settings import TypeEntity, YesNoChoice

        d = DataContainer.create(self)

        m = 'Use spectrum in structure calculation: %s/%s.' % (str(YES), str(NO))
        descr = 'If set to "%s", any existing peak assignment is discarded and peaks are re-assigned automatically. This option overrides the attribute "trust_assigned_peaks".' % str(NO)

        d['shifts'] = TypeEntity('ShiftData')
        d['peaks'] = TypeEntity('PeakData')
        d['enabled'] = YesNoChoice(description = m)
        d['use_assignments'] = YesNoChoice(description = descr)

        descr = \
'''If set to "%s", fully assigned peaks are always used for structure calculation and can never removed from the restraint list - even if the violation-analysis classifies it as violated. If set to "%s", fully assigned peaks are treated in the normal way, i.e. they can be removed from the restraint list when violated but its assignments are never changed.

Note: it is also possible to control that behaviour for every single cross-peak separately. Every cross-peak owns an attribute, "reliable" (cf. spectrum-xml). If set to "%s", the respective peak gets never rejected by the violation-analysis.''' % (str(YES), str(NO), str(YES))

        d['trust_assigned_peaks'] = YesNoChoice(description = descr)
        
        ## BARDIAUX 2.2
        descr = """Structural Rules are only valid for symmetric multimers assignments. If the 2 implicated atoms belong to the same secondary structure element and if they are separated by more than 5 residues (helix) or 4 residdues (beta strands), the assignement possibilty could be inter-molecular only. 
Secondary structure definition is taken from the "structure" attribute of each residue.

These filter is not use is the peak is set as reliable."""

        d['structural_rules_enabled'] = YesNoChoice(description = descr)

        ## BARDIAUX 2.2
        descr = """ If enabled, diagonal peaks will be removed."""

        d['filter_diagonal_peaks'] = YesNoChoice(description = descr)
        
        ## BARDIAUX 2.2
        descr = """ If set to %s, unassigned peaks will be removed.""" % (str(YES))

        d['filter_unassigned_peaks'] = YesNoChoice(description = descr)
        
        ## BARDIAUX rMat
        d['experiment_data'] = TypeEntity('ExperimentData')        

        return d
    
    def create_default_values(self):

        shifts = ShiftData()
        shifts.reset()

        peaks = PeakData()
        peaks.reset()

        ## BARDIAUX rMat
        experiment_data = ExperimentData()
        experiment_data.reset()
        
        return {'shifts': shifts,
                'peaks': peaks,
                'enabled': YES,
                'use_assignments': YES,
                'trust_assigned_peaks': NO,
                'structural_rules_enabled' : NO,
                'filter_diagonal_peaks' : NO,
                'filter_unassigned_peaks' : NO,                
                'experiment_data' : experiment_data} # BARDIAUX rMat

class TemplateData(SimpleDataContainer):

    known_formats = ('iupac', 'cns', 'dyana', 'ccpn')

    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_TEMPLATE_STRUCTURE)

    def create(self):

        from aria.Settings import ChoiceEntity

        d = SimpleDataContainer.create(self)

        d['filename'].mandatory(0)

        return d

    def create_default_values(self):
        d = SimpleDataContainer.create_default_values(self)

        q = {'format': 'iupac', 'enabled': YES}
        d.update(q)

        return d

class InitialStructureData(TemplateData):
    def __init__(self):
        SimpleDataContainer.__init__(self, DATA_INITIAL_STRUCTURE)
        
##         self.getEntity('filename')
##         e_filename.mandatory(0)

##         e_enabled = self.getEntity
        
    def create(self):
        d = TemplateData.create(self)

        ## Initial structure need not exist
        d['filename'].mandatory(0)

        descr = 'If enabled, the user-defined structure is used as initial structure for the simulated annealing protocol. If disabled, ARIA automatically creates an extended initial structure.'

        d['filename'].setDescription(descr)
        d['enabled'].setDescription(descr)

        return d

class AmbiguousParameters(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_ANNEALING_AMBIG)

    def create(self):

        from aria.Settings import NonNegativeInt, NonNegativeFloat

        keywords = {'first_iteration': NonNegativeInt(),
                    'k_hot': NonNegativeFloat(),
                    'k_cool1_initial': NonNegativeFloat(),
                    'k_cool1_final': NonNegativeFloat(),
                    'k_cool2': NonNegativeFloat()}

        return keywords
                
    def create_default_values(self):

        default_values = {'first_iteration': 0,
                          'k_hot': 10.0,
                          'k_cool1_initial': 10.0,
                          'k_cool1_final': 50.0,
                          'k_cool2': 50.0}

        return default_values

class UnambiguousParameters(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_ANNEALING_UNAMBIG)

    def create(self):
        
        from aria.Settings import NonNegativeInt, NonNegativeFloat

        keywords = {'first_iteration': NonNegativeInt(),
                    'k_hot': NonNegativeFloat(),
                    'k_cool1_initial': NonNegativeFloat(),
                    'k_cool1_final': NonNegativeFloat(),
                    'k_cool2': NonNegativeFloat()}

        return keywords

    def create_default_values(self):

        default_values = {'first_iteration': 0,
                          'k_hot': 10.0,
                          'k_cool1_initial': 10.0,
                          'k_cool1_final': 50.0,
                          'k_cool2': 50.0}

        return default_values

class HBondParameters(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_ANNEALING_HBOND)

    def create(self):

        from aria.Settings import NonNegativeInt, NonNegativeFloat

        keywords = {'first_iteration': NonNegativeInt(),
                    'k_hot': NonNegativeFloat(),
                    'k_cool1_initial': NonNegativeFloat(),
                    'k_cool1_final': NonNegativeFloat(),
                    'k_cool2': NonNegativeFloat()}

        return keywords

    def create_default_values(self):

        default_values = {'first_iteration': 0,
                          'k_hot': 10.0,
                          'k_cool1_initial': 10.0,
                          'k_cool1_final': 50.0,
                          'k_cool2': 50.0}

        return default_values

class DihedralParameters(DataContainer):

    def __init__(self):

        DataContainer.__init__(self, DATA_ANNEALING_DIHEDRAL)

    def create(self):
        
        from aria.Settings import NonNegativeFloat

        keywords = {'k_hot': NonNegativeFloat(),
                    'k_cool1': NonNegativeFloat(),
                    'k_cool2': NonNegativeFloat()}

        return keywords
                
    def create_default_values(self):

        default_values = {'k_hot': 5.0,
                          'k_cool1': 25.0,
                          'k_cool2': 200.0}

        return default_values

class KarplusParameters(DataContainer):

    def __init__(self):

        DataContainer.__init__(self, DATA_ANNEALING_KARPLUS)

    def create(self):

        from aria.Settings import ChoiceEntity, Float, NonNegativeFloat

        choices = [1, 2, 3, 4, 5]
        msg = 'Choices for class of Karplus restraint: ' + \
              ', '.join(map(str, choices))
        
        keywords = {'class': ChoiceEntity(choices, error_message = msg),
                    'a': Float(),
                    'b': Float(),
                    'c': Float(),
                    'd': Float(),
                    'k_hot': NonNegativeFloat(),
                    'k_cool1': NonNegativeFloat(),
                    'k_cool2': NonNegativeFloat()}

        return keywords

    def create_default_values(self):

        default_values = {'class': 1, 
                          'a': 6.98,
                          'b': -1.38,
                          'c': 1.72,
                          'd': -60.0,
                          'k_hot': 0.0,
                          'k_cool1': 0.2,
                          'k_cool2': 1.0}

        return default_values

class RDCParameters(DataContainer):

    def __init__(self):

        DataContainer.__init__(self, DATA_ANNEALING_RDC)

    def create(self):

        from aria.Settings import ChoiceEntity, Float, NonNegativeFloat, \
             NonNegativeInt

        keywords = {}

        choices = (1, 2, 3, 4, 5)
        msg = 'Choices for the class of RDC restraints: ' + \
              ', '.join(map(str, choices))
        descr = """Residual dipolar coupling (RDC) measurements can be grouped into 5 classes. Every class has its own parameter set which the user has to specify. It contains the specification of the alignment tensor and settings for the simulated annealing protocol."""
        
        keywords['class'] = ChoiceEntity(choices, description = descr,
                                         error_message = msg)
        
        choices = ('SANI', 'VANGLE')
        msg = 'Possible choice for RDC restraint: ' + \
              ', '.join(map(str, choices))
        
        descr = """Aria offers two approaches to use residual dipolar coupling data as restraints: SANI and VEAN. For SANI, you have to specify the rhombicity and magnitude of the alignment tensor. VEAN uses angular restraints which must be precalculated with a separate program."""
        
        keywords['method'] = ChoiceEntity(choices, description = descr,
                                          error_message = msg)

        descr = """Rhombicity of the alignment tensor."""
        keywords['r'] = Float(description = descr)
        descr = """Magnitude of the alignment tensor."""
        keywords['d'] = Float(description = descr) 
        keywords['first_iteration'] = NonNegativeInt()
        keywords['k_hot'] = NonNegativeFloat()
        keywords['k_cool1'] = NonNegativeFloat()
        keywords['k_cool2'] = NonNegativeFloat()
        keywords['border_hot_initial'] = NonNegativeFloat()
        keywords['border_cool1_initial'] = NonNegativeFloat()
        keywords['border_cool2_initial'] = NonNegativeFloat()
        keywords['border_hot_final'] = NonNegativeFloat()
        keywords['border_cool1_final'] = NonNegativeFloat()
        keywords['border_cool2_final'] = NonNegativeFloat()
        keywords['center_hot_initial'] = NonNegativeFloat()
        keywords['center_cool1_initial'] = NonNegativeFloat()
        keywords['center_cool2_initial'] = NonNegativeFloat()
        keywords['center_hot_final'] = NonNegativeFloat()
        keywords['center_cool1_final'] = NonNegativeFloat()
        keywords['center_cool2_final'] = NonNegativeFloat()

        return keywords
    
    def create_default_values(self):

        default_values = {'method': 'SANI',
                          'class': 1, 
                          'first_iteration': 0,
                          'k_cool1': 0.2,
                          'k_cool2': 1.0,
                          'k_hot': 0.0,
                          'd': 8.0,
                          'r': 0.4,
                          'border_hot_initial': 0.1,
                          'border_cool1_initial': 40.0,
                          'border_cool2_initial': 40.0,
                          'border_hot_final': 40.0,
                          'border_cool1_final': 40.0,
                          'border_cool2_final': 40.0,
                          'center_hot_initial': 0.1,
                          'center_cool1_initial': 10.0,
                          'center_cool2_initial': 10.0,
                          'center_hot_final': 0.1,
                          'center_cool1_final': 10.0,
                          'center_cool2_final': 10.0}

        return default_values

## BARDIAUX 2.2 SymmetryParameters
class SymmetryParameters(DataContainer):

    def __init__(self):

        DataContainer.__init__(self, DATA_ANNEALING_SYM)

    def create(self):
        
        from aria.Settings import NonNegativeFloat, YesNoChoice, NonNegativeInt

        ## TODO: check range and types

        keywords = {'k_packing_hot': NonNegativeFloat(),
                    'k_packing_cool1': NonNegativeFloat(),
                    'k_packing_cool2': NonNegativeFloat(),
                    'last_iteration_packing' : NonNegativeInt(),
                    'k_ncs': NonNegativeFloat()}

        return keywords

    def create_default_values(self):

        defaults = {'k_packing_hot': 15.,
                    'k_packing_cool1': 10.,
                    'k_packing_cool2': 5.,
                    'last_iteration_packing' : 8,
                    'k_ncs': 50.}

        return defaults


class FBHWParameters(DataContainer):

    def __init__(self):

        DataContainer.__init__(self, DATA_ANNEALING_FBHW)

    def create(self):
        
        from aria.Settings import NonNegativeFloat, Float

        ## TODO: check range and types

        keywords = {'m_rswitch_hot': NonNegativeFloat(),
                    'm_rswitch_cool1': NonNegativeFloat(),
                    'm_rswitch_cool2': NonNegativeFloat(),
                    'rswitch_hot': NonNegativeFloat(),
                    'rswitch_cool1': NonNegativeFloat(),
                    'rswitch_cool2': NonNegativeFloat(),
                    'm_asymptote_hot': Float(),
                    'm_asymptote_cool1': Float(),
                    'm_asymptote_cool2': Float(),
                    'asymptote_hot': NonNegativeFloat(),
                    'asymptote_cool1': NonNegativeFloat(),
                    'asymptote_cool2': NonNegativeFloat()}

        return keywords

    def create_default_values(self):

        defaults = {'m_rswitch_hot': 0.5,
                    'm_rswitch_cool1': 0.5, 
                    'm_rswitch_cool2': 0.5, 
                    'rswitch_hot': 0.5, 
                    'rswitch_cool1': 0.5, 
                    'rswitch_cool2': 0.5, 
                    'm_asymptote_hot': -1.0, 
                    'm_asymptote_cool1': -1.0, 
                    'm_asymptote_cool2': -0.1, 
                    'asymptote_hot': 1.0, 
                    'asymptote_cool1': 1.0, 
                    'asymptote_cool2': 0.1}

        return defaults
    
# BERNARD Aymeric 2.3 log-harmonic
class LogHarmonicParameters(DataContainer):

    def __init__(self):

        DataContainer.__init__(self, DATA_ANNEALING_LOGHARMONIC)

    def create(self):

        from aria.Settings import YesNoChoice, PositiveFloat

        ## TODO: check range and types

        keywords = {}

        #choices = ['soft' , 'lognormal' ]
        #help = ' Choice between '+' '.join(choices)
        #msg =  'error: '+help
        #keywords['potential_type'] = ChoiceEntity(choices, description = help, error_message = msg) #!+!lien

        descr = 'Use log-harmonic potential for SA second cooling stage ' + \
              ' / '.join((YES,NO))
        msg = "Choice between : " + ' / '.join((YES,NO))
        keywords['enabled'] = YesNoChoice(description = descr, error_message = msg)

        descr = 'Use automatic distance restraints weight for the log-harmonic potential ' + \
              ' / '.join((YES,NO))
        msg = "Choice between : " + ' / '.join((YES,NO))
        keywords['use_auto_weight'] = YesNoChoice(description = descr, error_message = msg)

        descr = 'Weight for unambiguous distance restraints with the log-harmonic potential\n  (if automatic weight is not enabled)'
        keywords["weight_unambig"] = PositiveFloat(description = descr)

        descr = 'Weight for ambiguous distance restraints with the log-harmonic potential\n  (if automatic weight is not enabled)'
        keywords["weight_ambig"] = PositiveFloat(description = descr)

        descr = 'Weight for hbond restraints with the log-harmonic potential\n  (if automatic weight is not enabled)'
        keywords["weight_hbond"] = PositiveFloat(description = descr)

        return keywords

    def create_default_values(self):

        defaults = {'enabled' : NO,
                    'use_auto_weight' : NO,
                    'weight_unambig' : 25.0,
                    'weight_ambig' : 25.0,
                    'weight_hbond' : 25.0 }

        return defaults

class MDParameters(DataContainer):

    def __init__(self):

        DataContainer.__init__(self, DATA_DYNAMICS)

    def create(self):
        
        ## TODO: check types and ranges, filesort?

        from aria.Settings import ChoiceEntity, PositiveInteger, NonNegativeInt, \
             PositiveFloat, NonNegativeFloat

        keywords = {'random_seed': PositiveInteger(),
                    'tad_temp_high': PositiveFloat(),
                    'cartesian_temp_high': PositiveFloat(),
                    'temp_cool1_final': NonNegativeFloat(),
                    'temp_cool2_final': NonNegativeFloat(),
                    'timestep': PositiveFloat(),
                    'tad_timestep_factor': PositiveFloat(),
                    'steps_high': PositiveInteger(),
                    'steps_refine': PositiveInteger(),
                    'steps_cool1': PositiveInteger(),
                    'steps_cool2': PositiveInteger(),
                    'cartesian_first_iteration': NonNegativeInt()}

        choices = ['cartesian', 'torsion']
        msg = 'MD-type for SA-protocol: ' + ' / '.join(choices)
        keywords['md_type'] = ChoiceEntity(choices, error_message = msg)

        return keywords

    def create_default_values(self):

        defaults = {'random_seed': 89764443, 
                    'tad_temp_high': 10000.,
                    'cartesian_temp_high': 2000.,
                    'temp_cool1_final': 1000.,
                    'temp_cool2_final': 50.,
                    'timestep': 0.003,
                    'tad_timestep_factor': 9.,
                    'steps_high': 10000,
                    'steps_refine': 4000,
                    'steps_cool1': 5000,
                    'steps_cool2': 4000,
                    'cartesian_first_iteration': 0, 
                    'md_type': 'torsion'}

        return defaults

class WaterRefinementParameters(DataContainer):

    def create(self):
        
        from aria.Settings import NonNegativeInt, ChoiceEntity
        from aria.Settings import YesNoChoice

        d = DataContainer.create(self)

        keywords = {'n_structures': NonNegativeInt()}

        choices = ['water', 'dmso']
        msg = 'Possible solvent for water refinement: ' + \
              ', '.join(choices)
        keywords['solvent'] = ChoiceEntity(choices, error_message = msg)

        choices = [YES, NO]
        msg = 'Water refinement for the last iteration? ' + \
              ' / '.join(choices)
        keywords['enabled'] = ChoiceEntity(choices, error_message = msg)
        
        descr = 'If enabled, PDB-files include solvent molecules.'
        keywords['write_solvent_molecules'] = YesNoChoice(description = descr)

        d.update(keywords)

        return d
    
    def create_default_values(self):

        defaults = {'n_structures': 10,
                    'solvent': 'water',
                    'enabled': YES,
                    'write_solvent_molecules': NO}

        return defaults

class AnnealingParameters(DataContainer):

    def __init__(self):
        DataContainer.__init__(self, DATA_ANNEALING)

        self[DATA_ANNEALING_KARPLUS] = {}
        self[DATA_ANNEALING_RDC] = {}
        
    def create(self):

        from aria.Settings import TypeEntity

        keywords = {DATA_ANNEALING_AMBIG:
                    TypeEntity('AmbiguousParameters'),

                    DATA_ANNEALING_UNAMBIG:
                    TypeEntity('UnambiguousParameters'),

                    DATA_ANNEALING_HBOND:
                    TypeEntity('HBondParameters'),

                    DATA_ANNEALING_DIHEDRAL:
                    TypeEntity('DihedralParameters'),
                    
                    DATA_ANNEALING_KARPLUS:
                    TypeEntity(DICT),

                    DATA_ANNEALING_RDC:
                    TypeEntity(DICT),
                    
                    DATA_ANNEALING_FBHW:
                    TypeEntity('FBHWParameters'),

                    DATA_ANNEALING_SYM:
                    TypeEntity('SymmetryParameters'),        # BARDIAUX 2.2
        
                    DATA_ANNEALING_LOGHARMONIC:
                    TypeEntity('LogHarmonicParameters')}     # BERNARD 2.3
        
        return keywords

    def addParameters(self, p):

        t = p.getType()

        if t in (DATA_ANNEALING_KARPLUS, DATA_ANNEALING_RDC):
            c = p['class']
            if self[t].has_key(c):
                self.warning('Parameters for restraint type "%s" ' % t + \
                             'must be unique for each class; values ' + \
                             'already specified for class "%d" ' % c + \
                             'will be kept.')
            else:
                self[t][c] = p
        else:
            self[t] = p

    def getParameters(self, type):

        val = self[type]

        if type in (DATA_ANNEALING_KARPLUS, DATA_ANNEALING_RDC):
            if not len(val):
                val = ()
            else:
                val = tuple(val.values())
        
        return val

# BARDIAUX 2.2
# previously in importFromCcpn
class CCPNData(Settings):

  def create(self):

    from aria.Settings import MultiTypeEntity, AbsolutePath, TypeEntity
    from aria.TypeChecking import TUPLE

    d = {}
    
    d['filename'] = AbsolutePath(exists=0)

    return d

  def create_default_values(self):
    
    d = {'filename': ''}
    
    return d

class CCPNDataXMLPickler(XMLBasePickler):

  order = 'filename',

  def create(self):
    return CCPNData()

  def _xml_state(self, x):

    from aria.xmlutils import XMLElement

    e = XMLElement()

    e.filename = x['filename']

    order = list(self.order)

    e.set_tag_order(order)

    return e

  def load_from_element(self, e):

    from aria.tools import as_tuple
    
    s = self.create()

    filename = str(e.filename).strip()

    E = s.getEntity('filename')

    E.reset()
    E.mandatory(filename <> '')

    if self.relaxed:
        E.mandatory(0)

    s['filename'] = filename

    if self.relaxed:
        E.mandatory(filename <> '')

    return s


class AnnealingParametersXMLPickler(XMLBasePickler):

    def _xml_state(self, x):

        order = ('unambiguous_restraints', 'ambiguous_restraints',
                 'hbond_restraints', 'dihedral_restraints',
                 'karplus_restraints', 'rdc_restraints',
                 'flat_bottom_harmonic_wall', 'symmetry_restraints',
                 'logharmonic_potential')

        e = XMLElement(tag_order = order)

        f = x.getParameters

        e.unambiguous_restraints = f(DATA_ANNEALING_UNAMBIG)
        e.ambiguous_restraints = f(DATA_ANNEALING_AMBIG)
        e.hbond_restraints = f(DATA_ANNEALING_HBOND)
        e.dihedral_restraints = f(DATA_ANNEALING_DIHEDRAL)
        e.flat_bottom_harmonic_wall = f(DATA_ANNEALING_FBHW)

        e.karplus_restraints = f(DATA_ANNEALING_KARPLUS)
        e.rdc_restraints = f(DATA_ANNEALING_RDC)
        # BARDIAUX 2.2
        e.symmetry_restraints = f(DATA_ANNEALING_SYM)
        # BERNARD 2.3
        e.logharmonic_potential = f(DATA_ANNEALING_LOGHARMONIC)
        
        return e

    def load_from_element(self, e):

        from aria.tools import as_tuple

        s = AnnealingParameters()

        f = lambda p, g = s.addParameters: g(p)

        [f(p) for p in as_tuple(e.unambiguous_restraints)]
        [f(p) for p in as_tuple(e.ambiguous_restraints)]
        [f(p) for p in as_tuple(e.hbond_restraints)]
        [f(p) for p in as_tuple(e.dihedral_restraints)]
        [f(p) for p in as_tuple(e.flat_bottom_harmonic_wall)]
        [f(p) for p in as_tuple(e.karplus_restraints)]
        [f(p) for p in as_tuple(e.rdc_restraints)]
        
        # BARDIAUX 2.2
        if hasattr(e, 'symmetry_restraints'):
            [f(p) for p in as_tuple(e.symmetry_restraints)]

        else:
            z = SymmetryParameters()
            z.reset()

            [f(p) for p in as_tuple(z)]
            
        # BERNARD 2.3
        if hasattr(e, 'logharmonic_potential'):
            [f(p) for p in as_tuple(e.logharmonic_potential)]
        else:
            z = LogHarmonicParameters()
            z.reset()
            [f(p) for p in as_tuple(z)]

        return s
   
class WaterRefinementXMLPickler(XMLBasePickler):

    order = ['solvent', 'n_structures', 'enabled',
             'write_solvent_molecules']

    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)

        e.solvent = x['solvent']
        e.n_structures = x['n_structures']
        e.enabled = x['enabled']
        e.write_solvent_molecules = x['write_solvent_molecules']

        return e

    def load_from_element(self, e):
        s = WaterRefinementParameters()

        s['solvent'] = str(e.solvent)
        s['n_structures'] = int(e.n_structures)
        s['enabled'] = str(e.enabled)
        s['write_solvent_molecules'] = str(e.write_solvent_molecules)

        return s

class DataContainerXMLPickler(XMLBasePickler):
    
    def create(self):
        return DataContainer()

    def _xml_state(self, x):
        return XMLElement()

    def load_from_element(self, e):
        return self.create()

class FileFormatXMLPickler(DataContainerXMLPickler):

    order = ('file', 'format', 'ccpn_id')
    
    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)

        e.file = x['filename']
        e.format = x['format']

        if 'ccpn_id' in x:
            e.ccpn_id = x['ccpn_id']
        else:
            e.ccpn_id = ''

        return e

    def create(self):
        return SimpleDataContainer()

    def load_from_element(self, e):
        s = self.create()

        ## relaxed loading?
        ## if so, do not check whether path exists

        if self.relaxed:
            entity = s.getEntity('filename')
            is_mandatory = entity.is_mandatory()
            entity.mandatory(0)

        s['filename'] = str(e.file)
        
        if self.relaxed:
            entity.mandatory(is_mandatory)
            
        s['format'] = str(e.format)

        if hasattr(e, 'ccpn_id'):
            s['ccpn_id'] = str(e.ccpn_id)
        else:
            s['ccpn_id'] = ''

        return s

class SimpleDCXMLPickler(FileFormatXMLPickler):

    order = FileFormatXMLPickler.order + ('enabled',)

    def _xml_state(self, x):
        
        e = FileFormatXMLPickler._xml_state(self, x)
        
        e.enabled = x['enabled']

        enabled = hasattr(e, 'enabled') and x['enabled'] == YES
            
        if enabled and e.format == 'ccpn' and e.ccpn_id == '':
            self.error(ValueError, 'Using the CCPN data model for data retrieval requires a valid CCPN id. Current id is "%s". Please check your project file.' % str(e.ccpn_id))

        if enabled and e.format <> 'ccpn' and e.file == '':
            self.error(ValueError, 'Format %s expects a filename.' % e.format)

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        s = FileFormatXMLPickler.load_from_element(self, e)

        s['enabled'] = str(e.enabled)

        ## consistency check.

        enabled = 'enabled' in s and s['enabled'] == YES

        if enabled and s['format'] == 'ccpn' and s['ccpn_id'] == '':
            
            self.error(ValueError, 'Using the CCPN data model for data retrieval requires a valid CCPN id.')

        if enabled and s['format'] <> 'ccpn' and s['filename'] == '':
            self.error(ValueError, 'Format "xml" expects a filename.')

        return s

class BoundCorrectionXMLPickler(XMLBasePickler):

    order = ['value', 'enabled']

    def __init__(self, _type = None):
        if not _type in (LowerBoundCorrection,
                         UpperBoundCorrection, None):
            s = 'Class LowerBoundCorrection or UpperBoundCorrection' + \
                'expected. %s given.'
            self.error(ValueError, s % str(_type))

        if _type is None:
            _type = BoundCorrection

        self.__type = _type

    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)
        e.value = x['value']
        e.enabled = x['enabled']

        return e

    def create(self):
        return self.__type()

    def load_from_element(self, e):

        s = self.create()
        s['value'] = float(e.value)
        s['enabled'] = str(e.enabled)

        return s

FFP = FileFormatXMLPickler

class SequenceDataXMLPickler(FFP):

    sub_order = ('name', 'filename')
    
    order = ('linkage_definition', 'parameter_definition',
             'topology_definition')

    def _xml_state_linkage(self, x):
        
        e = XMLElement(tag_order = self.sub_order)
        
        e.name = x['linkage_name']

        if e.name == 'user_defined':
            filename = x['linkage_filename']
        else:
            filename = ''

        e.filename = filename

        return e

    def _xml_state_topology(self, x):
        
        e = XMLElement(tag_order = self.sub_order)
        
        e.name = x['topology_name']

        if e.name == 'user_defined':
            filename = x['topology_filename']
        else:
            filename = ''

        e.filename = filename

        return e

    def _xml_state_parameter(self, x):
        
        e = XMLElement(tag_order = self.sub_order)
        
        e.name = x['parameter_name']

        if e.name == 'user_defined':
            filename = x['parameter_filename']
        else:
            filename = ''

        e.filename = filename

        return e

    def _xml_state(self, x):
        e = FFP._xml_state(self, x)
        e.set_tag_order(FFP.order + self.order)

        e.linkage_definition = self._xml_state_linkage(x)
        e.topology_definition = self._xml_state_topology(x)
        e.parameter_definition = self._xml_state_parameter(x)

        return e

    def load_from_element(self, e):
        x = FFP.load_from_element(self, e)

##         if x['format'] == 'ccpn' and x['ccpn_id'].count('|') <> 1:
##             self.error(ValueError, 'CCPN sequence identifier must have the following format: "molsystem_name|chain_code", e.g. "MS1|A" Current id is "%s". Please check for project file.' % x['ccpn_id'])

        if x['format'] == 'ccpn' and x['ccpn_id'].count('|') == 0:
            self.error(ValueError, 'CCPN sequence identifier must have the following format: "molsystem_name|chain_code1", e.g. "MS1|A" Current id is "%s". In case of of home-dimer, the format must be: "molsystem_name|first_chain_code|second_chain_code. Please check for project file.' % x['ccpn_id'])

        linkage = e.linkage_definition
        name = str(linkage.name)
        x['linkage_name'] = name

        if name == 'user_defined':
            x['linkage_filename'] = str(linkage.filename)

        topology = e.topology_definition
        name = str(topology.name)
        x['topology_name'] = name

        if name == 'user_defined':
             x['topology_filename'] = str(topology.filename)

        parameter = e.parameter_definition
        name = str(parameter.name)
        x['parameter_name'] = name

        if name == 'user_defined':
            x['parameter_filename'] = str(parameter.filename)

        return x
    
    def create(self):
        return SequenceData()

class ShiftDataXMLPickler(FFP):

    order = FFP.order + ('default_shift_error',)

    def _xml_state(self, x):

        e = FFP._xml_state(self, x)
        e.default_shift_error = x['default_shift_error']

        return e

    def load_from_element(self, e):

        x = FFP.load_from_element(self, e)
        x['default_shift_error'] = float(e.default_shift_error)

        return x
    
    def create(self):
        return ShiftData()

class PeakDataXMLPickler(FFP):

    order = FFP.order + ('peak_size',
                         'freq_window_proton1', 'freq_window_hetero1',
                         'freq_window_proton2', 'freq_window_hetero2',
                         'lower_bound_correction',
                         'upper_bound_correction')

    def create(self):
        return PeakData()

    def _xml_state(self, x):
        
        e = FFP._xml_state(self, x)
        e.freq_window_proton1 = x['proton1_shift_err']
        e.freq_window_hetero1 = x['hetero1_shift_err']
        e.freq_window_proton2 = x['proton2_shift_err']
        e.freq_window_hetero2 = x['hetero2_shift_err']
        e.peak_size = x['volume_or_intensity']
        e.lower_bound_correction = x['lower_bound_correction']
        e.upper_bound_correction = x['upper_bound_correction']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):
        
        s = FFP.load_from_element(self, e)
        
        if s['format'] == 'ccpn' and s['ccpn_id'].count('|') <> 3:
            
            self.error(ValueError, 'CCPN peak list identifier must have the following format: "experiment_name|data_source_name|serial_number", e.g. "noesy|xeasy|1". Current id is "%s". Please check you project file.' % s['ccpn_id'])
        
        s['proton1_shift_err'] = float(e.freq_window_proton1)
        s['hetero1_shift_err'] = float(e.freq_window_hetero1)
        s['proton2_shift_err'] = float(e.freq_window_proton2)
        s['hetero2_shift_err'] = float(e.freq_window_hetero2)
        s['volume_or_intensity'] = str(e.peak_size)
        s['lower_bound_correction'] = e.lower_bound_correction
        s['upper_bound_correction'] = e.upper_bound_correction

        return s

BASE = SimpleDCXMLPickler

class TemplateDataXMLPickler(BASE):
    def create(self):
        return TemplateData()

class InitialStructureDataXMLPickler(BASE):
    def create(self):
        return InitialStructureData()
    
## BARDIAUX 2.2
## class AmbiguousDistanceDataXMLPickler(BASE):

##     def create(self):
##         return AmbiguousDistanceData()

class DistanceDataXMLPickler(BASE):

    order = BASE.order + ('add_to_network','calibrate', 'run_network_anchoring', 'filter_contributions',)
    
    def create(self):
        return AmbiguousDistanceData()
    
    def _xml_state(self, x):
        
        e = BASE._xml_state(self, x)
        e.add_to_network = x['add_to_network']
        e.calibrate = x['calibrate']
        e.run_network_anchoring = x['run_network_anchoring']
        e.filter_contributions = x['filter_contributions']
        
        e.set_tag_order(self.order)
            
        if e.enabled == YES and e.format <> 'ccpn' :
            warn = 'you must use the CCPN format.'
            
            if  e.add_to_network == YES:
                self.error(ValueError, 'To add distance restraints for the Network-Anchoring, %s' % (s['filename'], warn))
            if  e.calibrate <> NO:
                self.error(ValueError, 'To calibrate distance restraints, %s' % (s['filename'], warn))
            if  e.run_network_anchoring <> NO:
                self.error(ValueError, 'To run the network-anchoring on restraints, %s' % (s['filename'], warn))                
            if  e.filter_contributions == YES:
                self.error(ValueError, 'To filter contributions of distance restraints, %s' % (s['filename'], warn))            
        return e
    
    def load_from_element(self, e):
        s = BASE.load_from_element(self, e)
        if hasattr(e, 'add_to_network'):
            s['add_to_network'] = str(e.add_to_network)
        else:
            s['add_to_network'] = NO
            
        if hasattr(e, 'calibrate'):
            s['calibrate'] = str(e.calibrate)
        else:
            s['calibrate'] = NO
            
        if hasattr(e, 'run_network_anchoring'):
            s['run_network_anchoring'] = str(e.run_network_anchoring)
        else:
            s['run_network_anchoring'] = NO
            
        if hasattr(e, 'filter_contributions'):
            s['filter_contributions'] = str(e.filter_contributions)
        else:
            s['filter_contributions'] = NO
            
        if s['enabled'] == YES and s['format'] <> 'ccpn' :
            warn = 'you must use the CCPN format. Please check for project file.'
            
            if  s['add_to_network'] == YES:
                self.error(ValueError, 'To add distance restraints for the Network-Anchoring, %s' % warn)
            if  s['calibrate'] <> NO:
                self.error(ValueError, 'To calibrate distance restraints, %s' % warn)
            if  s['run_network_anchoring'] == YES:
                self.error(ValueError, 'To run the network-anchoring on distance restraints, %s' % warn)       
            if  s['filter_contributions'] == YES:
                self.error(ValueError, 'To filter contributions of distance restraints, %s' % warn)
                
        return s
    
class AmbiguousDistanceDataXMLPickler(DistanceDataXMLPickler):

##    order = BASE.order + ('add_to_network','calibrate', 'filter_contributions',)
    
    def create(self):
        return AmbiguousDistanceData()

    def _xml_state(self, x):
        e = DistanceDataXMLPickler._xml_state(self, x)
        return e

    def load_from_element(self, e):
        s = DistanceDataXMLPickler.load_from_element(self, e)
        return s
        
##         e = BASE._xml_state(self, x)    
##     def _xml_state(self, x):
        
##         e = BASE._xml_state(self, x)
##         e.add_to_network = x['add_to_network']
##         e.calibrate = x['calibrate']
##         e.filter_contributions = x['filter_contributions']
        
##         e.set_tag_order(self.order)
            
##         if e.enabled == YES and e.format <> 'ccpn' :
##             warn = 'you must use the CCPN format.'
            
##             if  e.add_to_network == YES:
##                 self.error(ValueError, 'To add distance restraints for the Network-Anchoring, %s' % (s['filename'], warn))
##             if  e.calibrate <> NO:
##                 self.error(ValueError, 'To calibrate distance restraints, %s' % (s['filename'], warn))            
##             if  e.filter_contributions == YES:
##                 self.error(ValueError, 'To filter contributiosn of distance restraints, %s' % (s['filename'], warn))            
##         return e
    
##     def load_from_element(self, e):
##         s = BASE.load_from_element(self, e)
##         if hasattr(e, 'add_to_network'):
##             s['add_to_network'] = str(e.add_to_network)
##         else:
##             s['add_to_network'] = NO
            
##         if hasattr(e, 'calibrate'):
##             s['calibrate'] = str(e.calibrate)
##         else:
##             s['calibrate'] = NO

##         if hasattr(e, 'filter_contributions'):
##             s['filter_contributions'] = str(e.filter_contributions)
##         else:
##             s['filter_contributions'] = NO
            
##         if s['enabled'] == YES and s['format'] <> 'ccpn' :
##             warn = 'you must use the CCPN format. Please check for project file.'
            
##             if  s['add_to_network'] == YES:
##                 self.error(ValueError, 'To add distance restraints for the Network-Anchoring, %s' % (s['filename'], warn))
##             if  s['calibrate'] <> NO:
##                 self.error(ValueError, 'To calibrate distance restraints, %s' % (s['filename'], warn))            
##             if  s['filter_contributions'] == YES:
##                 self.error(ValueError, 'To filter contributiosn of distance restraints, %s' % (s['filename'], warn))
                
##         return s
    

## class UnambiguousDistanceDataXMLPickler(BASE):

##     def create(self):
##         return UnambiguousDistanceData()
    
## BARDIAUX 2.2
class UnambiguousDistanceDataXMLPickler(DistanceDataXMLPickler):

##    order = BASE.order + ('add_to_network',)
    
    def create(self):
        return UnambiguousDistanceData()
    
    def _xml_state(self, x):
        e = DistanceDataXMLPickler._xml_state(self, x)
        return e

    def load_from_element(self, e):
        s = DistanceDataXMLPickler.load_from_element(self, e)
        return s
    
##     def _xml_state(self, x):
        
##         e = BASE._xml_state(self, x)
##         e.add_to_network = x['add_to_network']
##         e.calibrate = x['calibrate']
##         e.filter_contributions = x['filter_contributions']
        
##         e.set_tag_order(self.order)
        
##         if e.enabled == YES and e.add_to_network ==YES and e.format <> 'ccpn':
##             self.error(ValueError, 'To add distance restraints to the Network-Anchoring, you must use the CCPN format.')
            
                       
##         return e
    
##     def load_from_element(self, e):
##         s = BASE.load_from_element(self, e)
##         if hasattr(e, 'add_to_network'):
##             s['add_to_network'] = str(e.add_to_network)
##         else:
##             s['add_to_network'] = NO

##         if s['enabled'] == YES and s['add_to_network'] == YES and s['format'] <> 'ccpn' :
##             self.error(ValueError, 'To add distance restraints to the Network-Anchoring, you must use the CCPN format. Please check for project file.' % s['filename'])
            
##         return s
   
class HBondDataXMLPickler(BASE):

    order = BASE.order + ('data_type',)
    
    def create(self):
        return HBondData()

    def _xml_state(self, x):
        
        e = BASE._xml_state(self, x)
        e.data_type = x['type']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):
        s = BASE.load_from_element(self, e)
        s['type'] = str(e.data_type)

        return s

class DihedralDataXMLPickler(HBondDataXMLPickler):
    def create(self):
        return DihedralData()

class OtherDataXMLPickler(HBondDataXMLPickler):
    def create(self):
        return OtherData()

class KarplusDataXMLPickler(BASE):

    order = BASE.order + ('parameter_class',)
    
    def create(self):
        return KarplusData()

    def _xml_state(self, x):
        
        e = BASE._xml_state(self, x)
        e.parameter_class = x['class']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):
        s = BASE.load_from_element(self, e)
        s['class'] = int(e.parameter_class)

        return s

class RDCDataXMLPickler(KarplusDataXMLPickler):
    def create(self):
        return RDCData()
    
## BARDIAUX 2.2
class CysPatchXMLPickler(DataContainerXMLPickler):

    order = ['residue', 'segid']

    def create(self):
        return CysPatch()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        e.residue = x['residue']
        e.segid = x['segid']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)
        x['residue'] = int(e.residue)
        x['segid'] = str(e.segid)

        return x

## BARDIAUX 2.2
class ZnPatchXMLPickler(DataContainerXMLPickler):

    order = ['type',
             'residue_zn', 'segid_zn',
             'residue1', 'segid1',
             'residue2', 'segid2',
             'residue3', 'segid3',
             'residue4', 'segid4',]

    def create(self):
        return ZnPatch()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        
        e.type = x['type']
        
        e.residue_zn = x['residue_zn']
        e.segid_zn = x['segid_zn']
        
        e.residue1 = x['residue1']
        e.segid1 = x['segid1']
        
        e.residue2 = x['residue2']
        e.segid2 = x['segid2']
        
        e.residue3 = x['residue3']
        e.segid3 = x['segid3']
        
        e.residue4 = x['residue4']
        e.segid4 = x['segid4']
        
        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)
        x['type'] = str(e.type)

        x['residue_zn'] = int(e.residue_zn)
        x['segid_zn'] = str(e.segid_zn)
        
        x['residue1'] = int(e.residue1)
        x['segid1'] = str(e.segid1)

        x['residue2'] = int(e.residue2)
        x['segid2'] = str(e.segid2)
        
        x['residue3'] = int(e.residue3)
        x['segid3'] = str(e.segid3)

        x['residue4'] = int(e.residue4)
        x['segid4'] = str(e.segid4)
        
        return x
    
## BARDIAUX 2.2 Disn
class SSBondDataXMLPickler(BASE):

    order = BASE.order + ('cyspatch',)
    
    def create(self):
        return  SSBondData()

    def _xml_state(self, x):
        
        e = BASE._xml_state(self, x)

        e.cyspatch = list(x['cyspatch'])

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):
        s = BASE.load_from_element(self, e)
        from aria.tools import as_tuple

        if hasattr(e, 'cyspatch'):
            [s['cyspatch'].append(x) for x in as_tuple(e.cyspatch)]
        else:
            s['cyspatch'] = []

        return s
## BARDIAUX rMat
class ExperimentDataXMLPickler(DataContainerXMLPickler):

    order = ('molecule_correlation_time', 'spectrum_mixing_time',
             'spectrometer_frequency', 'ambiguity_type')
    
    def create(self):
        return ExperimentData()

    def _xml_state(self, x):
        
        e = DataContainerXMLPickler._xml_state(self, x)
        e.molecule_correlation_time = x['molecule_correlation_time']
        e.spectrum_mixing_time = x['spectrum_mixing_time']
        e.spectrometer_frequency = x['spectrometer_frequency']     

        e.ambiguity_type = x['ambiguity_type']
        
        e.set_tag_order(self.order)

        return e
    
    def load_from_element(self, e):
        
        s = DataContainerXMLPickler.load_from_element(self, e)

        if hasattr(e,'ambiguity_type'):
            s['ambiguity_type'] = str(e.ambiguity_type)
        else:
            s['ambiguity_type'] = 'intra'

            
        for d in self.order[:-1]:

            v = getattr(e, d)
            if str(v) <> "":                
                s[d] = float(v)
            #else:
            #    s[d] = ''

        return s
    
class SpectrumDataXMLPickler(DataContainerXMLPickler):

    order = ('enabled', 'use_assignments', 'shifts', 'peaks',
             'trust_assigned_peaks', 'structural_rules', 'filter_diagonal_peaks', 'filter_unassigned_peaks', 'experiment_data') ## BARDIAUX 2.2, rMat
    
    def create(self):
        return SpectrumData()

    def _xml_state(self, x):
        
        e = DataContainerXMLPickler._xml_state(self, x)
        e.shifts = x['shifts']
        e.peaks = x['peaks']
        e.enabled = x['enabled']
        e.use_assignments = x['use_assignments']
        e.trust_assigned_peaks = x['trust_assigned_peaks']

        ## BARDIAUX 2.2
        e.structural_rules = x['structural_rules_enabled'] 
        e.filter_diagonal_peaks = x['filter_diagonal_peaks']
        e.filter_unassigned_peaks = x['filter_unassigned_peaks']
        
        ## BARDIAUX rMat
        e.experiment_data = x['experiment_data']
        
        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):
        s = DataContainerXMLPickler.load_from_element(self, e)
        
        s['shifts'] = e.shifts
        s['peaks'] = e.peaks

        ## consistency checks

        if s['shifts']['format'] <> s['peaks']['format']:
            self.error(ValueError, 'Peak list and chemical shift list must have the same format.')
            
        s['enabled'] = str(e.enabled)
        s['use_assignments'] = str(e.use_assignments)
        s['trust_assigned_peaks'] = str(e.trust_assigned_peaks)
        
        ## BARDIAUX 2.2
        if hasattr(e,'structural_rules'):
            s['structural_rules_enabled'] = str(e.structural_rules)
        else:
            s['structural_rules_enabled'] = NO

        ## BARDIAUX 2.2
        if hasattr(e,'filter_diagonal_peaks'):
            s['filter_diagonal_peaks'] = e.filter_diagonal_peaks
        else:
            s['filter_diagonal_peaks'] = NO

        ## BARDIAUX 2.2
        if hasattr(e,'filter_unassigned_peaks'):
            s['filter_unassigned_peaks'] = e.filter_unassigned_peaks
        else:
            s['filter_unassigned_peaks'] = NO
            
        ## BARDIAUX 2.2
        if hasattr(e,'experiment_data'):
            s['experiment_data'] = e.experiment_data
        else:
            exp = ExperimentData()
            exp.reset()
            s['experiment_data'] = exp

        

        return s

class SSBridgeXMLPickler(DataContainerXMLPickler):

    order = ['residue1', 'segid1', 'residue2', 'segid2']

    def create(self):
        return SSBridge()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        e.residue1 = x['residue1']
        e.segid1 = x['segid1']
        e.residue2 = x['residue2']
        e.segid2 = x['segid2']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)
        x['residue1'] = int(e.residue1)
        x['segid1'] = str(e.segid1)
        x['residue2'] = int(e.residue2)
        x['segid2'] = str(e.segid2)

        return x
    
class HisPatchXMLPickler(DataContainerXMLPickler):

    order = ['residue', 'segid', 'proton']

    def create(self):
        return HisPatch()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        e.residue = x['residue']
        e.segid = x['segid']
        e.proton = x['proton']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)
        x['residue'] = int(e.residue)
        x['segid'] = str(e.segid)
        x['proton'] = str(e.proton)

        return x
## BARDIAUX 2.2
class CisProPatchXMLPickler(DataContainerXMLPickler):

    order = ['residue', 'segid']

    def create(self):
        return CisProPatch()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        e.residue = x['residue']
        e.segid = x['segid']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)
        x['residue'] = int(e.residue)
        x['segid'] = str(e.segid)

        return x


## BARDIAUX : symmetry 2.2
class SymmetryXMLPickler(DataContainerXMLPickler):

    order = ['enabled', 'method', 'n_monomers','symmetry_type', 'ncs_enabled', 'packing_enabled']

    def create(self):
        return Symmetry()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        e.enabled = x['enabled']
        e.method = x['method']
        e.n_monomers = x['n_monomers']        
        e.symmetry_type = x['symmetry_type']
        e.ncs_enabled = x['ncs_enabled']
        e.packing_enabled = x['packing_enabled']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)
        x['enabled'] = str(e.enabled)
        x['n_monomers'] = int(e.n_monomers)
        x['symmetry_type'] = str(e.symmetry_type)
        x['ncs_enabled'] = str(e.ncs_enabled)
        x['packing_enabled'] = str(e.packing_enabled)

        if hasattr(e, 'method'):
            x['method'] = str(e.method)
        else:
            x['method'] = 'standard'

        return x

# BERNARD Aymeric 2.3 log-harmonic
class LogHarmonicParametersXMLPickler(DataContainerXMLPickler):

    order = ['enabled','use_auto_weight','weight_unambig','weight_ambig','weight_hbond']

    def create(self):
        return LogHarmonicParameters()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        [setattr(e, k, x[k]) for k in x.keys()]

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)

        x['enabled']         = str(e.enabled)
        x['use_auto_weight']  = str(e.use_auto_weight)
        x['weight_unambig']  = float(e.weight_unambig)
        x['weight_ambig']    = float(e.weight_ambig)
        x['weight_hbond']    = float(e.weight_hbond)

        return x

class MDParametersXMLPickler(DataContainerXMLPickler):

    order = ('dynamics', 'random_seed', 'tad_temp_high',
             'tad_timestep_factor', 'cartesian_temp_high',
             'cartesian_first_iteration', 'timestep',
             'temp_cool1_final', 'temp_cool2_final', 
             'steps_high', 'steps_refine',
             'steps_cool1', 'steps_cool2')
    
    def create(self):
        return MDParameters()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        e.dynamics = x['md_type']
        [setattr(e, k, x[k]) for k in x.keys() if k <> 'md_type']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)
        
        x['md_type'] = str(e.dynamics)
        x['random_seed'] = int(e.random_seed)
        x['tad_temp_high'] = float(e.tad_temp_high)
        x['cartesian_temp_high'] = float(e.cartesian_temp_high)
        x['temp_cool1_final'] = float(e.temp_cool1_final)
        x['temp_cool2_final'] = float(e.temp_cool2_final)
        x['timestep'] = float(e.timestep)
        x['tad_timestep_factor'] = float(e.tad_timestep_factor)
        x['steps_high'] = int(e.steps_high)
        x['steps_refine'] = int(e.steps_refine)
        x['steps_cool1'] = int(e.steps_cool1)
        x['steps_cool2'] = int(e.steps_cool2)
        x['cartesian_first_iteration'] = int(e.cartesian_first_iteration)

        return x

class FBHWParametersXMLPickler(DataContainerXMLPickler):

    order = ['m_rswitch_hot', 'm_rswitch_cool1', 'm_rswitch_cool2',
             'rswitch_hot', 'rswitch_cool1', 'rswitch_cool2',
             'm_asymptote_hot', 'm_asymptote_cool1', 'm_asymptote_cool2',
             'asymptote_hot', 'asymptote_cool1', 'asymptote_cool2']

    def create(self):
        return FBHWParameters()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        [setattr(e, k, x[k]) for k in x.keys()]

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)

        x['m_rswitch_hot'] = float(e.m_rswitch_hot)
        x['m_rswitch_cool1'] = float(e.m_rswitch_cool1)
        x['m_rswitch_cool2'] = float(e.m_rswitch_cool2)
        x['rswitch_hot'] = float(e.rswitch_hot)
        x['rswitch_cool1'] = float(e.rswitch_cool1)
        x['rswitch_cool2'] = float(e.rswitch_cool2)
        x['m_asymptote_hot'] = float(e.m_asymptote_hot)
        x['m_asymptote_cool1'] = float(e.m_asymptote_cool1)
        x['m_asymptote_cool2'] = float(e.m_asymptote_cool2)
        x['asymptote_hot'] = float(e.asymptote_hot)
        x['asymptote_cool1'] = float(e.asymptote_cool1)
        x['asymptote_cool2'] = float(e.asymptote_cool2)

        return x

class RDCParametersXMLPickler(DataContainerXMLPickler):

    order = ['parameter_class', 'method', 'first_iteration',
             'k_hot', 'k_cool1', 'k_cool2', 'r', 'd',
             'border_hot_initial', 'border_hot_final',
             'border_cool1_initial', 'border_cool1_final',
             'border_cool2_initial', 'border_cool2_final',
             'center_hot_initial', 'center_hot_final',
             'center_cool1_initial', 'center_cool1_final',
             'center_cool2_initial', 'center_cool2_final']
    
    def create(self):
        return RDCParameters()
    
    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        e.parameter_class = x['class']
        [setattr(e, k, x[k]) for k in x.keys() if k <> 'class']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)

        x['class'] = int(e.parameter_class)
        x['method'] = str(e.method)
        x['r'] = float(e.r)
        x['d'] = float(e.d) 
        x['first_iteration'] = int(e.first_iteration)
        x['k_hot'] = float(e.k_hot)
        x['k_cool1'] = float(e.k_cool1)
        x['k_cool2'] = float(e.k_cool2)
        x['border_hot_initial'] = float(e.border_hot_initial)
        x['border_cool1_initial'] = float(e.border_cool1_initial)
        x['border_cool2_initial'] = float(e.border_cool2_initial)
        x['border_hot_final'] = float(e.border_hot_final)
        x['border_cool1_final'] = float(e.border_cool1_final)
        x['border_cool2_final'] = float(e.border_cool2_final)
        x['center_hot_initial'] = float(e.center_hot_initial)
        x['center_cool1_initial'] = float(e.center_cool1_initial)
        x['center_cool2_initial'] = float(e.center_cool2_initial)
        x['center_hot_final'] = float(e.center_hot_final)
        x['center_cool1_final'] = float(e.center_cool1_final)
        x['center_cool2_final'] = float(e.center_cool2_final)

        return x
        
class KarplusParametersXMLPickler(DataContainerXMLPickler):

    order = ['parameter_class', 'a', 'b', 'c', 'd',
             'k_hot', 'k_cool1', 'k_cool2']
    
    def create(self):
        return KarplusParameters()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        e.parameter_class = x['class']
        [setattr(e, k, x[k]) for k in x.keys() if k <> 'class']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)

        x['class'] = int(e.parameter_class)
        x['a'] = float(e.a)
        x['b'] = float(e.b)
        x['c'] = float(e.c)
        x['d'] = float(e.d)
        x['k_hot'] = float(e.k_hot)
        x['k_cool1'] = float(e.k_cool1)
        x['k_cool2'] = float(e.k_cool2)

        return x

class DihedralParametersXMLPickler(DataContainerXMLPickler):

    order = ['k_hot', 'k_cool1', 'k_cool2']
    
    def create(self):
        return DihedralParameters()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        [setattr(e, k, x[k]) for k in x.keys()]

        e.set_tag_order(self.order)

        return e
    
    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)
        x['k_hot'] = float(e.k_hot)
        x['k_cool1'] = float(e.k_cool1)
        x['k_cool2'] = float(e.k_cool2)

        return x
    
## BARDIAUX 2.2: Symmetry specific restraints/potential
class SymmetryParametersXMLPickler(DataContainerXMLPickler):

    order = ['k_packing_hot', 'k_packing_cool1', 'k_packing_cool2', 'last_iteration_packing', 'k_ncs']
    
    def create(self):
        return SymmetryParameters()

    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        [setattr(e, k, x[k]) for k in x.keys()]

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):

        x = DataContainerXMLPickler.load_from_element(self, e)

        x['k_packing_hot'] = float(e.k_packing_hot)
        x['k_packing_cool1'] = float(e.k_packing_cool1)
        x['k_packing_cool2'] = float(e.k_packing_cool2)
        x['k_ncs'] = float(e.k_ncs)

        if hasattr(e, 'last_iteration_packing'):
            x['last_iteration_packing'] = int(e.last_iteration_packing)
        else:
            x['last_iteration_packing'] = 8
            
        return x
                
## CNS annealing parameters xml-pickler

class DRParametersXMLPickler(DataContainerXMLPickler):

    order = ['first_iteration', 'k_hot', 'k_cool1_initial',
             'k_cool1_final', 'k_cool2']
    
    def _xml_state(self, x):

        e = DataContainerXMLPickler._xml_state(self, x)
        e.first_iteration = x['first_iteration']
        e.k_hot = x['k_hot']
        e.k_cool1_initial = x['k_cool1_initial']
        e.k_cool1_final = x['k_cool1_final']
        e.k_cool2 = x['k_cool2']

        e.set_tag_order(self.order)

        return e

    def load_from_element(self, e):
        s = DataContainerXMLPickler.load_from_element(self, e)
        
        s['first_iteration'] = int(e.first_iteration)
        s['k_hot'] = float(e.k_hot)
        s['k_cool1_initial'] = float(e.k_cool1_initial)
        s['k_cool1_final'] = float(e.k_cool1_final)
        s['k_cool2'] = float(e.k_cool2)

        return s

class AmbiguousParametersXMLPickler(DRParametersXMLPickler):
    def create(self):
        return AmbiguousParameters()

class UnambiguousParametersXMLPickler(DRParametersXMLPickler):
    def create(self):
        return UnambiguousParameters()

class HBondParametersXMLPickler(DRParametersXMLPickler):
    def create(self):
        return HBondParameters()

## base class

DataContainer._xml_state = DataContainerXMLPickler()._xml_state
SimpleDataContainer._xml_state = SimpleDCXMLPickler()._xml_state

## data 

SequenceData._xml_state = SequenceDataXMLPickler()._xml_state
ShiftData._xml_state = ShiftDataXMLPickler()._xml_state
BoundCorrection._xml_state = BoundCorrectionXMLPickler()._xml_state
PeakData._xml_state = PeakDataXMLPickler()._xml_state
# BARDIAUX rMAt
ExperimentData._xml_state = ExperimentDataXMLPickler()._xml_state

TemplateData._xml_state = TemplateDataXMLPickler()._xml_state
InitialStructureData._xml_state = InitialStructureDataXMLPickler()._xml_state
HBondData._xml_state = HBondDataXMLPickler()._xml_state
DihedralData._xml_state = DihedralDataXMLPickler()._xml_state
KarplusData._xml_state = KarplusDataXMLPickler()._xml_state
RDCData._xml_state = RDCDataXMLPickler()._xml_state

CysPatch._xml_state = CysPatchXMLPickler()._xml_state
SSBondData._xml_state = SSBondDataXMLPickler()._xml_state


SSBridge._xml_state = SSBridgeXMLPickler()._xml_state
HisPatch._xml_state = HisPatchXMLPickler()._xml_state
# BARDIAUX2.2
Symmetry._xml_state = SymmetryXMLPickler()._xml_state
CCPNData._xml_state = CCPNDataXMLPickler()._xml_state
CisProPatch._xml_state = CisProPatchXMLPickler()._xml_state
ZnPatch._xml_state = ZnPatchXMLPickler()._xml_state
OtherData._xml_state = OtherDataXMLPickler()._xml_state

SpectrumData._xml_state = SpectrumDataXMLPickler()._xml_state
AmbiguousDistanceData._xml_state = AmbiguousDistanceDataXMLPickler()._xml_state
UnambiguousDistanceData._xml_state = UnambiguousDistanceDataXMLPickler()._xml_state

## misc

WaterRefinementParameters._xml_state = WaterRefinementXMLPickler()._xml_state

## CNS annealing parameters

AnnealingParameters._xml_state = AnnealingParametersXMLPickler()._xml_state
AmbiguousParameters._xml_state = AmbiguousParametersXMLPickler()._xml_state
UnambiguousParameters._xml_state = UnambiguousParametersXMLPickler()._xml_state
HBondParameters._xml_state = HBondParametersXMLPickler()._xml_state
MDParameters._xml_state = MDParametersXMLPickler()._xml_state
# BARDIAUX
SymmetryParameters._xml_state = SymmetryParametersXMLPickler()._xml_state

FBHWParameters._xml_state = FBHWParametersXMLPickler()._xml_state

# BERNARD2.3
LogHarmonicParameters._xml_state = LogHarmonicParametersXMLPickler()._xml_state

RDCParameters._xml_state = RDCParametersXMLPickler()._xml_state
KarplusParameters._xml_state = KarplusParametersXMLPickler()._xml_state
DihedralParameters._xml_state = DihedralParametersXMLPickler()._xml_state
