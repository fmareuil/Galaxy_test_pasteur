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
from aria.xmlutils import XMLElement, XMLBasePickler
from aria.Settings import Settings

EQUIV_METHYL = 'METHYL'
EQUIV_METHYLENE = 'METHYLENE'
EQUIV_ISOPROPYL = 'ISOPROPYL'
EQUIV_AROMATIC = 'AROMATIC'
EQUIV_NTERMINUS = 'NTERMINUS'

TERMINUS_N_STANDARD = 'N_STANDARD'
TERMINUS_N_AMINYL = 'N_AMINYL'
TERMINUS_C_STANDARD = 'C_STANDARD'
TERMINUS_C_AMIDO = 'C_AMIDO'
TERMINUS_C_CARBOXYL = 'C_CARBOXYL'
TERMINUS_C5_PRIME_PHOSPHATE = 'C5_PRIME_PHOSPHATE'
TERMINUS_C5_PRIME_HYDROXYL = 'C5_PRIME_HYDROXYL'
TERMINUS_C3_PRIME_HYDROXYL = 'C3_PRIME_HYDROXYL'

TYPE_AMINO_ACID = 'AMINO_ACID'
TYPE_DNA_BASE = 'DNA_BASE'
TYPE_RNA_BASE = 'RNA_BASE'
TYPE_NONBASE = 'NONBASE'

ATOM_TYPE_UNKNOWN = 'UNKNOWN'

TOPOLOGY_IUPAC = 'iupac.xml'

NOMENCLATURE_IUPAC = 'IUPAC'

ATOM_NOMENCLATURES = (NOMENCLATURE_IUPAC,)

class BaseAtomSettings(Settings):

    def create(self):

        from aria.Settings import ChoiceEntity, String

        d = {}

        types = ('H','O','C','N','S','P','ZN','Zn', ATOM_TYPE_UNKNOWN)
        d['type'] = ChoiceEntity(elements = types)
        d['hetero_atom_name'] = String(optional = 1)

        return d

class AtomGroupSettings(Settings):
    pass

class EquivalentGroupSettings(Settings):

    def create(self):

        from aria.Settings import ChoiceEntity, TypeEntity

        d = {}

        types = (EQUIV_METHYL, EQUIV_METHYLENE, EQUIV_ISOPROPYL,
                 EQUIV_AROMATIC, EQUIV_NTERMINUS)
        
        d['type'] = ChoiceEntity(elements = types)
        d['atom_names'] = TypeEntity(TUPLE)

        return d

class TerminusSettings(Settings):

    def create(self):

        from aria.Settings import ChoiceEntity, TypeEntity

        d = {}

        types = (TERMINUS_N_STANDARD, TERMINUS_N_AMINYL,
                 TERMINUS_C_STANDARD, TERMINUS_C_AMIDO,
                 TERMINUS_C_CARBOXYL, TERMINUS_C5_PRIME_PHOSPHATE,
                 TERMINUS_C5_PRIME_HYDROXYL, TERMINUS_C3_PRIME_HYDROXYL)
        
        d['type'] = ChoiceEntity(elements = types)
        d['atom_names'] = TypeEntity(TUPLE)

        return d

class BaseAtom(AriaBaseClass):

    def __init__(self, settings, name):
        check_type(settings, 'BaseAtomSettings')
        check_string(name)

        AriaBaseClass.__init__(self, settings)

        self.setName(name)
        self.setHeteroAtom(None)

        self.eq_groups = []

    def setName(self, name):
        check_string(name)
        self.name = name

    def getName(self):
        return self.name

    def getType(self):
        return self.getSettings()['type']

    def setType(self, t):
        self.getSettings()['type'] = t

    def _addEquivalentGroup(self, g):
        check_type(g, 'EquivalentGroup')

        if g in self.eq_groups:
            s = 'Equiv. group (type %s) does already exist for atom %s. '
            self.warning(s % (str(self.eq_groups), str(self)))
        else:
            self.eq_groups.append(g)

    def getEquivalentGroups(self):
        return self.eq_groups

    def setHeteroAtom(self, a):
        check_type(a, 'BaseAtom', NONE)

        if a is not None:
            hetero_atom = self.getSettings()['hetero_atom_name']

            if hetero_atom is None:
                self.error(ValueError, 'Atom cannot have hetero atom.')
            elif not a.getName() == hetero_atom:
                m = 'Atom "%s" can only have "%s" as hetero atom; ' + \
                    'atom "%s" passed.'
                self.error(ValueError, m % (self.getName(), a.getName(),
                                            hetero_atom))
                
        self.hetero = a

    def getHeteroAtom(self):
        return self.hetero

    def isProton(self):
        return self.getType()[0] == 'H'

    def equal_properties(self, other):
        """
        checks whether both atoms have the same

        - type
        - hetero_atom_name
        - name
        """

        check_type(other, 'BaseAtom')

        s = self.getSettings()
        other_s = other.getSettings()

        if s.as_dict() <> other_s.as_dict():
            return 0
        elif self.getName() <> other.getName():
            return 0

        return 1
    
    def __str__(self):

        s = self.getSettings()

        class_name = self.__class__.__name__
        
        return '%s(name=%s, type=%s, hetero_name=%s)' % \
               (class_name, self.getName(), s['type'],
                s['hetero_atom_name'])

    __repr__ = __str__

class BaseAtomXMLPickler(XMLBasePickler):

    order = ('name', 'atom_type', 'hetero')

    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)

        s = x.getSettings()
        e.name = x.getName()
        e.atom_type = s['type']
        e.hetero = s['hetero_atom_name']

        return e
    
    def load_from_element(self, e):

        s = BaseAtomSettings()
        s['type'] = str(e.atom_type)

        try:
            hetero = eval(e.hetero)
        except:
            hetero = str(e.hetero)
        
        s['hetero_atom_name'] = hetero

        x = BaseAtom(s, str(e.name))
        
        return x

class AtomGroup(AriaBaseClass):
    """
    we assume, that 'settings' have a entity 'type'
    """

    def __init__(self, settings):
        check_type(settings, 'Settings')
        AriaBaseClass.__init__(self, settings)

        self.setAtoms(())

    def setAtoms(self, a):
        check_type(a, LIST, TUPLE)
        check_elements(a, 'BaseAtom')

        self.atoms = a

    def getAtoms(self):
        return self.atoms

    def getType(self):
        return self.getSettings()['type']

    def setType(self, t):
        self.getSettings()['type'] = t

    def getAtomNames(self):
        return self.getSettings()['atom_names']

    def __str__(self):

        s = self.getSettings()

        class_name = self.__class__.__name__
        
        return '%s(type=%s, atom_names=%s)' % \
               (class_name, s['type'], str(s['atom_names']))

    __repr__ = __str__

class AtomGroupXMLPickler(XMLBasePickler):

    order = ['atom_name']

    def make_atom_list(self, names):

        l = []

        for name in names:
            e = XMLElement()
            e.name = name
            l.append(e)

        return l

    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)

        atom_names = list(x.getSettings()['atom_names'])
        atom_names.sort()
        
        e.atom_name = self.make_atom_list(atom_names)

        return e

    def create_settings(self):
        return AtomGroupSettings()

    def create(self, s):
        return AtomGroup(s)

    def load_from_element(self, e):
        from aria.tools import as_tuple
        
        s = self.create_settings()

        if hasattr(e, 'atom_name'):
            atom_names = [str(e.name) for e in as_tuple(e.atom_name)]
            s['atom_names'] = tuple(atom_names)
            
        x = self.create(s)

        return x

class EquivalentGroup(AtomGroup):
    def __init__(self, settings):
        check_type(settings, 'EquivalentGroupSettings')
        AtomGroup.__init__(self, settings)

    def setAtoms(self, atoms):
        AtomGroup.setAtoms(self, atoms)
        for a in atoms:
            a._addEquivalentGroup(self)

class EquivalentGroupXMLPickler(AtomGroupXMLPickler):

    order = AtomGroupXMLPickler.order + ['group_type']

    def _xml_state(self, x):

        e = AtomGroupXMLPickler._xml_state(self, x)
        e.group_type = x.getSettings()['type']

        e.set_tag_order(self.order)

        return e

    def create_settings(self):
        return EquivalentGroupSettings()

    def create(self, s):
        return EquivalentGroup(s)

    def load_from_element(self, e):

        x = AtomGroupXMLPickler.load_from_element(self, e)

        x.getSettings()['type'] = str(e.group_type)

        return x

class Terminus(AtomGroup):
    def __init__(self, settings):
        check_type(settings, 'TerminusSettings')
        AtomGroup.__init__(self, settings)

class TerminusXMLPickler(AtomGroupXMLPickler):

    order = AtomGroupXMLPickler.order + ['terminus_type']

    def _xml_state(self, x):

        e = AtomGroupXMLPickler._xml_state(self, x)
        e.terminus_type = x.getSettings()['type']

        e.set_tag_order(self.order)

        return e

    def create_settings(self):
        return TerminusSettings()

    def create(self, s):
        return Terminus(s)

    def load_from_element(self, e):

        x = AtomGroupXMLPickler.load_from_element(self, e)

        x.getSettings()['type'] = str(e.terminus_type)

        return x

class BaseResidueSettings(Settings):

    def create(self):

        from aria.Settings import ChoiceEntity

        d = {}

        types = (TYPE_AMINO_ACID, TYPE_DNA_BASE, TYPE_RNA_BASE, TYPE_NONBASE)
        d['type'] = ChoiceEntity(elements = types)

        return d

class BaseResidue(AriaBaseClass):

    def __init__(self, settings, name):

        check_type(settings, 'BaseResidueSettings')
        check_string(name)

        AriaBaseClass.__init__(self, settings)
        
        self.setName(name)
        self.atoms = {}
        self.termini = {}
        self.eq_groups = []
        
        self.setBackboneAtomNames(())
        self.setSidechainAtomNames(())

    def setName(self, name):
        check_string(name)
        self.name = name

    def getName(self):
        return self.name

    def getType(self):
        return self.getSettings()['type']

    def addAtom(self, x):
        check_type(x, 'BaseAtom')
        
        name = x.getName()

        if name in self.atoms:
            raise KeyError, 'Atom "%s" does already exist.' % name

        self.atoms[name] = x

    def getAtoms(self):
        return self.atoms.values()
    
    def addEquivalentGroup(self, x):
        check_type(x, 'EquivalentGroup')

        if x in self.eq_groups:
            s = 'Group (type %s) already defined in residue.'
            raise ValueError, s % x.getType()

        for atom_name in x.getAtomNames():
            if not atom_name in self.atoms:
                raise ValueError, 'Atom "%s" not known.' % atom_name

        self.eq_groups.append(x)

    def getEquivalentGroups(self):
        return self.eq_groups

    def addTerminus(self, x):
        check_type(x, 'Terminus')

        if x.getType() in self.termini:
            s = 'Terminus (type %s) already defined in residue.'
            raise ValueError, s % x.getType()

        for atom_name in x.getAtomNames():
            if not atom_name in self.atoms:
                raise ValueError, 'Atom "%s" not known.' % atom_name

        self.termini[x.getType()]= x

    def getTermini(self):
        return self.termini

    def setBackboneAtomNames(self, x):
        check_type(x, LIST, TUPLE)
        self.backbone_atom_names = x
        self.backbone_atoms = None

    def getBackboneAtomNames(self):
        return self.backbone_atom_names

    def getBackboneAtoms(self):
        return self.backbone_atoms

    def setBackboneAtoms(self, l):
        check_type(l, LIST, TUPLE)
        check_elements(l, 'BaseAtom')
        
        self.backbone_atoms = l
    
    def setSidechainAtomNames(self, x):
        check_type(x, LIST, TUPLE)
        self.sidechain_atom_names = x
        self.sidechain_atoms = None

    def getSidechainAtomNames(self):
        return self.sidechain_atom_names

    def getSidechainAtoms(self):
        return self.sidechain_atoms

    def setSidechainAtoms(self, l):
        check_type(l, LIST, TUPLE)
        check_elements(l, 'BaseAtom')
        
        self.sidechain_atoms = l

    def find(self, atom_name = None):
        check_type(atom_name, STRING, NONE)
        
        if atom_name is not None:
            if atom_name not in self.atoms:
                s = 'Residue %s: atom %s not known.'
                self.error(ValueError, s % (self.getName(), atom_name))
                
            return self.atoms[atom_name]

    def link(self):
        """
        converts atoms, referenced by its name, into
        BaseAtom instances.
        """

        ## backbone

        atoms = [self.find(atom_name = name) for name in \
                 self.getBackboneAtomNames()]

        self.setBackboneAtoms(atoms)

        ## sidechain

        atoms = [self.find(atom_name = name) for name in \
                 self.getSidechainAtomNames()]

        self.setSidechainAtoms(atoms)

        ## terminii

        for terminus in self.getTermini().values():
            atoms = [self.find(atom_name = name) for name in \
                     terminus.getAtomNames()]
            terminus.setAtoms(atoms)

        ## equivalent groups

        for group in self.getEquivalentGroups():
            atoms = [self.find(atom_name = name) for name in \
                     group.getAtomNames()]
            group.setAtoms(atoms)

        ## hetero atoms

        for atom in self.getAtoms():
            hetero_name = atom.getSettings()['hetero_atom_name']
            
            if hetero_name is None:
                continue

            try:
                hetero_atom = self.find(atom_name = hetero_name)
            except:
                print repr(hetero_name)

            atom.setHeteroAtom(hetero_atom)

    def __str__(self):

        name = self.__class__.__name__

        s = '%s(name=%s, type=%s, atoms=%s, equivalent_groups=%s, ' + \
            'termini=%s, ' + \
            'backbone_atom_names=%s, sidechain_atom_names=%s)'

        return s % (name, self.getName(), self.getType(), \
                    str(self.getAtoms()),
                    str(self.getEquivalentGroups()),
                    str(self.getTermini().values()),
                    str(self.getBackboneAtomNames()),
                    str(self.getSidechainAtomNames()))

    __repr__ = __str__

class BaseResidueXMLPickler(XMLBasePickler):

    order = ['name', 'residue_type', 'atom', 'backbone',
             'sidechain', 'terminus', 'equivalent_group']

    def make_atom_list(self, names):
        
        l = []

        for name in names:
            element = XMLElement()
            element.name = name
            l.append(element)
        
        return l
    
    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)
        e.name = x.getName()
        e.residue_type = x.getSettings()['type']
        
        atoms = x.getAtoms()
        atoms.sort(lambda a, b, cmp = cmp: cmp(a.getName(), b.getName()))
        e.atom = atoms
        
        names = list(x.getBackboneAtomNames())
        names.sort()

        bb = XMLElement()
        bb.atom_name = self.make_atom_list(names)
        e.backbone = bb

        names = list(x.getSidechainAtomNames())
        names.sort()

        sc = XMLElement()
        sc.atom_name = self.make_atom_list(names)
        e.sidechain = sc

        e.terminus = x.getTermini().values()
        
        eq_groups = x.getEquivalentGroups()
        
        if eq_groups:
            e.equivalent_group = eq_groups

        return e

    def get_atom_names(self, elements):
        from aria.tools import as_tuple
        
        return [str(e.name) for e in as_tuple(elements)]

    def load_from_element(self, e):

        from aria.tools import as_tuple

        s = BaseResidueSettings()
        s['type'] = str(e.residue_type)

        x = BaseResidue(s, str(e.name))
        
        [x.addAtom(a) for a in as_tuple(e.atom)]

        if hasattr(e.backbone, 'atom_name'):
            atom_names = self.get_atom_names(e.backbone.atom_name)
            x.setBackboneAtomNames(tuple(atom_names))

        if hasattr(e.sidechain, 'atom_name'):
            atom_names = self.get_atom_names(e.sidechain.atom_name)
            x.setSidechainAtomNames(tuple(atom_names))

        if hasattr(e, 'equivalent_group'):
            [x.addEquivalentGroup(a) for a in as_tuple(e.equivalent_group)]
            
        [x.addTerminus(a) for a in as_tuple(e.terminus)]

        return x

class Topology(AriaBaseClass):

    def __init__(self):
        AriaBaseClass.__init__(self)

        d = {}
        d[TYPE_AMINO_ACID] = {}
        d[TYPE_RNA_BASE] = {}
        d[TYPE_DNA_BASE] = {}
        d[TYPE_NONBASE] = {}
        
        self.residues = d
        
        self.setName('')
        self.setDescription('')

    def setName(self, s):
        check_string(s)
        self.name = s

    def getName(self):
        return self.name

    def setDescription(self, s):
        check_string(s)
        self.description = s

    def getDescription(self):
        return self.description

    def addResidue(self, r):
        check_type(r, 'BaseResidue')

        type = r.getType()
        name = r.getName()
        
        if name in self.residues[type]:
            s = 'Residue %s is already defined. Old definition' + \
                ' updated.'
            self.warning(s % name)

        self.residues[type][name] = r

        ## resolve atom names.
        
        r.link()

    def getAminoAcids(self):
        return self.residues[TYPE_AMINO_ACID].values()

    def getDNABases(self):
        return self.residues[TYPE_DNA_BASE].values()

    def getRNABases(self):
        return self.residues[TYPE_RNA_BASE].values()

    def getResidues(self):
        return self.getAminoAcids() + self.getRNABases() + \
               self.getDNABases()

    def __find(self, name, dict, err_msg):
        if name not in dict:
            self.error(KeyError, err_msg)
        else:
            return dict[name]

    def find(self, amino_acid = None, dna_base = None,
             rna_base = None):

        if amino_acid is not None:
            check_string(amino_acid)
            return self.__find(amino_acid, self.residues[TYPE_AMINO_ACID],
                               'Amino acid %s not known.' % amino_acid)

        if dna_base is not None:
            check_string(dna_base)
            return self.__find(dna_base, self.residues[TYPE_DNA_BASE],
                               'DNA base %s not known.' % dna_base)

        if rna_base is not None:
            check_string(rna_base)
            return self.__find(rna_base, self.residues[TYPE_RNA_BASE],
                               'RNA base %s not known.' % rna_base)

        return None

    def check(self, x):
        """
        checks whether the following attributes of a BaseResidue, 'x',
        comply with those defined by 'Topology':

        - name
        - atom names, atom type, hetero atom name

        returns None, if x is ok.
        """
        
        check_type(x, 'BaseResidue', 'Residue')

        ## check name

        if is_type(x, 'BaseResidue'):
            name = x.getName()
            type = x.getType()

        elif is_type(x, 'Residue'):
            import aria.Chain as C
            name = x.getType()
            chain_type = x.getChain().getType()
            if chain_type == C.TYPE_PROTEIN:
                type = TYPE_AMINO_ACID
            elif chain_type == C.TYPE_DNA:
                type = TYPE_DNA_BASE
            elif chain_type == C.TYPE_RNA:
                type = TYPE_RNA_BASE
            else:
                self.error(TypeError, 'Chain type "%s" not known.' \
                           % str(chain_type))

        if not name in self.residues[type]:
            return 'Residue %s not known.' % name

        ref = self.residues[type][name]

        ## check atom content

        missing = []

        for a_ref in ref.getAtoms():
            if not [a for a in x.getAtoms() if a.equal_properties(a_ref)]:
                missing.append(a_ref)

        additional = []

        for a in x.getAtoms():
            if not [a_ref for a_ref in ref.getAtoms()
                    if a_ref.equal_properties(a)]:
                
                additional.append(a)

        s = None

        if missing and not is_type(x, 'Residue'):
            s = ('Residue %s contains less atoms than its IUPAC ' + \
                ' definition: %s.') % (name, str(missing))

        if additional:
            if s is None:
                s = 'Residue %s' % x.getName()
            s += ' contains more atoms than its IUPAC definition: %s.' \
                 % str(additional)

        return s

    def __str__(self):

        name = self.__class__.__name__

        s = '%s(name=%s, description=%s, residues=%s)'

        return s % (name, self.getName(), self.getDescription(),
                    str(self.getResidues()))

    __repr__ = __str__

class TopologyXMLPickler(XMLBasePickler):

    order = ['name', 'description', 'residue']

    def __init__(self):

        self.topology = self
        self.residue = BaseResidueXMLPickler()
        self.atom = BaseAtomXMLPickler()
        self.equivalent_group = EquivalentGroupXMLPickler()
        self.terminus = TerminusXMLPickler()

    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)

        e.name = x.getName()
        e.description = x.getDescription()

        residues = x.getResidues()
        residues.sort(lambda a, b, cmp = cmp: cmp(a.getName(), b.getName()))
        
        e.residue = residues

        return e

    def load_from_element(self, e):

        from aria.tools import as_tuple

        x = Topology()
        x.setName(str(e.name))
        x.setDescription(str(e.description))

        if hasattr(e, 'residue'):
            [x.addResidue(a) for a in as_tuple(e.residue)]

        return x

def load_topology(name = TOPOLOGY_IUPAC):
    import os
    import aria.AriaXML as AriaXML

    filename = os.path.join(AriaBaseClass.data_path, name)

    p = AriaXML.AriaXMLPickler()

    return p.load(filename)
    
## needed for xml pickling (dump)
        
BaseAtom._xml_state = BaseAtomXMLPickler()._xml_state
AtomGroup._xml_state = AtomGroupXMLPickler()._xml_state
EquivalentGroup._xml_state = EquivalentGroupXMLPickler()._xml_state
Terminus._xml_state = TerminusXMLPickler()._xml_state
BaseResidue._xml_state = BaseResidueXMLPickler()._xml_state
Topology._xml_state = TopologyXMLPickler()._xml_state

