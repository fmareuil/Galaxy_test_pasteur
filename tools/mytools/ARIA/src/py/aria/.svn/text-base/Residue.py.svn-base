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



from aria.TypeChecking import *
from aria.xmlutils import XMLBasePickler
from aria.ariabase import AriaBaseClass

class Residue(AriaBaseClass):

    def __init__(self, number, residue_type = None, structure = ""):

        check_int(number)
        check_type(residue_type, STRING, NONE)

        AriaBaseClass.__init__(self)

        self.number = number
        self.type = residue_type

        self.eq_groups = []
        self.atoms = {}
        ## BARDIAUX 2.2
        self.structure = structure

        self.setChain(None)

    def __getitem__(self, atom_name):
        check_string(atom_name)

        if not self.hasAtom(atom_name):
            self.error(KeyError, 'Atom %s not in residue %s%d.' \
                       % (atom_name, self.getType(), self.getNumber()))

        return self.atoms[atom_name]

    def __len__(self):
        return len(self.atoms)

    def _setType(self, residue_type):
        check_type(residue_type, STRING, NONE)
        
        if not self.type is None:
            m = 'Cannot set residue-type "%s" in %s; type already set.'
            self.error(ValueError, m % (str(residue_type), str(self)))

        self.type = residue_type
        
    def getType(self):
        return self.type

    def getNumber(self):
        return self.number

    def setNumber(self, number):        
        check_int(number)

        chain = self.getChain()
        if chain is not None and chain.hasResidue(number):
            m = 'Renumbering not possible. Residue with number %d is ' + \
                'already exists in chain "%s".'                
            self.error(IndexError, m % (number, chain.getSegid()))

        previous = self.number
        self.number = number
        
        if chain is not None:
            chain.delResidue(previous)
            chain.addResidue(self)

    def getName(self):
        """
        returns TypeNumber, e.g. ARG10
        """
        number = self.getNumber()
        if number is None:
            return None

        return '%s%d' % (str(self.getType()), number)

    ## BARDIAUX 2.2
    def getStructure(self):
        return self.structure

    def hasAtom(self, atom_name):
        check_string(atom_name)

        return self.atoms.has_key(atom_name)
        
    def addAtom(self, atom):
        check_type(atom, 'Atom')

        name = atom.getName()

        if self.hasAtom(name):
            self.error(KeyError, 'Atom "%s" already in %s.' \
                       % (name, str(self)))
        
        self.atoms[name] = atom
        self.atoms[name].setResidue(self)

    def getAtoms(self):
        return self.atoms.values()

    def getAtom(self, atom_name):

        return self[atom_name]

    def setChain(self, chain):
        check_type(chain, 'Chain', NONE)

        self.chain = chain

    def getChain(self):
        return self.chain

    def addEquivalentGroup(self, x):
        check_type(x, 'EquivalentGroup')

        if x in self.eq_groups:
            s = 'Group (type %s) already defined in residue.'
            self.error(ValueError, s % x.getType())

        for atom_name in x.getAtomNames():
            if not atom_name in self.atoms:
                self.error(ValueError, 'Atom "%s" not known.' % atom_name)

        self.eq_groups.append(x)

    def getEquivalentGroups(self):
        return self.eq_groups

    def link(self):

        ## equivalent groups

        for group in self.getEquivalentGroups():
            atoms = [self[name] for name in group.getAtomNames()]
            group.setAtoms(atoms)

        ## hetero atoms

        for atom in self.getAtoms():
            hetero_name = atom.getSettings()['hetero_atom_name']
            
            if hetero_name is None:
                continue

            elif not self.hasAtom(hetero_name):
                m = 'Linking failed: hetero atom "%s" not known in %s.'
                self.error(ValueError, m % (hetero_name, str(self)))
                
            hetero_atom = self[hetero_name]
            atom.setHeteroAtom(hetero_atom)

    def __str__(self):
        atoms = self.getAtoms()
        atoms.sort(lambda a, b: cmp(a.getName(), b.getName()))

        class_name = self.__class__.__name__

        atoms = ('\n  ' + ' ' * len(class_name)).join(map(str, atoms))

        string = '%s(number=%s, type=%s,\n  ' + ' ' * len(class_name) + '%s)'

        return string % (class_name, self.getNumber(), self.getType(), atoms)

    __repr__ = __str__

class ResidueXMLPickler(XMLBasePickler):

## TODO
##    order = ('number', 'residue_type', 'atom', 'equivalent_group')

    def _xml_state(self, residue):
                   
        from aria.xmlutils import XMLElement

##        e = XMLElement(tag_order = self.order)
        e = XMLElement()

        e.residue_type = residue.getType()
        e.number = residue.getNumber()
        e.structure=residue.getStructure()
        atoms = residue.getAtoms()
        atoms.sort(lambda a, b: cmp(a.getName(), b.getName()))
        e.atom = atoms
        if len(residue.getEquivalentGroups()):
            e.equivalent_group = residue.getEquivalentGroups()

        return e

Residue._xml_state = ResidueXMLPickler()._xml_state
