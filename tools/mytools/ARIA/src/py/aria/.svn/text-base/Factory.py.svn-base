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



from aria.ariabase import AriaBaseClass
from aria.TypeChecking import *
from aria.Singleton import Singleton

class SpinPairFactory(AriaBaseClass, Singleton):

    __instances = {}
    
    def __init__(self):

        Singleton.__init__(self)
        AriaBaseClass.__init__(self)

    def __call__(self, atom1, atom2):

        from aria.SpinPair import SpinPair
        
        check_type(atom1, 'Atom')
        check_type(atom2, 'Atom')

        id1 = atom1.getId()
        id2 = atom2.getId()

        key = (min((id1, id2)), max((id1, id2)))

        instances = self.__class__.__instances
        
        if not instances.has_key(key):
            instances[key] = SpinPair(len(instances), atom1, atom2)

        return instances[key]

class SpinPairListFactory(AriaBaseClass, Singleton):

    __instances = {}
    
    def __init__(self):

        Singleton.__init__(self)
        AriaBaseClass.__init__(self)

    def __call__(self, spin_pairs):

        check_type(spin_pairs, LIST, TUPLE)

        ids = map(lambda sp: sp.getId(), spin_pairs)
        ids.sort()
        key = tuple(ids)

        instances = self.__class__.__instances
        if not instances.has_key(key):
            instances[key] = tuple(spin_pairs)

        return instances[key]

class AtomFactory(AriaBaseClass, Singleton):

    __atoms = {}
    
    def __init__(self):

        Singleton.__init__(self)

        self.unfreeze()

        AriaBaseClass.__init__(self)

    def reset(self):
        self.__class__.__atoms = {}

    def unfreeze(self):
        self.__frozen = 0

    def freeze(self):
        self.__frozen = 1

    def isfrozen(self):
        return self.__frozen

    def __get_atom_repr(self, atom_key):
        check_tuple(atom_key)

        s = '[segid=%s, residue_number=%d, name=%s]'
        return s % atom_key

    def createAtom(self, segid, residue_number, atom_name, atom_type = None,
                   hetero_atom_name = None):

        check_type(segid, STRING, NONE)
        check_int(residue_number)
        check_string(atom_name)
        check_type(atom_type, STRING, NONE)
        check_type(hetero_atom_name, STRING, NONE)

        from aria.Atom import Atom, AtomSettings

        if segid is not None:
            
            segid = '%4s' % segid

            if len(segid) > 4:

                m = 'SEGIDs are at most of length 4; "%s" given' % segid
                self.error(ValueError, m)

        atoms = self.__class__.__atoms

        key = (segid, residue_number, atom_name)

        if self.isfrozen() and not atoms.has_key(key):

            m = \
"""              
Could not create atom '%s': AtomFactory is frozen and atom has not been created yet. Usually, this error occurs if there is some inconsistency between molecule description (molecule XML file) and data-files (either spectra or chemical-shifts) or PDB files (is the segid correct?). Atoms might be missing in the molecule description but are referenced by the data-files. 
"""
            self.error(ValueError, m % self.__get_atom_repr(key))

        elif not atoms.has_key(key):

            atom = Atom(AtomSettings(), atom_name, len(atoms))
            atom._setSegid(segid)

            if atom_type is not None:
                atom.setType(atom_type)

            atom.getSettings()['hetero_atom_name'] = hetero_atom_name

            atoms[key] = atom

        return atoms[key]
        
