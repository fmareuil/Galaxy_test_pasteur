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

class SpinPair(AriaBaseClass):

    def __init__(self, id, atom1, atom2):

        check_int(id)
        check_type(atom1, 'Atom')
        check_type(atom2, 'Atom')

        if atom1 == atom2:
            self.error(ValueError, 'Atoms are equal')

        self.__atoms = (atom1, atom2)

        self.__id = id

    def getId(self):
        return self.__id

    def getAtoms(self):
        return self.__atoms

    def __getitem__(self, n):

        check_int(n)

        if n < 0 or n > 1:
            raise IndexError, 'index out of range. must be 0 or 1.'

        return self.getAtoms()[n]

    def __str__(self):

        return 'SpinPair(id=%s, %s, %s)' % (str(self.getId()),
                                            str(self[0]),
                                            str(self[1]))

    __repr__ = __str__

class SpinPairXMLPickler:

    def _xml_state(self, x):

        from aria.xmlutils import XMLElement

        e = XMLElement()

        e.atom = x.getAtoms()
        return e

    def load_from_element(self, e):

        from aria.Singleton import SpinPairFactory

        factory = SpinPairFactory()

        atom1, atom2 = e.atom
        return factory(atom1, atom2)

SpinPair._xml_state = SpinPairXMLPickler()._xml_state

    
