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


import aria.ariabase as ariabase
import aria.xmlutils as xmlutils

ASSIGNMENT_TYPE_MANUAL = 'manual'
ASSIGNMENT_TYPE_AUTOMATIC = 'automatic'

ASSIGNMENT_TYPES = (ASSIGNMENT_TYPE_MANUAL, ASSIGNMENT_TYPE_AUTOMATIC)

class Assignment(ariabase.AriaBaseClass):

    def __init__(self, atoms, assignment_type):

        ariabase.AriaBaseClass.__init__(self)

        ariabase.check_tuple(atoms)
        ariabase.check_elements(atoms, 'Atom')

        if not assignment_type in ASSIGNMENT_TYPES:
            e = 'Assignment type "%s" unknown; supported types %s'
            self.error(TypeError, e % (str(assignment_type),
                                       str(ASSIGNMENT_TYPES)))

        if not atoms:
            e = 'At least one atom expected.' 
            self.error(ValueError, e)
            
        self.atoms = atoms
        self.type = assignment_type

    def getType(self):

        return self.type

    def getAtoms(self):

        return self.atoms

    def __eq__(self, other):

        return self.type == other.type and \
               self.atoms == other.atoms

    def __str__(self):

        s = 'Assignment(atoms=%s, type=%s)'

        return s % (str(self.getAtoms()), self.getType())

    __repr__ = __str__

class AssignmentXMLPickler(xmlutils.XMLBasePickler):

    order = ('assignment_type', 'atom',)

    def __init__(self):
        xmlutils.XMLBasePickler.__init__(self)

    def _xml_state(self, a):

        e = xmlutils.XMLElement(tag_order = self.order)
        e.atom = a.getAtoms()
        e.assignment_type = a.getType()

        return e

    def load_from_element(self, e):

        import aria.tools as tools

        atoms = tools.as_tuple(e.atom)

        ## TODO: move to unpickler for version 1.1
        
        if hasattr(e, 'assignment_type'):
            type = str(e.assignment_type)
        else:
            type = ASSIGNMENT_TYPE_MANUAL

        return Assignment(atoms, type)

Assignment._xml_state = AssignmentXMLPickler()._xml_state
    
