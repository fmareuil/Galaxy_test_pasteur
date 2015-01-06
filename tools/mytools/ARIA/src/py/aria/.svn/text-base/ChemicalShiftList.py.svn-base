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
from aria.TypeChecking import check_type
from aria.xmlutils import XMLBasePickler

class ChemicalShiftList(AriaBaseClass):

    def __init__(self):

        AriaBaseClass.__init__(self)

        self.__shift_assignments = []

    def addShiftAssignment(self, assignment):
        
        check_type(assignment, 'ShiftAssignment')

        self.__shift_assignments.append(assignment)

    def getShiftAssignments(self):

        return tuple(self.__shift_assignments)

    def __len__(self):

        return len(self.__shift_assignments)

class ChemicalShiftListXMLPickler(XMLBasePickler):

    def __init__(self):

        from aria.ShiftAssignment import ShiftAssignmentXMLPickler
        from aria.ShiftAssignment import SpinSystemXMLPickler
        from aria.Atom import AtomXMLPickler
        from aria.Datum import ChemicalShiftXMLPickler

        XMLBasePickler.__init__(self)
        
        self.chemical_shift_list = self
        self.shift_assignment = ShiftAssignmentXMLPickler()
        self.spin_system = SpinSystemXMLPickler()
        self.atom = AtomXMLPickler()
        self.chemical_shift = ChemicalShiftXMLPickler()

    def _xml_state(self, shift_list):

        from aria.xmlutils import XMLElement

        e = XMLElement()
        e.shift_assignment = shift_list.getShiftAssignments()

        return e

    def load_from_element(self, e):
        
        x = ChemicalShiftList()

        [x.addShiftAssignment(a) for a in e.shift_assignment]

        return x

ChemicalShiftList._xml_state = ChemicalShiftListXMLPickler()._xml_state
