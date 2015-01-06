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


## TODO:
## - default values for chemical shift errors: None?
## - how to access assignments from the assigned atoms?

from aria.ariabase import *
from aria.xmlutils import XMLBasePickler

AVERAGING_METHOD_FAST = 'FAST'
AVERAGING_METHOD_SLOW = 'SLOW'
AVERAGING_METHOD_NONE = 'NONE'

AVERAGING_METHODS = (AVERAGING_METHOD_FAST,
                     AVERAGING_METHOD_SLOW,
                     AVERAGING_METHOD_NONE)

ASSIGNMENT_METHOD_STEREO_SPECIFIC = 'STEREO_SPECIFIC'
ASSIGNMENT_METHOD_EQUIVALENT = 'EQUIVALENT'
ASSIGNMENT_METHOD_FLOATING = 'FLOATING'

ASSIGNMENT_METHODS = (ASSIGNMENT_METHOD_STEREO_SPECIFIC,
                      ASSIGNMENT_METHOD_EQUIVALENT,
                      ASSIGNMENT_METHOD_FLOATING)

class SpinSystem(AriaBaseClass):
    
    def __init__(self, averaging_method):
        
        check_string(averaging_method)
        if not averaging_method in AVERAGING_METHODS:
            s = 'Averaging method "%s" not known. Valid methods are %s'
            self.error(ValueError, s % (averaging_method,
                                        str(AVERAGING_METHODS)))

        self.__averaging_method = averaging_method
        self.__atoms = ()
        self.__chemical_shifts = ()

    def getAveragingMethod(self):
        return self.__averaging_method

    def getChemicalShifts(self):
        return self.__chemical_shifts

    def setChemicalShifts(self, l):
        check_type(l, TUPLE)
        
        if not l:
            s = 'At least 1 chemical shift expected.'
            self.error(ValueError, s)
            
        check_elements(l, 'ChemicalShift')

        self.__chemical_shifts = l

    def setAtoms(self, l):
        check_tuple(l)
        check_elements(l, 'Atom')

        if not len(l):
            s = 'At least 1 atom expected.'
            self.error(ValueError, s)
            
        self.__atoms = l

    def getAtoms(self):
        return self.__atoms

    def __str__(self):

        name = self.__class__.__name__

        s = '%s(averaging_method=%s, shifts=%s, atoms=%s)'
        return s % (name,
                    self.getAveragingMethod(),
                    str(self.getChemicalShifts()),
                    str(self.getAtoms()))

    # BARDIAUX 2.2
    # introduce __eq__ for importFromCccpn.getAria2ConstraintList...
    # __hash__ required if we want SpinSyetem to be hashable
    def __hash__(self):
        return id(self)
    
    def __eq__(self, other):

        return self.__atoms == other.__atoms and \
               self.__averaging_method  == other.__averaging_method and \
               self.__chemical_shifts == other.__chemical_shifts
    

    __repr__ = __str__

class ShiftAssignment(AriaBaseClass):

    def __init__(self, assignment_method):
        """
        Creates an Assignment instance from a list of Atom instance
        and a list of ppm-values.
        """

        check_string(assignment_method)
        
        if not assignment_method in ASSIGNMENT_METHODS:
            s = 'Method %s not known. Valid methods are %s.'
            self.error(ValueError, s % (assignment_method,
                                        ASSIGNMENT_METHODS))
            
        AriaBaseClass.__init__(self)

        self.__method = assignment_method
        self.__spin_systems = ()

        self.is_valid(1)

    def setSpinSystems(self, l):
        check_type(l, TUPLE)
        check_elements(l, 'SpinSystem')

        if self.getMethod() <> ASSIGNMENT_METHOD_FLOATING and len(l) > 1:

            m = 'Assignment of type "%s" can only have one spin-system.'
            self.error(ValueError, m % self.getMethod())
            
        self.__spin_systems = l

    def getMethod(self):
        return self.__method
        
    def getSpinSystems(self):
        return self.__spin_systems

    def is_valid(self, valid = None):
        
        check_type(valid, INT, NONE)

        if valid is not None:
            self.__valid = valid
        else:
            return self.__valid == 1
    
    def __str__(self):

        class_name = self.__class__.__name__

        s = "%s(method=%s, spin_systems=%s)"
        
        return s % (class_name,
                    self.getMethod(),
                    str(self.getSpinSystems()))

    __repr__ = __str__

class ShiftAssignmentXMLPickler(XMLBasePickler):

    def _xml_state(self, a):
        from aria.xmlutils import XMLElement
        
        e = XMLElement()
        e.method = a.getMethod()
        e.spin_system = a.getSpinSystems()
        
        return e

    def load_from_element(self, e):
        from aria.tools import as_tuple

        method = str(e.method)
        ss = as_tuple(e.spin_system)
        
        sa = ShiftAssignment(method)
        sa.setSpinSystems(ss)

        return sa

class SpinSystemXMLPickler(XMLBasePickler):

    order = ['averaging_method', 'atom', 'chemical_shift']

    def _xml_state(self, a):

        from aria.xmlutils import XMLElement
        
        e = XMLElement(tag_order = self.order)

        e.atom = a.getAtoms()
        e.chemical_shift = a.getChemicalShifts()
        e.averaging_method = a.getAveragingMethod()

        return e

    def load_from_element(self, e):

        from aria.tools import as_tuple

        averaging_method = str(e.averaging_method)

        ss = SpinSystem(averaging_method)

        atoms = as_tuple(e.atom)
        shifts = as_tuple(e.chemical_shift)

        ss.setAtoms(atoms)
        ss.setChemicalShifts(shifts)

        return ss

SpinSystem._xml_state = SpinSystemXMLPickler()._xml_state
ShiftAssignment._xml_state = ShiftAssignmentXMLPickler()._xml_state



