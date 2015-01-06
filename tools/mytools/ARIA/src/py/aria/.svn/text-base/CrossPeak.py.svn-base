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
from aria.xmlutils import XMLBasePickler
from aria.tools import as_tuple

import aria.TypeChecking as TCheck
import aria.Assignment as Assignment

CROSSPEAK_ASSIGNMENT_TYPE_MANUAL = Assignment.ASSIGNMENT_TYPE_MANUAL
CROSSPEAK_ASSIGNMENT_TYPE_AUTOMATIC = Assignment.ASSIGNMENT_TYPE_AUTOMATIC
CROSSPEAK_ASSIGNMENT_TYPE_SEMIAUTOMATIC = 'semiautomatic'

class CrossPeak(AriaBaseClass):

    def __init__(self, number, volume, intensity):

        """
        Initialization via:
        - peak number
        - spectrum size, i.e. a tuple that consists of the peak's
          intensity / volume and its uncertainty
        """

        TCheck.check_int(number)
        TCheck.check_type(volume, 'Datum', TCheck.NONE)
        TCheck.check_type(intensity, 'Datum', TCheck.NONE)

        from aria.Datum import Datum

        if volume is None and intensity is None:
            self.error(ValueError, 'Either volume or intensity must ' + \
                       'be of type Datum')
        
        if number < 0:
            self.error(ValueError, 'Number must be >= 0')

        if volume is None:
            volume = Datum(None)

        if intensity is None:
            intensity = Datum(None)
            
        AriaBaseClass.__init__(self)

        self.__number = number
        self.__intensity = intensity
        self.__volume = volume
        self.__reliable = 0
        self.__updated = 0
        self.__ambiguity = None
        
        self.setDefaultValues()

    def setDefaultValues(self):

        from aria.Datum import ChemicalShift

        self.__proton1ChemicalShift = ChemicalShift(None)
        self.__proton2ChemicalShift = ChemicalShift(None)

        self.__hetero1ChemicalShift = ChemicalShift(None)
        self.__hetero2ChemicalShift = ChemicalShift(None)

        self.__spectrum = None

        self.__id = None

        self.__proton1_assignments = []
        self.__hetero1_assignments = []
        self.__proton2_assignments = []
        self.__hetero2_assignments = []

        self.__valid = 1

    def isReliable(self, e = None):
        """
        if argument is None:

        returns 1 if cross-peak will be used in calculation under
        any circumstances, 0 otherwise.

        if argumente is int (zero or one): crosspeaks 
        """

        if e is None:
            return self.__reliable

        TCheck.check_bool(e)

        self.__reliable = e > 0

    def isAssigned(self):
        """
        Returns 1 if the proton-dimensions have assignments
        """
        if self.getProton1Assignments() and self.getProton2Assignments():
            return 1
        else:
            return 0

    def getAssignmentType(self):

        types = [a.getType() for a in self.getProton1Assignments() + \
                 self.getProton2Assignments() + self.getHetero1Assignments()+ \
                 self.getHetero2Assignments()]

        if types.count(Assignment.ASSIGNMENT_TYPE_MANUAL) == len(types):

            return CROSSPEAK_ASSIGNMENT_TYPE_MANUAL

        elif types.count(Assignment.ASSIGNMENT_TYPE_AUTOMATIC) == len(types):

            return CROSSPEAK_ASSIGNMENT_TYPE_AUTOMATIC

        else:
            return CROSSPEAK_ASSIGNMENT_TYPE_SEMIAUTOMATIC

    def addProton1Assignment(self, a):

        TCheck.check_type(a, 'Assignment')
        self.__proton1_assignments.append(a)

    ## BARDIAUX 15/04/05
    def setProton1Assignments(self, assignments):

        TCheck.check_type(assignments, TCheck.LIST, TCheck.TUPLE)
        TCheck.check_elements(assignments, 'Assignment')

        self.__proton1_assignments = list(assignments)

    def setProton2Assignments(self, assignments):

        TCheck.check_type(assignments, TCheck.LIST, TCheck.TUPLE)
        TCheck.check_elements(assignments, 'Assignment')

        self.__proton2_assignments = list(assignments)


    def setHetero1Assignments(self, assignments):

        TCheck.check_type(assignments, TCheck.LIST, TCheck.TUPLE)
        TCheck.check_elements(assignments, 'Assignment')

        self.__hetero1_assignments = list(assignments)


    def setHetero2Assignments(self, assignments):

        TCheck.check_type(assignments, TCheck.LIST, TCheck.TUPLE)
        TCheck.check_elements(assignments, 'Assignment')

        self.__hetero2_assignments = list(assignments)

    ##
        

    def addHetero1Assignment(self, a):

        TCheck.check_type(a, 'Assignment')
        self.__hetero1_assignments.append(a)

    def addProton2Assignment(self, a):

        TCheck.check_type(a, 'Assignment')
        self.__proton2_assignments.append(a)

    def addHetero2Assignment(self, a):

        TCheck.check_type(a, 'Assignment')
        self.__hetero2_assignments.append(a)

    def getProton1Assignments(self):
        return self.__proton1_assignments

    def getHetero1Assignments(self):
        return self.__hetero1_assignments

    def getProton2Assignments(self):
        return self.__proton2_assignments

    def getHetero2Assignments(self):
        return self.__hetero2_assignments

    def setId(self, id):
        TCheck.check_int(id)
        self.__id = id

    def getId(self):
        return self.__id

    def setSpectrum(self, s):
        TCheck.check_type(s, 'NOESYSpectrum')
        self.__spectrum = s

    def getSpectrum(self):
        return self.__spectrum

    def getNumber(self):
        """
        Returns the peak's number.
        """
        return self.__number

    def getIntensity(self):
        """
        Returns the peak's intensity and its uncertainty.
        """
        return self.__intensity

    def getVolume(self):
        """
        Returns the peak's volume and its uncertainty.
        """
        return self.__volume

    def getAmbiguity(self):
        return self.__ambiguity

    def setAmbiguity(self, amb):
        self.__ambiguity = amb

    ## TODO: why set-methods for chemical shifts?
    ## presumably they will be set only once.

    def setProton1ChemicalShift(self, v):
        TCheck.check_type(v, 'ChemicalShift')
        self.__proton1ChemicalShift = v

    def getProton1ChemicalShift(self):
        return self.__proton1ChemicalShift

    def setHetero1ChemicalShift(self, v):
        TCheck.check_type(v, 'ChemicalShift')
        self.__hetero1ChemicalShift = v

    def getHetero1ChemicalShift(self):
        return self.__hetero1ChemicalShift

    def setProton2ChemicalShift(self, v):
        TCheck.check_type(v, 'ChemicalShift')
        self.__proton2ChemicalShift = v

    def getProton2ChemicalShift(self):
        return self.__proton2ChemicalShift

    def setHetero2ChemicalShift(self, v):
        TCheck.check_type(v, 'ChemicalShift')
        self.__hetero2ChemicalShift = v

    def getHetero2ChemicalShift(self):
        return self.__hetero2ChemicalShift

    def is_valid(self, valid = None):

        TCheck.check_type(valid, TCheck.NONE, TCheck.INT)
        
        if valid is not None:
            self.__valid = valid
        else:
            return self.__valid == 1
        
    def getDimensions(self, spSys1, spSys2):
        """
        Return dimensions (1 or 2) of input
        spinsystems
        """

        proton1, proton2 = spSys1.getAtoms()[0], spSys2.getAtoms()[0]

        assignments1 = [atom for assi in self.getProton1Assignments()
                        for atom in assi.getAtoms()]

        assignments2 = [atom for assi in self.getProton2Assignments()
                        for atom in assi.getAtoms()]

        dimensions1 = [proton1 in assignments1, proton1 in assignments2]
        dimensions2 = [proton2 in assignments1, proton2 in assignments2]    

        dims = [None, None]

        if dimensions1.count(1) == 1 and dimensions2.count(1) == 1:
            dims[0] = dimensions1.index(1) + 1
            dims[1] = dimensions2.index(1) + 1
            

        elif dimensions1.count(1) == 1 and dimensions2.count(1) == 2:

            dims[0] = dimensions1.index(1) + 1
            dims[1] = 1
            if dims[0] == 1:
                dims[1] = 2

        elif dimensions1.count(1) == 2 and dimensions2.count(1) == 1:
            dims[1] = dimensions2.index(1) + 1
            dims[0] = 1
            if dims[1] == 1:
                dims[0] = 2

        elif dimensions1.count(1) == 2 and dimensions2.count(1) == 2:
            dims = [1,2]
            h1 = proton1.getHeteroAtom()
            h2 = proton2.getHeteroAtom()
            
            h1assi = [atom for assi in self.getHetero1Assignments()
                      for atom in assi.getAtoms()]
            
            h2assi = [atom for assi in self.getHetero2Assignments()
                      for atom in assi.getAtoms()]

            if h1assi:
                if h2 in h1assi:
                    dims = [2, 1]
                    
            if h2assi:
                if h1 in h2assi:
                    dims = [2, 1]
            
        if None in dims:
            raise ValueError, 'spin pair could not be assigned to dimensions'

        else:
            return tuple(dims)
    
    ## BARDIAUX 25/04/05
    def getDimension(self,spinsystem):
        """
        Input : spinsystem
        Return : Dimension of the spinsytem in this Xpk
        """

        from numpy import array, sum

        p1chem = self.getProton1ChemicalShift().getValue()
        p2chem = self.getProton2ChemicalShift().getValue()
        cpchem = [p1chem,p2chem]

        sschem = [a.getValue() for a in spinsystem.getChemicalShifts() \
                  if a.getValue() is not None]

        d1 = [abs(s-p1chem) for s in sschem]
        d2 = [abs(s-p2chem) for s in sschem]
        d1.sort()
        d2.sort()

        if d1[0] < d2[0]:
            dim = 1
        else:
            dim = 2
            
##         sschem = sum(array(sschem)) / len(sschem)

##         if abs((sschem - p1chem)) <  abs((sschem - p2chem)):
##             dim = 1
##         else:
##             dim = 2

        return dim
    
    def __str__(self):
        
        p1a = self.getProton1Assignments()
        h1a = self.getHetero1Assignments()
        p2a = self.getProton2Assignments()
        h2a = self.getHetero2Assignments()

        string = 'CrossPeak(id=%s, number=%d, volume=%s, intensity=%s, ' + \
                 'reliable=%s, proton1ppm=%s, hetero1ppm=%s, ' + \
                 'proton2ppm=%s, hetero2ppm=%s, ' + \
                 'proton1assignments=%s, hetero1assignments=%s ' + \
                 'proton2assignments=%s, hetero2assignments=%s)'

        return string % (str(self.getId()),
                         self.getNumber(),
                         str(self.getVolume()), 
                         str(self.getIntensity()),
                         str(self.isReliable()),
                         str(self.getProton1ChemicalShift()),
                         str(self.getHetero1ChemicalShift()),
                         str(self.getProton2ChemicalShift()),
                         str(self.getHetero2ChemicalShift()),
                         str(p1a), str(h1a), str(p2a), str(h2a))

    __repr__ = __str__

class CrossPeakXMLPickler(XMLBasePickler):

    order = ('number', 'reliable', 'ambiguity', 'volume', 'intensity', 'proton1',
             'proton2','hetero1', 'hetero2')

    def __init__(self):
        XMLBasePickler.__init__(self)
        self._set_name('Cross peak')

    def _xml_state(self, peak):

        from aria.xmlutils import XMLElement
        from aria.Atom import AtomXMLPickler
        from aria.ariabase import YES, NO

        atom_xml_state = AtomXMLPickler()._xml_state
        
        e = XMLElement(tag_order = self.order)
        e.number = peak.getNumber()

        if peak.isReliable():
            e.reliable = YES
        else:
            e.reliable = NO
            
        e.volume = peak.getVolume()
        e.intensity = peak.getIntensity()

        e.ambiguity = peak.getAmbiguity()

        d = {'proton1': {'shift': peak.getProton1ChemicalShift,
                         'assignment': peak.getProton1Assignments},
             'proton2': {'shift': peak.getProton2ChemicalShift,
                         'assignment': peak.getProton2Assignments},
             'hetero1': {'shift': peak.getHetero1ChemicalShift,
                         'assignment': peak.getHetero1Assignments},
             'hetero2': {'shift': peak.getHetero2ChemicalShift,
                         'assignment': peak.getHetero2Assignments}}

        for dim, attr in d.items():

            dd = XMLElement()
            tag_order = ['shift']

            shift = attr['shift']()
            dd.shift = shift
                
            assignments = attr['assignment']()
            if assignments:
                dd.assignment = assignments

            setattr(e, dim, dd)
            
        return e

    def load_from_element(self, e):

        from aria.Datum import Datum
        from aria.ariabase import YES, NO
        
        number = int(e.number)
        
        value, error = e.volume
        
        if value is not None and value < 0.:
            value = abs(value)
            
        volume = Datum(value, error)

        value, error = e.intensity
        
        if value is not None and value < 0.:
            value = abs(value)
            
        intensity = Datum(value, error)

        p = CrossPeak(number, volume, intensity)
        reliable = e.reliable
        
        if hasattr(e, 'ambiguity'):
            a = str(e.ambiguity)
            if not(a): a = None
            p.setAmbiguity(a)
        else:
            p.setAmbiguity(None)

        if not reliable in (YES, NO):
            s = 'Crosspeak: attribute "reliable" must be either ' + \
                '"%s" or "%s". "%s" given.'
            self.error(ValueError, s % (str(YES), str(NO), reliable))

        p.isReliable(reliable == YES)

        d = {'proton1': {'shift': p.setProton1ChemicalShift,
                         'assignment': p.addProton1Assignment},
             'proton2': {'shift': p.setProton2ChemicalShift,
                         'assignment': p.addProton2Assignment},
             'hetero1': {'shift': p.setHetero1ChemicalShift,
                         'assignment': p.addHetero1Assignment},
             'hetero2': {'shift': p.setHetero2ChemicalShift,
                         'assignment': p.addHetero2Assignment}}

        for dim, attr in d.items():

            if not hasattr(e, dim): continue

            dd = getattr(e, dim)

            if hasattr(dd, 'shift'):
                attr['shift'](dd.shift)

            if hasattr(dd, 'assignment'):
                [attr['assignment'](a) for a in as_tuple(dd.assignment)]

        return p

CrossPeak._xml_state = CrossPeakXMLPickler()._xml_state
    
