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
from aria.xmlutils import XMLBasePickler

class NOESYSpectrum(AriaBaseClass):
    """
    Container class that holds a list of 'NOE's.
    """

    id = 0
    
    def __init__(self, name, noes, stype = "intra"):

        check_type(name, STRING, NONE)
        check_type(noes, LIST, TUPLE)

        check_elements(noes, 'CrossPeak')

        AriaBaseClass.__init__(self)

        self.__name = name
        # BARDIAUX 2.2
        self.__type = stype
        
        self.__noes = {}
        
        for noe in noes:
            self.__noes[noe.getNumber()] = noe
            noe.setSpectrum(self)

    def addPeak(self, peak):

        check_type(peak, 'CrossPeak')

        num = peak.getNumber()

        if num in self.__noes:

            self.error('Spectrum %s: Crosspeak with number %d already exists.' % (self.getName(), num))

        self.__noes[num] = peak
        peak.setSpectrum(self)
            
    def getPeaks(self):
        peaks = self.__noes.values()
        peaks.sort(lambda a, b: cmp(a.getNumber(), b.getNumber()))
        
        return peaks

    def findPeak(self, number):

        if self.__noes.has_key(number):
            return self.__noes[number]
        else:
            return None

    def getName(self):
        return self.__name

    def _setName(self, n):
        check_string(n)
        self.__name = n

    # BARDIAUX 2.2 type
    def getType(self):
        return self.__type

    def setType(self, t):
        check_string(t)
        self.__type = t

    def setDataSource(self, source):

        check_type(source, 'PeakData')
        self.__data_source = source

    def getDataSource(self):

        return self.__data_source
    
    # Malliavin/Bardiaux rMat
    def getExperimentData(self):
        return self.__experiment_data

    def setExperimentData(self, data):
        check_type(data, 'ExperimentData')        
        self.__experiment_data = data

    #
    def __getitem__(self, number):

        check_int(number)

        if not self.__noes.has_key(number):
            raise KeyError, "Peak no. %d not in spectrum." %number

        return self.__noes[number]
        
    def __len__(self):
        return len(self.__noes)

    def __str__(self):

        class_name = self.__class__.__name__

        return '%s(name=%s, n_peaks=%s)' % (class_name,
                                            self.getName(),
                                            str(len(self)))

    __repr__ = __str__

# BARDIAUX
class ConstraintList(NOESYSpectrum):

    def __init__(self, name, noes):

        NOESYSpectrum.__init__(self, name, noes)

    def setListSource(self, source):

        check_type(source, 'UnambiguousDistanceData', 'AmbiguousDistanceData')
        self.__list_source = source

    def getListSource(self):

        return self.__list_source

class NOESYSpectrumXMLPickler(XMLBasePickler):

    order = ('peak', 'name',)#, 'type', )## BARDIAUX 2.2 : type is intra|inter|all

    def __init__(self):

        from aria.Datum import DatumXMLPickler, ChemicalShiftXMLPickler
        from aria.Atom import AtomXMLPickler
        from aria.CrossPeak import CrossPeakXMLPickler
        from aria.Assignment import AssignmentXMLPickler
        
        self.spectrum = self
        self.volume = DatumXMLPickler()
        self.intensity = DatumXMLPickler()
        self.shift = ChemicalShiftXMLPickler()
        self.atom = AtomXMLPickler()
        self.peak = CrossPeakXMLPickler()
        self.assignment = AssignmentXMLPickler()
        
    def _xml_state(self, spectrum):

        from aria.xmlutils import XMLElement

        e = XMLElement(tag_order = self.order)
        e.name = spectrum.getName()
        #e.type = spectrum.getType()
        e.peak = spectrum.getPeaks()
        e.peak.sort(lambda a, b: cmp(a.getNumber(), b.getNumber()))
        
        return e

    def load_from_element(self, e):

        p = e.peak
        if not len(p):
            self.error(ValueError, 'Empty peak list.')

        try:
            name = eval(e.name)
        except:
            name = str(e.name)

        if hasattr(e, 'type'):
            stype = str(e.type)
        else:
            stype = "intra"

        s = NOESYSpectrum(name, e.peak, stype)
        

        return s


NOESYSpectrum._xml_state = NOESYSpectrumXMLPickler()._xml_state

    
