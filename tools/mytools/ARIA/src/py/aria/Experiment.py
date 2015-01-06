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

class Experiment(AriaBaseClass):

    def __init__(self, spectrum, shift_assignments):

        check_type(spectrum, 'NOESYSpectrum')
        check_type(shift_assignments, 'ChemicalShiftList')

        self.__spectrum = spectrum
        self.__shift_list = shift_assignments
        self.__name = spectrum.getName()

        self.setDefaultValues()

    def setDefaultValues(self):
        
        self.__data_source = None
        
    def getSpectrum(self):
        return self.__spectrum

    def getShiftList(self):
        return self.__shift_list

    def getName(self):
        return self.__name

    def setDataSource(self, ds):
        check_type(ds, 'SpectrumData')
        self.getSpectrum().setDataSource(ds['peaks'])
        self.__data_source = ds

        # bardiaux rMat
        self.getSpectrum().setExperimentData(ds['experiment_data'])

    def getDataSource(self):
        return self.__data_source

    def getFilteredShiftAssignments(self):
        shifts = self.getShiftList().getShiftAssignments()
        return [a for a in shifts if a.is_valid()]

    def getFilteredPeaks(self):
        return [p for p in self.getSpectrum().getPeaks() if p.is_valid()]

    def __str__(self):

        class_name = self.__class__.__name__

        return '%s(name=%s)' % (class_name, self.getName())

    __repr__ = __str__

