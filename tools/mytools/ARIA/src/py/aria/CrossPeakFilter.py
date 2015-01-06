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



DIAGONAL_PEAK = 'Diagonal peak.'
NO_PEAK_SIZE = 'No peak-size (volume / intensity) given.'
NEGATIVE_PEAK = 'Negative peak-size (volume / intensity).'
UNASSIGNED_PEAK = 'Unassigned peak.'

from aria.ariabase import AriaBaseClass
from aria.Settings import Settings

class CrossPeakFilterSettings(Settings):

    def create(self):

        from aria.Settings import NonNegativeFloat, PeakType, YesNoChoice

        keywords = {'proton1_shift_err': NonNegativeFloat(),
                    'proton2_shift_err': NonNegativeFloat(),
                    'volume_or_intensity': PeakType(),
                    'filter_diagonal_peaks' : YesNoChoice(),
                    'filter_unassigned_peaks' : YesNoChoice()}

        return keywords

class CrossPeakFilter(AriaBaseClass):

    def __init__(self, settings):

        from aria.TypeChecking import check_type

        check_type(settings, 'CrossPeakFilterSettings')
        
        AriaBaseClass.__init__(self)

        self.setSettings(settings)

    def __filter_proton_shifts(self, peak):

        from aria.TypeChecking import check_type, FLOAT

        check_type(peak, 'CrossPeak')

        from aria.ChemicalShiftFilter import ChemicalShiftFilter, \
             ChemicalShiftFilterSettings

        filter = ChemicalShiftFilter(ChemicalShiftFilterSettings())
        filter.getSettings()['ppm_type'] = FLOAT

        result = {}

        valid = filter(peak.getProton1ChemicalShift())
        result.update(filter.result)

        valid &= filter(peak.getProton2ChemicalShift())
        result.update(filter.result)
            
        return valid, result

    def __filter_diagonal_peak(self, peak):
        """
        Checks whether a cross-peak is diagonal, i.e. if the
        proton chemical-shift windows overlap.
        """

        #return 1, None

        ## TODO: treatment of hetero atoms!

        from aria.TypeChecking import check_type

        check_type(peak, 'CrossPeak')

        value1, error1 = peak.getProton1ChemicalShift()
        value2, error2 = peak.getProton2ChemicalShift()

        if error1 is None:
            error1 = self.getSettings()['proton1_shift_err']

        if error2 is None:
            error2 = self.getSettings()['proton2_shift_err']

        if value1 is not None and value2 is not None and \
           abs(value1-value2) < error1 + error2:
            return 0, DIAGONAL_PEAK

        else:
            return 1, None

    ## BARDIAUX 2.2 test
    def __filter_unassigned_peak(self, peak):
        """
        Check if a peak is unassigned and remove it from the
        peak list (only for filter_unassigned)
        """

        from aria.TypeChecking import check_type

        check_type(peak, 'CrossPeak')

        is_assigned = peak.isAssigned()
        
        if not is_assigned:
            return 0, UNASSIGNED_PEAK

        else:
            return 1, None
        
    def __filter_volume_and_intensity(self, peak):

        from aria.TypeChecking import check_type

        check_type(peak, 'CrossPeak')

        peak_type = self.getSettings()['volume_or_intensity']

        if peak_type == 'volume':
            peak_size = peak.getVolume()[0]

        elif peak_type == 'intensity':
            peak_size = peak.getIntensity()[0]

        if peak_size is None:
            return 0, NO_PEAK_SIZE

        ## peak_sizes are abs'd when loading the data.

        elif peak_size <= 0.0:
            return 0, NEGATIVE_PEAK

        else:
            return 1, None
        
    def __call__(self, peak):

        from aria.TypeChecking import check_type

        check_type(peak, 'CrossPeak')

        valid = 1
        result = {}

        ## first check wether the two proton chemical-shifts are
        ## valid

        v, r = self.__filter_proton_shifts(peak)
        valid &= v
        result.update(r)
            
        ## check wether peak is a diagonal peak
        from aria.ariabase import YES, NO
        if self.getSettings()['filter_diagonal_peaks'] == YES:
            v, r = self.__filter_diagonal_peak(peak)
        else:
            v, r = 1, None

        valid &= v
        result['diagonal'] = r

        ## BARDIAUX 2.2 test
        # check if a peak is assigned
        if self.getSettings()['filter_unassigned_peaks'] == YES:
            v, r = self.__filter_unassigned_peak(peak)
        else:
            v, r = 1, None

        valid &= v
        result['unassigned'] = r
        
        ## check wether a valid value for the volume or
        ## the intensity is given

        v, r = self.__filter_volume_and_intensity(peak)
        valid &= v
        result['size'] = r

        self.result = {peak: result}

        return valid
