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
from aria.Settings import Settings, PositiveFloat, ChoiceEntity, PeakType, YesNoChoice
from aria.xmlutils import XMLElement, XMLBasePickler

import aria.TypeChecking as TCheck

class DistanceCutoffEntity(PositiveFloat):

    def __init__(self):

        descr = ''
        err_msg = 'Distance-cutoff must be positive float or None'

        PositiveFloat.__init__(self, description = descr,
                               error_message = err_msg)

    def is_valid(self, value):

        v = PositiveFloat.is_valid(self, value)
        return v or (value is None)

class CalibratorSettings(Settings):

    def create(self):

        CE = ChoiceEntity

        keywords = {'distance_cutoff': DistanceCutoffEntity(),
                    'volume_or_intensity': PeakType(),
                    'estimator': CE(('ratio_of_averages',)),
                    'relaxation_matrix' : YesNoChoice(),
                    'error_estimator' : CE(('distance', 'intensity',))} # Malliavin/Bardiaux rMat
        
        return keywords

    def create_default_values(self):

        defaults = {}
        defaults['distance_cutoff'] = 6.0
        defaults['estimator'] = 'ratio_of_averages'

        return defaults

class Calibrator(AriaBaseClass):
    
    def __init__(self, settings, model):

        TCheck.check_type(settings, 'CalibratorSettings')
        TCheck.check_type(model, 'NOEModel')

        AriaBaseClass.__init__(self)

        self.setSettings(settings)
        self.setModel(model)

    def setModel(self, model):

        TCheck.check_type(model, 'NOEModel')

        self.__model = model

    def getModel(self):
        return self.__model

    def calculateEstimator(self, peaks, ensemble, store_analysis = 0,
                           use_cutoff = 1):
        """
        This method should be called if you want the Calibrator to
        calculate the calibration factor for the data-set 'spectrum'.
        'peaks': list of AriaPeaks. Every restaint in that list
        is only used to query the experimental peak-volume/intensity
        and (if 'store_analysis' is non-zero) to store the calculated
        peaksize in its analysis section.
        
        'store_analysis': if non-zero, calculated volumes are stored
        in the 'analysis' section of every AriaPeak.
        
        if 'use_cutoff' is zero, all distances are used to calculate
        the estimator.
        """

        TCheck.check_type(peaks, TCheck.LIST, TCheck.TUPLE)

        import numpy

        if not len(peaks):
            self.error(ValueError, 'No peaks specified.')
        
        #TCheck.check_elements(peaks, 'AriaPeak')
        TCheck.check_elements(peaks, 'AbstractPeak')
        TCheck.check_type(ensemble, 'StructureEnsemble')

        from aria.Datum import Datum
        
        settings = self.getSettings()

        if settings['volume_or_intensity'] == 'volume':
            exp_peak_sizes = [p.getReferencePeak().getVolume()[0] \
                              for p in peaks]
        else:
            exp_peak_sizes = [p.getReferencePeak().getIntensity()[0] \
                              for p in peaks]

        f = self.getModel().calculatePeaksize

        model_peak_sizes = numpy.array([f(p, ensemble) for p in peaks])

        ## consider only crosspeaks with volume/intensity
        ## larger than NOE_cutoff.
        
        if use_cutoff:
            NOE_cutoff = settings['distance_cutoff'] ** (-6.)
        else:
            NOE_cutoff = 0.

        strong_NOEs = numpy.greater_equal(model_peak_sizes, NOE_cutoff)

        sum_noe_model = numpy.sum(numpy.compress(strong_NOEs,
                                                     model_peak_sizes))
        sum_noe_exp = numpy.sum(numpy.compress(strong_NOEs,
                                                   exp_peak_sizes))

        ## if there are no NOEs larger than NOE_cutoff,
        ## return None.

        if sum_noe_model <= 1.e-30:
            return None

        ## calculate estimator
        if settings['estimator'] == 'ratio_of_averages':
            factor = sum_noe_exp / sum_noe_model

        ## store calculated peak-size

        if store_analysis:
            calculated_peak_sizes = model_peak_sizes * factor
            
            for i in range(len(peaks)):
                d = Datum(calculated_peak_sizes[i], None)
                peaks[i].analysis.setCalculatedPeaksize(d)

        return factor

class CalibratorXMLPickler(XMLBasePickler):

    def _xml_state(self, x):
        e = XMLElement()

        e.distance_cutoff = x['distance_cutoff']
        e.estimator = x['estimator']
        
        # Malliavin/Bardiaux rMat
        e.relaxation_matrix = x['relaxation_matrix']
        e.error_estimator = x['error_estimator']
        
        return e

    def load_from_element(self, e):
        s = CalibratorSettings()

        s['distance_cutoff'] = float(e.distance_cutoff)
        s['estimator'] = str(e.estimator)
        
        # Malliavin/Bardiaux rMat
        if hasattr(e, 'relaxation_matrix'):
            s['relaxation_matrix'] = str(e.relaxation_matrix)
        else:
            s['relaxation_matrix'] = 'no'

        if hasattr(e, 'error_estimator'):
            s['error_estimator'] = str(e.error_estimator)
        else:
            s['error_estimator'] = 'distance'        
        #
        
        return s

CalibratorSettings._xml_state = CalibratorXMLPickler()._xml_state

        
