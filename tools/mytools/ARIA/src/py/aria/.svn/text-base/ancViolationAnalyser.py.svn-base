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
from numpy import *
from aria.xmlutils import XMLElement, XMLBasePickler

class ViolationAnalyserSettings(Settings):

    def create(self):

        from aria.Settings import NonNegativeFloat, Weight, TypeEntity

        descr = 'Violation threshold.'
        err_msg = 'Violation threshold must be between 0. and 1.'

        vt = Weight(description = descr, error_message = err_msg)

        d = {'violation_tolerance': NonNegativeFloat(),
             'violation_threshold': vt,
             'lower_bound_correction': TypeEntity('LowerBoundCorrection'),
             'upper_bound_correction': TypeEntity('UpperBoundCorrection')}

        return d

    def create_default_values(self):

        import aria.DataContainer as DC

        d = {}
        d['violation_tolerance'] = 1.0
        d['violation_threshold'] = 0.5

        l_c = DC.LowerBoundCorrection()
        u_c = DC.UpperBoundCorrection()
        
        l_c.reset()
        u_c.reset()
        
        d['lower_bound_correction'] = l_c
        d['upper_bound_correction'] = u_c

        return d

class ViolationAnalyser(AriaBaseClass):

    def __init__(self, settings):

        check_type(settings, 'ViolationAnalyserSettings')

        from aria.Contribution import ContributionEvaluator

        AriaBaseClass.__init__(self)
        self.setSettings(settings)

        self.__evaluator = ContributionEvaluator()

        ## set callbacks to be notified when settings
        ## have changed

        s = self.getSettings()
        e = s.getEntity('violation_tolerance')
        e.set_callback(self.settings_changed)

        try:
            self.__tolerance = s['violation_tolerance']
        except:
            pass

    def settings_changed(self, entity):
        self.__tolerance = self.getSettings()['violation_tolerance']

    def get_tolerance(self):
        return self.__tolerance

##     def check(self, d_calc, peak):

##         ## confidence level

##         L = 0.5

##         length, a, b = confidenceInterval(d_calc, L)

##         tol = self.getSettings()['violation_tolerance']

##         l = peak.getLowerBound() - tol
##         u = peak.getUpperBound() + tol

##         if l >= a:
##             if u < b:
##                 p = 1.
##             else:
##                 p = (b - l) * L + (u - b) * (1. - L)

##         else:
##             if u < a:
##                 p = 0.
##             else:
##                 p = (l - a) * (1. - L) + (u - a) * L

##         return p

    def analysePeak(self, ensemble, peak, store_analysis = 0,
                    lower_correction = None, upper_correction = None):

        """
        """

        check_type(ensemble, 'StructureEnsemble')
        check_type(peak, 'AriaPeak')
        check_int(store_analysis)
        check_type(lower_correction, FLOAT, NONE)
        check_type(upper_correction, FLOAT, NONE)

        from aria.mathutils import _average as average

        self.__evaluator.setStructureEnsemble(ensemble)

        ## for every structure: calculate effective contributon-distance
        ## d_avg is a [n_c x n_s] dim. array
        ## n_c: number of contributions
        ## n_s: number of structures in ensemble

        d_avg = [self.__evaluator.effective_distances(c) \
                 for c in peak.getContributions()]

        ## For each structure: Calculate effective distance,
        ## from theoretical NOE.:
        ##
        ## d_eff = (\sum_i d_i^{-6})^{-1/6}
        ##
        ## d_i: distance of atoms comprising i-th contribution
        ##
        ## d_avg has length #structures.
        
        d_avg = power(sum(power(d_avg, -6), axis = 0), -1./6)

        ## Effective lower/upper bounds

        tol = self.get_tolerance()

        if lower_correction is not None:
            lower = lower_correction
        else:
            lower = peak.getLowerBound()

        if upper_correction is not None:
            upper = upper_correction
        else:
            upper = peak.getUpperBound()

        ## calculate fraction of violated distances
        ## 1: distance is violated, 0: distance lies within bounds

        violated_lower = less(d_avg, lower - tol)
        violated_upper = greater(d_avg, upper + tol)
        
        violated = logical_or(violated_lower, violated_upper)
        
        n_lower = sum(violated_lower)
        n_upper = sum(violated_upper)

        R_viol = float(sum(violated)) / float(len(violated))

        ## If we shall store some itermediate results,
        ## i.e. if 'result' is non-None, make it so.
        
        if store_analysis:

            from aria.mathutils import standardDeviation
            from aria.Datum import Datum

            ## TODO: what else can we store? shall we use median
            ## instead of average?

            result = peak.analysis

            ## Average distance

            if len(d_avg) > 1:
                sd = standardDeviation(d_avg)
            else:
                sd = None

            result.setAverageDistance(Datum(average(d_avg), sd))

            ## We calculating bound-violations, we do not consider
            ## the violation tolerance! The violation tolerance is
            ## used only to evaluate when a restraint is violated.
            ## In other words, lower/upper bound violations relect
            ## the structures 'as they are' and do not depend on
            ## parameter settings of a particular iteration.

            violated_lower = less(d_avg, lower)
            violated_upper = greater(d_avg, upper)

            n_lower = sum(violated_lower)
            n_upper = sum(violated_upper)

            ## Lower bound violation

            if n_lower:
                d_viol = compress(violated_lower, d_avg)
                d = lower - average(d_viol)

                if n_lower > 1:
                    sd = standardDeviation(d_viol)
                else:
                    sd = None
            else:
                d = sd = 0.
                
            result.setLowerBoundViolation(Datum(d, sd))

            ## Upper bound violation

            if n_upper:
                d_viol = compress(violated_upper, d_avg)
                d = average(d_viol) - upper
                
                if n_upper > 1:
                    sd = standardDeviation(d_viol)
                else:
                    sd = None
            else:
                d = sd = 0.

            result.setUpperBoundViolation(Datum(d, sd))

            ## Violation percentage
            
            result.setDegreeOfViolation(R_viol)

        return R_viol

class ViolationAnalyserXMLPickler(XMLBasePickler):

    def _xml_state(self, x):
        e = XMLElement()

        e.violation_tolerance = x['violation_tolerance']
        e.violation_threshold = x['violation_threshold']
        
        return e

    def load_from_element(self, e):
        s = ViolationAnalyserSettings()

        s['violation_tolerance'] = float(e.violation_tolerance)
        s['violation_threshold'] = float(e.violation_threshold)
        
        return s

ViolationAnalyserSettings._xml_state = ViolationAnalyserXMLPickler()._xml_state
