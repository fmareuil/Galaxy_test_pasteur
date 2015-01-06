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
from aria.xmlutils import XMLElement, XMLBasePickler

import aria.Settings as Settings
import aria.TypeChecking as TCheck

class WeightCutoff(Settings.Weight):

    def is_valid(self, value):
        return value is None or Settings.Weight.is_valid(self, value)

class ContributionAssignerSettings(Settings.Settings):

    def create(self):

        return {'weight_cutoff': WeightCutoff(),
                'max_contributions': Settings.PositiveInteger()}

    def create_default_values(self):

        return {'weight_cutoff': 1.0,
                'max_contributions': 20}

class ContributionAssigner(AriaBaseClass):

    def __init__(self, settings):

        TCheck.check_type(settings, 'ContributionAssignerSettings')

        from aria.Contribution import ContributionEvaluator
        
        AriaBaseClass.__init__(self)

        self.setSettings(settings)
        self.__evaluator = ContributionEvaluator()

    def filter_weights(self, weights, cutoff, max_n):
        """
        Let I be the index-list of weights whose
        sum is >= cutoff. The function returns indices
        in range(0, len(weights)) which are not in I
        """

        TCheck.check_array(weights)
        TCheck.check_float(cutoff)
        TCheck.check_int(max_n)

        import numpy

        ## sort weights in descending order

        indices = numpy.argsort(weights)
        indices = numpy.take(indices, numpy.arange(len(indices)-1,-1,-1))
        s_weights = numpy.take(weights, indices)

        x = numpy.add.accumulate(s_weights)
        
        try:
            index = numpy.flatnonzero(numpy.greater(x, cutoff))[1]
        except:
            index = len(indices)

        ## we limit the number of contributing
        ## weights to max_n.

        ## BARDIAUX
        # test maxn remove peak
        #index = min(index, max_n)

        ## Return set of large and small weights.
            
        return indices[:index], indices[index:]

    def setDefaultAssignments(self, restraint_list):
        
        """
        standard procedure of assigning weights to contributions.
        w_i 1. / n_contribs
        """

        for restraint in restraint_list:
            contribs = restraint.getContributions()
            w = 1. / len(contribs)
            
            [c.setWeight(w) for c in contribs]

    def assign(self, restraint_list, ensemble, filter_contributions = 1):

        TCheck.check_type(restraint_list, TCheck.LIST, TCheck.TUPLE)
        #TCheck.check_elements(restraint_list, 'AriaPeak')
        TCheck.check_elements(restraint_list, 'AbstractPeak')
        TCheck.check_type(ensemble, 'StructureEnsemble', TCheck.NONE)

        from aria.mathutils import standardDeviation
        from aria.mathutils import _average as average
        from aria.Datum import Datum

        import numpy

        if ensemble is None:
            self.setDefaultAssignments(restraint_list)
            return

        self.__evaluator.setStructureEnsemble(ensemble)

        all_contributions = []
        weights = []

        for restraint in restraint_list:

            distances = []

            contributions = restraint.getContributions()
            all_contributions.append(contributions)

            for contribution in contributions:

                ## for every structure: get effective distance
                ## for 'contribution'

                d = self.__evaluator.effective_distances(contribution)

                d_avg = average(d)
                distances.append(d_avg)
            
                if len(d) > 1:
                    sd = standardDeviation(d, avg = d_avg)
                else:
                    sd = None

                contribution.setAverageDistance(Datum(d_avg, sd))

            ## calculate partial NOE wrt to ensemble-averaged
            ## distance. The partial NOE serves as weight which
            ## subsequently will be normalized to 1.

            w = numpy.power(distances, -6.)
            
            ## normalize weights and store weights

            w /= numpy.sum(w)

            weights.append(w)
            
        settings = self.getSettings()
        cutoff = settings['weight_cutoff']
    
        #if cutoff is not None:
        if cutoff is not None and filter_contributions:

            ## 1. disable all contributions according to
            ##    the partial-assignment scheme.
            ## 2. allow at most 'max_contributions' contributions

            max_n = settings['max_contributions']
            
            for i in range(len(weights)):

                w = weights[i]
                c = all_contributions[i]
                
                on, off = self.filter_weights(w, cutoff, max_n) 

                ## active / deactive contributions
            
                [c[index].setWeight(0.) for index in off]
                [c[index].setWeight(w[index]) for index in on]
            
        else:

            ## if cutoff is not set, enable all contribution.
            
            ## note: setting 'max_contributions' does not
            ## apply in that case since we have no rule 
            ## how to select contributions which remain
            ## active.

            for i in range(len(weights)):
                
                contributions = all_contributions[i]
                w = weights[i]

                for j in range(len(w)):
                    contributions[j].setWeight(w[j])
                    
    def assign_before(self, restraint_list, ensemble):

        TCheck.check_type(restraint_list, TCheck.LIST, TCheck.TUPLE)
        TCheck.check_elements(restraint_list, 'AriaPeak')
        TCheck.check_type(ensemble, 'StructureEnsemble', TCheck.NONE)

        from aria.mathutils import standardDeviation
        from aria.mathutils import _average as average
        from aria.Datum import Datum

        import numpy

        if ensemble is None:
            self.setDefaultAssignments(restraint_list)
            return

        self.__evaluator.setStructureEnsemble(ensemble)

        all_contributions = []
        weights = []

        for restraint in restraint_list:

            distances = []

            contributions = restraint.getContributions()
            all_contributions.append(contributions)

            for contribution in contributions:

                ## for every structure: get effective distance
                ## for 'contribution'

                d = self.__evaluator.effective_distances(contribution)

                d_avg = average(d)
                distances.append(d_avg)
            
                if len(d) > 1:
                    sd = standardDeviation(d, avg = d_avg)
                else:
                    sd = None

                contribution.setAverageDistance(Datum(d_avg, sd))

            ## calculate partial NOE wrt to ensemble-averaged
            ## distance. The partial NOE serves as weight which
            ## subsequently will be normalized to 1.

            w = numpy.power(distances, -6.)
            
            ## normalize weights and store weights

            w /= numpy.sum(w)

            weights.append(w)
            
        settings = self.getSettings()
        cutoff = settings['weight_cutoff']


        ## if cutoff is not set, enable all contribution.
        
        ## note: setting 'max_contributions' does not
        ## apply in that case since we have no rule 
        ## how to select contributions which remain
        ## active.

        for i in range(len(weights)):
                
            contributions = all_contributions[i]
            w = weights[i]

            for j in range(len(w)):
                contributions[j].setWeight(w[j])
                

    def filter_contributions(self, restraint_list):

        TCheck.check_type(restraint_list, TCheck.LIST, TCheck.TUPLE)
        #TCheck.check_elements(restraint_list, 'AriaPeak')
        TCheck.check_elements(restraint_list, 'AbstractPeak')

        from aria.mathutils import standardDeviation
        from aria.mathutils import _average as average
        from aria.Datum import Datum

        import numpy

        all_contributions = []
        weights = []

        for restraint in restraint_list:

            distances = []

            contributions = restraint.getContributions()
            all_contributions.append(contributions)

            w = numpy.array([c.getWeight() for c in contributions])

            weights.append(w)
            
        settings = self.getSettings()
        cutoff = settings['weight_cutoff']
    
        if cutoff is not None:

            ## 1. disable all contributions according to
            ##    the partial-assignment scheme.
            ## 2. allow at most 'max_contributions' contributions

            max_n = settings['max_contributions']
            
            for i in range(len(weights)):

                w = weights[i]
                c = all_contributions[i]
                
                on, off = self.filter_weights(w, cutoff, max_n) 

                ## active / deactive contributions
            
                [c[index].setWeight(0.) for index in off]
                [c[index].setWeight(w[index]) for index in on]
            
        else:

            ## if cutoff is not set, enable all contribution.
            
            ## note: setting 'max_contributions' does not
            ## apply in that case since we have no rule 
            ## how to select contributions which remain
            ## active.

            for i in range(len(weights)):
                
                contributions = all_contributions[i]
                w = weights[i]

                for j in range(len(w)):
                    contributions[j].setWeight(w[j])

class ContributionAssignerXMLPickler(XMLBasePickler):

    order = ['weight_threshold', 'max_contributions']

    def _xml_state(self, x):
        e = XMLElement(tag_order = self.order)

        e.weight_threshold = x['weight_cutoff']
        e.max_contributions = x['max_contributions']
        
        return e

    def load_from_element(self, e):
        s = ContributionAssignerSettings()

        s['weight_cutoff'] = float(e.weight_threshold)
        s['max_contributions'] = int(e.max_contributions)
        
        return s

ContributionAssignerSettings._xml_state = ContributionAssignerXMLPickler()._xml_state
