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

class AssignmentFilter(AriaBaseClass):

    def __init__(self):

        self.n_residues = None

    def filter_weights(self, weights, cutoff, max_n):

        check_array(weights)
        check_float(cutoff)
        check_int(max_n)

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
            
        index = min(index, max_n)

        ## Return set of large and small weights.
            
        return indices[:index], indices[index:]

    def getResidueNumbers(self, contributions):
        
        res_numbers = []

        for c in contributions:

            sp = c.getSpinPairs()[0]

            n1 = sp[0].getResidue().getNumber() - 1
            n2 = sp[1].getResidue().getNumber() - 1

            res_numbers.append((n1, n2))

        return res_numbers

    def getNResidues(self, restraints):

        n_max = 0

        for p in restraints:

            for c in p.getContributions():

                sp = c.getSpinPairs()[0]
                n1 = sp[0].getResidue().getNumber()
                n2 = sp[1].getResidue().getNumber()

                if max(n1, n2) > n_max:
                    n_max = max(n1, n2)

        return n_max

    def buildMatrix(self, restraints, n_residues, weight_matrix = None):

        import numpy

        z = numpy.zeros((self.n_residues, self.n_residues), numpy.float)

        for restraint in restraints:

            contribs = restraint.getContributions()
            n_contribs = len(contribs)

            res_numbers = self.getResidueNumbers(contribs)

            if weight_matrix is None:
                weights = numpy.ones(n_contribs, numpy.Float)
            else:
                weights = numpy.array([weight_matrix[i[0], i[1]] \
                                         for i in res_numbers])

                if numpy.sum(weights) < 1.e-10:
                    weights = numpy.ones(n_contribs, numpy.Float)

            weights = weights / numpy.sum(weights)

            for i in range(len(contribs)):

                n1, n2 = res_numbers[i]
                w = weights[i]

                z[n1, n2] += w

                if n1 <> n2:
                    z[n2, n1] += w

        return z

    def buildContactMatrix(self, restraint_list, n_iterations):

        from pystartup import Dump

        if self.n_residues is None:
            self.n_residues = self.getNResidues(restraint_list)

        contact_matrix = None

        for i in range(n_iterations):

            contact_matrix = self.buildMatrix(restraint_list, \
                                              self.n_residues, \
                                              contact_matrix)

            Dump(contact_matrix, '/tmp/cm%d' % (i+1))

        return contact_matrix

    def filterContributions(self, restraint_list, n_iterations = 5,\
                            cutoff = 0.8, max_n = 10):

        import numpy

        cm = self.buildContactMatrix(restraint_list, n_iterations)

        cutoff = numpy.sort(numpy.ravel(cm))[-1000]

        cm = cm * numpy.greater(cm, cutoff)

        for restraint in restraint_list:

            contribs = restraint.getContributions()
            res_numbers = self.getResidueNumbers(contribs)

            weights = numpy.array([cm[i[0], i[1]] for i in res_numbers])

            if numpy.sum(weights) < 1.e-10:
                restraint.isActive(0)
                continue
            
            weights = weights / numpy.sum(weights)

            on, off = self.filter_weights(weights, 1.1, max_n)

            ## active / deactive contributions

            [contribs[i].setWeight(0.) for i in off]
            [contribs[i].setWeight(weights[i]) for i in on]
