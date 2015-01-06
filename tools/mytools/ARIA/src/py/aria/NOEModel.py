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
#import numpy as N
from numpy import *
from aria.Settings import Settings

class NOEModel(AriaBaseClass):

    """
    The main class for calculating NOE spectra from structures.

    Update: Malliavin/Bardiaux
    Becomes an abstract class to didtinct ISPA and RelaxationMatrix
    to include the spin-duffusion correction of distances 
    """
    
    def __init__(self):
        AriaBaseClass.__init__(self)
        
        from aria.Contribution import ContributionEvaluator
        self.__evaluator = ContributionEvaluator()

        self.is_spin_diff = None

class ISPA(NOEModel):

    def __init__(self):
        NOEModel.__init__(self)

        from aria.Contribution import ContributionEvaluator
        self.__evaluator = ContributionEvaluator()


    def calculatePeaksize(self, peak, ensemble):
        """
        for the given peak (AriaPeak) this method computes
        the intensity of a simulated NOE wrt the instance' ensemble
        of structures.

        n_c: number of contributions

        c_i: i-th contribution, n contributions
        <c_i>: ensemble-average for i-th contribution.
        
        NOE = \sum_{i=0}^n_c <c_i>^{-6}

        i.e. it is summed over ensemble-averaged contributions.
        """

        check_type(peak, 'AriaPeak')
        check_type(ensemble, 'StructureEnsemble')

        if not peak:
            self.error(ValueError, 'No contributions in xpk: %d' %
                       peak.getId())

        from aria.mathutils import average

        self.__evaluator.setStructureEnsemble(ensemble)

        ## for each structure: calculate effective distance
        ## for contribution, i.e. distances between atoms
        ## of every spinpair are averaged according to the
        ## type of the given contribution.
        
        f = self.__evaluator.effective_distances
        avg_distances = [f(c) for c in peak.getContributions()]

        ## for each contribution: calculate ensemble-average
        ## TODO: average -> _average, probably faster
        avg_distances = average(avg_distances, axis = 1)

        ## calculate NOEs
        d = power(avg_distances, -6.)

        ## NOE is sum over partial NOEs

        return sum(d)


class SpinDiffusionCorrection(NOEModel):

    def __init__(self):
        
        NOEModel.__init__(self)

        from aria.Contribution import ContributionEvaluator
        self.__evaluator = ContributionEvaluator()

        self.__intensity_matrix = {}

        
    def prepare(self, molecule, ensemble):

        from aria.Relaxation import Relaxation

        self.relaxation = Relaxation()
        self.relaxation.initialize(molecule, ensemble)
        self._spin_first_atom = self.relaxation.getNonEquivalentSpinList()
        self.spin_ids = self.relaxation.spin_list_id
        self._spin_multiplicity = self.relaxation.getSpinMultiplicity()


    def setIntensityMatrix(self, spectrum):
        m = self.relaxation.calculateIntensityMatrix(spectrum)
        spectrum_name = spectrum.getName()
        self.__intensity_matrix[spectrum_name] = m
    
    def getIntensityMatrix(self, name):
        return self.__intensity_matrix[name]
        
    def calculatePeaksize(self, peak, ensemble):

        ## Malliavin 2005/2006
        
        """
        for the given peak (AriaPeak) this method computes
        the intensity of a simulated NOE wrt the instance' ensemble
        of structures.

        n_c: number of contributions

        c_i: i-th contribution, n contributions
        <c_i>: ensemble-average for i-th contribution.
        
        NOE = \sum_{i=0}^n_c <c_i>^{-6}

        i.e. it is summed over ensemble-averaged contributions.
        """

        check_type(peak, 'AriaPeak')
        check_type(ensemble, 'StructureEnsemble')

        if not peak:
            self.error(ValueError, 'No contributions in xpk: %d' %
                       peak.getId())

        from aria.mathutils import average
        from time import clock
        
        self.__evaluator.setStructureEnsemble(ensemble)

        # Modification Therese Malliavin, December 16, 2005
        
        spectrum = peak.getReferencePeak().getSpectrum()
        spectrum_name = spectrum.getName()
        
        intensities = self.getIntensityMatrix(spectrum_name)
        
        
        atoms = [tuple(sp.getAtoms()) for c in peak.getContributions() for sp in c.getSpinPairs()]
        lstintens = []

        spsys = [c.getSpinSystems() for c in peak.getContributions()]

        atoms = [(s[0].getAtoms()[0], s[1].getAtoms()[0]) for s in spsys]

        for a1, a2 in atoms:
            sp1 = self.spin_ids[a1.getId()]
            sp2 = self.spin_ids[a2.getId()]
            lstintens.append(intensities[sp1,sp2])
            
##         for a1, a2 in atoms:

##            #uu = [sp.count(a1) for sp in SpinFirstAtom]
##            #sp1 = uu.index(1)
##            #sp1 = self._get_spin_first_atom_id(a1)
##            sp1 = self.spin_ids[a1.getId()]
        
##            if self._spin_multiplicity[sp1] > 1 and a1 != self._spin_first_atom[sp1][0]:
##               sp1 = 0
           
##            #sp2 = self._get_spin_first_atom_id(a2)
##            sp2 = self.spin_ids[a2.getId()]
##            #uu = [sp.count(a2) for sp in SpinFirstAtom]
##            #sp2 = uu.index(1)
        
##            if self._spin_multiplicity[sp2] > 1 and a2 != self._spin_first_atom[sp2][0]:
##               sp2 = 0

##            if sp1 != 0 and sp2 != 0:
##               lstintens.append(intensities[sp1,sp2])

##         for a1, a2 in atoms:
##             sp1 = self.spin_ids[a1.getId()]
##             sp2 = self.spin_ids[a2.getId()]
##             lstintens.append(intensities[sp1,sp2])
            
        int_aria_pk = sum(lstintens)

        peak.setTheoricVolume(int_aria_pk)

        ## TEST ISPA
        ispa = []
        for a1, a2 in atoms:
            sp1 = self.spin_ids[a1.getId()]
            sp2 = self.spin_ids[a2.getId()]
            ispa.append(self.relaxation.distance_matrix[sp1,sp2])

        
        peak.setIspa(sum(ispa))

        return int_aria_pk

    def _get_spin_first_atom_id(self, a):

        for i in range(len(self._spin_first_atom)):
            if a in self._spin_first_atom[i]: return i
            


