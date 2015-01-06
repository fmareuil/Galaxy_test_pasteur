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
from numpy import *


class Relaxation(AriaBaseClass):

    def __init__(self):

        AriaBaseClass.__init__(self, name = 'Relaxation')

        self.spin_list = []
        self.spin_list_id = []
        self.distance_matrix = None
        self.cut_distance_matrix = None
        self._spin_multiplicity = None
        self._sq_spin_multiplicity = None

        self.relaxation_matrix = {}
        self.intensity_matrix = {}


    def initialize(self, molecule, ensemble):

        """
        spectrum independant
        """
        
        ## List of non equivelent protons
        self.spin_list = self.setNonEquivalentSpinList(molecule)
        
        
        ## spin multiplicity
        self._spin_multiplicity = self.setSpinMultiplicity()
        self._sq_spin_multiplicity = self.setSquareSpinMultiplicity()

        ## distance matrix
        from time import clock
        t =clock()
        
        self.message('Calculating distances matrix ...')
        self.distance_matrix = self.setDistanceMatrix(ensemble)
        self.message('Done ...')
        self.debug('Time: %ss' % str(clock() - t))
        
    def setSpinMultiplicity(self):
        return [len(a) for a in self.spin_list]

    def getSpinMultiplicity(self):
        return self._spin_multiplicity
    
    def setSquareSpinMultiplicity(self):        
        return sqrt(self._spin_multiplicity)

    def getSquareSpinMultiplicity(self):        
        return self._sq_spin_multiplicity
    
    def setDistanceMatrix(self, ensemble):
        # Bardiaux rMat
               
        """
        Calculate the distances**-6 matrix between non equivalent spins.
        """
        from aria.mathutils import _average

        # disable cache to save memory
        AriaBaseClass.cache = 0
        
        spin_list = self.spin_list
        nbspin = len(spin_list)

        distances = zeros((nbspin, nbspin), 'd')
        # half-matrix + diagonal

        for i in range(nbspin):
            for j in range(i, nbspin):
                distances[i,j] = self._dist_spin(spin_list[i], spin_list[j], ensemble)

        # symmetrization        
        tdistances = distances + transpose(distances)
        putmask(tdistances, identity(nbspin), distances)

        # restore cache
        AriaBaseClass.cache = 1
        
        #return self.apply_cutoff(tdistances)
        self.cut_distance_matrix = self.apply_cutoff(tdistances)
        
        return tdistances
    
    def getDistanceMatrix(self):
        return self.distance_matrix
    
    def getCutDistanceMatrix(self):
        return self.cut_distance_matrix
        
    def setNonEquivalentSpinList(self, molecule):

        """
        set a list of list of unequivalent protons.
        METHYL atoms are equivalent.
        NTERMINUS atoms are equivalent.
        METHYLS from ISOPROPYL are equivalent.

        METHYLENE are not equivalent.

        Others (H, HA..) are not equivalent.
        What about AROMATIC.

        original implementation : Malliavin
        refactoring : Bardiaux
        """

        spin_list = []

        n_at = 0
        
        for c in molecule.get_chains():
            for r in c.getResidues():
                
                n_at += len(r.getAtoms())
                
                protons = [a for a in r.getAtoms() if a.isProton()]
                protons.sort(lambda a, b: cmp(a.getId(), b.getId()))

                it = iter(protons)
                while True:
                    try:
                        p = it.next()
                    except StopIteration:
                        break
                    
                    g = p.getEquivalentGroups()
                    if not g:
                        spin_list.append([p])
                        
                    else:
                        g = g[0]
                        if g.getType() == 'METHYLENE':
                            spin_list.append([p])

                        elif g.getType() in ['METHYL', 'NTERMINUS']:
                            at = g.getAtoms()
                            spin_list.append(list(at))
                            [it.next() for j in range(len(at)-1)]

                        elif g.getType() == 'ISOPROPYL':
                            at = g.getAtoms()
                            iso1, iso2 = at[0:3], at[3:6]
                            spin_list.append(list(iso1))
                            spin_list.append(list(iso2))
                            x = [it.next() for j in range(len(at)-1)]

                        elif g.getType() == 'AROMATIC':
                            if p.getName() == 'HZ':
                                spin_list.append([p])
                            else:    
                                at = g.getAtoms()
                                spin_list.append(list(at))
                                [it.next() for j in range(len(at)-1)]

        # NEW:
        self.spin_list_id = zeros(n_at)
        for a in range(len(spin_list)):
            atoms = spin_list[a]
            for i in atoms:
                self.spin_list_id[i.getId()] = a
        
        
        return spin_list

    def getNonEquivalentSpinList(self):
        return self.spin_list
        
    def calculateIntensityMatrix(self, spectrum):

        ## Malliavin 2005/2006

        # spin multiplicity
        # calculate the approximated intensity matrix

        rmat = self.calculateRelaxationMatrix(spectrum)

        matrixmultiply = dot
        # mixing_time
        #mixing_time = spectrum.getMixingTime()
        mixing_time = spectrum.getExperimentData()['spectrum_mixing_time']
        
        Kpower = 10
        delta = float(mixing_time/(1000. * 2**Kpower))
        

        nbat = shape(rmat)[0]    
        intens0 = zeros((nbat,nbat),'d')
        putmask(intens0, identity(nbat), self._sq_spin_multiplicity)
        
        i_Rdt = array(identity(nbat)-(rmat*delta), 'd')

        for i in range(Kpower):
            i_Rdt = matrixmultiply(i_Rdt, i_Rdt)
            intens = matrixmultiply(i_Rdt,intens0)

        intens = matrixmultiply(i_Rdt, intens0)
         
        intens1 = zeros((nbat,nbat),'d')
        putmask(intens1, identity(nbat), self._sq_spin_multiplicity)
        intens = matrixmultiply(intens1,intens)
        
        return intens
    
    def getIntensityMatrix(self, spectrum):
        return self.intensity_matrix[spectrum]
    
    
    def calculateRelaxationMatrix(self, spectrum):

        ## Malliavin 2005/2006
        
        """
        Calculate the relaxation matrix.
        """
        
        
        ## pi = 4*arctan(1)

        # constants
        # http://physics.nist.gov/cgi-bin/cuu/Value?gammap

        # calculate the relaxation matrix
        mgyr = 2.6751965
        mhbar = 1.0545919
        sord = 1.0

        # correction factors
        EOM=9                     #GHz-->Hz
        ETAU=-9                   #ns-->s
        EDIST=-10                 #A -->m
        EMU=-14                   #Mu0**2/16pi**2
        #      PARAMETER (MGYR = 2.6751965D0, EGYR=8)
        #     PARAMETER (MHBAR = 1.0545919D0, EHBAR=-34)
        crat=1.0                  #10**(4*EGYR+2*EHBAR+ETAU-6*EDIST+EMU-1)

        #magnet_field = spectrum.getMagneticField()
        magnet_field = spectrum.getExperimentData()['spectrometer_frequency']
        
        omega = magnet_field/1000.0
        #tauc = self._correlation_time
        tauc = spectrum.getExperimentData()['molecule_correlation_time']
        
        # calculs intermediaires
        con0 = (mgyr**4)*(mhbar**2)
        con1 = 3.0*con0
        con2 = 6.0*con0
        #omsq = 2.0*pi*omega*1.0e-9
        omsq = 2.0*pi*omega
        omsq = omsq*omsq

        # get spectral densities
        #        tau = tauc*10e-9
        tau = tauc
        tsws = tau*tau*omsq
        
        ts2ws = 4*tsws
        s0 = sord*con0*tauc
        s1 = sord*con1*tauc/(1+tsws)
        s2 = sord*con2*tauc/(1+ts2ws)

        # replace d by the inverse sixth power of distances 
        d = self.getCutDistanceMatrix()
        nbat = shape(d)[0]

        rmat = zeros((nbat,nbat),'d')
        
        # spin multiplicity
        matrixmultiply = dot
        
        # calculate the transition probability
        dtri = (s1+s2+s0)*matrixmultiply(self._spin_multiplicity, d)

        # pour recuperer grp utiliser getAtoms ds AtomGroup ds Topology.py
        # calculate the relaxation matrix non-diagonal elements
        self._sq_spin_multiplicity = sqrt(self._spin_multiplicity)
        
        uu = multiply.outer(self._sq_spin_multiplicity, self._sq_spin_multiplicity)
        rmat = (s2-s0)*d*uu

        # calculate the relaxation matrix diagonal elements
        vv = 2.0*(self._spin_multiplicity - ones((nbat),'d'))*(((s1/2.0)+s2)*ones((nbat),'d'))

        vv = vv + dtri 

        # add z leakage and transition probabilities to diagonal terms
        # zleak = 0.0 in CNS
        #        vv = vv + dtri + (zleak*ones((nbat),'d'))
        #        vv = vv + dtri 

        putmask(rmat, identity(nbat), vv)
        
        return rmat

##     def getRelaxationMatrix(self, spectrum):
##         return self.relaxation_matrix[spectrum]
        
##     def integrate(self, rmat, mixing_time):

##         ## Malliavin 2005/2006

##         # spin multiplicity
##         # calculate the approximated intensity matrix

##         Kpower = 10
##         delta = float(mixing_time/(1000. * 2**Kpower))
        

##         nbat = shape(rmat)[0]    
##         intens0 = zeros((nbat,nbat),'d')
##         putmask(intens0, identity(nbat), self._sq_spin_multiplicity)

##         i_Rdt = array(identity(nbat)-(rmat*delta), 'd')

##         for i in range(Kpower):
##             i_Rdt = matrixmultiply(i_Rdt, i_Rdt)
##             #intens = matrixmultiply(val,intens)

##         intens = matrixmultiply(i_Rdt, intens0)
        
##         intens1 = zeros((nbat,nbat),'d')
##         putmask(intens1, identity(nbat), self._sq_spin_multiplicity)
##         intens = matrixmultiply(intens1,intens)

##         return intens

    

    
    def apply_cutoff(self, dmat):

        # cutoff hard-coded = 6 A
        cutoff = power(5., -6.)
        size = shape(dmat)

        from copy import copy
        cdmat = copy(dmat)
        #sel = where(dmat>cutoff, zeros(size[0],size[1],Float),ones(size[0],size[1],Float))

        sel = where(cdmat < cutoff, zeros(size,'d'), ones(size,'d'))

        for i in range(0,size[0]-1):
           list1 = sel[i,:]
           for j in range(0,i):
              list2 = sel[:,j]

              if dot(list1,list2) == 0.: 
                 cdmat[i,j] = 0.0
                 cdmat[j,i] = 0.0

        return cdmat
    
    # Bardiaux rMat
    def _dist_spin(self, spin1, spin2, ensemble):
        
        """
        average dist**-6 between spin1 and spin2 regarding
        ensemble.
        zero distances return 0 volume.
        """
        from aria.mathutils import _average

        distances = [ensemble.getDistances(at1,at2) for at1 in spin1 for at2 in spin2 if at1 <> at2]

        if distances:
            d = _average(distances)
        else:
            return 0.

        return _average(power(d, -6.))
    
        #if d == 0.:
        #    return 0.
        #else:
        #    return power(d, -6.)
        
        #[distances.extend(ensemble.getDistances(at1,at2)) for at1 in spin1 for at2 in spin2]
        
        #if spin1 == spin2:
        #    
        #    distances = array(distances)
        #    non_zero_d = nonzero(distances)
        #    put(distances, non_zero_d, power(take(distances, non_zero_d), -6.))
        #    return _average(distances)
        
        #else:
        #    return _average(power(distances, -6.))

            
