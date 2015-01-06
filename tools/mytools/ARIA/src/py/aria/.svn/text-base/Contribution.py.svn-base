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

import aria.TypeChecking as TCheck

## TODO: which types are actually supported?
## To we need types anyway?

CONTRIBUTION_TYPE_TEST = '__test__'
CONTRIBUTION_TYPE_NORMAL = 'normal'
CONTRIBUTION_TYPE_FAST_EXCHANGE = 'fast_exchange'
CONTRIBUTION_TYPE_FLOATING_ASSIGNMENT = 'floating_assignment'

class ContributionEvaluator(AriaBaseClass):

    """
    Helper class which allows the calculation of the
    average-distance for a given contribution (depend
    on its type) and structure-ensemble.
    """

##     avg_funcs = {CONTRIBUTION_TYPE_NORMAL:
##                  lambda d: average(d, exponent = -6),
                 
##                  CONTRIBUTION_TYPE_TEST:
##                  lambda d: average(d, exponent = -6),

##                  CONTRIBUTION_TYPE_FAST_EXCHANGE:
##                  lambda d: average(d, exponent = -6)}

    def __init__(self):
        AriaBaseClass.__init__(self)
        self.__ensemble = None

    def setStructureEnsemble(self, ensemble):
        TCheck.check_type(ensemble, 'StructureEnsemble')
        self.__ensemble = ensemble

    def getStructureEnsemble(self):
        return self.__ensemble

    def effective_distances(self, contribution):

        from numpy import sum, power

        if self.__ensemble is None:
            s = 'No structure ensemble set.'
            self.error(StandardError, s)

        ## get distance for all spin-pairs

        f = self.__ensemble.getDistances

        d = [f(*sp.getAtoms()) for sp in contribution]

        ## for each structure: calculate partial volume

        volumes = sum(power(d, -6.), axis = 0)

        return power(volumes, -1./6)

class Contribution(AriaBaseClass):
    
    KNOWN_TYPES = (CONTRIBUTION_TYPE_TEST,
                   CONTRIBUTION_TYPE_NORMAL,
                   CONTRIBUTION_TYPE_FAST_EXCHANGE,
                   CONTRIBUTION_TYPE_FLOATING_ASSIGNMENT)

    __counter = 0 # BARDIAUX 2.2
    
    ## TODO: forget spin_pairs, introduce spin_systems!

    def __init__(self, id, type, spin_pairs = None, spin_systems=None):
        """
        type must be a string.
        spin_pairs: list or tuples, elements must be SpinPair-instances
        """
        
        TCheck.check_int(id)
        TCheck.check_string(type)
        TCheck.check_type(spin_pairs, TCheck.LIST, TCheck.TUPLE, TCheck.NONE)

        ## TODO: type checking for spin_systems

        if spin_pairs is not None:
            TCheck.check_elements(spin_pairs, 'SpinPair')

        from aria.Singleton import SpinPairListFactory
        
        if type not in self.KNOWN_TYPES:
            raise 'Type "%s" not known' % type

        if type == CONTRIBUTION_TYPE_NORMAL and spin_pairs and \
               len(spin_pairs) > 1:
            s = 'For type "%s", exactly one spin_pair expected.'
            raise ValueError, s % type
        
        AriaBaseClass.__init__(self)

        #self.__id = id # BARDIAUX 2.2
        self.__id = self.__class__.__counter
        self.__class__.__counter += 1
        
        self.__type = type

        if spin_pairs is not None:
            factory = SpinPairListFactory()
            self.__spin_pairs = factory(spin_pairs)

        else:
            self.__spin_pairs = []

        ## TODO: get rid of spin pairs

        ## TODO: introduce factory for spin_systems 

        if spin_systems is not None:
            self.setSpinSystems(spin_systems)
        else:
            self.__spin_systems = []

        self.setDefaultValues()

        ## BARDIAUX 2.2 Networl anchoring
        self.__score = 0.
        self.__network_score = 0.

    def setDefaultValues(self):

        from aria.Datum import Datum

        self.__figure_of_merit = None
        self.__average_distance = Datum(None)
        self.__weight = None
        self.__degree_of_violation = None

    def setSpinSystems(self, s):
        TCheck.check_tuple(s)
        if len(s) <> 2:
            s = '2 spin systems expected. %d given.'
            self.error(ValueError, s % len(s))

        self.__spin_systems = s

    def getSpinSystems(self):
        return self.__spin_systems
    
    def getId(self):
        return self.__id

    def getWeight(self):
        return self.__weight

    def setWeight(self, w):
        TCheck.check_float(w)

        self.__weight = w

    def getSpinPairs(self):
        return self.__spin_pairs

    def getType(self):
        return self.__type

    def setFigureOfMerit(self, f):
        TCheck.check_type(f, TCheck.FLOAT, TCheck.NONE)
        self.__figure_of_merit = f

    def getFigureOfMerit(self):
        return self.__figure_of_merit

    def setAverageDistance(self, d):
        TCheck.check_type(d, 'Datum')
        self.__average_distance = d

    def getAverageDistance(self):
        return self.__average_distance

    ## BARDIAUX 2.2 for Network Anxchoring
    def setScore(self, s):
        self.__score = s

    def getScore(self):
        return self.__score

    def setNetworkScore(self, ns):
        self.__network_score = ns

    def getNetworkScore(self):
        return self.__network_score    

    def isInter(self):
        a1 = self.__spin_systems[0].getAtoms()[0].getSegid()
        a2 = self.__spin_systems[1].getAtoms()[0].getSegid()
        if a1 <> a2:
            return 1
        else:
            return 0

    ##

    def __getitem__(self, index):
        return self.__spin_pairs[index]

    def __len__(self):

        return len(self.__spin_pairs)

    def __str__(self):
        s = 'Contribution(id=%s, type="%s", weight=%s, ' + \
            'figure_of_merit=%s, average_distance=%s, ' + \
            'spin_pairs=%s)'

        return s % (str(self.getId()),
                    str(self.getType()),
                    str(self.getWeight()),
                    str(self.getFigureOfMerit()),
                    str(self.getAverageDistance()),
                    str(self.getSpinPairs()))

    __repr__ = __str__

class ContributionXMLPickler:

    order = ['id', 'weight', 'figure_of_merit', 'spin_system',
             'average_distance']

    def _xml_state(self, x):

        from aria.xmlutils import XMLElement

        e = XMLElement(tag_order = self.order)
        e.id = x.getId()
        e.type = x.getType()
        e.figure_of_merit = x.getFigureOfMerit()
        e.weight = x.getWeight()
        e.average_distance = x.getAverageDistance()
        e.spin_system = x.getSpinSystems()

        ss = x.getSpinSystems()
        
        if not ss:
            self.error(ValueError, 'Spin Systems not set.')
        
        e.spin_system = ss
        
        return e

    def load_from_element(self, e):

        from aria.tools import as_tuple

        try:
            t = str(e.type)
        except:
            t = CONTRIBUTION_TYPE_NORMAL

        sp = as_tuple(e.spin_system)
            
##         c = Contribution(int(e.id), str(e.type), as_tuple(e.spin_pair))
        c = Contribution(int(e.id), t, spin_systems=sp)

        c.setWeight(float(e.weight))
        c.setFigureOfMerit(e.figure_of_merit)
        c.setAverageDistance(e.average_distance)
        c.setSpinSystems(as_tuple(e.spin_system))

        return c

Contribution._xml_state = ContributionXMLPickler()._xml_state

    
