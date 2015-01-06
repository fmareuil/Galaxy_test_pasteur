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
from aria.xmlutils import XMLBasePickler, XMLElement
from aria.Settings import Settings, Entity

class Molecule(AriaBaseClass):

    def __init__(self, name):
        check_string(name)

        AriaBaseClass.__init__(self)

        self.__name = name
        self.__chains = []
        self.data_source = None

        ## to keep track of the molecule in a ccpn project

        self.code = None
        
    def getName(self):
        return self.__name

    def add_chain(self, c):
        check_type(c, 'Chain')

        self.__chains.append(c)

        c.setMolecule(self)

    def get_chains(self):
        return self.__chains
        
    def setDataSource(self, s):
        check_type(s, 'SequenceData')
        self.data_source = s

    def getDataSource(self):
        return self.data_source

    def getChain(self, segid):

        for chain in self.get_chains():
            if chain.getSegid() == segid:
                return chain
        self.error(IndexError, 'Chain with segid "%s" does not exist.' % segid)

    def getType(self):
        ## TODO: hack
        molecule_type = [0, None]

        for chain in self.get_chains():
            if molecule_type[0] < len(chain.getResidues()):
                molecule_type[1] = chain.getType()
                
        return molecule_type[1]
    
    ## BARDIAUX 2.2 extend
    def getTypes(self):
        """
        return list of chain types
        """
        molecule_types = []

        for chain in self.get_chains():
            if len(chain.getResidues()) > 0:
                if chain.getType() not in molecule_types:
                    molecule_types.append(chain.getType())
                
        return molecule_types

class MoleculeXMLPickler(XMLBasePickler):

    order = ('name', 'chain')

    def __init__(self):
        from aria.Chain import ChainXMLPickler
        from aria.Residue import ResidueXMLPickler
        from aria.Topology import EquivalentGroupXMLPickler
        from aria.Atom import AtomXMLPickler

        self.molecule = self
        self.chain = ChainXMLPickler()
        self.residue = ResidueXMLPickler()
        self.atom = AtomXMLPickler()
        self.equivalent_group = EquivalentGroupXMLPickler()

    def _xml_state(self, m):

        e = XMLElement(tag_order = self.order)
        e.name = m.getName()
        e.chain = m.get_chains()
##        e.atom_nomenclature = Topology.NOMENCLATURE_IUPAC 

        return e

    def load_from_element(self, e):

        from aria.tools import as_tuple

        name = str(e.name)
        m = Molecule(name)

        chains = as_tuple(e.chain)

        [m.add_chain(c) for c in chains]

        return m

Molecule._xml_state = MoleculeXMLPickler()._xml_state
