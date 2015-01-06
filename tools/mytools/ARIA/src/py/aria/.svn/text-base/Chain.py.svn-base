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
from aria.xmlutils import XMLBasePickler
from aria.Settings import Settings, Entity

TYPE_PROTEIN = 'PROTEIN'
TYPE_DNA = 'DNA'
TYPE_RNA = 'RNA'
TYPE_NONPOLYMER = 'NONPOLYMER'

class AtomNomenclatureEntity(Entity):

    def __init__(self):

        descr = 'Atom-nomenclature follows the IUPAC standard.'
        err_msg = 'Atom-nomenclature must be "IUPAC"'
        
        Entity.__init__(self, descr, err_msg)
    
    def is_valid(self, v):

        from aria.Topology import ATOM_NOMENCLATURES
        
        return type(v) == type('') and v in ATOM_NOMENCLATURES

class ChainSettings(Settings):

    def create(self):

        from aria.Settings import ChoiceEntity, TypeEntity

        chain_types = [TYPE_PROTEIN, TYPE_DNA, TYPE_RNA, TYPE_NONPOLYMER]
        
        return {'type': ChoiceEntity(elements = chain_types),
                'nomenclature': TypeEntity('Topology'),
                'atom_nomenclature': AtomNomenclatureEntity()}

    def create_default_values(self):

        import aria.Topology as Topology

        d = {'nomenclature': Topology.load_topology()}
        d['atom_nomenclature'] = Topology.NOMENCLATURE_IUPAC

        return d

class Chain(AriaBaseClass):
    
    def __init__(self, settings, segid):

        check_type(settings, 'ChainSettings')
        check_type(segid, STRING)

        from aria.tools import string_to_segid

        if len(segid) > 4:
            self.error('Max. length of segid is 4')
                
        AriaBaseClass.__init__(self, settings)

        self.__segid = string_to_segid('%4s' % segid)
            
        self.residues = {}
        self.data_source = None

        self.__molecule = None

    def __getitem__(self, number):

        check_int(number)

        if not self.hasResidue(number):

            self.error(KeyError, 'Residue no. %d not in chain.' % number)

        return self.residues[number]

    def setMolecule(self, molecule):

        check_type(molecule, 'Molecule')

        self.__molecule = molecule

    def getMolecule(self):

        return self.__molecule
    
    def getType(self):
        return self.getSettings()['type']

    def getSegid(self):
        return self.__segid

    def hasResidue(self, number):
        check_int(number)

        return self.residues.has_key(number)

    def delResidue(self, number):

        check_int(number)

        if not self.hasResidue(number):
            m = 'Residue no. %d is not contained in chain "%s".'
            self.error(IndexError, m % (number, self.getSegid()))

        self.residues[number].setChain(None)        
        del self.residues[number]

    def addResidue(self, residue):

        check_type(residue, 'Residue')

        number = residue.getNumber()

        if self.hasResidue(number):
            self.error(KeyError, 'Residue no. %d already in chain.' % number)
            
        self.residues[number] = residue
        self.residues[number].setChain(self)

    def getResidues(self):
        residues = self.residues.values()
        residues.sort(lambda a, b: cmp(a.getNumber(), b.getNumber()))

        return residues

    def linkResidues(self):
        [r.link() for r in self.getResidues()]

    def checkNomenclature(self):
        n = self.getSettings()['nomenclature']
        s = [n.check(r) for r in self.getResidues()]

        [self.warning(x) for x in s if x is not None]

    def setDataSource(self, s):
        check_type(s, 'SequenceData')
        self.data_source = s

    def getDataSource(self):
        return self.data_source

    def __len__(self):
        return len(self.residues)

class ChainXMLPickler(XMLBasePickler):

    order = ['chain_type', 'segid', 'residue']

    def _xml_state(self, chain):

        from aria.xmlutils import XMLElement

        s = chain.getSettings()

        e = XMLElement(tag_order = self.order)
        e.segid = chain.getSegid()
        e.chain_type = s['type']
        e.residue = tuple(chain.getResidues())

        return e

    def load_from_element(self, e):

        from aria.Singleton import AtomFactory
        from aria.Residue import Residue
        from aria.tools import as_tuple
        
        factory = AtomFactory()

        segid = str(e.segid)

        settings = ChainSettings()
        settings.reset()
        settings['type'] = str(e.chain_type)
            
        chain = Chain(settings, segid)
        segid = chain.getSegid()

        for r in as_tuple(e.residue):

            number = int(r.number)
            code = str(r.residue_type)
            ## BARDIAUX 2.2
            if hasattr(r, 'structure'):
                structure = str(r.structure)
            else:
                structure = ""

            if chain.hasResidue(number):

                m = 'Chain segid=%s: multiple definitions for residue' + \
                    ' "%d" found in xml file.'
                self.error(ValueError, m % (str(segid), number))

            residue = Residue(number, code, structure)
            chain.addResidue(residue)

            for a in as_tuple(r.atom):

                atom_name = str(a.name)
                atom_type = str(a.atom_type)

                if a.hetero_name == '':
                    hetero = None
                else:
                    hetero = str(a.hetero_name)

                atom = factory.createAtom(segid, number, atom_name, atom_type,
                                          hetero)

                residue.addAtom(atom)

            if hasattr(r, 'equivalent_group'):
                [residue.addEquivalentGroup(g) for g in
                 as_tuple(r.equivalent_group)]

        chain.linkResidues()
        if settings['type'] in [TYPE_PROTEIN, TYPE_DNA, TYPE_RNA]:
            chain.checkNomenclature()
        
        return chain
        
Chain._xml_state = ChainXMLPickler()._xml_state
