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



from aria.xmlutils import XMLBasePickler, XMLElement
from aria.Topology import BaseAtom as _Base
from aria.Topology import BaseAtomSettings as _BaseSettings
from aria.TypeChecking import *
from aria.Settings import NOT_INIT

class AtomSettings(_BaseSettings):
    pass

class Atom(_Base):

    def __init__(self, settings, name, id):

        _Base.__init__(self, settings, name)
        
        check_type(settings, 'AtomSettings')
        check_int(id)
        
        self.id = id
        self.__segid = None
        self.setResidue(None)
        # BARDIAUX 2.2
        self.homologous = []

    def setId(self, id):
        check_int(id)
        self.id = id

    def getId(self):
        return self.id

    def _setSegid(self, id):
        check_type(id, STRING)
        if len(id) <> 4:
            self.error('Length of segid (%s) must be 4.' % id)

        if self.__segid is not None: print self

        self.__segid = id

    def getSegid(self):
        return self.__segid

    def getResidue(self):
        return self.residue

    def setResidue(self, residue):
        check_type(residue, 'Residue', NONE)
        self.residue = residue

    # BARDIAUX 2.2
    def addHomologous(self, homo):
        check_type(homo, 'Atom', NONE)
        self.homologous.append(homo)

    def setHomologous(self, homo):
        check_type(homo, LIST, NONE)
        self.homologous = homo

    def getHomologous(self):
        return self.homologous

    def __str__(self):

        s = '%s(id=%s, name=%s, type=%s, hetero_name=%s, res=%s, segid="%s")'

        res = self.getResidue()
        
        if res is not None:

            try:
                r_type = res.getType()
            except:
                r_type = NOT_INIT
            
            r_name = str(r_type) + str(res.getNumber())

            chain = res.getChain()

            if chain is not None:
                segid = chain.getSegid()
            else:
                segid = NOT_INIT
            
        else:
            r_name = str(None)
            segid = str(None)

        settings = self.getSettings()

        return s %  (self.__class__.__name__, str(self.getId()),
                     self.getName(),
                     settings.str('type'),
                     settings.str('hetero_atom_name'), r_name, segid)

    __repr__ = __str__
    
class AtomXMLPickler(XMLBasePickler):

    order_ppm = ('segid', 'residue', 'name')
    order_chain = ('name', 'atom_type', 'hetero_name')

    def __init__(self):

        import aria.AriaXML as AX

        XMLBasePickler.__init__(self)

        d = {AX.DTD_CHEMICAL_SHIFT_LIST: self._xml_state_ppm,
             AX.DTD_MOLECULE: self._xml_state_chain,
             AX.DTD_NOESY_SPECTRUM: self._xml_state_ppm,
             AX.DTD_NOE_RESTRAINT: self._xml_state_ppm}
        
        self.dtd_writer_map = d

        d = {AX.DTD_CHEMICAL_SHIFT_LIST: self._ppm_reader,
             AX.DTD_MOLECULE: self._chain_reader,
             AX.DTD_NOESY_SPECTRUM: self._ppm_reader,
             AX.DTD_NOE_RESTRAINT: self._ppm_reader}
        
        self.dtd_reader_map = d

    def _xml_state_ppm(self, atom):
        e = XMLElement(tag_order = self.order_ppm)

        residue = atom.getResidue()

        e.name = atom.getName()
        e.segid = atom.getSegid()
        e.residue = residue.getNumber()

        return e

    def _xml_state_chain(self, atom):
        e = XMLElement(tag_order = self.order_chain)

        s = atom.getSettings()
        e.name = atom.getName()
        e.atom_type = s['type']
        e.hetero_name = s['hetero_atom_name']

        return e
    
    def __find_writer(self):
        
        ## determine DTD

        dtd_name = self.get_dtd_name()

        if dtd_name is not None:

            for name, func in self.dtd_writer_map.items():
                if dtd_name.find(name) >= 0:
                    return func

        return None

    def __find_reader(self):
        
        ## determine DTD

        dtd_name = self.get_dtd_name()

        if dtd_name is None:
            return None

        for name, func in self.dtd_reader_map.items():
            if dtd_name.find(name) >= 0:
                return func

        return None

    def _xml_state(self, atom):

        ## determine DTD-dependent pickler

        pickler = self.__find_writer()

        if pickler is None:
            class_name = self.__class__.__name__
            s = '%s: No XML writer found for DTD "%s".'
            self.error(Exception, s % (class_name, str(self.get_dtd_name())))

        return pickler(atom)

    def _ppm_reader(self, e):

        from aria.Singleton import AtomFactory
        from aria.tools import string_to_segid
        
        factory = AtomFactory()
        
        segid = string_to_segid(str(e.segid))

        number = int(e.residue)

        atom = factory.createAtom(segid, number, str(e.name))

        return atom

    def _chain_reader(self, e):
        return e

    def load_from_element(self, e):
        reader = self.__find_reader()

        if reader is None:
            class_name = self.__class__.__name__
            s = '%s: No XML reader found for DTD "%s".'
            self.error(Exception, s % (class_name, str(self.get_dtd_name())))

        return reader(e)

Atom._xml_state = AtomXMLPickler()._xml_state
