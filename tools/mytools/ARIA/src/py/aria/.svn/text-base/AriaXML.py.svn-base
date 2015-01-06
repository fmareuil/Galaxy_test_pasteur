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
from aria.xmlutils import *
from exceptions import Exception
from aria.tools import last_traceback

## DTD base names

DTD_CHEMICAL_SHIFT_LIST = 'chemical_shift_list'
DTD_MOLECULE = 'molecule'
DTD_NOESY_SPECTRUM = 'noesy_spectrum'
DTD_PROJECT = 'project'
DTD_NOE_RESTRAINT = 'noe_restraint'

## DOCTYPEs

DOCTYPE_NOE_RESTRAINTS = 'noe_restraint_list'

class DOCManager(AriaBaseClass):

    VERSION = 1.0

    def __init__(self):

        from aria.Project import Project
        from aria.AriaPeak import AriaPeak
        from aria.Topology import Topology
        from aria.Molecule import Molecule
        from aria.NOESYSpectrum import NOESYSpectrum
        from aria.ChemicalShiftList import ChemicalShiftList
        from aria.Project import ProjectSingleton
        from aria.conversion import ConverterSettings

        AriaBaseClass.__init__(self)

        ## DTD base names
        
        d = {Project: DTD_PROJECT,
             AriaPeak: DTD_NOE_RESTRAINT,
             Topology: 'topology',
             Molecule: DTD_MOLECULE,
             NOESYSpectrum: DTD_NOESY_SPECTRUM,
             ChemicalShiftList: DTD_CHEMICAL_SHIFT_LIST,
             ProjectSingleton: DTD_PROJECT,
             ConverterSettings: 'conversion'}

        for key, value in d.items():
            d[key] = self._compile_name(value)

        self.__dtd_base_names = d

        ## DOCTYPEs

        d = {Project: 'project',
             AriaPeak: DOCTYPE_NOE_RESTRAINTS,
             Topology: 'topology',
             Molecule: 'molecule',
             NOESYSpectrum: 'spectrum',
             ChemicalShiftList: 'chemical_shift_list',
             ProjectSingleton: 'project',
             ConverterSettings: 'conversion'}

        self.__doc_types = d

    def _compile_name(self, base):
        return '%s%s.dtd' % (base, self.VERSION)

    def getDOCTypes(self):
        return self.__doc_types.values()

    def getDTDs(self):
        return self.__dtd_base_names.values()

    def is_validDOCType(self, name):
        return name in self.getDOCTypes()

    def is_validDTD(self, name):
        return name in self.getDTDs()

    def getDTD(self, c, is_class=0):
        """
        returns DTD for class 'c'
        """

        if not is_class:

            if not type(c) == type(self):
                s = 'Instance expected, "%s" given.' % str(type(c))
                self.error(s)

            c = c.__class__

        if not c in self.__dtd_base_names:
            s = 'DTD for class "%s" unknown'
            self.error(ValueError, s % c.__name__)

        dtd_name = self.__dtd_base_names[c]
        
        return dtd_name

    def getDOCType(self, c, is_class=0):

        if not is_class:

            if not type(c) == type(self):
                s = 'Instance expected, "%s" given.' % str(type(c))
                self.error(s)

            c = c.__class__

        if not c in self.__doc_types:
            s = 'DOCTYPE for class "%s" unknown'
            self.error(ValueError, s % c.__class__.__name__)
            
        return self.__doc_types[c]

class AriaXMLContentHandler(XMLContentHandler):
    
    def __init__(self):

        from aria.Project import ProjectXMLPickler
        from aria.AriaPeak import AriaPeakXMLPickler
        from aria.Topology import TopologyXMLPickler
        from aria.Molecule import MoleculeXMLPickler
        from aria.NOESYSpectrum import NOESYSpectrumXMLPickler
        from aria.ChemicalShiftList import ChemicalShiftListXMLPickler
        from aria.conversion import ConverterSettingsXMLPickler
	## < Mareuil
        from JobManager import JobManagerXMLPickler
        ## Mareuil >

        XMLContentHandler.__init__(self)
        
        pickler = {}

        ## register pickler for different DOCTYPEs
        
        pickler['project'] = ProjectXMLPickler()
        pickler[DOCTYPE_NOE_RESTRAINTS] = AriaPeakXMLPickler()
        pickler['topology'] = TopologyXMLPickler()
        pickler['molecule'] = MoleculeXMLPickler()
        pickler['spectrum'] = NOESYSpectrumXMLPickler()
        pickler['chemical_shift_list'] = ChemicalShiftListXMLPickler()
        pickler['conversion'] = ConverterSettingsXMLPickler()
        ## < Mareuil 0.1
        pickler['jobmanager'] = JobManagerXMLPickler()
        ## Mareuil 0.1 >

        self.pickler = pickler

    def select(self, doc_type):
        if doc_type in self.pickler.keys():
            self.sub_pickler = self.pickler[doc_type]
        else:
            self.sub_pickler = None

    def load_from_element(self, name, e):

        try:
            pickler = getattr(self.sub_pickler, name)
        except:
            return e

        return pickler.load_from_element(e)

class AriaXMLPickler(XMLPickler, AriaBaseClass):

    def __init__(self):

        AriaBaseClass.__init__(self)

        handler = AriaXMLContentHandler()
        XMLPickler.__init__(self, handler)

        self.doc_manager = DOCManager()

    def createParser(self):
        parser = XMLPickler.createParser(self)

        if hasattr(parser, 'is_self_made'):
            msg = 'Expat XML parser unavailable, using Python version.'
            self.message(msg)

        return parser

    def check_document(self, doctype, dtd):
        if not self.doc_manager.is_validDOCType(doctype):
            s = 'Unknown DOCTYPE: "%s". Check your XML file.'
            self.error(TypeError, s % doctype)

        if not self.doc_manager.is_validDTD(dtd):
            s = 'Unknown DTD: "%s". Check your XML file.'
            self.error(TypeError, s % dtd)

    def parse_doc_type(self, f):

        doc_type, dtd = XMLPickler.parse_doc_type(self, f)

        self.check_document(doc_type, dtd)

        handler = self.getContentHandler()
        handler.select(doc_type)

        return doc_type, dtd

    def __load(self, filename, gzip):

        try:
            doc = XMLPickler.load(self, filename, gzip)

        except XMLTagError, msg:
            print last_traceback()
            msg = 'XML file %s: %s' % (filename, msg)
            msg += '\nMake sure that all tags/attributes are ' + \
                   'spelled correctly and that no mandatory tags' + \
                   '/attributes are missing.'
            self.error(XMLReaderError, msg )
            
        except Exception, msg:
            print last_traceback()
            self.error(Exception, '%s: %s' % (filename, msg))

        ## TODO: hacked

        return getattr(doc, doc.get_doc_type().split(':')[-1])

    def load(self, filename, gzip = 0):
        XMLBasePickler.relaxed = 0
        return self.__load(filename, gzip)

    def load_relaxed(self, filename, gzip = 0):
        XMLBasePickler.relaxed = 1
        x = self.__load(filename, gzip)
        XMLBasePickler.relaxed = 0

        return x

    def create_document(self, object):
        
        dtd = self.doc_manager.getDTD(object)
        doc_type = self.doc_manager.getDOCType(object)

        doc = XMLDocument(doc_type, dtd)
        setattr(doc, doc_type, object)

        return dtd, doc
                
    def dumps(self, object):

        if is_type(object, LIST) or is_type(object, TUPLE):

            try:
                check_elements(object, 'AriaPeak')
                is_peak = 1
            except:
                is_peak = 0

            if is_peak:

                if object:
                    dtd = self.doc_manager.getDTD(object[0])
                    doc_type = self.doc_manager.getDOCType(object[0])

                else:

                    from aria.AriaPeak import AriaPeak
                    
                    dtd = self.doc_manager.getDTD(AriaPeak, is_class=1)
                    
                    doc_type = self.doc_manager.getDOCType(AriaPeak,
                                                           is_class=1)
                    
                doc = XMLDocument(doc_type, dtd)
                
                new_object = XMLElement()
                new_object.peak = object
                object = new_object

                setattr(doc, doc_type, object)

            ## TODO: what happens if is_peak == 0 ?
                
        else:
            dtd, doc = self.create_document(object)

        handler = self.getContentHandler()

        handler.set_dtd_name(doc.get_dtd())

        xml = XMLPickler.dumps(self, doc)

        handler.release()

        return xml

if __name__ == '__main__':

    pass
