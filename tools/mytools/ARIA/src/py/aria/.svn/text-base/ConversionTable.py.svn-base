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

CONVERSION_DATA = 'atomnames.xml'

IUPAC_CONVENTION = 'iupac'
DYANA_CONVENTION = 'dyana'
CNS_CONVENTION = 'cns'

NAMING_CONVENTIONS = (IUPAC_CONVENTION, DYANA_CONVENTION, CNS_CONVENTION)

## TODO: implement conversion error class

class AtomnameConversionError(Exception):
    pass

class ConversionTable(AriaBaseClass):

    formats = NAMING_CONVENTIONS

    def __init__(self):

        from os.path import join

        AriaBaseClass.__init__(self)
        
        path = join(AriaBaseClass.data_path, CONVERSION_DATA)
        self.load_from_xml(path)

    def has_format(self, format):

        return format in self.formats    

    def load_from_xml(self, xml_file):

        import aria.xmlutils as xmlutils

        content_handler = xmlutils.XMLContentHandler()
        pickler = xmlutils.XMLPickler(content_handler)
        
        conversion_table = pickler.load(xml_file).conversion_table

        table = {}

        for residue in conversion_table.residue:

            key = str(residue.residue_type)

            if not table.has_key(key):

                table[key] = {}

                for format in self.formats: table[key][format] = {}

            r_iupac = str(residue.iupac_name)
            r_dyana = str(residue.dyana_name)
            r_cns = str(residue.cns_name)

            table[key][IUPAC_CONVENTION][r_iupac] = {DYANA_CONVENTION: r_dyana,
                                                     CNS_CONVENTION: r_cns}
            
            table[key][DYANA_CONVENTION][r_dyana] = {IUPAC_CONVENTION: r_iupac,
                                                     CNS_CONVENTION: r_cns}
            
            table[key][CNS_CONVENTION][r_cns] = {DYANA_CONVENTION: r_dyana,
                                                 IUPAC_CONVENTION: r_iupac}

            # BARDIAUX 2.3
            # if residue has only one atom
            atoms = residue.atom
            if type(residue.atom) <> type([]):
                atoms = [residue.atom]
                
            for atom in atoms:

                a_iupac = str(atom.iupac_name)
                a_dyana = str(atom.dyana_name)
                a_cns = str(atom.cns_name)

                table[key][IUPAC_CONVENTION][r_iupac][a_iupac] = \
                        {DYANA_CONVENTION: a_dyana, CNS_CONVENTION: a_cns}
                
                table[key][DYANA_CONVENTION][r_dyana][a_dyana] = \
                        {IUPAC_CONVENTION: a_iupac, CNS_CONVENTION: a_cns}
                
                table[key][CNS_CONVENTION][r_cns][a_cns] = \
                        {IUPAC_CONVENTION: a_iupac, DYANA_CONVENTION: a_dyana}

        self.table = table

    def error(self, exception = None, error = '', msg = None, raise_error = 1):
        if raise_error <> 1: return

        AriaBaseClass.error(self, exception, error, msg)
                
    def convert_residue(self, name, format, target_format, type,
                        raise_error = 1):

        from aria.TypeChecking import check_string

        check_string(name)
        check_string(type)

        if not self.has_format(format):
            self.error(ValueError, 'Residue format "%s" not known.' \
                       % str(format), raise_error = raise_error)
            return None
        
        if not self.has_format(target_format):
            self.error(ValueError, 'Residue format "%s" not known.' \
                       % str(target_format), raise_error = raise_error)
            return None
        
        table = self.table

        if not table.has_key(type):
            
            self.error(ValueError, 'Residue of type "%s" not supported.' \
                       % type, raise_error = raise_error)
            return None

        if not table[type][format].has_key(name):

            m = 'Residue "%s" of type "%s" and format "%s" not known.' 
            self.error(ValueError, m % (name, type, format),
                       raise_error = raise_error)
            return None
            
        if target_format == format:
            return name
        else:
            return table[type][format][name][target_format]

    def convert_atom(self, residue, atom, format, target_format, type,
                     raise_error = 1):
        """
        Known bug: Conversion of atom O from dyana-format to cns / iupac
        will always yield O' for all amino acids.
        """
        from aria.TypeChecking import check_string

        check_string(residue)
        check_string(atom)
        check_string(type)

        if not self.has_format(format):
            self.error(ValueError, 'Residue format "%s" not known.' \
                       % format, raise_error = raise_error)
            return None
            
        if not self.has_format(target_format):
            self.error(ValueError, 'Residue format "%s" not known.' \
                       % target_format, raise_error = raise_error)
            return None

        table = self.table

        if not table.has_key(type):            
            self.error(ValueError, 'Residue of type "%s" not supported.' \
                       % type, raise_error = raise_error)
            return None

        if not table[type][format].has_key(residue):
            self.error(ValueError, 'Residue "%s" of type "%s" not known.' \
                       % (residue, type), raise_error = raise_error)
            return None
            
        if not table[type][format][residue].has_key(atom):
            self.error(ValueError, 'Atom "%s" not known in residue "%s".' \
                       % (atom, residue), raise_error = raise_error)
            return None

        if target_format == format:
            return atom
        else:
            return table[type][format][residue][atom][target_format]

    def convert_atoms(self, residue, atoms, format, target_format, type,
                      raise_error = 1):

        from aria.TypeChecking import check_string, check_elements, check_type
        from aria.TypeChecking import LIST, TUPLE, STRING
        from aria.Topology import TYPE_AMINO_ACID
        
        check_string(residue)
        check_type(atoms, LIST, TUPLE)
        check_elements(atoms, STRING)
        check_string(type)

        if not self.has_format(format):
            self.error(ValueError, 'Residue format "%s" not known.' \
                       % format, raise_error = raise_error)
            return None
            
        if not self.has_format(target_format):
            self.error(ValueError, 'Residue format "%s" not known.' \
                       % target_format, raise_error = raise_error)
            return None
        
        table = self.table

        if not table.has_key(type):            
            self.error(ValueError, 'Compound of type "%s" not supported.' \
                       % type, raise_error = raise_error)
            return None

        if not table[type][format].has_key(residue):
            self.error(ValueError, 'Residue "%s" of type "%s" not known.' \
                       % (residue, type), raise_error = raise_error)
            return None
            
        if target_format == format: return tuple(atoms)

        converted_atoms = []

        for atom in atoms:

            if not table[type][format][residue].has_key(atom):                
                self.error(ValueError, 'Atom "%s" not known in residue "%s".' \
                           % (atom, residue), raise_error = raise_error)
                return None
                
            converted_atom = table[type][format][residue][atom][target_format]
            
            if format == 'dyana' and type == TYPE_AMINO_ACID and \
               atom == 'O' and not 'OXT' in atoms:

                converted_atom = atom

            converted_atoms.append(converted_atom)

        return tuple(converted_atoms)
