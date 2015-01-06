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



## TODO: since grepping of energies is hard-coded. changes in
## pdb-format could be implemented in a sub-class.

from aria.ariabase import AriaBaseClass

from aria.Chain import TYPE_PROTEIN, TYPE_RNA, TYPE_DNA, TYPE_NONPOLYMER
from aria.Topology import TYPE_AMINO_ACID, TYPE_DNA_BASE, TYPE_RNA_BASE, \
     TYPE_NONBASE
from aria.ConversionTable import IUPAC_CONVENTION, DYANA_CONVENTION, \
     CNS_CONVENTION

BASE_TYPES = {TYPE_PROTEIN: TYPE_AMINO_ACID,
              TYPE_RNA: TYPE_RNA_BASE,
              TYPE_DNA: TYPE_DNA_BASE,
              TYPE_NONPOLYMER: TYPE_NONBASE}

class PDBReader(AriaBaseClass):

    def __init__(self, table = None):
        
        from aria.TypeChecking import check_type, NONE

        check_type(table, 'ConversionTable', NONE)

        AriaBaseClass.__init__(self)

        from aria.ConversionTable import ConversionTable
        #from Scientific.IO.PDB import PDBFile
        from aria.scientific.PDB import PDBFile

        if table is None:

            from aria.ConversionTable import ConversionTable

            table = ConversionTable()

        self.table = table
        self.PDB = PDBFile

        self.base_types = BASE_TYPES

    def get_energy(self, filename, type = 'total_energy'):

        import os

        file = self.PDB(os.path.expanduser(filename))
        eof = 0

        # BARDIAUX 2.2 = add 'restraint_energy', 'noe_violations', 'restraint_violations'
        criterion = {'total_energy' : 0.,
                     'restraint_energy' : 0.,
                     'noe_violations' : 0.,
                     'restraint_violations' : 0.}

        ln = 1
        while not eof:
            try:
                ttype, data = file.readLine()
                ln += 1
                
            except Exception, msg:
                import aria.tools as tools
                self.warning(tools.last_traceback())
                self.error(Exception, error = 'Could not read PDB-file %s, line %d' % (filename, ln), msg = msg)

            if ttype == 'END': eof = 1

            elif ttype == 'REMARK':
                
                if data.find('energies') >= 0:
                    
                    ## assumed here: special cns energy-header
                    criterion['total_energy'] = data[11:].split(',')[0]
                    rest_e = data[11:].split(',')[-5:]
                    criterion['restraint_energy'] = sum([float(e) for e in rest_e])
                    #return float(energy[type])
                
                if data.find('violations') >= 0:
                    criterion['noe_violations'] = data[13:].split(',')[0]
                    criterion['restraint_violations'] = sum([float(e) for e in data[13:].split(',')])

                return float(criterion[type])

        return None

    def read(self, filename, chain_types = None, format = CNS_CONVENTION):
        """
        Reads all (het)atoms in a PDB-file and returns a dict
        of dicts of the following hierachy:
          * segment level (keys are SEGIDs)
          * residue level (keys are residue numbers)
          * atom level (keys are atom names)
        Residue type and atom names in the returned dict will
        follow the IUPAC convention.
        """

        from aria.TypeChecking import check_string
        from aria.tools import string_to_segid
        from aria.OrderedDict import OrderedDict

        import os
        
        check_string(filename)

        if not self.table.has_format(format):
            self.error(ValueError, 'Atom naming "%s" not supported.' \
                       % str(format))

        if chain_types is not None:
            for segid, chain_type in chain_types.items():
                if not self.base_types.has_key(chain_type):
                    m = 'Type "%s" for chain "%s" not supported.'
                    self.error(TypeError, m % (str(chain_type), str(segid)))

        file = self.PDB(os.path.expanduser(filename))

        eof = 0

        dict = OrderedDict()
        dict['total_energy'] = None
        # BARDIAUX 2.2
        dict['restraint_energy'] = None
        dict['noe_violations'] = None
        dict['restraint_violations'] = None
        
        ln = 1
        while not eof:
            try:
                type, data = file.readLine()
                ln += 1
                
            except Exception, msg:
                import aria.tools as tools
                self.warning(tools.last_traceback())
                self.error(Exception, error = 'Could not read PDB-file %s, line %d' % (filename, ln), msg = msg)

            ## PDB file terminated

            if type == 'END': eof = 1

            ## try to grep energies

            elif type == 'REMARK':

                if data.find('energies') >= 0:

                    ## assumed here: special cns energy-header

                    energy = data[11:].split(',')[0]
                    dict['total_energy'] = float(energy)
                    # BARDIAUX 2.2
                    rest_e = data[11:].split(',')[-5:]
                    dict['restraint_energy'] = float(sum([float(e) for e in rest_e]))

                if data.find('violations') >= 0:
                    dict['noe_violations'] = float(data[13:].split(',')[0])
                    dict['restraint_violations'] = sum([float(e) for e in data[13:].split(',')])
                    
            ## read atom / hetatom information 

            elif type == 'ATOM' or type == 'HETATM':

                residue_type = data['residue_name'].strip()

                atom_name = data['name'].strip()
                atom_coor = data['position']#.array

                residue_number = data['residue_number']
                
                segid = string_to_segid(data['segment_id'])

                if not dict.has_key(segid):

                    dict[segid] = {}

                if not dict[segid].has_key(residue_number):

                    dict[segid][residue_number] = {}
                    dict[segid][residue_number]['residue_type'] = residue_type

                if dict[segid][residue_number].has_key(atom_name):

                    m = 'Atom "%s" ocurred twice in residue "%s%d".' \
                        % (atom_name, residue_type, residue_number)

                    self.error(ValueError, m)

                dict[segid][residue_number][atom_name] = atom_coor

        ## convert residue and atom names to IUPAC format

        if chain_types is None:
            chain_types = {}

        for segid in dict.keys():

            if not len(segid) == 4: continue

            if not chain_types.has_key(segid):
                chain_types[segid] = None
            else:
                base_type = self.base_types[chain_types[segid]]

            for number, residue in dict[segid].items():

                residue_type = residue['residue_type']

                if chain_types[segid] is None:
                    for chain_type in [TYPE_PROTEIN, TYPE_DNA, TYPE_RNA]:
                        base_type = self.base_types[chain_type]
                        args = (residue_type, format, IUPAC_CONVENTION,
                                base_type, 0)
                        if self.table.convert_residue(*args):
                            chain_types[segid] = chain_type
                            break
                    else:
                        chain_types[segid] = TYPE_NONPOLYMER
                        base_type = TYPE_NONBASE

                if base_type == TYPE_NONBASE: break

                atom_names = residue.keys()
                atom_names.remove('residue_type')
                atom_coors = [residue[x] for x in atom_names]
                for x in atom_names: del residue[x]

                args = (residue_type, atom_names, format, IUPAC_CONVENTION,
                        base_type, 0)
                new_names = self.table.convert_atoms(*args)
                if new_names is not None: atom_names = new_names
                
                for x, y in zip(atom_names, atom_coors): residue[x] = y

                args = (residue_type, format, IUPAC_CONVENTION, base_type, 0)

                new_type = self.table.convert_residue(*args)
                if new_type is not None: residue_type = new_type
                
                residue['residue_type'] = residue_type

            dict[segid]['chain_type'] = chain_types[segid]
                
        return dict

    def write(self, dict, filename, input_format = IUPAC_CONVENTION,
              output_format = IUPAC_CONVENTION):
        """
        Writes a dict of dicts of the hierarchy
        {<seg_id>:
         {<residue_number>:
          {'residue_type': <three-letter-code>,
           <atom-name>: <coordinates>}}}
        to a file in PDB format. 
        """

        table = self.table

        if not table.has_format(input_format):
            self.error(ValueError, 'Atom naming "%s" not supported.' \
                       % str(input_format))

        if not table.has_format(output_format):
            self.error(ValueError, 'Output format "%s" not supported.' \
                       % str(output_format))

        from aria.tools import string_to_segid

        import os

        file = self.PDB(os.path.expanduser(filename), 'w')

        serial_number = 1

        for segid in dict.keys():

            if len(segid) > 4: continue

            chain_type = dict[segid]['chain_type']
            if not self.base_types.has_key(chain_type):
                m = 'Chain "%s" has type "%s", this type is not supported.'
                self.error(TypeError, m % (str(segid), str(base_type)))

            base_type = self.base_types[dict[segid]['chain_type']]
            
            data = {}

            data['segment_id'] = string_to_segid(segid)

            residues = dict[segid]
            residue_numbers = residues.keys()
            residue_numbers.remove('chain_type')
            residue_numbers.sort()

            for residue_number in residue_numbers:

                residue = residues[residue_number]

                residue_type = residue['residue_type']

                if base_type is not TYPE_NONBASE:
                    args = (residue_type, input_format, output_format,
                            base_type, 0)
                    new_type = table.convert_residue(*args)
                    if new_type is not None: residue_type = new_type
                
                data['residue_name'] = residue_type
                data['residue_number'] = residue_number

                atom_names = residue.keys()
                atom_names.remove('residue_type')

                atom_coors = [residue[x] for x in atom_names]

                if base_type is not TYPE_NONBASE:
                    args = (residue['residue_type'], atom_names, input_format,
                            output_format, base_type, 0)
                    new_names = table.convert_atoms(*args)

                    if new_names is not None: atom_names = new_names
                
                for atom_name, atom_coor in zip(atom_names, atom_coors):

                    ## enforce correct alignment of atom names

                    atom_name = ' ' * (len(atom_name) < 4) + atom_name

                    data['name'] = atom_name
                    data['position'] = atom_coor
                    data['serial_number'] = serial_number

                    serial_number += 1

                    file.writeLine('ATOM', data)

        file.close()

