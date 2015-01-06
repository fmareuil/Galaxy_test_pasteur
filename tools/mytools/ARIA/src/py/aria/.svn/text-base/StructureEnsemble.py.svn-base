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
from aria.Settings import NonNegativeInt as _NonNegativeInt
from aria.Settings import Settings as _Settings
from aria.TypeChecking import *
from aria.ConversionTable import CNS_CONVENTION, IUPAC_CONVENTION, \
     DYANA_CONVENTION

class NumberOfBestStructures(_NonNegativeInt):

    def __init__(self, description = None, error_message = None):

        if error_message is None:
            error_message = 'Number of best structures must be a ' + \
                            'non-negative integer or "all"; %s given.'

        _NonNegativeInt.__init__(self, description, error_message)

    def getErrorMessage(self, value):
        return _NonNegativeInt.getErrorMessage(self, value) % str(value)
    
    def is_valid(self, value):
        return _NonNegativeInt.is_valid(self, value) or value == 'all'

class StructureEnsembleSettings(_Settings):

    ## (DONE (BARDIAUX)) TODO: allow 'restraint_energy' to be sorting criterion 
    
    def create(self):

        from aria.Settings import ChoiceEntity

        choices = ('total_energy', 'restraint_energy', 'noe_violations', 'restraint_violations',) # BARDIAUX 2.2
        msg = 'Structures can be sorted according to: ' + \
              ' \ '.join(choices)
        keywords = {'sort_criterion': ChoiceEntity(choices, msg),
                    'number_of_best_structures': NumberOfBestStructures()}

        return keywords

class StructureEnsemble(AriaBaseClass):

    def __init__(self, settings):

        check_type(settings, 'StructureEnsembleSettings')

        from aria.PDBReader import PDBReader

        AriaBaseClass.__init__(self)

        self.setSettings(settings)

        self.__X = None
        self.__unknown_atoms = []
        self.__info = None

        self.reader = PDBReader()

        names = ('sort_criterion', 'number_of_best_structures')
        
        for name in names:
            entity = settings.getEntity(name)
            entity.set_callback(lambda entity, s = self, n = name: \
                                s.entity_has_changed(entity, n))
        
        self.__invalidate_cache()

    def entity_has_changed(self, entity, name):

        if name in ('sort_criterion', 'number_of_best_structures'):
            self.__invalidate_cache()

            if name == 'sort_criterion':

                self.debug('Resorting...')
                self.sort()

    def __invalidate_cache(self):
        self.__cache = {'hit': 1, 'total':1}

    def settingsChanged(self):
        self.sort()

    def sort(self):

        from numpy import argsort, take

        criterion = self.getSettings().getEntity('sort_criterion')

        if self.__info is None or not criterion.is_initialized():
            return

        energies = [d[criterion.get()] for d in self.__info[:,1]]
        indices = argsort(energies)

        self.__info = take(self.__info, indices, 0)
        self.__X = take(self.__X, indices, 0)

    def read(self, files, molecule, format = CNS_CONVENTION,
             float_files = None):
        """
        Reads a list of PDB-files and stores the coordinates in an
        array where the outer index corresponds to the atom-id.
        'files' may also be a dict. If so, keys are filenames,
        values are format-stringa.
        """

        check_type(files, LIST, TUPLE, DICT)
        check_elements(files, STRING)        
        check_type(molecule, 'Molecule')
        check_type(float_files, NONE, LIST, TUPLE)
        check_string(format)

        if float_files and len(files) <> len(float_files):
            m = 'Number of structure- and float-files must be equal'
            self.error(ValueError, m)

        from aria.Topology import EQUIV_METHYL
        from aria.Chain import TYPE_NONPOLYMER
        from aria.FloatFile import FloatFile
        from numpy import zeros, array
        
        if type(files) == type({}):
            format = files.values()
            files = files.keys()
        else:
            format = [format] * len(files)

        self.message('Reading PDB files ...')

        ## store all atoms in some convenient way

        atoms = [a for c in molecule.get_chains() for r in c.getResidues()
                 for a in r.getAtoms()]
        atoms.sort(lambda a, b: cmp(a.getId(), b.getId()))

        atom_dict = {}

        for atom in atoms:

            segid = atom.getSegid()
            name = atom.getName()
            number = atom.getResidue().getNumber()

            atom_dict[(segid, name, number)] = atom

        n = len(atoms)
        m = len(files)

        ## array for storing the cartesian coordinates

        X = zeros((m, n, 3), 'd')
        info = zeros((m, 2), 'O')

        ## for storing the atoms that have no coordinates in the PDB-file
        
        unknown_atoms = []

        chain_types = {}
        for chain in molecule.get_chains():
            chain_types[chain.getSegid()] = chain.getSettings()['type']

        for j in range(m):

            pdb_dict = self.reader.read(files[j], chain_types, format[j])

            if float_files:
                swapped_atoms = float_files[j]
            else:
                swapped_atoms = {}

            ## store structure-specific data

            d = {'total_energy': pdb_dict['total_energy'],
                 'restraint_energy': pdb_dict['restraint_energy'],
                 'noe_violations': pdb_dict['noe_violations'],
                 'restraint_violations': pdb_dict['restraint_violations'],
                 'index': j}

            info[j] = files[j], d

            for i in range(n):

                atom = atoms[i]

                name = atom.getName()
                num = atom.getResidue().getNumber()
                segid = atom.getSegid()

                try:
                    X[j][i][:] = pdb_dict[segid][num][name][:]

                except:

                    if i not in unknown_atoms: unknown_atoms.append(i)

                    s = 'Atom not found in structure %d (%s): [segid=%s' + \
                        ', residue_number=%s, name=%s].'

                    self.warning(s %(j, info[j][0], segid, num, name))

            ## swap atoms according to float file

            for key, value in swapped_atoms.items():

                segid1, number1, name1 = key
                segid2, number2, name2 = value

                chain1 = molecule.getChain(segid1)
                chain2 = molecule.getChain(segid2)

                if chain1.getType() == TYPE_NONPOLYMER or \
                   chain2.getType() == TYPE_NONPOLYMER:
                    continue

                residue1 = chain1[number1]
                residue2 = chain2[number2]

                base_type1 = self.reader.base_types[chain1.getType()]
                base_type2 = self.reader.base_types[chain2.getType()]

                ## convert atom names

                args = (residue1.getType(), name1, CNS_CONVENTION, \
                        IUPAC_CONVENTION, base_type1)

                name1 = self.reader.table.convert_atom(*args)

                args = (residue2.getType(), name2, CNS_CONVENTION, \
                        IUPAC_CONVENTION, base_type2)

                name2 = self.reader.table.convert_atom(*args)

                ## handle residues with floating methyl groups

                ## first residue

                if residue1.getType() in ['VAL', 'LEU', 'ILE']:

                    methyls = [g for g in residue1.getEquivalentGroups()
                               if g.getType() == EQUIV_METHYL]

                    methyl = [g for g in methyls if name1 in g.getAtomNames()]

                    if len(methyl) == 0:

                        names1 = (name1,)  

                    elif len(methyl) == 1:
                        
                        methyl1 = methyl[0]

                        names1 = list(methyl1.getAtomNames())
                        names1.sort()
                        names1 = tuple(names1)
                    
                    elif len(methyl) > 1:

                        m = 'Inconsistency in float file: could not ' + \
                            'identify methyl group for atom "%s" in ' + \
                            'residue "%s".'

                        self.warning(m % (name1, residue1.getName()))

                        continue

                else:

                    names1 = (name1,)

                ## second residue

                if residue2.getType() in ['VAL', 'LEU', 'ILE']:

                    methyls = [g for g in residue2.getEquivalentGroups()
                               if g.getType() == EQUIV_METHYL]

                    methyl = [g for g in methyls if name2 in g.getAtomNames()]

                    if len(methyl) == 0:

                        names2 = (name2,)

                    elif len(methyl) == 1:
                              
                        methyl2 = methyl[0]

                        names2 = list(methyl2.getAtomNames())
                        names2.sort()
                        names2 = tuple(names2)

                    elif len(methyl) > 1:

                        m = 'Inconsistency in float file: could not ' + \
                            'identify methyl group for atom "%s" in ' + \
                            'residue "%s".'

                        self.warning(m % (name2, residue2.getName()))

                        continue

                else:

                    names2 = (name2,)

                ## should not occur

                if names1 == names2:
                    
                    names1 = (name1,)
                    names2 = (name2,)

                for name1 in names1:
                    for name2 in names2:

                        atom1 = residue1[name1].getId()
                        atom2 = residue2[name2].getId()

                        if atom1 in unknown_atoms or atom2 in unknown_atoms:

                            m = 'Try to swap unknown atoms %s and %s.' 

                            self.warning(m % (str(residue1[name1]),
                                              str(residue2[name2])))

                            continue

                        X[j,atom1], X[j,atom2] = array(X[j,atom2]), \
                                                 array(X[j,atom1])

        unknown_atoms.sort()

        self.__X = X
        self.__info = info
        self.__unknown_atoms = unknown_atoms

        self.sort()

        ## invalidate cache

        self.__invalidate_cache()

        self.message('PDB files read.')

    def getFiles(self):
        """
        getFiles(self)
        returns a tuple containing the names of pdb-files
        currently loaded into the structure-ensemble.
        files are sorted according to the settings
        'sort_criterion'
        """
        
        if self.__info is None:
            return None

        return tuple(self.__info[:,0])

    def getDistances(self, atom1, atom2):

        check_type(atom1, "Atom")
        check_type(atom2, "Atom")

        if AriaBaseClass.cache:
            key = [id(atom1), id(atom2)]
            key.sort()
            key = tuple(key)
            if key in self.__cache:
                self.__cache['hit'] += 1
                self.__cache['total'] += 1
                return self.__cache[key]

        from numpy import sqrt, sum, power

        id1 = atom1.getId()
        id2 = atom2.getId()

        if self.__X is None:        
            self.error(ValueError, 'No coordinates have been read in.')

        ## check if the coordinates for the atoms are stored
        ## in the ensemble

        unknown = self.__unknown_atoms

        m = 'Atom "%s" not found in residue "%s".'

        if id1 in unknown:
            
            self.error(ValueError, m % (atom1.getName(),
                                        atom1.getResidue().getName()))

        elif id2 in unknown:

            self.error(ValueError, m % (atom2.getName(),
                                        atom2.getResidue().getName()))
            
        ## calculate distances

        x = self.__X[:,id1,:]
        y = self.__X[:,id2,:]

        d = sqrt(sum(power(x - y, 2), 1))

        n = self.getSettings()['number_of_best_structures']

        if n <> 'all':
            d = d[:n]
            
        if AriaBaseClass.cache:
            self.__cache[key] = d
            self.__cache['total'] += 1

        return d

    def __len__(self):

        if not self.__X:
            return 0
        else:
            return self.__X.shape[0]
