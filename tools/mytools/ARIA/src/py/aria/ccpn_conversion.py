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
from aria.Topology import TYPE_AMINO_ACID, TYPE_DNA_BASE, TYPE_RNA_BASE, \
     TERMINUS_N_STANDARD, TERMINUS_C_STANDARD, TERMINUS_N_AMINYL, \
     TERMINUS_C_CARBOXYL, TERMINUS_C5_PRIME_PHOSPHATE, \
     TERMINUS_C5_PRIME_HYDROXYL, TERMINUS_C3_PRIME_HYDROXYL, \
     TYPE_NONBASE
from aria.Chain import TYPE_PROTEIN, TYPE_RNA, TYPE_DNA, TYPE_NONPOLYMER


import sys, os

ccpn_path = os.path.join(AriaBaseClass.install_path, '3rd_party')
if ccpn_path not in sys.path: sys.path.insert(0, ccpn_path)

from ccpn.Global import false

import ccpn.Api as api

## conversion class

class CCPNConverter(AriaBaseClass):

    mol_types = {TYPE_AMINO_ACID: 'protein', TYPE_RNA_BASE: 'RNA',
                 TYPE_DNA_BASE: 'DNA', TYPE_PROTEIN: 'protein',
                 TYPE_RNA: 'RNA', TYPE_DNA: 'DNA', TYPE_NONBASE: 'nonpolymer',
                 TYPE_NONPOLYMER: 'nonpolymer'}

    termini = {'middle':
               {TYPE_AMINO_ACID: (TERMINUS_N_STANDARD, TERMINUS_C_STANDARD),
                TYPE_RNA_BASE: (TERMINUS_C5_PRIME_PHOSPHATE, None),
                TYPE_DNA_BASE: (TERMINUS_C5_PRIME_PHOSPHATE, None)},

               'start':
               {TYPE_AMINO_ACID: (TERMINUS_N_AMINYL, TERMINUS_C_STANDARD),
                TYPE_RNA_BASE: (TERMINUS_C5_PRIME_HYDROXYL, None),
                TYPE_DNA_BASE: (TERMINUS_C5_PRIME_HYDROXYL, None)},

               'end':
               {TYPE_AMINO_ACID: (TERMINUS_N_STANDARD, TERMINUS_C_CARBOXYL),
                TYPE_RNA_BASE: (TERMINUS_C5_PRIME_PHOSPHATE,
                                TERMINUS_C3_PRIME_HYDROXYL),
                TYPE_DNA_BASE: (TERMINUS_C5_PRIME_PHOSPHATE,
                                TERMINUS_C3_PRIME_HYDROXYL)},

               'any':
               {TYPE_NONBASE: (None, None)}} 

    def _create_chem_comp_locations(self, project, base_residue, location):

        type = self.mol_types[base_residue.getType()]
        code = base_residue.getName()

        if not CCPNConverter.termini[location].has_key(base_residue.getType()):
            return

        chem_comp_loc = api.ChemCompLoc(project ,location = location,
                                        molType = type, code3Letter = code,
                                        code = code)

        ## create chem atom locations

        base_atoms = list(base_residue.getBackboneAtoms())
        base_atoms += list(base_residue.getSidechainAtoms())

        n, c = CCPNConverter.termini[location][base_residue.getType()]

        termini = base_residue.getTermini()

        if n in termini:
            base_atoms += list(termini[n].getAtoms())
        if c in termini:
            base_atoms += list(termini[c].getAtoms())

        for base_atom in base_atoms:

            name = base_atom.getName()
            symbol = base_atom.getType()
            if not symbol: symbol = 'None'
            api.ChemAtomLoc(chem_comp_loc, chemAtomName = name,
                            name = name, elementSymbol = symbol)

    def _create_chem_comp(self, project, base_residue):

        type = self.mol_types[base_residue.getType()]
        code = base_residue.getName()

        ## TODO: hack, we will just load CCPN reference data later...

        chem_comp = api.NonStdChemComp(project, molType = type,
                                       code3Letter = code, code = code)

        for atom in base_residue.getAtoms():

            name = atom.getName()
            symbol = atom.getType()

            if not symbol: symbol = 'None'
            
            api.ChemAtom(chem_comp, name = name, elementSymbol = symbol)
            
        ## create bonds

        for atom in base_residue.getAtoms():

            hetero = atom.getHeteroAtom()
            if hetero is None: continue
            
            atom = chem_comp.findFirstChemAtom(name = atom.getName())
            hetero = chem_comp.findFirstChemAtom(name = hetero.getName())

            api.ChemBond(chem_comp, chemAtoms = (atom, hetero))

        ## create chem comp locations

        for location in ('start', 'middle', 'end', 'any'):
            self._create_chem_comp_locations(project, base_residue, location)

    def convert_topology(self, topology, project = None):

        if project is None:
            project = api.Project(name = topology.getName())
            
        [self._create_chem_comp(project, r) for r in topology.getResidues()]

        return project

    def convert_molecule(self, molecule, project):

        for chain in molecule.get_chains():
            self.convert_chain(chain, project)

    def convert_chain(self, chain, project):

        from aria.Topology import EQUIV_NTERMINUS

        polymer = api.Polymer(project, name = chain.getSegid())

        mol_type = self.mol_types[chain.getType()]

        for residue in chain.getResidues():

            location = 'middle'
            groups = residue.getEquivalentGroups()
            if EQUIV_NTERMINUS in [g.getType() for g in groups]:
                location = 'start'

            if project.findFirstChemComp(molType = mol_type,
                                         code3Letter = residue.getType()) \
                is None:
                
                from aria.Topology import BaseResidue, BaseResidueSettings, \
                     BaseAtom, BaseAtomSettings, TYPE_NONBASE

                base_residue = BaseResidue(BaseResidueSettings(),
                                           residue.getType())
                base_residue.getSettings()['type'] = TYPE_NONBASE

                base_atoms = [BaseAtom(BaseAtomSettings(), a.getName())
                              for a in residue.getAtoms()]
                [b.getSettings().update(a.getSettings())
                 for a, b in zip(residue.getAtoms(), base_atoms)]
                base_residue.setBackboneAtoms(base_atoms)
                base_residue.setSidechainAtoms([])
                self._create_chem_comp(project, base_residue)

                location = 'any'

            ## TODO: get chain-location right

            api.MolResidue(polymer, seqID = residue.getNumber(),
                           molType = mol_type, code3Letter = residue.getType(),
                           chainLocation = location)

        ## create a chain and check atom content

        molsystem = api.MolSystem(project, code = chain.getSegid())
        ccpn_chain = api.NmrChain(molsystem, code = chain.getSegid(),
                                  moleculeName = chain.getSegid())
                                           
        for aria_residue in chain.getResidues():

            number = aria_residue.getNumber()
            residue = ccpn_chain.findFirstResidue(seqID = number)

            for atom in residue.getAtoms():
                if not aria_residue.hasAtom(atom.getName()):
                    atom.delete()

    def convert_non_floating_assignment(self, shift_assignment,
                                        shift_list, atom_mapping):

        project = shift_list.getProject()

        spin_system = shift_assignment.getSpinSystems()[0]

        atoms = spin_system.getAtoms()
        chemical_shift = spin_system.getChemicalShifts()[0]

        if chemical_shift[0] == None: return

        resonance = api.Resonance(project)
        shift = api.Shift(shift_list, value = chemical_shift[0],
                          resonance = resonance)

        ccpn_atoms = [atom_mapping[atom] for atom in atoms]
        atom_set = api.AtomSet(project, atoms = ccpn_atoms)

        api.ResonanceSet(project, atomSets = (atom_set, ),
                         resonances = (resonance, ))
        
    def convert_floating_assignment(self, shift_assignment, shift_list,
                                    atom_mapping):

        project = shift_list.getProject()

        ## create atom-sets

        atom_sets = []
        for system in shift_assignment.getSpinSystems():
            atoms = [atom_mapping[a] for a in system.getAtoms()]
            atom_sets.append(api.AtomSet(project, atoms = atoms))

        ## create resonances and shifts

        resonances = []
        for s in shift_assignment.getSpinSystems()[0].getChemicalShifts():

            if s[0] is None: continue
            resonance = api.Resonance(project)
            shift = api.Shift(shift_list, value = s[0], resonance = resonance)

            resonances.append(resonance)

        if len(resonances) == 0:
            m = 'Floating assignment without any chemical-shift: ' + \
                str(shift_assignment)
            self.error(ValueError, m)
            
        api.ResonanceSet(project, atomSets = atom_sets,
                         resonances = resonances)

    def convert_shiftlist(self, aria_shift_list, project):

        from aria.ShiftAssignment import ASSIGNMENT_METHOD_FLOATING, \
             ASSIGNMENT_METHOD_STEREO_SPECIFIC, ASSIGNMENT_METHOD_EQUIVALENT

        shift_list = api.ShiftList(project)

        ## TODO: sort of clumsy

        molecule = aria_shift_list.getShiftAssignments()[0].\
                   getSpinSystems()[0].getAtoms()[0].getResidue().\
                   getChain().getMolecule()
        
        atom_mapping = map_aria_atoms(molecule, project)
        
        for assignment in aria_shift_list.getShiftAssignments():

            if assignment.getMethod() == ASSIGNMENT_METHOD_STEREO_SPECIFIC \
            or assignment.getMethod() == ASSIGNMENT_METHOD_EQUIVALENT :

                self.convert_non_floating_assignment(assignment, shift_list,
                                                     atom_mapping)

            elif assignment.getMethod() == ASSIGNMENT_METHOD_FLOATING:

                self.convert_floating_assignment(assignment, shift_list,
                                                 atom_mapping)

    # get method if it exists, otherwise create

    def getProjectMethod(self, project, methodType, methodName):

        if (not project):
            return None

        procedure = 'Aria_CCPN_' + methodType
        parameters = 'method=' + methodName

        method = project.findFirstMethod(procedure = procedure,
                                         parameters = parameters)
        if (not method):
            method = api.Method(project, procedure=procedure,
                                 parameters=parameters)
    
        return method

    def convert_noesy_spectrum(self, noesy_spectrum, project):

        from numpy import sum, flatnonzero

        # TODO: this should be done only once

        heightMethod = self.getProjectMethod(project, methodType='height',
                                             methodName='manual')
        volumeMethod = self.getProjectMethod(project, methodType='volume',
                                             methodName='manual')

        mapping = map_dimensions(noesy_spectrum)

        n_dimensions = sum(mapping)
        spectrum_name = noesy_spectrum.getName()
        if spectrum_name is None:
            spectrum_name = 'noesy'

        experiment = api.Experiment(project, name = spectrum_name,
                                    numDim = n_dimensions)

        datasource = api.DataSource(experiment, name = spectrum_name,
                                    numDim = n_dimensions,
                                    dataType = 'processed')

        ## create dimensions of the spectrum

        data_dim_refs = []

        for i in range(1, n_dimensions + 1):

            exp_dim = experiment.findFirstExpDim(dim = i)

            ## sf is arbitrary, need later for conversion of
            ## peak dim positions to and from points and ppm
            exp_dim_ref = api.ExpDimRef(exp_dim, sf = 1.0, unit = 'ppm')

            ## TODO: dummy values for numPointsOrig, etc.

            ## numPointsOrig, valuePerPoint and numPoints are
            ## arbitrary, see comment about sf above            
            freq_dim = api.FreqDataDim(datasource, numPointsOrig = 512,
                                       valuePerPoint = 1.0, dim = i,
                                       numPoints = 512, isComplex = false,
                                       expDim = exp_dim)
            
            data_dim_ref = api.DataDimRef(freq_dim, expDimRef = exp_dim_ref)

            data_dim_refs.append(data_dim_ref)

        ## create peak-list

        peaklist = api.PeakList(datasource)

        get_shifts = {0: 'Proton1', 1: 'Proton2', 2: 'Hetero1', 3: 'Hetero2'}

        for cross_peak in noesy_spectrum.getPeaks():

            shifts = []
            assignments = []
            
            for i in flatnonzero(mapping):

                shift = getattr(cross_peak, 'get%sChemicalShift' \
                                % get_shifts[i])()
                shifts.append(shift)

                assignment = getattr(cross_peak, 'get%sAssignments' \
                                     % get_shifts[i])()
                assignments.append(assignment)

            ## skip if not all positions are defined

            if None in [s[0] for s in shifts]: continue

            peak = api.Peak(peaklist)

            volume = api.PeakIntensity(peak, type = 'volume',
                                       methodSerial = volumeMethod.serial)
            
            if cross_peak.getVolume()[0]:
                volume.value = cross_peak.getVolume()[0]
                
            if cross_peak.getVolume()[1]:
                volume.error = cross_peak.getVolume()[1]

            intensity = api.PeakIntensity(peak, type = 'height',
                                          methodSerial = heightMethod.serial)
                                              
            if cross_peak.getIntensity()[0]:
                intensity.value = cross_peak.getIntensity()[0]

            if cross_peak.getIntensity()[1]:
                intensity.error = cross_peak.getIntensity()[1]

            ## set peak position

            for i in range(len(shifts)):

                peakdim = peak.findFirstPeakDim(dim = i + 1)

                if shifts[i][0] is not None:
                    peakdim.value = shifts[i][0]

                if shifts[i][1] is not None:
                    peakdim.value = shifts[i][1]

            ## store peak assignments

            for i in range(len(assignments)):

                if assignments[i] is None: continue

                peakdim = peak.findFirstPeakDim(dim = i + 1)

                atom_sets = []

                for atoms in assignments[i]:

                    ccpn_atoms = [find_ccpn_atom(atom, project) for atom in
                                  atoms]
                    atom_set = self.find_atom_set(ccpn_atoms)
                    if atom_set is None:
                       atom_set = api.AtomSet(project, atoms = ccpn_atoms)
                    atom_sets.append(atom_set)

                    resonance = self.find_resonance(atom_set)
                    if resonance is None:
                        resonance = self.create_resonance(project, atom_set)

                    # Note that in CCPN API a PeakContrib will automatically
                    # be created which these PeakDimContribs will link to.
                    # In future might want several PeakContribs...
                    api.PeakDimContrib(peakdim, resonance = resonance) 

    def find_resonance(self, atom_set):

        ## TODO: should use dict perhaps keyed off Aria atom sets or ...
        return None

    def create_resonance(self, project, atom_set):

        ## TODO: should do properly
        
        resonance = api.Resonance(project)
        resonance_set = api.ResonanceSet(project, resonances=(resonance,),
                                         atomSets=(atom_set,))

        return resonance
    
    def find_atom_set(self, ccpn_atoms):

        if not ccpn_atoms: return None 

        atom_set = ccpn_atoms[0].atomSet

        if atom_set is None: return None

        if self.have_atom_set_match(atom_set, ccpn_atoms):
            return atom_set

        return None

    def have_atom_set_match(self, atom_set, ccpn_atoms):

        if len(atom_set.atoms) != len(ccpn_atoms):
            return 0

        for atom in atom_set.atoms:
            if atom not in ccpn_atoms:
                return 0

        return 1
    
    def convert_experiment(self, aria_exp, project):

        noesy_spectrum = aria_exp.getSpectrum()
        shift_list = aria_exp.getShiftList()

        self.convert_noesy_spectrum(noesy_spectrum, project)
        self.convert_shiftlist(shift_list, project)

## auxiliary functions
        
def map_aria_atoms(molecule, project):
    """
    Create a mapping between all aria atoms in a molecule and the
    chem atoms found in the CCPN project.
    """

    atom_mapping = {}

    for chain in molecule.get_chains():

        molsystem = project.findFirstMolSystem(code = chain.getSegid())
        ccpn_chain = molsystem.findFirstChain(code = chain.getSegid())

        for residue in chain.getResidues():

            code = residue.getType()
            number = residue.getNumber()

            ccpn_residue = ccpn_chain.findFirstResidue(seqID = number)

            if not ccpn_residue.getMolResidue().getCode3Letter() == code:

                raise ValueError, 'Inconsistent chain.'

            for atom in residue.getAtoms():

                ccpn_atom = ccpn_residue.findFirstAtom(name = atom.getName())

                atom_mapping[atom] = ccpn_atom

    return atom_mapping

def find_ccpn_atom(atom, project):

    residue = atom.getResidue()

    chain = residue.getChain()
    molsystem = project.findFirstMolSystem(code = chain.getSegid())
    ccpn_chain = molsystem.findFirstChain(code = chain.getSegid())
    residue = ccpn_chain.findFirstResidue(seqID = residue.getNumber())

    if residue is None:
        raise ValueError, 'Cannot find residue in CCPN chain.'

    return residue.findFirstAtom(name = atom.getName())

def map_dimensions(spectrum):

    from numpy import array, compress, not_equal, sum, equal

    dimensions = [(peak.getProton1ChemicalShift()[0],
                   peak.getProton2ChemicalShift()[0],
                   peak.getHetero1ChemicalShift()[0],
                   peak.getHetero2ChemicalShift()[0])
                  for peak in spectrum.getPeaks()]

    dimensions = array(dimensions)

    ## throw-out invalid measurements

    dimensions = compress(not_equal(dimensions[:,0], None), dimensions, 0)
    dimensions = compress(not_equal(dimensions[:,1], None), dimensions, 0)

    mapping = [1, 1, 1, 1]

    if sum(equal(dimensions[:,2], None)) == len(dimensions):
        mapping[2] = 0
    
    if sum(equal(dimensions[:,3], None)) == len(dimensions):
        mapping[3] = 0
    
    return mapping

