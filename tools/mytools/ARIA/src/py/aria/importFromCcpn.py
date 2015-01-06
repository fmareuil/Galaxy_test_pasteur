from aria.Molecule import Molecule
from aria.Chain import TYPE_PROTEIN, TYPE_DNA, TYPE_RNA, TYPE_NONPOLYMER, Chain, ChainSettings
from aria.Residue import Residue
from aria.Atom import Atom, AtomSettings
from aria.ccpn2top import setupEquivalentGroup


from aria.Assignment import Assignment, ASSIGNMENT_TYPE_MANUAL
from aria.NOESYSpectrum import NOESYSpectrum
from aria.CrossPeak import CrossPeak
from aria.Datum import ChemicalShift, Datum
from aria.ChemicalShiftList import ChemicalShiftList

from aria.ShiftAssignment import AVERAGING_METHOD_FAST, AVERAGING_METHOD_SLOW, AVERAGING_METHOD_NONE, \
                            ASSIGNMENT_METHOD_STEREO_SPECIFIC, ASSIGNMENT_METHOD_EQUIVALENT, \
                            ASSIGNMENT_METHOD_FLOATING, SpinSystem, ShiftAssignment

from aria.tools import string_to_segid
from ccpnmr.format.general.Constants import ccpNmr_kw, isAriaInput_kw

from memops.api import Implementation
from memops.universal.Constants import True

from aria.Settings import Settings
from aria.xmlutils import XMLBasePickler

chainTypeMapping = {'protein'   : TYPE_PROTEIN,
                    'DNA'       : TYPE_DNA,
                    'RNA'       : TYPE_RNA,
                    'other': TYPE_NONPOLYMER}

## BARDIAUX BB log
from aria.ariabase import AriaBaseClass
messager = AriaBaseClass()
messager._name = 'CCPN import'
##
## BARDIAUX BB SecStruct
prevRes = None
ssDist = {}
##

def getKeysFromString(word, delimiter='|'):
    """Descrn: Get a list of CCPN data model keys given an object identifier string.
       Inputs: String
       Output: List of Keys (Words or Ints)
    """

    items = word.split(delimiter)

    keys = []
    for item in items:

        try:
            key = int(item)
        except:
            key = item

        keys.append(key)

    return keys


def getObjectKeys(object):
    """Descrn: Get a list of CCPN data model keys for a CCPN data model object.
       Inputs: CCPN data model object
       Output: List of Keys (Words or Ints or a further list of keys)
    """

    return object.getExpandedKey()


def getObjectKeyString(object, delimiter='|'):
    """Descrn: Make an object identifier string for a CCPN data model object
       Inputs: CCPN data model object
       Output: String
    """

    keys = object.getExpandedKey()

    for i in range(len(keys)):
        key = keys[i]

        keyType = type(key)
        if keyType is type([]):
            keys[i] = delimiter.join([str(k) for k in key])
        elif keyType is not type(''):
            keys[i] = str(key)

    return delimiter.join(keys)

# BARDIAUX 2.2
def getObjectsKeyString(o):

    s = []
    for a in o:
        s.append(getKeysFromString(getObjectKeyString(a)))

    base = s[0]

    for i in zip(*s):
        for a in i[1:]:
            if i[0] <> a:
                base.append(a)

    return "|".join(base)

def getShiftLists(project):
    """Descrn: Get a list of shift lists from a CCPN project
       Inputs: Implementation.Project
       Output: List of ccp.nmr.Nmr.ShiftLists
    """

    nmrProject = project.currentNmrProject

    shiftLists = []
    for measurementList in nmrProject.sortedMeasurementLists():
        if measurementList.className == 'ShiftList':
            shiftLists.append(measurementList)

    return shiftLists

# BARDIAUX
def getModels(project, ccpChains):
    """Descrn: Get a list of models from a CCPN project that contain ccpChains
       Inputs: Implementation.Project, list of ccp.api.molecule.MolStructure.Chain
       Output: List of ccp.api.molecule.MolStructure.Model 
    """
    
    nmrProject = project.currentNmrProject
    molSystem = ccpChains[0].molSystem
    
    models = []

    ccpChains.sort()
    
    for structureEnsemble in molSystem.sortedStructureEnsembles():
        chains = [molChain.chain for molChain in structureEnsemble.coordChains]
        chains.sort()
        if ccpChains == chains:
            for model in structureEnsemble.sortedModels():
                models.append(model)

    return models

# BARDIAUX
def getStructureEnsembles(project, ccpChains):
    """Descrn: Get a list of from a CCPN project
       Inputs: Implementation.Project, list of ccp.api.molecule.MolStructure.Chain
       Output: List of ccp.api.molecule.MolStructure.StructureEnsemble
       """

    nmrProject = project.currentNmrProject
    molSystem = ccpChains[0].molSystem

    ensembles =  []

    ccpChains.sort()    

    for structureEnsemble in molSystem.sortedStructureEnsembles():
        chains = [molChain.chain for molChain in structureEnsemble.coordChains]
        chains.sort()
        if ccpChains == chains:
            if structureEnsemble.models:
                ensembles.append(structureEnsemble)
            
    return ensembles

def getNoesyPeakLists(project, molSystem=None):
    """Descrn: Get the NOE peak lists from a CCPN project. Can filter if appropriate to a given
               molecular system if passed in.
       Inputs: Implementation.Project, ccp.molecule.MolSystem.MolSystem
       Output: List of ccp.nmr.Nmr.Peaks
    """

    peakLists = []
    for experiment in project.currentNmrProject.experiments:
        if experiment.molSystems:
            if molSystem and (molSystem not in experiment.molSystems):
                continue

        if experiment.refExperiment:
            if experiment.refExperiment.nmrExpPrototype.name.find('NOESY') < 0:
                continue
        else:
            transfer = experiment.findFirstExpTransfer(transferType='NOESY')
            if not transfer:
                continue

        for spectrum in experiment.dataSources:
            if (spectrum.dataType == 'processed') and (spectrum.numDim > 1):

                isotopes = []
                for dataDim in spectrum.dataDims:
                    for expDimRef in dataDim.expDim.expDimRefs:
                        if expDimRef.measurementType in ('shift','Shift'):
                            isotope = ','.join(expDimRef.isotopeCodes)
                            isotopes.append(isotope)
                            break

                if isotopes.count('1H') > 1:
                    for peakList in spectrum.peakLists:
                        if peakList.findFirstPeak():
                            peakLists.append(peakList)


    return peakLists


def getCcpnPeakList(ccpnProject, peakListKey):
    """Descrn: Fetch a CCPN peak list from a CCPN project using a list of object keys
       Inputs: Implementation.Project, List of keys (Words or Ints)
       Output: ccp.nmr.Nmr.PeakList
    """

    nmrProjectName, exptSerial, specSerial, plSerial = peakListKey

    nmrProject = ccpnProject.findFirstNmrProject(name=nmrProjectName)
    if not nmrProject:
        ValueError('No NMR project with name %s' % nmrProjectName)

    experiment = nmrProject.findFirstExperiment(serial=exptSerial)
    if not experiment:
        ValueError('Could not find experiment %d in CCPN project.' % (exptSerial))
        return

    spectrum = experiment.findFirstDataSource(serial=specSerial)
    if not spectrum:
        ValueError('Could not find spectrum %d in experiment %d.' % (specSerial, exptSerial))
        return

    peakList = spectrum.findFirstPeakList(serial=plSerial)
    if not peakList:
        ValueError('Could not find peak list serial %d in spectrum %d:%d' % (plSerial, exptSerial, specSerial))

    return peakList


def getCcpnShiftList(ccpnProject, shiftListKey):
    """Descrn: Fetch a CCPN shift list from a CCPN project using a list of object keys
       Inputs: Implementation.Project, List of keys (Words or Ints)
       Output: ccp.nmr.Nmr.ShiftList
    """

    nmrProjectName, shiftListSerial = shiftListKey

    nmrProject = ccpnProject.findFirstNmrProject(name=nmrProjectName)
    if not nmrProject:
        ValueError('No NMR project with name %s' % nmrProjectName)

    shiftList = nmrProject.findFirstMeasurementList(serial=shiftListSerial)
    if not shiftList:
        ValueError('Could not find shift list %d in CCPN project' % (int(shiftListSerial)))

    return shiftList


def getCcpnExperimentShiftList(ccpnExperiment):
    """Descrn: Fetch the CCPN shift list that corresponds to a CCPN NMR experiment
       Inputs: ccp.nmr.Nmr.Experiment
       Output: ccp.nmr.Nmr.ShiftList
    """

    shiftList = ccpnExperiment.shiftList

    if not shiftList:
        for spectrum in ccpnExperiment.dataSources:
            for peakList in spectrum.peakLists:
                for peak in peakList.peaks:
                    for peakDim in peak.peakDims:
                        for contrib in peakDim.peakDimContribs:
                            shift = contrib.resonance.findFirstShift()
                            if shift:
                                return shift.parentList

    return shiftList


def getCcpnChain(ccpProject, chainKey):
    """Descrn: Fetch a CCPN chain from a CCPN project using a list of object keys
       Inputs: Implementation.Project, List of keys (Words or Ints)
       Output: ccp.molecule.MolSystem.Chan
    """

    molSystemCode, chainCode = chainKey

    molSystem = ccpProject.findFirstMolSystem(code=molSystemCode)
    if not molSystem:
        ValueError('No molecular system with code "%s" in CCPN project' % molSystemCode)

    chain = molSystem.findFirstChain(code=chainCode)
    if not chain:
        ValueError('No chain found with code "%s" in molecular system "%s"' % (chainCode, molSystemCode))

    return chain

# BARDIAUX
def getCcpnChains(ccpProject, chainKey):
    """Descrn: Fetch a CCPN chain from a CCPN project using a list of object keys
       Inputs: Implementation.Project, List of keys (Words or Ints)
       Output: list of ccp.molecule.MolSystem.Chains
    """

    molSystemCode, chainCodes = chainKey[0], chainKey[1:]

    molSystem = ccpProject.findFirstMolSystem(code=molSystemCode)
    if not molSystem:
        ValueError('No molecular system with code "%s" in CCPN project' % molSystemCode)

    chains = []
    for chainCode in chainCodes:
        chain = molSystem.findFirstChain(code=chainCode)
        if not chain:
            ValueError('No chain found with code "%s" in molecular system "%s"' % (chainCode, molSystemCode))

        chains.append(chain)

    return chains

def getCcpnPeakAndShiftLists(ccpProject, molSystem, peakShiftListKeys):
    """Descrn: Fetch a list of peak list and corresponding shift list pairs for a CCPN project
               given pairs of corresponding object key lists
       Inputs: Implementation.Project, List of 2-List of keys (Words or Ints)
       Output: List of (ccp.nmr.Nmr.PeakList, ccp.nmr.Nmr.ShiftList)
    """

    peakAndShiftLists = []

    for peakListKey, shiftListKey in peakShiftListKeys:
        shiftlist = None
        peakList  = getCcpnPeakList(ccpProject, peakListKey)

        if peakList:
            experiment = peakList.dataSource.experiment
            if (experiment.molSystems) and (molSystem not in experiment.molSystems):
                data = (molSystem.code, experiment.serial, experiment.name)
                messager.warning('WARNING: Selected molecular system "%s" does not match experiment %d (%s) molecular systems' % data)

            shiftList  = getCcpnExperimentShiftList(experiment)

        shiftList0 = getCcpnShiftList(ccpProject, shiftListKey)

        if shiftList0:
            shiftList = shiftList0
            if shiftList and (shiftList0 is not shiftList):
                messager.warning('WARNING: Selected CCPN shift list does not match experiment %d (%s) shift list' % (experiment.serial, experiment.name))

        if peakList and shiftList:
            peakAndShiftLists.append([peakList, shiftList])

    if not peakAndShiftLists and peakShiftListKeys:
        messager.message('No selectable CCPN peak lists')
        return []

    return peakAndShiftLists

# BARDIAUX
def getCcpnModels(ccpnProject, model_names):
    """Descrn: Fetch CCPN models from a CCPN project using list
               of object keys
       Inputs: Implementation.Project, list of ccpn_id
       Output: List if ccp.api.molecule.MolStructure.Model  
    """

    
    ccpn_models = []
    for model_key in model_names:
        keys = getKeysFromString(model_key)
        model = getCcpnModel(ccpnProject, keys)
        if model:
            ccpn_models.append(model)

    return ccpn_models

# BARDIAUX
def getCcpnModel(ccpnProject, keys):
    """Descrn: Fetch a CCPN model from a CCPN project using list
               of object keys
       Inputs: Implementation.Project, key
       Output: List if ccp.api.molecule.MolStructure.Model  
    """
    model = None
    
    molSystemCode, structureEnsembleSerial, modelSerial = keys
    
    molSystem = ccpnProject.findFirstMolSystem(code=molSystemCode)
    if not molSystem:
        raise ValueError('No MolSystesm "%s" in CCPN project' % molSystem)
    else:
        sE = molSystem.findFirstStructureEnsemble(ensembleId=structureEnsembleSerial)

    if not sE:
        raise ValueError('No StructureEnsemble with serial "%d" in CCPN project' % structureEnsembleSerial)
    else:
        model = sE.findFirstModel(serial=modelSerial)

    if model is None:
        raise ValueError('No Model with serial "%d" in StructureEnsemble %d' % (modelSerial, structureEnsembleSerial))

    return model
#

def getCcpnConstraintLists(ccpnProject, restraintsNames):
    """Descrn: Fetch a CCPN constraint list from a CCPN project using list
               of restraint type and corresponding object keys
       Inputs: Implementation.Project, ARIA restraint names (what object is this?)
       Output: List if ccp.nmr.NmrConstraint.AbstractConstraintLists
       restraintsNames = dict with ARIA_DATA_types as keys, and list
       of ccpn_id as items.
    """

    constraintLists = {}
    for restraintType, constraintListKeys in restraintsNames.items():
        constraintLists[restraintType] = []

        for constraintListKey in constraintListKeys:
            keys = getKeysFromString(constraintListKey)
            constraintList = getCcpnConstraintList(ccpnProject, keys)

            if constraintList:
                constraintLists[restraintType].append(constraintList)


    return constraintLists


def getCcpnConstraintList(ccpnProject, keys):
    """Descrn: Fetch a CCPN constraint list from a CCPN project using list of object keys
       Inputs: Implementation.Project, List of keys (Words or Ints)
       Output: ccp.nmr.NmrConstraint.AbstractConstraintList
    """

    constraintList  = None
    #nmrProjectKey, storeSerial, serial = keys
    storeSerial, serial = keys

##   nmrProject = ccpnProject.findFirstNmrProject(name=nmrProjectKey)
##   if not nmrProject:
##     raise ValueError('No NMR project with name %s' % nmrProjectKey)

##  constraintStore = nmrProject.findFirstNmrConstraintStore(serial=storeSerial)
    constraintStore = ccpnProject.findFirstNmrConstraintStore(serial=storeSerial)
    if not constraintStore:
        raise ValueError('No NMR constraint store with serial "%d" in CCPN project' % storeSerial)
    else:
        constraintList = constraintStore.findFirstConstraintList(serial=serial)

    if constraintList is None:
        raise ValueError('No constraint list with serial "%d" in store "%d"' % (serial,storeSerial))

    return constraintList


def makeAriaMolecule(ccpMolSystem, ccpChain=None):
    """Descrn: Make an ARIA Molecule given a CCPN MolSystem object
       Inputs: ccp.molecule.MolSystem.MolSystem
       Output: ARIA Molecule
    """

    aria_molecule = Molecule(name=ccpMolSystem.code)

    ## BARDIAUX 2.2 extend
    ## sort chains according to residue.seqCode
    ccpChains = [[c,c.sortedResidues()[0].seqCode] for c in ccpChain]
    ccpChains.sort(lambda a, b: cmp(a[1], b[1]))

#  for chain in ccpMolSystem.chains:
    for chain, junk in ccpChains:
        if not ccpChain or chain in ccpChain:
            aria_chain = makeAriaChain(chain)
            aria_molecule.add_chain(aria_chain)

    return aria_molecule


def makeAriaChain(ccpChain):
    """Descrn: Make an ARIA Chain given a CCPN Chain object
       Inputs: ccp.molecule.MolSystem.Chain
       Output: ARIA Chain
    """


    # Does below work for DNA/RNA?

    aria_settings = ChainSettings()
    aria_settings['type'] = chainTypeMapping[ccpChain.molecule.molType]

    aria_chain = Chain(settings=aria_settings, segid=string_to_segid(ccpChain.code))

    prevRes = None
    ssDist = {}
    for residue in ccpChain.sortedResidues():

        aria_residue = makeAriaResidue(residue)
        aria_chain.addResidue(aria_residue)
        # BARDIAUX: add sec. struc
        prevRes, ssDist = addSecStrucToAriaRes(aria_residue, residue,  prevRes, ssDist)

    return aria_chain


def makeAriaResidue(ccpResidue):
    """Descrn: Make an ARIA Residue given a CCPN Residue object
       Inputs: ccp.molecule.MolSystem.Residue
       Output: ARIA Residue
    """


    ## TODO: have to define chemComp ARIA sysName for residues
    ## that match, throw error if not recognized (have to define
    ## own topology file for CNS for now)

    ## BARDIAUX : DIANA sysName for DNA 'cause 3Letter not ok
    if ccpResidue.molType == 'DNA':
        resCode = ccpResidue.molResidue.chemComp.\
                  stdChemComp.findFirstNamingSystem(name ='XPLOR').mainChemCompSysName.sysName
    else:
        resCode = ccpResidue.molResidue.chemComp.code3Letter

    try:
        aria_residue = Residue(number=ccpResidue.seqCode, residue_type=resCode)

    except:
        message = 'Residue %s not recognized in ARIA - define topology!'
        raise Exception(message %  ccpResidue.ccpCode)

    ariaAtomDict = {}

    idNum = 0
    for atom in ccpResidue.sortedAtoms():
        aria_atom = makeAriaAtom(atom, idNum)
        aria_residue.addAtom(aria_atom)
        ariaAtomDict[atom.name] = aria_atom
        idNum += 1

    setupEquivalentGroup(ccpResidue.chemCompVar, aria_residue, ariaAtomDict)

    return aria_residue


def makeAriaAtom(ccpAtom, idNum, heteroElements=('N', 'C')):
    """Descrn: Make an ARIA Atom given a CCPN Atom object
               Also sets up the hetero atom name.
       Inputs: ccp.molecule.MolSystem.Atom
       Output: ARIA Atom
    """

    chemAtom = ccpAtom.chemAtom
    elementSymbol = chemAtom.elementSymbol

    heteroName = None
    if elementSymbol == 'H':
        for bound in chemAtom.findFirstChemBond().chemAtoms:
            if bound.elementSymbol in heteroElements:
                heteroName = bound.name
                break

    aria_settings = AtomSettings()
    aria_settings['type'] = elementSymbol
    aria_settings['hetero_atom_name' ] = heteroName

    aria_atom = Atom(settings=aria_settings, name=ccpAtom.name, id=idNum)
    aria_atom._setSegid(string_to_segid(ccpAtom.residue.chain.code))

    return aria_atom


def makeAriaChemicalShift(ccpShift):
    """Descrn: Make an ARIA ChemicalShift given a CCPN Shift object
       Inputs: ccp.nmr.Nmr.Shift
       Output: ARIA ChemicalShift
    """

    if ccpShift is None:
        val = None
        err = None
    else:
        val = ccpShift.value
        err = ccpShift.error

    return ChemicalShift(val, err)


def makeAriaShiftList(ccpShiftList, ccpChains, ariaMolecule):
    """Descrn: Make a populated ARIA ChemicalShiftList for a given molecule given a CCPN ShiftList object
       Inputs: ccp.nmr.Nmr.Shift, ccp.molecule.MolSystem.Chains, ARIA Molecule
       Output: ARIA ChemicalShiftList
    """

    aria_shiftList = ChemicalShiftList()

    ccpMolSystem = ccpChains[0].molSystem
    
    ss = []

    from aria.OrderedDict import OrderedDict

    shiftLookup = {}
    for shift in ccpShiftList.sortedMeasurements():
        shiftLookup[shift.resonance] = shift

    #for chain in ccpMolSystem.sortedChains():
    for chain in ccpChains:
        for residue in chain.sortedResidues():
            atomSetDict = OrderedDict()

            for atom in residue.sortedAtoms():
                atomSet = atom.atomSet
                if atomSet and atomSet.resonanceSets:
                    atomSetDict[atomSet] = True

            for atomSet in atomSetDict.keys():
                if atomSetDict[atomSet] is False:
                    # TJS: This is very unlikely to happen, but it is not impossible
                    # in the CCPN model to have atomSets defined twice for same atoms
                    continue


                resonanceSets = atomSet.sortedResonanceSets()

                if len(resonanceSets) > 1:
                    messager.message('CCPN atom set %d%s %s has multiple resonance assignments' % (residue.seqCode, residue.ccpCode, atomSet.name))
                    nonstereo = None
                    stereo    = None
                    for resonanceSet0 in resonanceSets:
                        if len(resonanceSet0.atomSets) > 1:
                            nonstereo = resonanceSet0
                        else:
                            stereo = resonanceSet0

                    if stereo:
                        if nonstereo:
                            atomSets2 = nonstereo.sortedAtomSets()
                            atomSets2.remove(atomSet)
                            resonanceSet = nonstereo
                            for resonanceSet0 in atomSets2[0].sortedResonanceSets(): # resonance set of other atomSet in pair
                                if len(resonanceSet0.atomSets) == 1: # if other atom set stereospecifically assigned
                                    resonanceSet = stereo
                        else:
                            resonanceSet = stereo
                    else:
                        resonanceSet = nonstereo
                else:
                    resonanceSet = resonanceSets[0]

                nA = len(resonanceSet.atomSets)
                nR = len(resonanceSet.resonances)

                if (nA==1) and (nR==1):
                    # TJS mod
                    resonance = resonanceSet.findFirstResonance()
                    ariaAtoms = getAriaAtomsFromResonance(resonance, ariaMolecule)
                    nAtoms =  len(ariaAtoms)

                    shift = shiftLookup.get(resonance)

                    if shift:
                        aria_shift = makeAriaChemicalShift(shift)
                        atomSetDict[atomSet] = False

                        if (nAtoms == 2) and (resonance.isotopeCode == '13C'):
                            # TJS horrid fix for carbons that are equivalent in CCPN
                            # but not in ARIA e.g. CD1/CD2 of PHE
                            # Makes two ARIA spin systems rather than one equivalent
                            # shift is repeated

                            # First carbon
                            aria_spinsystem = SpinSystem(AVERAGING_METHOD_NONE)
                            aria_spinsystem.setAtoms((ariaAtoms[0],))
                            aria_spinsystem.setChemicalShifts((aria_shift,))

                            aria_shiftassignment = ShiftAssignment(ASSIGNMENT_METHOD_STEREO_SPECIFIC)
                            aria_shiftassignment.setSpinSystems((aria_spinsystem,))
                            aria_shiftList.addShiftAssignment(aria_shiftassignment)

                            # Second carbon
                            aria_spinsystem = SpinSystem(AVERAGING_METHOD_NONE)
                            aria_spinsystem.setAtoms((ariaAtoms[1],))
                            aria_spinsystem.setChemicalShifts((aria_shift,))

                            aria_shiftassignment = ShiftAssignment(ASSIGNMENT_METHOD_STEREO_SPECIFIC)
                            aria_shiftassignment.setSpinSystems((aria_spinsystem,))
                            aria_shiftList.addShiftAssignment(aria_shiftassignment)

                        else:
                            if nAtoms == 1:
                                method = AVERAGING_METHOD_NONE
                                assignment_type = ASSIGNMENT_METHOD_STEREO_SPECIFIC

                            else:
                                method = AVERAGING_METHOD_FAST
                                assignment_type = ASSIGNMENT_METHOD_EQUIVALENT

                            aria_spinsystem = SpinSystem(method)
                            aria_spinsystem.setAtoms(tuple(ariaAtoms))
                            aria_spinsystem.setChemicalShifts((aria_shift,))

                            if (aria_spinsystem,) in ss:
                                continue

                            ss.append((aria_spinsystem,))

                            aria_shiftassignment = ShiftAssignment(assignment_type)
                            aria_shiftassignment.setSpinSystems((aria_spinsystem,))
                            aria_shiftList.addShiftAssignment(aria_shiftassignment)


                elif (nA==2) and (nR<3):

                    ariaShifts = []
                    for resonance in resonanceSet.sortedResonances():
                        shift = shiftLookup.get(resonance)
                        if shift:
                            ariaShifts.append(makeAriaChemicalShift(shift))

                    while len(ariaShifts) < nA:
                        ariaShifts.append(ChemicalShift(None))

                    aria_spinsystems = []

                    # TJS mod
                    for resonance in resonanceSet.sortedResonances():
                        ariaAtoms = getAriaAtomsFromResonance(resonance, ariaMolecule)

                        if len(ariaAtoms) == 1:
                            method = AVERAGING_METHOD_NONE
                        else:
                            method = AVERAGING_METHOD_FAST

                        aria_spinsystem = SpinSystem(method)
                        aria_spinsystem.setAtoms(tuple(ariaAtoms))
                        aria_spinsystem.setChemicalShifts(tuple(ariaShifts))
                        aria_spinsystems.append(aria_spinsystem)

                        atomSetDict[atomSet] = False


                    if aria_spinsystems in ss:
                        # TJS: This ensures prochirals will not be put in twice
                        # Both are put in on first run-through - resonance
                        # link to both atom groups
                        continue

                    ss.append(aria_spinsystems)

                    aria_shiftassignment = ShiftAssignment(ASSIGNMENT_METHOD_FLOATING)
                    aria_shiftassignment.setSpinSystems(tuple(aria_spinsystems))
                    aria_shiftList.addShiftAssignment(aria_shiftassignment)


    return aria_shiftList

def getAcqRefExpDimRef(refExperiment):
    """ RefExpDimRef that corresponds to acquisition dimension
    """

    # get acquisition measurement
    expGraph = refExperiment.nmrExpPrototype.findFirstExpGraph()

    # even if there are several the
    # acquisition dimension should be common.
    expSteps = [(es.stepNumber, es) for es in expGraph.expSteps]
    expSteps.sort()

    if refExperiment.isReversed:
        acqMeasurement = expSteps[0][1].expMeasurement
    else:
        acqMeasurement = expSteps[-1][1].expMeasurement

    # get RefExpDimRef that fits measurement
    acqRefs = []
    for refExpDim in refExperiment.sortedRefExpDims():
        for refExpDimRef in refExpDim.sortedRefExpDimRefs():
            if refExpDimRef.expMeasurement is acqMeasurement:
                acqRefs.append(refExpDimRef)

    if len(acqRefs) == 1:
        return acqRefs[0]

    else:
        msg = "%s has no unambiguous RefExpDimRef for acqMeasurement (%s)"
        messager.warning(msg % (refExperiment, acqMeasurement))


def checkExpTransfers(experiment):
    """Descrn: Check whether a CCPN Experiment's transfers are setup.
               Links any missing refReperiment information.
       Inputs: ccp.nmr.Nmr.Experiment
       Output: None
    """

    refExperiment = experiment.refExperiment
    if not refExperiment:
        return

    refExpDims = refExperiment.sortedRefExpDims()
    if not refExpDims:
        # Something is wrong with the reference data
        return

    # TJS: best not overwrite; this is exp probably curated
    # Guessing dim mappings is probably more error prone...
    if experiment.expTransfers:
        return

    # Clear everything out, just in case it is wrong
    for expDim in experiment.expDims:
        if expDim.refExpDim:
            for expDimRef in expDim.expDimRefs:
                if expDimRef.refExpDimRef:
                    expDimRef.setRefExpDimRef(None)
            expDim.setRefExpDim(None)

    acqRefExpDim = getAcqRefExpDimRef(refExperiment).refExpDim
    acqExpDim = experiment.findFirstExpDim(isAcquisition=True)

    for expDim in experiment.sortedExpDims():
        expData = []
        
        for expDimRef in expDim.expDimRefs:
            isotopes = frozenset(expDimRef.isotopeCodes)
            
            if isotopes:
                mType = expDimRef.measurementType.lower()
                expData.append((mType, isotopes))
        
        if not expData:
            continue

        for refExpDim in refExpDims:
            refData = [] 
            
            for refExpDimRef in refExpDim.refExpDimRefs:
                expMeasurement = refExpDimRef.expMeasurement
                isotopes = frozenset([x.isotopeCode for x in expMeasurement.atomSites])
                mType = expMeasurement.measurementType.lower()
                refData.append((mType, isotopes))
       
            if expData == refData:
                expDim.setRefExpDim(refExpDim)
                refExpDims.remove(refExpDim)
                
                if refExpDim is acqRefExpDim:
                    if not acqExpDim:
                        expDim.isAcquisition = True
                        acqExpDim = expDim
            
                break

    for expDim in experiment.sortedExpDims():
        if not expDim.expDimRefs:
            continue

        if not expDim.refExpDim:
            expDim.setRefExpDim(refExpDims.pop(0))

        # set reference data comparison list
        refExpDimRefs = list(expDim.refExpDim.refExpDimRefs)
        refData = []
        for refExpDimRef in refExpDimRefs:
            expMeasurement = refExpDimRef.expMeasurement
            atomSites = expMeasurement.atomSites
            refData.append((frozenset(x.isotopeCode for x in atomSites),
                           expMeasurement.measurementType.lower(),
                           frozenset(x.name for x in atomSites),
                           refExpDimRef))

        # set experiment data comparison list
        inData = []
        for expDimRef in expDim.expDimRefs:
            inData.append((frozenset(expDimRef.isotopeCodes),
                          expDimRef.measurementType.lower(),
                          frozenset(((expDimRef.displayName or expDimRef.name),)),
                          expDimRef))

        # match expDimRef to refExpDimRef. comparing isotopeCodes,
        # if equal measurementTypes, if equal name/displayname
        for end in (-1,-2,-3):
            for ii in range(len(inData)-1, -1, -1):
                for jj in range(len(refData)-1, -1, -1):
                    if inData[ii][:end] == refData[jj][:end]:
                        inData[ii][-1].setRefExpDimRef(refData[jj][-1])
                        del inData[ii]
                        del refData[jj]
                        break

    for expTransfer in experiment.expTransfers:
        expTransfer.delete()
        
    transferDict = {}
    for expDim in experiment.expDims:
        for expDimRef in expDim.expDimRefs:
            if not expDimRef.refExpDimRef:
                continue

            measurement = expDimRef.refExpDimRef.expMeasurement
            if measurement.measurementType in ('Shift','shift'):
                atomSite = measurement.findFirstAtomSite()

                for expTransfer in atomSite.expTransfers:
                    if transferDict.get(expTransfer) is None:
                        transferDict[expTransfer] = []
                    transferDict[expTransfer].append(expDimRef)

    newTransfer = experiment.newExpTransfer

    for refTransfer in transferDict.keys():
        expDimRefs = frozenset(transferDict[refTransfer])

        if len(expDimRefs) == 2:
            transferType = refTransfer.transferType
            newTransfer(transferType=transferType, expDimRefs=expDimRefs)


# TJS: Addef option to filter out peaks rejected on CCPN

def makeAriaSpectrum(peakList, ariaMolecule, filterRejected=True):
    """Descrn: Make an ARIA NOESY Spectrum (and assign it)
               given a CCPN peak list and ARIA Molecule
       Inputs: ccp.nmr.Nmr.PeakList, ARIA Molecule
       Output: ARIA NOESYSpectrum
    """

    ariaDimNames = ('Proton1','Proton2','Hetero1','Hetero2')

    spectrum   = peakList.dataSource
    experiment = spectrum.experiment
    shiftList  = experiment.shiftList or experiment.nmrProject.findFirstMeasurementList(className='ShiftList')

    # TJS: Added checks for experiment not fully linked
    # to its refExperiment
    checkExpTransfers(experiment)

    transfer = experiment.findFirstExpTransfer(transferType='NOESY')

    if not transfer:
        raise Exception('%s is not a NOESY.' % experiment.name)

    light_dim1 = None
    light_dim2 = None
    heavy_dim1 = None
    heavy_dim2 = None

    expDimRefDict = {}

    for expDimRef in transfer.sortedExpDimRefs():
        if expDimRef.isotopeCodes != ('1H',):
            raise Exception('Not an H-H NOESY')

        onebondTransfer = expDimRef.findFirstExpTransfer(transferType='onebond')

        if onebondTransfer:
            expDimRefs = list(onebondTransfer.expDimRefs)
            expDimRefs.remove(expDimRef)
            expDimRefX = expDimRefs[0]

            if light_dim1 is None:
                heavy_dim1 = findCcpnDataDim(spectrum, expDimRefX.expDim)
                light_dim1 = findCcpnDataDim(spectrum, expDimRef.expDim)

            else:
                heavy_dim2 = findCcpnDataDim(spectrum, expDimRefX.expDim)
                light_dim2 = findCcpnDataDim(spectrum, expDimRef.expDim)

        else:

            if light_dim1 is None:
                heavy_dim1 = None
                light_dim1 = findCcpnDataDim(spectrum, expDimRef.expDim)

            else:
                heavy_dim2 = None
                light_dim2 = findCcpnDataDim(spectrum, expDimRef.expDim)


    cross_peaks = []
    for peak in peakList.peaks:
        # TJS added to remove rejected peaks
        if filterRejected and peak.figOfMerit == 0.0:
            continue

        volume = peak.findFirstPeakIntensity(intensityType='volume')

        if volume:
            value = volume.value
            err   = volume.error
        else:
            value = None
            err   = None

        volume = Datum(value, err)

        height = peak.findFirstPeakIntensity(intensityType='height')
        if height:
            value = height.value
            err   = height.error
        else:
            value = None
            err   = None

        intensity = Datum(value, err)

        cross_peak = CrossPeak(number=peak.serial, volume=volume, intensity=intensity)

        peakDims = []
        for dim in (light_dim1, light_dim2, heavy_dim1 ,heavy_dim2):
            peakDims.append(peak.findFirstPeakDim(dim=dim))


        for i in range(len(ariaDimNames)):
            peakDim = peakDims[i]
            if peakDim:
                setShiftFunc =  getattr(cross_peak, 'set%sChemicalShift' % ariaDimNames[i])
                setShiftFunc(ChemicalShift(peakDim.value, peakDim.valueError))

        assignments = []
        for peakDim in peakDims:
            assignments0 = []

            if peakDim:
                for contrib in peakDim.peakDimContribs:
                    resonance = contrib.resonance
                    # TJS mod
                    ariaAtoms = getAriaAtomsFromResonance(resonance, ariaMolecule)
                    # BB when no ariaAtoms
                    if not ariaAtoms:
                        continue
                    assi = Assignment(tuple(ariaAtoms), assignment_type=ASSIGNMENT_TYPE_MANUAL)
                    if assi not in assignments0:
                        assignments0.append(assi)

            assignments.append(assignments0)

        for i in range(len(ariaDimNames)):
            addAssignFunc =  getattr(cross_peak, 'add%sAssignment' % ariaDimNames[i])
            for assignment in assignments[i]:
                addAssignFunc(assignment)


        cross_peaks.append(cross_peak)

    spectrum = NOESYSpectrum(name=experiment.name, noes=cross_peaks)

    return spectrum

def getAriaAtomsFromResonance(resonance, ariaMolecule, cache={}):
    """Descrn: Get the corresponding Aria2 atom in a molecule
               given the input CCPN Resonance object
       Inputs: Nmr.Resonance, Aria2 Molecule object
       Output: List of List of Aria2 Atom objects
    """

    from aria.tools import string_to_segid

    ariaAtoms = cache.get(resonance)
    if ariaAtoms:
        return ariaAtoms
    else:
        ariaAtoms = []

    resonanceSet = resonance.resonanceSet
    if not resonanceSet:
        messager.warning('Attempt to get atoms from unassigned CCPN resonance (%d)' % resonance.serial)
        return []

    atomSets = resonanceSet.sortedAtomSets()

    residue = atomSets[0].findFirstAtom().residue

    ariaChain = ariaMolecule.getChain(string_to_segid(residue.chain.code)) # Func warns if failure

    ariaResidue = ariaChain.residues.get(residue.seqCode)
    if not ariaResidue:
        messager.warning('Could not find ARIA Residue for CCPN residue %d%s' % (residue.seqCode,residue.ccpCode))
        return []


    #for atomSet in atomSets:
    # - Removed this loop and added below.
    # - Loop was giving two ARIA atoms for prochirals: now using only one
    # because the restraint re-export was giving two resonances where
    # it couldn't tell which was which: both linked to the same two atoms

    # TJS Added
    index   = resonanceSet.sortedResonances().index(resonance)
    atomSet = atomSets[min(index,len(atomSets)-1)]

    # TJS modify to return just a list of atoms, rather than a list of list
    ariaAtoms = []
    for atom in atomSet.atoms:
        ariaAtom = ariaResidue.atoms.get(atom.name)
        if not ariaAtom:
            messager.warning('Could not find ARIA Atom for CCPN atom %d%s %s' % (residue.seqCode,residue.ccpCode, atom.name))
            return []

        ariaAtoms.append(ariaAtom)

    cache[resonance] = ariaAtoms

    return ariaAtoms


def findCcpnDataDim(spectrum, expDim):
    """Descrn: Get the data dimension number that corresponds
               to a given experimental dimension of a spectrum
       Inputs: ccp.nmr.Nmr.DataSource, ccp.nmr.Nmr.ExpDim
       Output: Int
    """

    dataDim = None
    if expDim:
        dataDim = spectrum.findFirstDataDim(expDim=expDim)

    if dataDim:
        return dataDim.dim

# BARDIAUX
def dumpModel(model, molSystem, data_dir):
    
    from ccpnmr.format.converters.CnsFormat import CnsFormat
    import aria.DataContainer as DC
    import os

    ccpProject = model.root

    cns_format = CnsFormat(ccpProject)

    cns_format.molSystem = molSystem

    fileName = getObjectKeyString(model, delimiter='_').replace(" ", "")
    fileName = "%s.pdb" % fileName

    dst = os.path.join(data_dir, fileName)

    kw = {'minimalPrompts'     : 1,
          'forceNamingSystemName'  : 'XPLOR',
          'forceExportChainId' : " ",
          'verbose' : False}

    cns_format.writeCoordinates(dst, structures=[model],**kw)

    return dst
    
# BARDIAUX
# TJS: Added molSystem
def dumpRestraintsList(constraintList, molSystem, data_dir, cType):
    """
    Dump CCPN RestraintList to CNS Format TBL file.
    """

    import aria.DataContainer as DC
    import os
    from ccpnmr.format.converters.CnsFormat import CnsFormat

    ccpProject = constraintList.root#nmrConstraintHead.project

    cns_format = CnsFormat(ccpProject)
    # TJS: Also included molSystem
    cns_format.molSystem = molSystem
    
    f = {DC.DATA_UNAMBIGUOUS : cns_format.writeDistanceConstraints,
         DC.DATA_AMBIGUOUS   : cns_format.writeDistanceConstraints,
         DC.DATA_HBONDS      : cns_format.writeHBondConstraints,
         DC.DATA_DIHEDRALS   : cns_format.writeDihedralConstraints,
         DC.DATA_RDCS        : cns_format.writeRdcConstraints,
         DC.DATA_KARPLUS     : cns_format.writeJCouplingConstraints,
         DC.DATA_SSBONDS     : cns_format.writeDistanceConstraints}


    kw = {'compressResonances' : 0,
          'minimalPrompts'     : 1,
          'forceNamingSystemName'  : 'XPLOR',
          'resetMapping' : 1,
          'verbose' : False}

    fileName = getObjectKeyString(constraintList, delimiter='_').replace(" ", "")
    fileName = "%s.tbl" % fileName

    dst = os.path.join(data_dir, fileName)

    # assumes that Resonnances are already linked
    # overwrite existing file without warning

    # attribute 'usePeakInfo' has default True now, but not all ccpn projects has
    # this attribute and it's not necessary to ARIA 2 any longer
    if cType in (DC.DATA_AMBIGUOUS, DC.DATA_UNAMBIGUOUS): #by AWSS
        kw['usePeakInfo'] = False  #by AWSS

    f[cType](dst, constraintList=constraintList,**kw)

    return dst

# BARDIAUX
def getAriaDistanceRestraintsList(constraint_list, constraint_type, aria_mol):
    """
    Convert CCPN DistanceContraintList
    to list of ARIA DistanceRestraint
    """

    from aria.ShiftAssignment import SpinSystem, AVERAGING_METHOD_NONE
    import aria.Contribution as C
    from aria.Singleton import SpinPairFactory
    from aria.AriaPeak import DistanceRestraint

    # BARDIAUX 2.2 NEW
    from aria.CrossPeak import CrossPeak
    from aria.NOESYSpectrum import ConstraintList
    from numpy import power
    from aria.DataContainer import PeakData#, LowerBoundCorrection, UpperBoundCorrection
    from aria.Datum import Datum

    ccpn_project = constraint_list.root


    SpinPair = SpinPairFactory()

    restraints = []
    peaks = []

    #source = getObjectKeyString(constraint_list)
    source = "%s %s" % (constraint_list.name, getObjectKeyString(constraint_list))

    msg = "%s ignored."

    for distConstr in constraint_list.constraints:

        spin_systems_pairs = []

        # BARDIAUX 2.2 NEW
        if distConstr.origData:
            vol = distConstr.origData
        else:
            vol = power(distConstr.targetValue, -6)

        vol = Datum(vol, 0)
        xpk = CrossPeak(distConstr.serial, vol, vol)

        contributions = []

        restraint = DistanceRestraint()

        target, upper, lower = distConstr.targetValue, distConstr.upperLimit, distConstr.lowerLimit
        weight = distConstr.weight

        restraint.setDistance(target)
        restraint.setLowerBound(lower)
        restraint.setUpperBound(upper)
        restraint.setWeight(weight)


        for constrItem in distConstr.items:

            reso1, reso2 = constrItem.resonances

            # TJS fix for mapping prochirals
            # always use real resonnances where possible
            if reso1.resonance:
              reso1 = reso1.resonance
            if reso2.resonance:
              reso2 = reso2.resonance

            ss1 = SpinSystem(AVERAGING_METHOD_NONE)
            atoms = getAriaAtomsFromResonance(reso1, aria_mol)

            if not atoms:
                messager.warning(msg % (constrItem))
                continue
            # TJS remove; now only get a simple list of ARIA atoms
            #elif len(atoms)> 1:
            #  atoms = map(lambda x: x[0], atoms)
            #else:
            #  atoms = atoms[0]

            ss1.setAtoms(tuple(atoms))

            ss2 = SpinSystem(AVERAGING_METHOD_NONE)
            atoms = getAriaAtomsFromResonance(reso2, aria_mol)

            if not atoms:
                messager.warning(msg % (constrItem))
                continue
            # TJS remove; now only get a simple list of ARIA atoms
            #elif len(atoms)> 1:
            #  atoms = map(lambda x: x[0], atoms)
            #else:
            #  atoms = atoms[0]

            ss2.setAtoms(tuple(atoms))

            spin_pairs = []
            for a1 in ss1.getAtoms():
                for a2 in ss2.getAtoms():

                    if a1 == a2:
                        continue

                    sp = SpinPair(a1, a2)
                    spin_pairs.append(sp)

            if not len(spin_pairs):
                continue

            if (ss1,ss2) in spin_systems_pairs or \
                   (ss2,ss1) in spin_systems_pairs :
                continue

            spin_systems_pairs.append((ss1,ss2))

            c = C.Contribution(0,
                               C.CONTRIBUTION_TYPE_FAST_EXCHANGE,
                               spin_pairs, spin_systems = (ss1, ss2))



            ## TODO: hack
            c.setSpinSystems((ss1, ss2))

            contributions.append(c)

            ## Need a weight if no further Contribution evaluator
            c.setWeight(1.)

        if len(contributions) == 0:
            messager.warning(msg % (distConstr))
            continue

        restraint.setContributions(contributions)

        peaks.append(xpk)
        restraint.setReferencePeak(xpk)

        restraints.append(restraint)

    # assign a ARIA ConstraintList as Spectrum
    s = ConstraintList(source, peaks)
    pd = PeakData()
    pd.reset()
    pd['ccpn_id'] =  getObjectKeyString(constraint_list)
    s.setDataSource(pd)

    return restraints, s
    # return restraints

# BARDIAUX 2.2
def getCcpnExperimentData(peakList):

    experiment = peakList.dataSource.experiment
    frequency = mixing = None

    # Spec freq
    spec = experiment.spectrometer
    if spec:
        frequency = spec.protonFreq

    # Mixing time
    # CCPN return seconds, ARIA needs ms
    transfer = experiment.findFirstExpTransfer(transferType='NOESY')
    if transfer:
        mixing = transfer.mixingTime
        if mixing:
            mixing = mixing *1e3

    # Correlation Time ?
    return frequency, mixing


# BARDIAUX : add sec. struc. info
def addSecStrucToAriaRes(ariaRes, ccpnRes, prevRes, ssDict):

    ss = [rg.secStrucCode for rg in ccpnRes.resonanceGroups if rg.secStrucCode is not None]
    if not ss:
        return prevRes, ssDict

    ss = ss[0]
    if prevRes is not None:
        old_ss, old_nb = prevRes
        if old_ss == ss:
            k = prevRes
        else:
            old_nb = ssDict.get(ss)
            if old_nb:
                k = (ss, old_nb+1)
                ssDict[ss] +=1
            else:
                k = (ss, 1)
                ssDict[ss] = 1
    else:
        k = (ss, 1)
        ssDict[ss] = 1

    prevRes = k

    val = '%s%d' % k
    ariaRes.structure = val
    return prevRes, ssDict

