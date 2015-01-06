from ccpnmr.analysis.core.ConstraintBasic import makeNmrConstraintStore, makeStructureGeneration, getFixedResonance
from ccpnmr.analysis.core.ExperimentBasic import getOnebondDataDims
from ccpnmr.analysis.core.AssignmentBasic import assignAtomsToRes, assignResToDim
from ccpnmr.analysis.core.PeakBasic       import pickPeak, setManualPeakIntensity
from ccpnmr.analysis.core.MoleculeBasic   import DEFAULT_ISOTOPES
from ccpnmr.analysis.core.StructureBasic  import getBestChemComp

from ccp.util.Molecule import addMolResidues
from aria.importFromCcpn import getCcpnPeakList, getKeysFromString

from aria.Chain import TYPE_PROTEIN

## BARDIAUX BB log
from aria.ariabase import AriaBaseClass
messager = AriaBaseClass()
messager._name = 'CCPN export'
##

def getAria2ObjectsFromPickle(fileName):
  """Descrn: Given an ARIA2 pickle file (and an istallation) get the ARIA2
             data model objects.
     Inputs: Word (pickle file name)

     Output: List of Aria2 Restraint objects
  """
  import cPickle, os

  print "Reading Pickle file... "

  fileName = os.path.expanduser(fileName)
  if fileName[-3:].lower() == '.gz':
    import gzip
    file = gzip.open(fileName)
  else:
    file = open(fileName)

  objects = None

  while 1:
    try:
      objects= cPickle.load(file)
    except EOFError:
      break

  file.close()

  print " ...done."

  return objects


def importAria2RunData(dirName, constraintSet=None, project=None, loadStructures=False):
  """Descrn: Import ARIA2 structures plus restraints and violations from a Pickle file
             given an iteration director of an ARIA2 run.
     Inputs: Word (run directory name), Nmr.NmrConstraintStore, Implementation.MemopsRoot
     Output: NmrConstraint.ConstraintList, NmrConstraint.ViolationList

  """

  assert constraintSet or project
  import aria.ariabase as ariabase
  import os, sys

  if constraintSet:
    nmrProject = constraintSet.nmrConstraintStore.nmrProject
  else:
    nmrProject = project.currentNmrProject

  ariaPath = os.environ.get(ariabase.ARIA_ENV)
  if not ariaPath:
    print 'Warning','Could not find ARIA2 environmant variable to locate the ARIA installation'
    return

  ariaPath = os.path.join(ariaPath, 'src/py')
  sys.path.insert(0, ariaPath)

  print "Getting constraint set"
  constraintSet = makeNmrConstraintStore(nmrProject) # This function is present locally

  pickleFile  = None
  fileNames   = os.listdir(dirName)
  for fileName in fileNames:
    if fileName[-7:] == '.pickle':
      pickleFile = fileName
      break
    elif fileName[-10:] == '.pickle.gz':
      pickleFile = fileName
      break

  if not pickleFile:
    print 'Failure','Found no .pickle file in directory %s' % dirName
    return

  print "Loading PICKLE file %s" % pickleFile
  pickleFile = os.path.join(dirName,pickleFile)
  objects = getAria2ObjectsFromPickle(pickleFile)
  if not objects:
    print 'Failure','Unable to find any Aria2 restraint objects'
    sys.path.remove(ariaPath)
    return

  print "Getting chain"
  # BARDIAUX : Update return list of chains
  chains = getChainsFromAria2(objects, project)
  if not chains:
    print 'Failure','No molecular system chain found for ARIA2 data'
    return

  structures = None
  if loadStructures:
    print "loading structures"
    # BARDIAUX : give chains list
    structures = getStructuresFromAria2Dir(dirName, chains)
    if not structures:
      print 'Failure','Unable to load any structures from Aria2 iteration directory'
      return

    structureGeneration = makeStructureGeneration(structures, constraintSet)
    structureGeneration.name = 'ARIA2 import'

  sys.path.remove(ariaPath)

  print "Interpreting ARIA2 restraints"
  # BARDIAUX : Update  needs chainS
  return getConstraintsFromAria2(objects, chains, constraintSet, structures=structures)


def getChainFromAria2(objects, project, aria_chain=None):
  """Descrn: Get matching (or new) Ccp.MolSystem Chain from Aria2 restraint objects
     Inputs: List of Aria2 restreaint objects, Implementation.MemopsRoot
     Output: MolSystem.chain
  """

  from ccpnmr.analysis.core.MoleculeBasic   import findMatchingChain

  try:
    from memops.gui.DataEntry import askString
  except:
    from memops.universal.DataEntry import askString

  # TBD: deal with multiple chains
  # BARDIAUX: Done. See getChainsFromAria2

  if aria_chain is None:
    ariaChain = objects[0].getReferencePeak().getProton1Assignments()[0].getAtoms()[0].residue.chain
  else:
    ariaChain = aria_chain

  sequence = [r.type for r in ariaChain.getResidues()]


  chains = []

  #by AWSS begin
  # issue: ccpn project with several molecular systems, user want use 3rd molSystem.
  # molSystem should be the one user defined in the aria project.
  # My solution:
  molSystemName = aria_chain.getMolecule().getName()
  molSystem = project.findFirstMolSystem(code=molSystemName)
  #chain0 = findMatchingChain(molSystem, sequence, doWarning=True)[0] #by AWSS
  chain0 = molSystem.findFirstChain(code = ariaChain.getSegid().strip()) # BB only works with if initial CCPN project
  if chain0:
    chains.append(chain0)

  ## it was:
  #for molSystem in project.molSystems:
  #  chain0 = findMatchingChain(molSystem, sequence, doWarning=True)[0] #by AWSS
  #  if chain0:
  #    chains.append(chain0)
  # by AWSS end

  chain = None
  if len(chains) == 1:
    chain = chains[0]


  elif len(chains) > 1:
    texts = ['%s:%s' % (ch.molSystem.code, ch.code) for ch in chains]
    while chain is None:
      response = askString('Enter Mol System and chain codes (Available: %s)' % ' '.join(texts),texts[0])
      if not response:
        continue

      codes = response.split(':')
      if len(codes) != 2:
        continue

      molSystem = project.findFirstMolSystem(code=codes[0])
      if not molSystem:
        continue

      chain = molSystem.findFirstChain(code=codes[1])


  # TBD: consider making molecule, chain and molSystem

  return chain

# BARDIAUX
def getChainsFromAria2(objects, project, aria_chain=None):
  """Descrn: Get matching (or new) Ccp.MolSystem Chains from Aria2 restraint objects
     Inputs: List of Aria2 restreaint objects, Implementation.MemopsRoot
     Output: List of MolSystem.chain
  """

  from ccpnmr.analysis.core.MoleculeBasic import findMatchingChains

  if aria_chain is not None:
    ariaChain = aria_chain
  else:
    ariaChain = objects[0].getReferencePeak().getProton1Assignments()[0].getAtoms()[0].residue.chain


  chains = []
  molSystemName = aria_chain.getMolecule().getName()
  molSystem = project.findFirstMolSystem(code=molSystemName)

  # TJS adjust for when not using CCPN to start with
  if molSystem:
    for aria_chain in ariaChain.getMolecule().get_chains():
      segid = aria_chain.getSegid().strip()
      if segid == '':
        segid = ' '
      chain0 = molSystem.findFirstChain(code = segid) # BB only works with if initial CCPN project
      if chain0:
        chains.append(chain0)

 
  if not chains:
    # TJS: Have to do sequence alignment
    for aria_chain in ariaChain.getMolecule().get_chains():
       sequence = getCcpnSeqfromAriaChain(aria_chain)

       for molSystem in project.molSystems:
         for molChain in molSystem.chains:
           mChains, mappings = findMatchingChains(molSystem, sequence,
                                                  excludeChains=chains, doWarning=True)
           if mChains:
             chains.append(mChains[0])

    if not chains:
      # TJS: Nothing matches - make new chains
      
      i = 1
      code = 'MS%d' % i
      while project.findFirstMolSystem(code=code):
        i += 1
        code = 'MS%d' % i
      
      molSystem = project.newMolSystem(code=code)
      
      for aria_chain in ariaChain.getMolecule().get_chains():
        sequence = getCcpnSeqfromAriaChain(aria_chain)
        
        i = 1
        name = aria_chain.getMolecule().getName()
        while project.findFirstMolecule(name=name):
          name = '%s_%d' % (aria_chain.getMolecule().getName(),i)
          i += 1
        
        molecule = project.newMolecule(name=name)
        molType = aria_chain.getType()
        
        if molType == TYPE_PROTEIN: # ARIA uses uppercase
          molType = 'protein'
          
        addMolResidues(molecule, molType, sequence)
        chainCode = aria_chain.getSegid().strip() or ' '
        chain = molSystem.newChain(code=chainCode, molecule=molecule)
        chains.append(chain)
        
  return chains

def getCcpnSeqfromAriaChain(ariaChain, project):
  """Descrn: Get "ccpCodes" for CCPN residues/chemComps based
             upon three-letter code sequence from an aria chain
     Inputs: araia.Chain, Implementation.MemopsRoot
     Output: List of Words (ChemComp.CcpCodes)
  """

  sequence = []
  
  for residue in ariaChain.getResidues():
    atomNames = [a.name for a in residue.getAtoms()]
    chemComp = getBestChemComp(project, residue.type, atomNames)
    sequence.append(chemComp.ccpCode)
  
  return sequence

# BARDIAUX 2.3 : change reading of ARIA PDB from Analysis to FormatConverter
# some issues with getStructureFromFile from Analysis:
# wrong HB2/3 and O atoms mapping, molType other connect to ebi.ac.uk..
# no ANI or TIP3 filetring allowed

def getStructuresFromAria2Dir_Analysis(dirName, chains):
  """Descrn: Load structures of the input chains, from the PDB style files in an Aria2 iteration directory
     Inputs: Word - Aria2 run directory path. List of MolSystem.Chains
     Output: MolStrcuture.StructureEnsemble
  """

  from ccpnmr.analysis.core.StructureBasic  import getStructureFromFile, makeEnsemble
  import os

  structures = []
  fileNames  = os.listdir(dirName)
  for fileName in fileNames:
    if fileName[-4:] == '.pdb':
      fileName = os.path.join(dirName,fileName)
      fileName = os.path.expanduser(fileName)
      structures.append( getStructureFromFile(chains[0].molSystem, fileName, doWarnings=False) )

  return makeEnsemble(structures)

# BARDIAUX 2.3
def getStructuresFromAria2Dir_Format(dirName, chains):
  """Descrn: Load structures of the input chains, from the PDB style files in an Aria2 iteration directory
     Inputs: Word - Aria2 run directory path. List of MolSystem.Chain
     Output: List of MolSystem.MolStructures
  """

  from ccpnmr.format.converters.CnsFormat import CnsFormat
  import os

  project = chains[0].root
  molSystem = chains[0].molSystem

  cnsFormat = CnsFormat(project, verbose = 0)

  cnsFormat.namingSystemName = 'XPLOR'

  num = 1
  while project.findFirstStructureEnsemble(molSystem=molSystem, ensembleId=num):
    num += 1
  ensemble = project.newStructureEnsemble(molSystem=molSystem, ensembleId=num)

  files = os.listdir(dirName)
  fileNames  = []
  for fileName in files:
    if fileName[-4:] == '.pdb':
      fileName = os.path.join(dirName,fileName)
      fileName = os.path.expanduser(fileName)
      fileNames.append(fileName)

  kw = {'minimalPrompts' : True,
        'linkAtoms' : False,
        'ignoreResNames' : ['ANI','TIP3'],
        'ignoreUnknownChemComps' : True,
        'verbose' : False}

  #if len(chains) > 1:
  mappings = _getMappings(chains)

  kw.update({'resetMapping' : True,
             'forceChainMappings' : mappings})

  cnsFormat.readCoordinates(fileNames, molSystem=molSystem, structureEnsemble=ensemble,  **kw)

  return ensemble

# BARDIAUX 2.3
getStructuresFromAria2Dir = getStructuresFromAria2Dir_Format

# BARDIAUX 2.3
def _getMappings(chains):

  mappings = []
  for chain in chains:
    seqCode = chain.sortedResidues()[0].seqCode
    code = chain.code
    mappings.append([code, code, 1, seqCode-1])
  return mappings

def getConstraintFromAria2(ariaRestraint, constraintList):
  """Descrn: Make a distance constraint from an Aria2 restraint object
     Inputs: Aria2 restraint object, NmrConstraint.DistanceConstraintList
     Output: NmrRestraints.DistanceConstraint
  """

  nmrProject = constraintList.topObject.nmrProject
  if not hasattr(constraintList, 'peakLookup'):
    constraintList.peakLookup = {}

  weight = ariaRestraint.getWeight()
  upperB = ariaRestraint.getUpperBound()
  lowerB = ariaRestraint.getLowerBound()
  rfPeak = ariaRestraint.getReferencePeak()
  active = ariaRestraint.isActive()
  target = ariaRestraint.getDistance()
  contrb = ariaRestraint.getContributions()

  # TBD deal with weight - exclude some constraints given weight threshold
  peak = None
  peaksDict = constraintList.peakLookup


  spectrumName = rfPeak.getSpectrum().getDataSource().get('ccpn_id')

  # BARDIAUX 2.2
  from aria.ariabase import is_type
  # create the CCPN DistanceConstraint first

  origData = rfPeak.getIntensity().getValue()
  if not origData:
    origData = rfPeak.getVolume().getValue()

  # BARDIAUX 2.2 backCalcVolume
  backCalcVolume = ariaRestraint.analysis.getCalculatedPeaksize()[0]

  constraint  = constraintList.newDistanceConstraint(weight=weight, error=abs(upperB-lowerB),
                                                     origData=origData, targetValue=target,
                                                     upperLimit=upperB, lowerLimit=lowerB,
                                                     backCalcVolume = backCalcVolume)


  # BARDIAUX 2.2 now deals with ConstraintPeakContrib #typo fixed by AWSS
  if spectrumName:

    ccpnNames = getKeysFromString(spectrumName) # TJS

    ## BARDIAUX 2.2
    if is_type(rfPeak.getSpectrum(), 'ConstraintList'):

      # TJS adjust below for new CCPN model:
      # NmrConstraintStores are nor childeren of memopsRoot
      # Constraint comes from a CCPN constraintList, not a spectrum
      # first, retrieve the original constraintList
      storeSerial, serial = ccpnNames

      constraintStore = constraintList.root.findFirstNmrConstraintStore(serial=storeSerial)
      if not constraintStore:
        raise ValueError('No NMR constraint store with serial "%d" in CCPN project' % storeSerial)
      else:
        origConstraintList = constraintStore.findFirstConstraintList(serial=serial)

      if origConstraintList is None:
          raise ValueError('No constraint list with serial "%d" in store "%d"' % (serial,storeSerial))

      # second, find the original DistanceConstraint
      constraintSerial = rfPeak.getNumber()
      origConstraint = origConstraintList.findFirstConstraint(serial=constraintSerial)
      # track the original constraint
      constraint.details = 'Orig Id %d' % constraintSerial

      # now, get the original ConstraintPeakContribs and transfer them to the
      # new DistanceConstraint
      origConstraintPeakContribs = origConstraint.peakContribs

      if origConstraintPeakContribs:

        for origConstraintPeakContrib in origConstraintPeakContribs:

          experimentSerial  = origConstraintPeakContrib.experimentSerial
          dataSourceSerial =  origConstraintPeakContrib.dataSourceSerial
          peakListSerial =  origConstraintPeakContrib.peakListSerial
          peakSerial =  origConstraintPeakContrib.peakSerial

          peakContrib = constraint.newConstraintPeakContrib(experimentSerial= experimentSerial,
                                                            dataSourceSerial= dataSourceSerial,
                                                            peakListSerial= peakListSerial,
                                                            peakSerial=peakSerial)


    else:
      # BARDIAUX 2.2
      # Restraints from a Spectrum

      if len(ccpnNames) == 4:
        #Fix TJS: NmrProject is part of the key
        #keys are serial num for experiment and spectrum not name
        nName, eSerial, sSerial, plSerial = ccpnNames

        experiment = nmrProject.findFirstExperiment(serial=eSerial)
        if experiment:

          spectrum = experiment.findFirstDataSource(serial=sSerial)
          if spectrum:

            peakList = spectrum.findFirstPeakList(serial=plSerial)
            if peakList:
              if peaksDict.get(peakList) is None:
                peaksDict[peakList] = {}
                for peak in peakList.peaks:
                  peaksDict[peakList][peak.serial] = peak

              peak = peaksDict[peakList].get(rfPeak.getNumber())



  if peak:
    peakContrib = constraint.newConstraintPeakContrib(experimentSerial=experiment.serial,
                                                      dataSourceSerial=spectrum.serial,
                                                      peakListSerial=peakList.serial,
                                                      peakSerial=peak.serial)

  return constraint

def getViolationFromAria2(ariaViolation, constraint, violationList):
  """Descrn: Make a constraint's violation object in the input violation list given an Aria2
             violation analysis object.
     Inputs: Aria2 violation analysis object, NmrConstraint.DistanceConstraint, NmrConstraint.ViolationList
     Output: NmrConstraint.Violation
  """

  violation  = None
  if ariaViolation.isViolated():
    calcDist  = ariaViolation.getAverageDistance().getValue()
    calcError = ariaViolation.getAverageDistance().getError()
    fracViols = max(0.0, min(ariaViolation.getDegreeOfViolation(), 1.0))
    violValue = ariaViolation.getUpperBoundViolation().getValue()
    violation = violationList.newViolation(violation=violValue,calcValue=calcDist,
                                           calcValueError=calcError,constraint=constraint,
                                           fractionViolated=fracViols)

  return violation

def getAria2SpinSystems(ariaRestraint):
  """Descrn: Get the pairs of ARIA spin systems that correspond to the contributions
             to an ARIA restraint
     Inputs: Aria2 restraint object
     Output: List of 2-List of ARIA SpinSystem objects
  """

  spinSystemPairs = []
  for contribution in ariaRestraint.getContributions():
    ariaSpinSystems = contribution.getSpinSystems()
    if len(ariaSpinSystems) == 2:
      spinSystemPairs.append(list(ariaSpinSystems))

  return spinSystemPairs


def getAria2AtomSetPairs(ariaRestraint):
  """Descrn: Get a list of pairs of aria atom objects from an Aria2 restraint
             Which correspond to the atom sets of a DistanceConstraint.
     Inputs: Aria2 restraint object
     Output: List of 2-List of Aria2 atom objects
  """

  ariaAtomPairs = []
  for contribution in ariaRestraint.getContributions():
    if contribution.getWeight() == 0:
      continue

    ariaAtomSetPair = []

    for ariaSpinSystem in contribution.getSpinSystems():
      ariaSpinSystemAtoms = ariaSpinSystem.getAtoms()
      if len(ariaSpinSystemAtoms) != 2:
        ariaAtomSetPair.append([ariaSpinSystemAtoms[0],])

      else:
        for ariaSpinSystemAtom in ariaSpinSystemAtoms:
          if ariaSpinSystemAtom.getName()[-1] == '1':
            ariaAtomSetPair.append([ariaSpinSystemAtom,])
            break

        else:
          ariaAtomSetPair.append(ariaSpinSystemAtoms)

    if len(ariaAtomSetPair) == 2:
      ariaAtomPairs.append(ariaAtomSetPair)

  return ariaAtomPairs

def getConstraintsFromAria2(objects, chains, constraintSet, structures=None, details='ARIA2 import'):
  """Descrn: Make a constraint and violation list from Aria2 restraint objects in a constraint set.
             Constraint asignments will be made relative to the input chain.
             The violation analysis will be linked to the input structures.
     Inputs: List of Aria2 restraint objects, List of MolSystem.chain,
             MolSystem.MolStructures, Nmr.NmrConstraintStore, String
     Output: NmrRestraints.DistanceConstraintList, NmrConstraint.ViolationList
  """
  # BARDIAUX: Update for multiple chains.

  resonanceDict  = {}
  constraintList = constraintSet.newDistanceConstraintList(details=details)

  if structures:
    models = structures.models
  else:
    #models = None
    models = []
    
  # TJS simplify and adjust for fewer violation lists
  violationList = constraintSet.findFirstViolationList(details=details, molStructures=models)
  if not violationList:
    violationList  = constraintSet.newViolationList(details=details, molStructures=models)

  # TJS collect inactive restraints for inspection in CCPN
  # Only make a DistanceConstraintList for these if required
  rejectConstraintList = None

  # TJS collect experiments to add to new restraint lists
  rejectExperimentDict = {}
  usedExperimentDict = {}

  for ariaRestraint in objects:

    ariaAtomPairs = getAria2AtomSetPairs(ariaRestraint)
    if not ariaAtomPairs:
      msg = "No valid atom pairs for ARIA restraint %s" % ariaRestraint.getId()
      messager.warning(msg)
      continue

#    print "Restr", ariaRestraint.getWeight(), ariaRestraint.getDistance(),  ariaRestraint.getUpperBound(), ariaRestraint.getLowerBound()

    for contribution in ariaRestraint.getContributions():
      if contribution.getWeight() == 0:
        continue

#      print "  Contrib", contribution.getId(), contribution.getWeight(), contribution.getAverageDistance()

    # TJS adjust to collect rejects
    if ariaRestraint.isActive():
      constraint = getConstraintFromAria2(ariaRestraint, constraintList)
      # TJS add
      for peakContrib in constraint.peakContribs:
        rejectExperimentDict[peakContrib.experimentSerial] = True

    else:
      # TJS we want the restraints unmerged if one is active
      # (because both have valid CCPN peaks)
      equivPeak = ariaRestraint.getEquivalentPeak()
      if ariaRestraint.isMerged() and equivPeak and equivPeak.isActive():
        constraint = getConstraintFromAria2(ariaRestraint, constraintList)
        # TJS add
        for peakContrib in constraint.peakContribs:
          usedExperimentDict[peakContrib.experimentSerial] = True

      else:
        if not rejectConstraintList:
          rejectConstraintList = constraintSet.newDistanceConstraintList(details=details + ' *REJECTS*')

        constraint = getConstraintFromAria2(ariaRestraint, rejectConstraintList)
        # TJS add
        for peakContrib in constraint.peakContribs:
          usedExperimentDict[peakContrib.experimentSerial] = True

    if not constraint:
      continue

    violation = getViolationFromAria2(ariaRestraint.analysis, constraint, violationList)

    doneItems = {}
    for ariaAtomSet1, ariaAtomSet2 in ariaAtomPairs:
      for ariaAtom1 in ariaAtomSet1:

        chain1 = getChainFromAria2Atom(ariaAtom1, chains)

        fixedResonance1 = resonanceDict.get(ariaAtom1)
        if fixedResonance1 is None:
          fixedResonance1 = getFixedResonanceFromAria2Atom(ariaAtom1, chain1, constraintSet)

        for ariaAtom2 in ariaAtomSet2:

          chain2 = getChainFromAria2Atom(ariaAtom2, chains)

          fixedResonance2 = resonanceDict.get(ariaAtom2)
          if fixedResonance2 is None:
            fixedResonance2 = getFixedResonanceFromAria2Atom(ariaAtom2, chain2, constraintSet)

          # BARDIAUX 2.2
          # check if fixedResonance1 too
          if fixedResonance2 and fixedResonance1 and ( fixedResonance2 is not fixedResonance1 ):
            serials = [fixedResonance1.serial,fixedResonance2.serial]
            serials.sort()
            serials = tuple(serials)

            if doneItems.get(serials) is None:
              item = constraint.newDistanceConstraintItem(resonances=[fixedResonance1,fixedResonance2])
              doneItems[serials] = True

    if not constraint.items:
      constraint.delete()
      #print "Constraint %d has no items" % constraint.serial
      #print "Aria restraint atoms", ariaAtomPairs

  # TJS substituted old method
  constraintList.setExperimentSerials(usedExperimentDict.keys())
  if rejectConstraintList:
    rejectConstraintList.setExperimentSerials(rejectExperimentDict.keys())

  for chain in chains:
    if hasattr(chain,'residueLookup'):
      del chain.residueLookup

  return constraintList, rejectConstraintList, violationList


def getPeakAssignmentsFromAria2(project, ariaRestraints, namesDict=None,
                                aria_chain=None):
  """Descrn: Assign a CCPN peak list (new or existing) according to ARIA restraint assignments
     Inputs: Implementation.MemopsRoot, List of ARIA Restraints, Boolean
     Output: ccp.nmr.Nmr.PeakList
  """

  # BARDIAUX: Update for multiple chains
  chains = getChainsFromAria2(ariaRestraints, project, aria_chain=aria_chain)
  # Really need to deal with multiple chains

  ariaDimDict  = {}
  peakListDict = {}
  peaksDict    = {}
  bondedDimsDict = {}

  boundResonances = {}

  for ariaRestraint in ariaRestraints:
    ariaAtomPairs = getAria2AtomSetPairs(ariaRestraint)
    ariaSpinSystemsPairs = getAria2SpinSystems(ariaRestraint)

    resonanceDimPairs = []

    refPeak = ariaRestraint.getReferencePeak()

    for i in range(len(ariaAtomPairs)):

      atomPair = ariaAtomPairs[i]
      spinSystemPair = ariaSpinSystemsPairs[i]

      atomPair = [a[0] for a in atomPair]

      # BARDIAUX: Update for multiple chains
      chain0, chain1 = getChainFromAria2Atom(atomPair[0], chains), getChainFromAria2Atom(atomPair[1], chains)

      resonance1 = getResonanceFromAria2Atom(chain0, atomPair[0])
      resonance2 = getResonanceFromAria2Atom(chain1, atomPair[1])

##       dim1 = (refPeak.getDimension(spinSystemPair[0]) -1) *2
##       dim2 = (refPeak.getDimension(spinSystemPair[1]) -1) *2
      # BB new method to retrieve spinSystem dimension (less error prone)
      dim1, dim2 = refPeak.getDimensions(*spinSystemPair)
      dim1 = (dim1-1)*2
      dim2 = (dim2-1)*2

      if dim1 == dim2:
        messager.warning("WARNING: Problem resolving dimensions of peak assignments for CCPN")
        if dim2 == 2:
          dim2 = 0
        else:
          dim2 = 2

      resonanceDimPairs.append((dim1,resonance1,dim2,resonance2))

    peakListKeys = tuple(getKeysFromString(refPeak.getSpectrum().getDataSource().get('ccpn_id')))

    if not peakListKeys:
      messager.warning("Aria Restraint ", ariaRestraint, " has no associated CCPN spectrum/peak list")
      continue

    peakList = peakListDict.get(peakListKeys)

    if peakList is None:

      peakList = getCcpnPeakList(project, peakListKeys)
      spectrum = peakList.dataSource

      if namesDict:

        new_name = namesDict[refPeak.getSpectrum().getDataSource().get('ccpn_id')]

        peakList = spectrum.newPeakList(name=new_name)
        peakList.details = new_name

      peakListDict[peakListKeys] = peakList

    else:
      spectrum = peakList.dataSource


    if peaksDict.get(peakList) is None:
      peaksDict[peakList] = {}
      for peak0 in peakList.peaks:
        peaksDict[peakList][peak0.serial] = peak0

    onebondDims = bondedDimsDict.get(spectrum)
    if onebondDims is None:
      onebondDims = {}
      for dataDim1, dataDim2 in  getOnebondDataDims(spectrum):
        onebondDims[dataDim1.dim] = dataDim2.dim
        onebondDims[dataDim2.dim] = dataDim1.dim

      bondedDimsDict[spectrum] = onebondDims

    ppmH1 = refPeak.getProton1ChemicalShift().getValue()
    ppmX1 = refPeak.getHetero1ChemicalShift().getValue()
    ppmH2 = refPeak.getProton2ChemicalShift().getValue()
    ppmX2 = refPeak.getHetero2ChemicalShift().getValue()
    ppms = [ppmH1, ppmX1, ppmH2, ppmX2]

    ariaDims = ariaDimDict.get(spectrum)
    if not ariaDims:
      ariaDims = [] #[0,1,2]

      dataDims = spectrum.sortedDataDims()
      if len(dataDims) == 3:
        for dataDim in dataDims:
          expDimRef = dataDim.findFirstDataDimRef().expDimRef
          if '1H' in expDimRef.isotopeCodes: # 0 or 2
            if onebondDims.get(dataDim.dim):
              if ppmX1 is None:
                ariaDims.append(2)
              else:
                ariaDims.append(0)

            else:
              if ppmX1 is None:
                ariaDims.append(0)
              else:
                ariaDims.append(2)

          else: # 1 or 3
            if ppmX1 is None:
              ariaDims.append(3)
            else:
              ariaDims.append(1)

      else:
        transfer = spectrum.experiment.findFirstExpTransfer(transferType='NOESY')

        for dataDim in dataDims:
          expDimRefs = [dataDimRef.expDimRef for dataDimRef in dataDim.dataDimRefs]
          i = 0
          for expDimRef in transfer.sortedExpDimRefs():
            if expDimRef in expDimRefs:
              ariaDims.append(i)
              boundDim = onebondDims.get(dataDim.dim)
              if boundDim:
                ariaDims.append(i+1)

            i += 2

      ariaDimDict[spectrum] = ariaDims

    if namesDict:
      # Fix TJS # # # # # # # # # # # # # # # #
      position = [ppms[dim] for dim in ariaDims]

      peak = pickPeak(peakList, position, unit='ppm', doFit=False)

      setManualPeakIntensity(peak, refPeak.getVolume().getValue(), intensityType='volume')

    else:
      peak = peaksDict[peakList].get(refPeak.getNumber())

    if not peak:
      data = [spectrum.experiment.name,
              spectrum.name,
              peakList.serial,
              peak.serial]
      messager.warning('Cannot find CCPN Peak %s %s %d number %d' % data)
      continue # TJS should not try to assign

    peak.details = "Orig. #: %d" % refPeak.getNumber()
    
    # TJS mark rejected peaks - would not be used next time by default
    if not ariaRestraint.isActive():
      peak.figOfMerit = 0.0
      peak.details += ', ARIA2 REJECT'

    # TJS rearranges this to put in checks for isotope mismatch
    for dim1, resonance1, dim2, resonance2 in resonanceDimPairs:
      peakAssignments = []

      ccpnDim = 0 # Tim's
      for dim in ariaDims:
        ccpnDim += 1 # Fix TJS # # # # # # # # # # # # # # # #
        peakDim = peak.findFirstPeakDim(dim=ccpnDim)

        if dim1 == dim:
          peakAssignments.append((peakDim, resonance1))
          dimX = onebondDims.get(ccpnDim)
          if dimX:
            peakDimX = peak.findFirstPeakDim(dim=dimX)
            if boundResonances.has_key(resonance1):
              bound = boundResonances[resonance1]
            else:
              bound = getBoundResonance(resonance1)
              if not bound:
                bound = resonance1.findFirstCovalentlyBound()
              boundResonances[resonance1] = bound
            ## BB due to wrong behavior of findFirstCovalentlyBound
            ## I replaced it by getBoundResonance
            #bound = resonance1.findFirstCovalentlyBound()
            if bound:
              peakAssignments.append((peakDimX, bound))

        elif dim2 == dim:
          peakAssignments.append((peakDim, resonance2))
          dimX = onebondDims.get(ccpnDim)
          if dimX:
            peakDimX = peak.findFirstPeakDim(dim=dimX)
            if boundResonances.has_key(resonance2):
              bound = boundResonances[resonance2]
            else:
              bound = getBoundResonance(resonance2)
              if not bound:
                bound = resonance2.findFirstCovalentlyBound()
              boundResonances[resonance2] = bound
            #bound = resonance2.findFirstCovalentlyBound()
            if bound:
              peakAssignments.append((peakDimX, bound))

      # TJS the actual isotope mismatch check
      for peakDim, resonance in peakAssignments:
        if resonance.isotopeCode not in peakDim.dataDimRef.expDimRef.isotopeCodes:
          messager.warning('WARNING: Spectrum %s : Attempt to assign %s resonance to peak dim %s' % \
                           (refPeak.getSpectrum().getName(),resonance.isotopeCode,peakDim))
          break

      else:
        # TJS do the peak assignment if all OK
        for peakDim, resonance in peakAssignments:
          if resonance.isotopeCode == '1H':
            tolerance=0.5
          else:
            tolerance=2.0

          assignResToDim(peakDim, resonance, tolerance=tolerance, doWarning=False)

  return peakListDict.keys()

def getResonanceFromAria2Atom(chain, ariaAtom):
  """Descrn: Get the corresponding (or new) resonance for the input Aria2 atom
             given the input chain
     Inputs: MolSystem.Chain, Aria2 atom object
     Output: Nmr.Resonance
  """

  project = chain.root
  if hasattr(project, 'ariaAtomCache'):
    cache = project.ariaAtomCache
  else:
    project.ariaAtomCache = cache = {}

  resonance = cache.get(ariaAtom)
  if resonance:
    return resonance

  if hasattr(chain,'residueLookup'):
    residueDict = chain.residueLookup
  else:
    # Dictionary for quick residue lookup
    residueDict = {}
    for residue in chain.residues:
      residueDict[residue.seqCode] = residue
    chain.residueLookup = residueDict

  tlc   = ariaAtom.getResidue().getType()
  seqId = ariaAtom.getResidue().getNumber()
  residue = residueDict.get(seqId)
  if (not residue) or (_get3LetterCcpCode(residue) != tlc):
    messager.warning('Cannot find residue %s %s' % (seqId,tlc))
    return

  atom = residue.findFirstAtom(name=ariaAtom.name)
  if not atom:
    messager.warning('Cannot find atom %s %s %s' % (seqId,tlc,ariaAtom.name))
    return

  resonance = None

  atomSet = atom.atomSet
  if not atomSet:
    messager.warning('Missing atom set %s %s %s' % (seqId,tlc,atom.name))
    # BARDIAUX if no atom set, return FixedAtomSet
    atomSet = atom.findFirstFixedAtomSet()
    if not atomSet:
      return resonance
    if not list(atomSet.resonanceSets):
      return resonance

  resonanceSets = list(atomSet.resonanceSets)
  if resonanceSets:
    for resonanceSet in resonanceSets:
      numSets = len(resonanceSet.atomSets)
      if numSets == 1:
        resonance = resonanceSet.findFirstResonance()

    if not resonance:
      resonanceSet = resonanceSets[0]
      atomSets     = resonanceSet.sortedAtomSets()
      #index        = atomSets.index(atom.atomSet)
      index        = atomSets.index(atomSet)
      resonances   = resonanceSet.sortedResonances()
      resonance    = resonances[min(index,len(resonances)-1)]

  else:
    # Make new resonance if atom not assigned
    isotope = DEFAULT_ISOTOPES.get(atom.chemAtom.chemElement.symbol, '1H')
    resonance = chain.root.currentNmrProject.newResonance(isotopeCode=isotope)    
    assignAtomsToRes([atomSet,],resonance)

  cache[ariaAtom] = resonance
  return resonance


def getFixedResonanceFromAria2Atom(ariaAtom, chain, constraintSet):
  """Descrn: Get the corresponding (or new) FixedResonance for the input Aria2 atom
             given the input chain
     Inputs: Aria2 atom object, MolSystem.Chain, Nmr.NmrConstraintHead
     Output: NmrRestraints.FixedResonance
  """

  resonance = getResonanceFromAria2Atom(chain, ariaAtom)

  # BARDIAUX 2.2
  # in case ccpnAtom has no atomSet
  # and FixedResonance has no resonance
  if not resonance:
    return resonance

  if resonance.className == 'FixedResonance':
    if not resonance.resonance:
      return copyFixedResonance(resonance, constraintSet)


  return getFixedResonance(constraintSet,resonance)

def _get3LetterCcpCode(ccpResidue):
  if ccpResidue.molType == 'DNA':
    resCode = ccpResidue.molResidue.chemComp.\
              stdChemComp.findFirstNamingSystem(name ='XPLOR').mainChemCompSysName.sysName
  else:
    resCode = ccpResidue.molResidue.chemComp.code3Letter

  return resCode


# BARDIAUX
def getChainFromAria2Atom(ariaAtom, chains):
  """Descrn: Get the corresponding MolSystem.chain for the input Aria2 atom
     Inputs: Aria2 atom object, List of MolSystem.Chain
     Output:  MolSystem.Chain
  """

  from aria.tools import string_to_segid

  project = chains[0].root
  if hasattr(project, 'ariaAtomCache'):
    cache = project.ariaAtomCache
  else:
    project.ariaAtomCache = cache = {}

  for chain in chains:

    if ariaAtom.getSegid() <> string_to_segid(chain.code):
      continue

    if hasattr(chain,'residueLookup'):
      residueDict = chain.residueLookup
    else:
      # Dictionary for quick residue lookup
      residueDict = {}
      for residue in chain.residues:
        residueDict[residue.seqCode] = residue
      chain.residueLookup = residueDict

    tlc   = ariaAtom.getResidue().getType()
    seqId = ariaAtom.getResidue().getNumber()
    residue = residueDict.get(seqId)

    if (not residue) or (_get3LetterCcpCode(residue) != tlc):
      messager.warning('Cannot find residue %s %s' % (seqId,tlc))
      continue

    atom = residue.findFirstAtom(name=ariaAtom.name)

    if not atom:
      messager.warning('Cannot find atom %s %s %s' % (seqId,tlc,ariaAtom.name))
      continue
    else:
      return chain



## Tim Stevens code
def copyFixedResonance(fixedResonance, nmrConstraintStore):
  """Descrn: Make an equivalent fixed resonance in the input constrsint store
     Inputs: NmrConstraint.FixedResonance, Nmr.NmrConstraintStore
     Output: NmrConstraint.FixedResonance
  """
  # BARDIAUX 2.2
  # this function is used in case of fixedResonance have no Nmr.Resonance
  fixedResonance2 = None

  if fixedResonance.resonanceSerial:
    fixedResonance2 = nmrConstraintStore.findFirstFixedResonance(resonanceSerial=fixedResonance.resonanceSerial)

  if fixedResonance2:
    return fixedResonance2

  serial       = fixedResonance.resonanceSerial
  name         = fixedResonance.name
  details      = fixedResonance.details
  isotope      = fixedResonance.isotopeCode
  resonanceSet = fixedResonance.resonanceSet

  if resonanceSet:
    atomSets2 = [getFixedAtomSet(nmrConstraintStore, aS.atoms) for aS in resonanceSet.atomSets]
    resonanceSet2 = nmrConstraintStore.findFirstFixedResonanceSet(atomSets=atomSets2)

    if resonanceSet2:
      index = list(resonanceSet.resonances).index(fixedResonance)
      fixedResonance2 = list(resonanceSet2.resonances)[index]


    else:
      fixedResonance2 =  nmrConstraintStore.newFixedResonance(name=name,
                                                              details=details,
                                                              isotopeCode=isotope,
                                                              resonanceSerial=serial)
      fixedResonances2 = [fixedResonance2]

      for fixedResonance1 in resonanceSet.resonances:
        if fixedResonance1 is not fixedResonance:
          fixedResonance3 = nmrConstraintStore.newFixedResonance(name=fixedResonance1.name,
                                                                 details=fixedResonance1.details,
                                                                 isotopeCode=fixedResonance1.isotopeCode,
                                                                 resonanceSerial=fixedResonance1.resonanceSerial)
          fixedResonances2.append(fixedResonance3)

      resonanceSet2 = nmrConstraintStore.newFixedResonanceSet(atomSets=atomSets2,resonances=fixedResonances2)

  else:
    fixedResonance2 =  nmrConstraintStore.newFixedResonance(name=name,
                                                            details=details,
                                                            isotopeCode=isotope,
                                                            resonanceSerial=serial)

  return fixedResonance2


def getFixedAtomSet(nmrConstraintStore, atoms):
  """Descrn: Finds or creates a fixed set of atoms that is used in an NMR constraint head object (equivalent to one
             NmrConstraint file). Creating fixed atom sets allows assignments to change but old constraints to be preserved.
     Inputs: Nmr.NmrConstraintStore, List of MolSystem.Atoms
     Output: NmrConstraint.FixedAtomSet
  """
  atoms = list(atoms) # AWSS atoms is frozenset
  if not hasattr(nmrConstraintStore, 'quickAtomSets'):
    nmrConstraintStore.quickAtomSets = {}

  fixedAtomSet = nmrConstraintStore.quickAtomSets.get(atoms)

  if not fixedAtomSet:
    fixedAtomSet = nmrConstraintStore.findFirstFixedAtomSet(atoms=atoms)

  if not fixedAtomSet:
    if atoms[0].atomSet:
      atomSet = atoms[0].atomSet
      fixedAtomSet = nmrConstraintStore.newFixedAtomSet(atoms=atomSet.atoms,name=atomSet.name)

  # BARDIAUX 2.2
  # if atom have no atomSet
  if not fixedAtomSet:
    fixedAtomSet2 = atoms[0].findFirstFixedAtomSet(atoms=atoms)
    if fixedAtomSet2:
      fixedAtomSet = nmrConstraintStore.newFixedAtomSet(atoms=fixedAtomSet2.atoms,name=fixedAtomSet2.name)

  nmrConstraintStore.quickAtomSets[atoms] = fixedAtomSet

  return fixedAtomSet


## BARDIAUX
def getBoundResonance(resonance):
  """
  return one-bound resonance (structure only)
  required since it seems that resonance.covalentlyBound is broken
  taken from core.analysis.AssignmentBasic, exp. bonds part removed
  """

  from ccpnmr.analysis.core.MoleculeBasic import getBoundAtoms

  atomResonances = {} # Linked by bound atoms irrespective of spectra
  resonanceSet   = resonance.resonanceSet

  if resonanceSet:
    #residue  = resonanceSet.findFirstAtomSet().findFirstAtom().residue
    atomSets = resonanceSet.atomSets

    for atomSet in atomSets:
      #for atom in atomSet.atoms:
      atom = atomSet.findFirstAtom()

      for atom2 in getBoundAtoms(atom):
        atomSet2 = atom2.atomSet

        if atomSet2 and atomSet2.resonanceSets:
          if len(atomSets) > 1:
            chemAtomSet = atom2.chemAtom.chemAtomSet
            if chemAtomSet:
              if chemAtomSet.isProchiral:
                continue
              if chemAtomSet.chemAtomSet and chemAtomSet.chemAtomSet.isProchiral:
                continue # LEU CD-HD*, VAL CG-HG*

          for resonanceSet2 in atomSet2.resonanceSets:
            for resonance2 in resonanceSet2.resonances:
              atomResonances[resonance2] = True

  resonances = atomResonances.keys()
  if resonances:
    return resonances[0]

  return None

# TJS: Store ARIA settings in a CCPN project
def getAriaRun(molSystem, appName='ARIA'):
  """Descrn: Retrieve/setup an NmrCalc.Run object to store details of the ARIA run in CCPN.
             Stores input and output peakLists, structures etc. as a group.
     Inputs: MolSystem.MolSystem,
     Output: NmrCalc.Run
  """

  project = molSystem.root

  if not hasattr(project, 'nmrCalcStores'):
    # CCPN installation too old
    return

  nmrProject = project.currentNmrProject or \
               project.findFirstNmrProject() or \
               project.newNmrProject(name=project.name)

  calcStore = project.findFirstNmrCalcStore(name=appName) or \
              project.newNmrCalcStore(name=appName, nmrProject=nmrProject)

  runs = []
  for run in calcStore.sortedRuns():
    if run.inputs and (run.status not in ('completed', 'failed')):
      dataObj = run.findFirstData(className='MolSystemData',
                                  ioRole='input', chainCodes=())

      if dataObj and dataObj.molSystem is molSystem:
        runs.append(run)

  if runs:
    run = runs[-1]


  else:
    # Make a new one
    run = calcStore.newRun(status='provisional')
    run.newMolSystemData(molSystemCode=molSystem.code,
                         ioRole='input', chainCodes=())


  return run
