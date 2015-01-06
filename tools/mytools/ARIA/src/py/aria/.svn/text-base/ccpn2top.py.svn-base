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

Updated and almost completely recoded by Tim Stevens on 1st Nov 2006

$Author: bardiaux $
$Revision: 1.1.1.1 $
$Date: 2010/03/23 15:27:24 $
"""

import os, sys

from aria.Topology import Topology, BaseResidueSettings, BaseResidue, \
                     BaseAtomSettings, BaseAtom, \
                     TYPE_AMINO_ACID, TYPE_DNA_BASE, TYPE_RNA_BASE, \
                     EQUIV_NTERMINUS, EQUIV_METHYL, EQUIV_AROMATIC, \
                     EQUIV_METHYLENE, EQUIV_ISOPROPYL, \
                     EquivalentGroupSettings, EquivalentGroup

from ccp.general.Io import getChemComp

from memops.api     import Implementation


ariaPath = os.environ['ARIA2'] + '/src/py/'
if not ariaPath in sys.path:
  sys.path.insert(0, ariaPath)

 
ccpCodes = {'protein':('ALA', 'ARG', 'ASN', 'ASP',                                    
                       'CYS', 'GLN', 'GLU', 'GLY',                                    
                       'HIS', 'ILE', 'LEU', 'LYS',                                    
                       'MET', 'PHE', 'PRO', 'SER',                                    
                       'THR', 'TRP', 'TYR', 'VAL'),
            'RNA':('A', 'C', 'G', 'U'),
            'DNA':('A', 'C', 'G', 'T')
            }

#PROTEIN_BACKBONE = ('N', 'CA', 'C')
#DNA_BACKBONE = ("C1'", "H1'", "C2'", "H2'", "C3'", "O3'", "H3'", "C4'",
#                "O4'", "H4'", "O5'", "C5'", "H5'", "H5''", "H2''")
#RNA_BACKBONE = ("C1'", "H1'", "C2'", "H2'", "C3'", "O3'", "H3'", "C4'",
#                "O4'", "H4'", "O5'", "C5'", "H5'", "H5''", "O2'", "HO2'")


ariaResidueCode = {}
ariaResidueCode['A'] = 'ADE'
ariaResidueCode['C'] = 'CYT'
ariaResidueCode['G'] = 'GUA'
ariaResidueCode['T'] = 'THY'
ariaResidueCode['U'] = 'URI'
for ccpCode in ccpCodes['protein']:
  ariaResidueCode[ccpCode] = ccpCode

residueTypeDict = {'protein':TYPE_AMINO_ACID,
                   'RNA': TYPE_RNA_BASE,
                   'DNA':TYPE_DNA_BASE}


def createAriaTopology():
  """
  Extract ARIA topology from CCPN project.
  """
  
  ariaTopology = Topology()
  ccpProject   = Implementation.Project(name='tmp')

  for molType in ccpCodes.keys():
    for ccpCode in ccpCodes[molType]:
      ariaResidue = createAriaResidue(ccpnProject,ccpCode,molType)
      ariaTopology.addResidue(ariaResidue)
      
  return ariaTopology



def createAriaResidue(ariaResidueType, ccpProject, ccpCode,
                      ccpMoltype, backbone):


  chemComp = getChemComp(ccpProject, ccpMoltype, ccpCode)
  
  ariaSettings = BaseResidueSettings()
  ariaSettings['type'] = residueTypeDict[ccpMoltype]
  ariaResidue = BaseResidue(ariaSettings, name=ariaResidueCode[ccpCode])
  
  ariaAtomDict = {}
  doneChemAtomSets = {}
  for chemCompVar in chemComp.chemCompVars:
    if chemCompVar.linking in ('start', 'middle', 'end'):
      for chemAtom in chemCompVar.chemAtoms:
        if chemAtom.className == 'LinkAtom':
          continue
 
        name = chemAtom.name
        if ariaAtomDict.get(name):
          continue
 
        ariaAtom = createAriaAtom(chemAtom)
        if ariaAtom is None:
          print 'WARNING: Could not create atom "%s" for residue type "%s' % (name, chemComp.ccpCode)
          continue
          
        ariaAtomDict[name] = ariaAtom
    
      setupEquivalentGroup(chemCompVar, ariaResidue, ariaAtomDict, doneChemAtomSets)


  return aria_residue


def createAriaAtom(ccpChemAtom, heteroElements=('N', 'C')):
  
  elementSymbol = ccpChemAtom.elementSymbol

  heteroName = None
  if elementSymbol == 'H':
    for bound in ccpChemAtom.findFirstChemBond().chemAtoms:
      if bound.elementSymbol in heteroElements:
        heteroName = bound.name
        break

  ariaSettings  = BaseAtomSettings()
  ariaSettings['type'] = elementSymbol
  ariaSettings['hetero_atom_name' ] = heteroName
  
  ariaAtom = BaseAtom(ariaSettings, name=ccpChemAtom.name)
 
  return ariaAtom
 

def setupEquivalentGroup(chemCompVar, ariaResidue, ariaAtomDict, done=None):

  if done is None:
    done = {}

  for chemAtomSet in chemCompVar.chemAtomSets:
    if done.get(chemAtomSet):
      continue
    done[chemAtomSet] = True

    chemAtom = chemAtomSet.findFirstChemAtom()

    if chemAtom is None or chemAtom.elementSymbol <> 'H':
      continue

    ariaAtoms = []
    chemAtomSets = chemAtomSet.chemAtomSets or [chemAtomSet,]
    for chemAtomSet0 in chemAtomSets:
      ariaAtoms += [ariaAtomDict[chemAtom.name] for chemAtom in chemAtomSet.chemAtoms]
    
    if chemAtomSet.isEquivalent is True:
      if (chemCompVar.linking == 'start') and (chemAtomSet.name == 'H*'):
        equivType = EQUIV_NTERMINUS
      else:
        equivType = EQUIV_METHYL
      
    elif chemAtomSet.isEquivalent is None:
      equivType = EQUIV_AROMATIC
    
    else:#if chemAtomSet.isEquivalent is False:
    
      if chemAtomSet.chemAtomSets:
        equivType = EQUIV_ISOPROPYL
    
      elif chemAtomSet.isProchiral or len(chemAtomSet.chemAtoms) == 2:
        equivType = EQUIV_METHYLENE
      
      else: # Shouldn't really get here
        equivType = EQUIV_AROMATIC
    
    settings = EquivalentGroupSettings()
    settings['atom_names'] = tuple([a.getName() for a in ariaAtoms])
    settings['type'] = equivType
    
    equivalenceGroup = EquivalentGroup(settings)
    equivalenceGroup.setAtoms(ariaAtoms)
    
    ariaResidue.addEquivalentGroup(equivalenceGroup)


if __name__ == '__main__':

  from aria.AriaXML import AriaXMLPickler
  
  ariaPickler = AriaXMLPickler()

  ariaTopology = createAriaTopology()

  ariaPickler.dump(ariaTopology, '/tmp/t')
