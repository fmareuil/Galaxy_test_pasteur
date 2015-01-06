"""
A module to convert atomnames to IUPAC nomenclature

authors: Jens Linge, EMBL & Institut Pasteur
         Jurgen Doreleijers, BMRB
"""
__author__   = "$Author: bardiaux $"
__revision__ = "$Revision: 1.1.1.1 $"
__date__     = "$Date: 2010/03/23 15:27:24 $"

import re, sys, os, string
import AminoAcid

## Initializing some variables
## Use unix notation as a standard
## Dir with aria distribution 
## Jens, can this be derived from this files' location or so???
BaseDir                 = '/usr/home/jurgen/aria/aria1.0/Aria'

## Filename with path to atom name library having XPLOR and IUPAC 
AtomLib_Xplor_FileName  = BaseDir + '/Nomenclature/AtomLIB-xplor'

def ConvertCnsProtonNames(residueName, atomName):
    """
    convert an atomname from XPLOR/CNS to IUPAC nomenclature (or vice versa)
    residueName: a string which contains 1- or 3-letter code, e.g. 'A' or 'ALA'
                 only the 20 common aminoacids are supported!
    atomName:    a string containing the atomnames, e.g. 'HG12'

    returns a string with the new atomname (all characters are uppercase)
    If the atom name doesn't have to be changed, it will return the input
    atom name (stripped and uppercase)
    """
    #I. get a clean three-letter code and strip & uppercase the atomName
    threeLetter = AminoAcid.AminoAcid(residueName)[1]
    if threeLetter[2] == '':
        print 'WARNING: residue name', residueName, 'not understood'
        return atomName
    atomName = string.upper(string.strip(atomName))
    
    #II. methylenes
    #1. GLY HA:
    if threeLetter == 'GLY' and atomName == 'HA1':
        atomName = 'HA2'
    elif threeLetter == 'GLY' and atomName == 'HA2':
        atomName = 'HA1'
        
    #2. ARG, ASN, ASP, CYS, GLN, GLU, HIS, LEU, LYS, MET, PHE, PRO, SER, TRP, TYR HB%:
    elif threeLetter in ('ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'LEU', 'LYS',\
                         'MET', 'PHE', 'PRO', 'SER', 'TRP', 'TYR') and \
                         atomName == 'HB3':
        atomName = 'HB1'
    elif threeLetter in ('ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'LEU', 'LYS',\
                         'MET', 'PHE', 'PRO', 'SER', 'TRP', 'TYR') and \
                         atomName == 'HB1':
        atomName = 'HB3'

    #3. ARG, GLN, GLU, LYS, MET, PRO HG%:
    elif threeLetter in ('ARG', 'GLN', 'GLU', 'LYS', 'MET', 'PRO') and\
         atomName == 'HG1':
        atomName = 'HG3'
    elif threeLetter in ('ARG', 'GLN', 'GLU', 'LYS', 'MET', 'PRO') and\
         atomName == 'HG3':
        atomName = 'HG1'
    #4. ILE HG1%:
    elif threeLetter == 'ILE' and atomName == 'HG13':
        atomName = 'HG11'
    elif threeLetter == 'ILE' and atomName == 'HG11':
        atomName = 'HG13' 
    #5. ARG, ASN, LYS, PRO HD:
    elif threeLetter in ('ARG', 'ASN', 'LYS', 'PRO') and atomName == 'HD1':
        atomName = 'HD3'
    elif threeLetter in ('ARG', 'ASN', 'LYS', 'PRO') and atomName == 'HD3':
        atomName = 'HD1'
    #6. LYS HE:
    elif threeLetter == 'LYS' and atomName == 'HE3':
        atomName = 'HE1'
    elif threeLetter == 'LYS' and atomName == 'HE1':
        atomName = 'HE3'
        
    #III. methyls:
    #1. ALA beta:
    elif threeLetter == 'ALA' and atomName == 'HB2':
        atomName = 'HB1'
    elif threeLetter == 'ALA' and atomName == 'HB1':
        atomName = 'HB2'
    #2. VAL gamma1:
    elif threeLetter == 'VAL' and atomName == 'HG11':
        atomName = 'HG12'
    elif threeLetter == 'VAL' and atomName == 'HG12':
        atomName = 'HG11'
    #3. ILE, VAL gamma2:
    elif threeLetter in ('ILE', 'VAL') and atomName == 'HG21':
        atomName = 'HG22'
    elif threeLetter in ('ILE', 'VAL') and atomName == 'HG22':
        atomName = 'HG21'
    #4. ILE, LEU delta1:
    elif threeLetter in ('ILE', 'LEU') and atomName == 'HD11':
        atomName = 'HD12'
    elif threeLetter in ('ILE', 'LEU') and atomName == 'HD12':
        atomName = 'HD11'    
    #5. LEU delta2:
    elif threeLetter == 'LEU' and atomName == 'HD21':
        atomName = 'HD22'
    elif threeLetter == 'LEU' and atomName == 'HD22':
        atomName = 'HD21'    
    #6. MET epsilon:
    elif threeLetter == 'MET' and atomName == 'HE1':
        atomName = 'HE2'
    elif threeLetter == 'MET' and atomName == 'HE2':
        atomName = 'HE1'
    #7. zeta:
    elif atomName == 'HZ1':
        atomName = 'HZ2'
    elif atomName == 'HZ2':
        atomName = 'HZ1'     
        
    #IV. ARG NHs:
    elif threeLetter == 'ARG' and atomName == 'HH11':
        atomName = 'HH12'
    elif threeLetter == 'ARG' and atomName == 'HH12':
        atomName = 'HH11'
    elif threeLetter == 'ARG' and atomName == 'HH21':
        atomName = 'HH22'
    elif threeLetter == 'ARG' and atomName == 'HH22':
        atomName = 'HH21'    

    return atomName


def ConvertCnsPseudoAtomNames(residueName, atomName):
    """
    convert an pseudo atomname from XPLOR/CNS to IUPAC nomenclature (or vice versa)
    residueName: a string which contains full code, e.g. 'A' is not 'ALA'
                 but is allowed for DADE and RADE (this needs work)
                 only the 20 common aminoacids are supported!
    atomName:    a string containing the atomnames, e.g. 'HG1#'

    returns a string with the new atomname (all characters are uppercase)
    If the atom name doesn't have to be changed, it will return the input
    atom name (stripped and uppercase)
    
    e.g. 	ALA HB#  -> MB
            LEU HD1# -> MD1
            VAL HG#  -> QG
            
    This last example might be more than expected since HG# doesn't expand
    to HG11 and the 5 other protons.
    """
    #I. get a clean three-letter code and strip & uppercase the atomName
    residueCode = AminoAcid.AminoAcid(residueName)[1]
    atomName = string.upper(string.strip(atomName))
    
    #II. Set up the mapping table, this should of course be done only once
    ## and the results stored in a global parameter
    ## Jens, do you know how to set that up? I really prefer to have the 
    ## info outside the source code.
    iupac_name = {}    
    
    ## Regular expression pattern matching just the definitions from the
    ## AQUA atom name library
    pattern = re.compile(r"""
        ^def      \s+      # Start with def
         (\w+)    \s+      # Residue name
         \*       \s+      # Actual '*'
         (\w+\#)  \s+      # CNS pseudo atom name with trailing '#'
         (\w+)    \s*$     # IUPAC pseudo atom name
             """, re.IGNORECASE | re.MULTILINE | re.VERBOSE )
        
    ## Read file contents, leave file handle open later on?
    file_content = open(AtomLib_Xplor_FileName, 'r').read()
    match_list = pattern.findall( file_content ) 
    if ( match_list ):
        print "Read %d definitions from file: %s"         % ( len(match_list), AtomLib_Xplor_FileName )
    else:
        print "ERROR: No definitions read from file: %s"    % AtomLib_Xplor_FileName
        sys.exit(1)
    for match in match_list:
        iupac_name[ match[0] ] = { match[1]: match[2] }
    
    if ( iupac_name[ residueCode ].has_key( atomName ) ):
        print "Changed atom name: %s to: " % atomName
        atomName = iupac_name[ residueCode ][atomName]
        
    return atomName

if __name__ == "__main__":
    print ConvertCnsPseudoAtomNames( "ALA", "HB#" )
