"""
This module provides some classes for the conversion of TALOS output files to
phi/psi- and hbond restraints.

When invoked from the command line it prints to specified filenames in
the current directory otherwise it needs filenames/parameters as argument.

Author: Jens P. Linge, Institut Pasteur, July 2001
        Alexandre Bonvin, Utrecht University
"""

__author__   = "$Author: bardiaux $"
__revision__ = "$Revision: 1.1.1.1 $"
__date__     = "$Date: 2010/03/23 15:27:24 $"

import os, string, StringIO


#some utility functions, taken from Niklas Blomberg's Csi2Aria.py file:
def makePhiRestraint(aa,angle,error):
    return '! Talos derived phi restraint:\nassign (resid %3i and name C)\n       (resid %3i and name N)\n       (resid %3i and name CA)\n       (resid %3i and name C)\n       1.0 %3i %3i 2\n\n' % (aa-1,aa,aa,aa,angle,error)

def makePsiRestraint(aa,angle,error):
    return '! Talos derived psi restraint:\nassign (resid %3i and name N)\n       (resid %3i and name CA)\n       (resid %3i and name C)\n       (resid %3i and name N)\n       1.0 %3i %3i 2\n\n' %  (aa,aa,aa,aa+1,angle,error)

def makeBondRestraint(atom,parameters):
        return '! Talos derived psi restraint:\nassign (resid %3i and name %s) (resid %3i:%3i and name %s) %3.1f %3.1f %3.1f \n\n' % (atom,parameters[1],atom-parameters[0][0]-parameters[0][1],atom-parameters[0][0]+parameters[0][2],parameters[2],parameters[3],parameters[4],parameters[5])


class ReadTalos:
    """
    a parser for Talos files

    public methods:
      parseTalos        reads a talos file
      parseTalosString  reads a string containing a talos file
      parseTalosStream  reads a stream containing a talos file
      writeCNS          writes a CNS formatted restraint file, e.g. for ARIA

    attributes:
      residues          a list of residue instances
    """
    def __init__(self):
        self.residues = []
        self.talosFN = None
        self.talosString = None
        
    def parseTalos(self, talosFileName):
        """
        opens a file, feeds one big String as a stream to parseTalosString
        and returns the list from parseTalosString
        """
        self.talosFN = talosFileName
        if os.path.exists(talosFileName) == 0:
            print 'WARNING:', talosFileName, 'does not exist.'
            print '         ParseTalosmethod aborted!'
            return
        print 'reading', talosFileName
        fileStream = open(talosFileName, 'r')
        return self.parseTalosStream(fileStream)

    def parseTalosString(self, inputString):
        """
        accepts a string, converts it to a stream and uses parseTalosStream
        """
        self.talosString = inputString
        SIO = StringIO(inputString)
        return self.parseTalosStream(SIO)

    def parseTalosStream(self, inputStream):
        """
        parseTalosStream returns a list of lists which contains:
        0 residue number
        1 one-letter code
        2 phi
        3 psi
        4 delta phi
        5 delta psi
        """
        self.residues = [] #removes old residues!
        varsList = []
        formatList = []
        for line in inputStream.readlines():
            linelist=string.split(line)
            if len(linelist) < 5 or linelist[0] == 'REMARK' or \
               linelist[0] == 'DATA':
                continue
            elif linelist[0] == 'VARS':
                varsList = linelist
                #for safety reasons, do a simple check:
                expectedLine = "VARS   RESID RESNAME PHI PSI DPHI DPSI DIST COUNT CLASS"
                if (len(varsList) < 10) or \
                   varsList[1] != "RESID" or \
                   varsList[2] != "RESNAME" or \
                   varsList[3] != "PHI" or \
                   varsList[4] != "PSI" or \
                   varsList[5] != "DPHI" or \
                   varsList[6] != "DPSI" or \
                   varsList[7] != "DIST" or \
                   varsList[8] != "COUNT" or \
                   varsList[9] != "CLASS":
                    print "WARNING: The Talos format VARS description is corrupted in your file!"
                    print "         expected:"
                    print expectedLine
                    print "         found in your file:"
                    print line
                    print "WARNING: Please check your Talos input file!"
                continue
            elif linelist[0] == 'FORMAT':
                formatList = linelist
                weAreInTheHeader = 0
                continue

            #start parsing the actual data:
            RES = residue(linelist[0],\
                          linelist[1],\
                          linelist[2],\
                          linelist[3],\
                          linelist[4],\
                          linelist[5],\
                          linelist[6],\
                          linelist[7],\
                          linelist[8])
            self.residues.append(RES)

    def writeCNS(self, fileName=None, errorFactor=2.0, phiError=None, psiError=None, onlyPhiPsi=None):
        """
        writes a CNS formatted phi/psi restraint file for CNS & ARIA

        INPUT:
          fileName of the CNS dihedral restraint file
          errorFactor is multiplied on the given phi/psi Talos error ranges
          if phiError is specified, it is used as minimum error
          if psiError is specified, it is used as minimum error
          if onlyPhiPsi='phi', writes only phi restraints
          if onlyPhiPsi='psi', writes only psi restraints
        """
        try:
            file=open(fileName,'w')
        except:
            print 'WARNING: could not open ' + fileName
            print '         method writeCNS aborted!'
            return

        #1. phi & psi restraints:
        phiRestraints = ''
        psiRestraints = ''
        for iii in self.residues:
            if iii.className == 'Good':
                #1. phi:
                if phiError:
                    phiRestraints = phiRestraints +\
                                    makePhiRestraint(string.atoi(iii.resid), string.atof(iii.phi), \
                                    max(string.atof(iii.dphi)*string.atof(errorFactor),string.atof(phiError)))
                else:
                    phiRestraints = phiRestraints +\
                                    makePhiRestraint(string.atoi(iii.resid), string.atof(iii.phi), \
                                    string.atof(iii.dphi)*string.atof(errorFactor))
                #2. psi:
                if psiError:
                    psiRestraints = psiRestraints +\
                                    makePsiRestraint(string.atoi(iii.resid), string.atof(iii.psi), \
                                    max(string.atof(iii.dpsi)*string.atof(errorFactor),string.atof(psiError)))
                else:
                    psiRestraints = psiRestraints +\
                                    makePsiRestraint(string.atoi(iii.resid), string.atof(iii.psi), \
                                    string.atof(iii.dpsi)*string.atof(errorFactor))

        headerString = """! phi and psi dihedral restraint file generated by Talos2Aria.py
!
! TALOS filename:
! %s
!
! settings: min phiError=%s, min psiError=%s, errorFactor=%s
!

""" % (self.talosFN, phiError, psiError, errorFactor)
        file.write(headerString)
        if onlyPhiPsi != 'psi':
            file.write(phiRestraints)
        if onlyPhiPsi != 'phi':
            file.write(psiRestraints)
        
        file.close()
        return




class residue:
    """
    a data container for:
    RESID RESNAME PHI PSI DPHI DPSI DIST COUNT CLASS
    """
    def __init__(self, resid, resname, phi, psi, dphi, dpsi, dist, count, className):
        self.resid = resid
        self.resname = resname
        self.phi = phi
        self.psi = psi
        self.dphi = dphi
        self.dpsi = dpsi
        self.dist = dist
        self.count = count
        self.className = className



###############################################################################
if __name__ == "__main__":
    import sys, os

    args = [None, './talos_phi_psi.tbl', None, None]

    if len(sys.argv) == 1:
        print "Usage: Talos2Aria talos_file [cns_file] [phi_error] [psi_error]"
        sys.exit(1)
        
    else:
        args[:len(sys.argv) - 1] = sys.argv[1:]
        
    RT = ReadTalos()
    RT.parseTalos(args[0])
    print 'writing CNS/ARIA dihedral restraint file: %s' % args[1]
    RT.writeCNS(args[1],  errorFactor=2.0, phiError = args[2], psiError = args[3])
        
