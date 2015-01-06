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


UPPER_BOUND_EXT = 'upl'
LOWER_BOUND_EXT = 'lol'

from aria.ariabase import AriaBaseClass
from aria.TypeChecking import check_type, check_elements, check_string, LIST, TUPLE

class MolMolPickler(AriaBaseClass):

    line_format = '%(number1)3i %(residue1)-4s %(atom1)-5s ' + \
                  '%(number2)3i %(residue2)-4s %(atom2)-7s ' + \
                  '%(bound)-5.2f %(weight)-5.2f ' 

    def writePeak(self, peak):

        check_type(peak, 'AriaPeak')

        ## create lower and upper bound entries
        
        bound = {'bound': None}

        lower = peak.getLowerBound()
        upper = peak.getUpperBound()

        lower_bounds = []
        upper_bounds = []
        
        for contribution in peak.getContributions():

            weight = contribution.getWeight()

            if weight == 0.0: continue

            bound['weight'] = weight

            for pair in contribution.getSpinPairs():

                a1, a2 = pair.getAtoms()

                bound['residue1'] = a1.getResidue().getType()
                bound['number1'] = a1.getResidue().getNumber()
                bound['atom1'] = a1.getName()
                
                bound['residue2'] = a2.getResidue().getType()
                bound['number2'] = a2.getResidue().getNumber()
                bound['atom2'] = a2.getName()

            bound['bound'] = lower
            lower_bounds.append(self.line_format % bound)

            bound['upper'] = upper
            upper_bounds.append(self.line_format % bound)

        return lower_bounds, upper_bounds

    def writeBounds(self, peaks, lower_bound_file, upper_bound_file, gzip = 0):
        
        import os

        check_type(peaks, LIST, TUPLE)
        check_elements(peaks, 'AriaPeak')
        check_string(lower_bound_file)
        check_string(upper_bound_file)

        lower_bounds = []
        upper_bounds = []

        ## create entries for each peak

        for peak in peaks:

            l, u = self.writePeak(peak)
            
            lower_bounds += l
            upper_bounds += u

        ## write restraints to .lol and .upl files

        if os.path.exists(lower_bound_file):
            w = 'Lower bound file "%s" and will be overwritten.'
            self.warning(w % lower_bound_file)

        if os.path.exists(upper_bound_file):
            w = 'Upper bound file "%s" and will be overwritten.'
            self.warning(w % upper_bound_file)
            
        if gzip:
            from aria.tools import gzip_open as open_func
        else:
            open_func = open

        file = open_func(os.path.expanduser(lower_bound_file), 'w')
        file.write('\n'.join(lower_bounds) + '\n')
        file.close()

        file = open_func(os.path.expanduser(upper_bound_file), 'w')
        file.write('\n'.join(upper_bounds) + '\n')
        file.close()

def write_noe_restraints(peaks, dir, gzip = 0):

    check_type(peaks, LIST, TUPLE)
    check_elements(peaks, 'AriaPeak')

    import os

    pickler = MolMolPickler()

    ## sort peaks according to the spectrum of the reference peak

    spectra = {}

    for peak in peaks:

        spectrum = peak.getReferencePeak().getSpectrum().getName()
        if not spectra.has_key(spectrum): spectra[spectrum] = []
        spectra[spectrum].append(peak)

    ## write .lol and .upl files for each spectrum

    for spectrum, peaks in spectra.items():

        root = os.path.join(dir, spectrum)

        ## unambiguous peaks

        file = root + '_unambiguous'

        unambiguous = [p for p in peaks if not p.isAmbiguous()]

        accepted = [p for p in unambiguous if not p.analysis.isViolated()]
        lol_file = file + '_acc.' +  LOWER_BOUND_EXT
        upl_file = file + '_acc.' +  UPPER_BOUND_EXT
        pickler.writeBounds(accepted, lol_file, upl_file, gzip)
                            
        rejected = [p for p in unambiguous if p.analysis.isViolated()]
        lol_file = file + '_rej.' +  LOWER_BOUND_EXT
        upl_file = file + '_rej.' +  UPPER_BOUND_EXT
        pickler.writeBounds(rejected, lol_file, upl_file, gzip)

        ## ambiguous peaks

        file = root + '_ambiguous'

        ambiguous = [p for p in peaks if p.isAmbiguous()]

        accepted = [p for p in ambiguous if not p.analysis.isViolated()]
        lol_file = file + '_acc.' +  LOWER_BOUND_EXT
        upl_file = file + '_acc.' +  UPPER_BOUND_EXT
        pickler.writeBounds(accepted, lol_file, upl_file, gzip)
                            
        rejected = [p for p in ambiguous if p.analysis.isViolated()]
        lol_file = file + '_rej.' +  LOWER_BOUND_EXT
        upl_file = file + '_rej.' +  UPPER_BOUND_EXT
        pickler.writeBounds(rejected, lol_file, upl_file, gzip)

        
