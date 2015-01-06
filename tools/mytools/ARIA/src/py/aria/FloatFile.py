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


from aria.ariabase import AriaBaseClass as _AriaBaseClass

class FloatFile(_AriaBaseClass):

    def parse(self, file):

        from aria.tools import string_to_segid

        import re

        atom = 'segid "(?P<segid%(i)d>.*)" and ' + \
               'resid (?P<residue%(i)d>[0-9]+).*and ' + \
               'name (?P<atom%(i)d>H[A-Z0-9]+)'

        line = 'REVE.*\(\(.*%s.*\).*OR.*\(.*%s.*\)\)' \
               % (atom % {'i': 1}, atom % {'i': 2})
        
        regex = re.compile(line)

        table = regex.findall(open(file).read())

        swapped_atoms = {}

        for row in table:

            row = [x.strip() for x in row]
            row = [f(x) for f, x in zip([str, int, str, str, int, str], row)]

            row[0] = string_to_segid(row[0])
            row[3] = string_to_segid(row[3])

            key = tuple(row[:3])
            value = tuple(row[3:])

            if swapped_atoms.has_key(key):                
                m = 'Inconsistency: atom "%s" swapped twice.' % str(key)
                self.error(ValueError, m)

            swapped_atoms[key] = value
            
        return swapped_atoms
    
