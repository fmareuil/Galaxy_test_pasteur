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
from aria.Settings import Settings
from aria.TypeChecking import check_type

INVALID_SHIFT_VALUE = 'Invalid chemical shift value'
NEGATIVE_SHIFT_ERROR = 'Negative chemical shift error'


class ChemicalShiftFilterSettings(Settings):

    def create(self):

        from aria.Settings import ChoiceEntity
        from aria.TypeChecking import NONE, FLOAT

        keywords = {'ppm_type': ChoiceEntity((NONE, FLOAT, (NONE, FLOAT),
                                              (FLOAT, NONE)))}
        return keywords

    def create_default_values(self):

        from aria.TypeChecking import NONE, FLOAT

        return {'ppm_type': (NONE, FLOAT)}

class ChemicalShiftFilter(AriaBaseClass):

    def __init__(self, settings):

        check_type(settings, 'ChemicalShiftFilterSettings')

        AriaBaseClass.__init__(self)

        self.setSettings(settings)

    def __call__(self, shift):

        check_type(shift, 'ChemicalShift')

        from aria.tools import as_tuple
        from aria.TypeChecking import is_type
        
        value, error = shift

        messages = []
        valid = 1

        ppm_types = list(as_tuple(self.getSettings()['ppm_type']))

        if not 1 in [is_type(value, t) for t in ppm_types]:
            messages.append(INVALID_SHIFT_VALUE)
            valid = 0

        if error is not None and error < 0.0:
            messages.append(NEGATIVE_SHIFT_ERROR)
            valid = 0

        self.result = {shift: messages}
        
        return valid
