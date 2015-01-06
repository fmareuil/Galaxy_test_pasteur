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


TOO_MANY_SHIFTS = 'Maximum number of chemical shifts exceeded.'
INVALID_SHIFT_VALUES = 'At least one chemical shift has to be not None.'
EMPTY_SPIN_SYSTEM = 'No atoms in spin system.'

from aria.ariabase import AriaBaseClass
from aria.Settings import Settings, Entity
from aria.TypeChecking import check_type

class SpinSystemFilterSettings(Settings):

    def create(self):

        from aria.Settings import PositiveInteger, ChoiceEntity
        from aria.TypeChecking import FLOAT, NONE

        ppm_type = ChoiceEntity((FLOAT, NONE, (FLOAT, NONE), (NONE, FLOAT)))

        keywords = {'ppm_type': ppm_type,
                    'max_no_shifts': PositiveInteger()}

        return keywords

class ShiftAssignmentFilterSettings(Settings):

    def create(self):

        from aria.Settings import PositiveInteger, ChoiceEntity
        from aria.TypeChecking import FLOAT, NONE

        ppm_type = ChoiceEntity((FLOAT, NONE, (FLOAT, NONE), (NONE, FLOAT)))
        
        keywords = {'ppm_type': ppm_type, 
                    'max_no_shifts': PositiveInteger()}

        return keywords

class ChemicalShiftListFilterSettings(Settings):

    def create(self):

        from aria.Settings import PositiveInteger, ChoiceEntity
        from aria.TypeChecking import FLOAT, NONE
        
        ppm_type = ChoiceEntity((FLOAT, NONE, (FLOAT, NONE), (NONE, FLOAT)))
        
        keywords = {'ppm_type': ppm_type,
                    'max_no_shifts': PositiveInteger()}

        return keywords

class SpinSystemFilter(AriaBaseClass):

    def __init__(self, settings):

        check_type(settings, 'SpinSystemFilterSettings')

        AriaBaseClass.__init__(self)

        self.setSettings(settings)

    def __filter_chemical_shifts(self, spin_system):

        check_type(spin_system, 'SpinSystem')

        result = {}

        ## filter each chemical-shift in the shift-list

        from aria.ChemicalShiftFilter import ChemicalShiftFilter, \
             ChemicalShiftFilterSettings

        settings = ChemicalShiftFilterSettings(default_settings =
                                               self.getSettings())
        filter = ChemicalShiftFilter(settings)

        valid = 1

        for shift in spin_system.getChemicalShifts():
            valid &= filter(shift)
            result.update(filter.result)

        return valid, result

    def __call__(self, spin_system):

        check_type(spin_system, 'SpinSystem')

        import aria.ShiftAssignment as SA

        valid = 1
        result = {}

        ## first check chemical-shift list

        shifts = spin_system.getChemicalShifts()
        max_no = self.getSettings()['max_no_shifts']

        if not 0 < len(shifts) <= max_no:
            valid = 0
            result['n_shifts'] = TOO_MANY_SHIFTS
        else:
            result['n_shifts'] = None

        ## works because shifts are sorted according to their ppm-values

        if len(shifts) == 1 and shifts[0][0] is None:
            valid = 0
            result['ppm_values'] = INVALID_SHIFT_VALUES
        elif len(shifts) == 2 and shifts[0][0] is None and shifts[1][0] is None:
            valid = 0
            result['ppm_values'] = INVALID_SHIFT_VALUES
        else:
            result['ppm_values'] = None

##         if not shifts[0][0]:
##             valid = 0
##             result['ppm_values'] = INVALID_SHIFT_VALUES
##         else:
##             result['ppm_values'] = None

        v, r  = self.__filter_chemical_shifts(spin_system)
        valid &= v
        result.update(r)

        ## filter by no. of atoms

        result['spin_system'] = []

        atoms = spin_system.getAtoms()

        if not len(atoms):
            valid = 0
            result['spin_system'].append(EMPTY_SPIN_SYSTEM)

        elif len(atoms) == 1:

            method = spin_system.getAveragingMethod()

            if not method == SA.AVERAGING_METHOD_NONE:

                valid = 0
                m = ('Averaging "%s" not possible if spin system ' + \
                     'contains only one atom.') % method
                result['spin_system'].append(m)

        ## check if atoms belong to a group of (potentially)
        ## equivalent atoms

        if len(atoms) > 1:

            method = spin_system.getAveragingMethod()

            if not method == SA.AVERAGING_METHOD_FAST and \
               not method == SA.AVERAGING_METHOD_SLOW:

                valid = 0
                m = ('Averaging "%s" not possible for spin system ' + \
                     'consisting of %d atoms.') % (method, len(atoms))
                result['spin_system'].append(m)

            atom_names = [a.getName() for a in atoms]
            atom_names.sort()

            for atom in atoms:

                has_equivalent_group = 0

                for group in atom.getEquivalentGroups():

                    equivalent_atoms = list(group.getAtomNames())
                    equivalent_atoms.sort()

                    if atom_names == equivalent_atoms:
                        has_equivalent_group = 1
                        break

                if not has_equivalent_group:
                    
                    valid = 0
                    m = 'Atoms "%s" do not define an equivalent group.'
                    result['spin_system'].append(m % str(atoms))
                    break

        self.result = {spin_system: result}

        return valid
            
class ShiftAssignmentFilter(AriaBaseClass):

    def __init__(self, settings):

        check_type(settings, 'ShiftAssignmentFilterSettings')

        AriaBaseClass.__init__(self)

        self.setSettings(settings)

    def __filter_spin_systems(self, assignment):

        settings = SpinSystemFilterSettings()
        settings.update(self.getSettings())
        filter = SpinSystemFilter(settings)

        valid = 1
        result = {}

        for spin_system in assignment.getSpinSystems():
            valid &= filter(spin_system)
            result.update(filter.result)

        return valid, result

    def __check_consistency(self, assignment):

        import aria.ShiftAssignment as SA

        method = assignment.getMethod()

        spin_systems = assignment.getSpinSystems()

        valid = 1
        result = []

        if method == SA.ASSIGNMENT_METHOD_STEREO_SPECIFIC:

            if not len(spin_systems) == 1:

                m = 'Multiple spin systems in stereo-specific assignment.'
                result.append(m)
                valid = 0

            elif not len(spin_systems[0].getAtoms()) == 1:

                m = 'Multiple atoms in stereo-specifically assigned ' + \
                    'spin system.'
                result.append(m)
                valid = 0

        elif method == SA.ASSIGNMENT_METHOD_EQUIVALENT:

            if not len(spin_systems) == 1:

                m = 'Multiple spin systems found in degenerate atom ' + \
                    'assignment.'
                result.append(m)
                valid = 0

            else:

                method = spin_systems[0].getAveragingMethod()

                if not method == SA.AVERAGING_METHOD_FAST:

                    m = 'Averaging "%s" specified but only "%s" ' + \
                        'averaging possible.'
                    result.append(m % (method, SA.AVERAGING_METHOD_FAST))
                    valid = 0

        elif method == SA.ASSIGNMENT_METHOD_FLOATING:

            if not len(spin_systems) == 2:
                m = 'Ambiguous assignment only possible for two ' + \
                    'spin systems.'
                result.append(m)
                valid = 0

            elif [xx[0] for xx in spin_systems[0].getChemicalShifts()] <> \
                 [xx[0] for xx in spin_systems[1].getChemicalShifts()]:
                m = 'Chemical shift lists of floating chemical shift assignments must be equal'
                result.append(m)
                valid = 0

        return valid, result

    def __call__(self, assignment):

        check_type(assignment, 'ShiftAssignment')

        result = []
        valid = 1
        
        if not len(assignment.getSpinSystems()):

            valid = 0
            result.append('Empty shift assignment.')

        v, r = self.__check_consistency(assignment)

        if not v:

            valid = 0
            result.append(r)

        v, r = self.__filter_spin_systems(assignment)

        if not v:

            valid = 0
            [result.append('%s: %s' % (str(k), r[k])) for k in r.keys()]

        self.result = {assignment: result}

        return valid

class ChemicalShiftListFilter(AriaBaseClass):

    def __init__(self, settings):

        check_type(settings, 'ChemicalShiftListFilterSettings')

        AriaBaseClass.__init__(self)

        self.setSettings(settings)

    def __call__(self, shifts):

        check_type(shifts, 'ChemicalShiftList')

        assignment_filter = ShiftAssignmentFilter(\
            ShiftAssignmentFilterSettings())

        assignment_filter.getSettings().update(self.getSettings())

        result = {}

        for assignment in shifts.getShiftAssignments():

            valid = assignment_filter(assignment)
            assignment.is_valid(valid)

            if not valid: result.update(assignment_filter.result)

        self.result = result
        
        return [a for a in shifts.getShiftAssignments() if a.is_valid()]

