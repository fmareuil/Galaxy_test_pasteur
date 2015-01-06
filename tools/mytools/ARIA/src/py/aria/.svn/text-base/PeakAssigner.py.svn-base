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



from aria.ariabase import *
from aria.Settings import Settings
from numpy import *
from aria.xmlutils import XMLElement, XMLBasePickler
from aria.ShiftAssignment import *

def peak_shape(deltas, errors):

    exponent = - 0.5 * (deltas / errors) ** 2

    return exp(clip(exponent, -709., 709.))

class PeakAssignerSettings(Settings):

    def create(self):

        from aria.Settings import NonNegativeFloat
        from aria.Settings import YesNoChoice, FourLetterString, String

        p = 'Proton-shift error must be a non-negative floating point number.'
        h = 'Hetero-atom-shift error must be a non-negative ' + \
            'floating point number.'
        a = 'Assignment-shift error must be a non-negative' + \
            ' floating point number.'

        f = lambda s, c = NonNegativeFloat: c(error_message = s)

        ## TODO: xxxx_shift_err is actually misleading.
        ## should be sth like proton2_freq_window_size 

        keywords = {'proton1_shift_err': f(p),
                    'hetero1_shift_err': f(h),
                    'proton2_shift_err': f(p),
                    'hetero2_shift_err': f(h),
                    'default_shift_error': f(a),
                    'use_assignments': YesNoChoice(),
                    'ref_segid' : FourLetterString(),
                    'spec_type' : String(),
                    'sym_type' : String(),
                    'structural_rules_enabled' : YesNoChoice()}

        return keywords

    def create_default_values(self):
        d = {'use_assignments': YES,
             'ref_segid' : "    ",
             'spec_type' : "intra",
             'sym_type' : "None",
             'structural_rules_enabled' : NO}

        return d

    def __str__(self):

        s =  """Proton1 freq. window size [ppm]: %.3f,
Proton2 freq. window size [ppm]: %.3f,
Hetero1 freq. window size [ppm]: %.3f,
Hetero2 freq. window size [ppm]: %.3f,
Default shift error [ppm]:       %.3f,
Use existing assignments:        %s"""

        return s % (self['proton1_shift_err'], self['proton2_shift_err'],
                    self['hetero1_shift_err'], self['hetero2_shift_err'],
                    self['default_shift_error'], self['use_assignments'])

class PeakAssigner(AriaBaseClass):

    contribution_counter = 0

    def __init__(self, settings):
        """
        shift_assignments is expected to be a list/tuple
        of ShiftAssignments.
        """

        check_type(settings, 'PeakAssignerSettings')
        AriaBaseClass.__init__(self, settings = settings)

    def setShiftAssignments(self, shift_assignments):
        
        check_type(shift_assignments, LIST, TUPLE)
        check_elements(shift_assignments, 'ShiftAssignment')

        self.decompose_shift_assignments(shift_assignments)

    def decompose_shift_assignments(self, assignment_list):

        ## compile list of shifts and spin_systems
        
        ## in case of a floating assignment, we assign
        ## chemical shifts according to the following rule:

        ## spin-system1:   shifts: 1, 2    atom: A
        ## spin-system2:   shifts: 1, 2    atom: B

        ## assignments: A: 1 (the first assignment)
        ##              B: 2 (the other one.)
        ##                   if the other shift is None, we
        ##                   forget atom B.

        spin_systems = []
        shifts = []

        for assignment in assignment_list:

            method = assignment.getMethod()

            ## special treatment for floating chirality

            if method == ASSIGNMENT_METHOD_FLOATING:

                ## assumed: exactly 2 spin-systems
                ## ensured by the assignment filter
                
                s1, s2 = assignment.getSpinSystems()

                shift1 = s1.getChemicalShifts()[0]
                shift2 = s2.getChemicalShifts()[1]

                if shift1[0] is not None:
                    spin_systems.append(s1)
                    shifts.append(shift1)

                if shift2[0] is not None:
                    spin_systems.append(s2)
                    shifts.append(shift2)

            else:

                ## other methods (averaging, unique)
                ## are treated normally, since we have an
                ## one-to-one correspondance between shift
                ## and spin-system.

                s = assignment.getSpinSystems()[0]
                shift = s.getChemicalShifts()[0]

                shifts.append(shift)
                spin_systems.append(s)

        ppms = [p[0] for p in shifts]
        errs = [p[1] for p in shifts]

        default_shift_err = self.getSettings()['default_shift_error']

        ## If we don't have a shift-error, we use a default value

        for i in range(len(errs)):
            if errs[i] is None:
                errs[i] = default_shift_err

        ## complile index-list for all protons

        atom_lists = [s.getAtoms() for s in spin_systems]

        proton_indices = [i for i in range(len(atom_lists))
                          if atom_lists[i][0].isProton()] 

        ## get chemical-shifts/errs for all protons
        
        self.proton_shift_ppms = take(ppms, proton_indices)
        self.proton_shift_errs = take(errs, proton_indices)

        ## list of lists of spin-systems which
        ## contain protons only.

        self.proton_spin_systems = [spin_systems[i] for i in proton_indices]
        ## get chemical-shifts of hetero atoms corresponding
        ## to the protons

        values = []
        errors = []

        ## for every proton-list: collect shift/shift-errors of
        ## corresponding hetero atom.
        ## if proton has no hetero atom, we set a dummy value/error
        ## to -999./0.

        proton_lists = [atom_lists[i] for i in proton_indices]

        for proton_list in proton_lists:

            ## we assume here, that all protons are attached
            ## to the same hetero atom.

            hetero = proton_list[0].getHeteroAtom()

            value = -999.
            error = 0.

            if hetero is not None:
                
                for i in range(len(atom_lists)):
                    if hetero in atom_lists[i]:
                        value = ppms[i]
                        error = errs[i]
                        break

            values.append(value)
            errors.append(error)

        self.hetero_shift_ppms = array(values)
        self.hetero_shift_errs = array(errors)

    def _get_proton_spin_systems(self, proton_shift, hetero_shift,
                                 default_proton_err, default_hetero_err):

        proton_ppm, proton_err = proton_shift
        hetero_ppm, hetero_err = hetero_shift

        if proton_err is None:
            proton_err = default_proton_err
            
        if hetero_err is None:
            hetero_err = default_hetero_err

        ## get index-mask for protons whose shift S lies in the
        ## freq. window
        ## proton_shift - proton_shift_err <= S <=
        ## proton_shift + proton_shift_err

        mask = less_equal(abs(self.proton_shift_ppms - proton_ppm),
                          self.proton_shift_errs + proton_err)

        ## same for the hetero atom

        if hetero_ppm is not None:
            mask2 = less_equal(abs(self.hetero_shift_ppms - hetero_ppm),
                               self.hetero_shift_errs + hetero_err)
            
            mask = logical_and(mask, mask2)

        if self.use_restraint_weights and sum(mask):

            delta = compress(mask, self.proton_shift_ppms - proton_ppm)
            err = compress(mask, self.proton_shift_errs + proton_err)
            
            weights = peak_shape(delta, err)

            if hetero_ppm is not None:
                           
                delta = compress(mask, self.hetero_shift_ppms - hetero_ppm)
                err = compress(mask, self.hetero_shift_errs + hetero_err)
            
                weights *= peak_shape(delta, err)
            
            self.restraint_weights = weights

        else:
            self.restraint_weights = [1.]

        indices = flatnonzero(mask)

        ## return proton spin-systems which meet all requirements

        return [self.proton_spin_systems[i] for i in indices]

    ## BARDIAUX 2.2
    def check_structural_rule(self, a1, a2):
        """
        Structural rule is mainly use for dimers.
        If the 2 atoms belong to the same structure
        element and if they're separated by more than x
        residues, the contribution can be inter-monomer only.
        One can extend them to a monomer, ie interaction i-i+6
        within an helix are impossible ?
        """

        if self.getSettings()['structural_rules_enabled'] == NO:
            return 1
            
        res1 =int(a1.getResidue().getNumber())
        seg1 = a1.getSegid()
        str1 = a1.getResidue().getStructure()
        
        
        res2 = int(a2.getResidue().getNumber())
        seg2 = a2.getSegid()
        str2 = a2.getResidue().getStructure()
        
        if str1 == "" or str2 == "":
            return 1
                
        if seg1 <> seg2:
            return 1

        # here seg1 == seg2
        distance = abs(res1 - res2)

        both_H = str1 == str2 and str1[0] == 'H'
        both_B = str1 == str2 and str1[0] == 'B'

        if both_H:
            if distance <= 5:
                return 1
            else:
                return 0

        elif both_B:
            if distance <= 4:
                return 1
            else:
                return 0

        else:
            return 1
            

    def build_spin_system(self, atom_list, shift):
        ss = SpinSystem(AVERAGING_METHOD_NONE)
        ss.setAtoms(tuple(atom_list))
        ss.setChemicalShifts((shift,))

        return ss

    def build_contributions(self, xpk, list1, list2):
        """
        list1/2: list/tuple of SpinSystems
        """

        import aria.Contribution as C
        from aria.Singleton import SpinPairFactory
        
        SpinPair = SpinPairFactory()

        contributions = []

        # BARDAIUX 2.2
        # some contributiosn are doubled when using CCPN
##         spin_systems_pairs = []
        
        ## TODO: here, one could check the averaging
        ##       method of both spin_systems in order
        ##       to get the correct rule for creating
        ##       'contributions'. how a contribution
        ##       is then encoded/converted into a
        ##       restraint, should be up to the
        ##       structure engine!

        
        for i in range(len(list1)):

            spin_system1 = list1[i]

            for j in range(len(list2)):

                spin_system2 = list2[j]

                spin_pairs = []

                ## create spin-pairs for
                ## all possible atom-combinations

                for a1 in spin_system1.getAtoms():
                    for a2 in spin_system2.getAtoms():

                        if a1 == a2:
                            continue

                        ## BARDIAUX 2.2 : diagonal peaks !!!
                        ## what do we do with diagonal inter-paks
                        ## => remove them
                        if a1 == a2.getHomologous():
                            continue

                        if self.check_structural_rule(a1, a2):
                            
                            sp = SpinPair(a1, a2)
                            spin_pairs.append(sp)

                if not len(spin_pairs):
                    continue
                
                c = C.Contribution(self.__class__.contribution_counter,
                                   C.CONTRIBUTION_TYPE_FAST_EXCHANGE,
                                   spin_pairs, (spin_system1, spin_system2))
                

                c.spectral_weight = 1.
                
                if self.use_restraint_weights:

                    if xpk.proton1_weights is not None and \
                           xpk.proton2_weights is not None:

                        c.spectral_weight = xpk.proton1_weights[i] * \
                                            xpk.proton2_weights[j]

                self.__class__.contribution_counter += 1

                ## TODO: hack
                c.setSpinSystems((spin_system1, spin_system2))

                ## multimers
                inter = spin_system1.getAtoms()[0].getSegid() <> spin_system2.getAtoms()[0].getSegid()
                spec_amb = xpk.getSpectrum().getExperimentData()['ambiguity_type']
                peak_amb = xpk.getAmbiguity()

                amb = peak_amb
                if amb is None:
                    amb = spec_amb

                if (amb == 'inter' and inter) or (amb == 'intra' and not inter) or amb == 'all':
                    contributions.append(c)

        return contributions

    def assign(self, xpk):

        import aria.Assignment as Assignment
        
        check_type(xpk, 'CrossPeak')

        ## collect all shifts which are equal to
        ## the crosspeak's proton1 shift

        ## if we shall use manual assignments, get these first.

        use_assignments = self.getSettings()['use_assignments'] == YES

        proton1_spin_systems = None
        proton2_spin_systems = None

        if use_assignments:
            
            proton1_assignments = xpk.getProton1Assignments()
            proton2_assignments = xpk.getProton2Assignments()

            ## build list of proton spin-systems on-the-fly

            if proton1_assignments:
                shift = xpk.getProton1ChemicalShift()
                proton1_spin_systems = [self.build_spin_system(a.getAtoms(),
                                                               shift) \
                                        for a in proton1_assignments] 
                
            if proton2_assignments:
                shift = xpk.getProton2ChemicalShift()
                proton2_spin_systems = [self.build_spin_system(a.getAtoms(),
                                                               shift) \
                                        for a in proton2_assignments]

        if proton1_spin_systems is None:

            proton_shift = xpk.getProton1ChemicalShift()
            hetero_shift = xpk.getHetero1ChemicalShift()

            settings = self.getSettings()
            default_proton_err = settings['proton1_shift_err']
            default_hetero_err = settings['hetero1_shift_err']

            args = [proton_shift, hetero_shift, default_proton_err,
                    default_hetero_err]

            proton1_spin_systems = self._get_proton_spin_systems(*args)

            ## BARDIAUX 2.2
            ## For symmetry : if spec is inter, proton 1 must refer to only one
            ## segid (first one !)
            amb_type_spec = settings['spec_type']
            amb_type_xpk = xpk.getAmbiguity()
            is_inter =  amb_type_xpk <> 'intra' or amb_type_spec in ['inter', 'all']
            
            #if settings['spec_type'] in ['inter', 'all'] and settings['sym_type'] in ["C2","C3","D2"]:
            if is_inter and settings['sym_type'] in ["C2","C3","D2","C5"]:
            
                ref_seg = settings['ref_segid']
                proton1_spin_systems = [p for p in proton1_spin_systems if p.getAtoms()[0].getSegid() == ref_seg]
                                

            self.build_assignments(xpk, proton1_spin_systems,
                                   Assignment.ASSIGNMENT_TYPE_AUTOMATIC, 1)

            if self.use_restraint_weights:
                xpk.proton1_weights = self.restraint_weights
            
        else:
            if self.use_restraint_weights:
                xpk.proton1_weights = None
                
        ## collect all shifts which are equal to
        ## the crosspeak's proton2 shift

        if proton2_spin_systems is None:

            proton_shift = xpk.getProton2ChemicalShift()
            hetero_shift = xpk.getHetero2ChemicalShift()

            settings = self.getSettings()
            default_proton_err = settings['proton2_shift_err']
            default_hetero_err = settings['hetero2_shift_err']

            args = [proton_shift, hetero_shift, default_proton_err,
                    default_hetero_err]

            proton2_spin_systems = self._get_proton_spin_systems(*args)
            
            ## BARDIAUX 2.2
            ## For symmetry : if spec is inter, proton 2 must refer the other segid
            ## segid (first one !)
            amb_type_spec = settings['spec_type']
            amb_type_xpk = xpk.getAmbiguity()
            is_inter =  amb_type_xpk == 'inter' or amb_type_spec == 'inter'

            #if settings['spec_type'] == 'inter' and settings['sym_type'] in ["C2", "C3", "D2"]:
            if is_inter and settings['sym_type'] in ["C2","C3","D2","C5"]:
                ref_seg = settings['ref_segid']
                proton2_spin_systems = [p for p in proton2_spin_systems if p.getAtoms()[0].getSegid() <> ref_seg]
                            
            self.build_assignments(xpk, proton2_spin_systems, 
                                   Assignment.ASSIGNMENT_TYPE_AUTOMATIC, 2)

            if self.use_restraint_weights:
                xpk.proton2_weights = self.restraint_weights

        else:
            if self.use_restraint_weights:
                xpk.proton2_weights = None
            
        ## TODO: use hetero assignments if proton assignments are missing

        return self.build_contributions(xpk, proton1_spin_systems,
                                        proton2_spin_systems)

    def build_assignments(self, xpk, proton_spin_systems, assignment_type,
                          dimension):

        from aria.Assignment import Assignment

        d = {}
        d[1] = (xpk.addProton1Assignment, xpk.getHetero1ChemicalShift,
                xpk.addHetero1Assignment)

        d[2] = (xpk.addProton2Assignment, xpk.getHetero2ChemicalShift,
                xpk.addHetero2Assignment)

        addP, getH, addH = d[dimension]

        for spin_system in proton_spin_systems:
            
            protons = spin_system.getAtoms()
            proton_assignment = Assignment(protons, assignment_type)
            addP(proton_assignment)
            
            if not getH()[0]: continue
            
            heteros = [p.getHeteroAtom() for p in protons]

##             if not heteros.count(heteros[0]) == len(heteros):
##                 e = 'Spin-system refers to different hetero atoms.'
##                 self.error(ValueError, e)
                
##             hetero_assignment = Assignment(tuple(heteros[:1]), assignment_type)
##             addH(hetero_assignment)

            for hetero in heteros:

                hetero_assignment = Assignment((hetero,), assignment_type)
                addH(hetero_assignment)
    
class PeakAssignerTextPickler(AriaBaseClass):
    
    line = '%6s%12s%10s%10s%10s%10s'

    def __init__(self):
        AriaBaseClass.__init__(self, name = 'Assigner.TP')
        self.reset()

    def reset(self):
        self.peak_lists = []
        self.info = None

    def dumps_header(self):
        return self.line % ('## no.', 'spectrum', 'proton1', 'proton2',
                            'hetero1', 'hetero2') + '\n'

    def dumps_info(self):
        if self.info is None:
            return []

        info = []

        for line in self.info:
            info.append('## %s' % line)

        return info
    
    def set_info(self, info):
        check_type(info, LIST, TUPLE)
        self.info = info

    def add_peaks(self, peaks):
        check_type(peaks, LIST, TUPLE)
        check_elements(peaks, 'CrossPeak')
        
        self.peak_lists.append(peaks)

    def dump(self, filename):

        import os
        
        check_string(filename)

        report = self.dumps_info()
        report.append(self.dumps_header())

        for l in self.peak_lists:
            report += self.__encode(l)

        if os.path.exists(filename):
            s = 'PeakAssigner report file %s already exists and ' + \
                'will be overwritten.'
            self.warning(s % filename)

        f = open(filename, 'w')
        f.write('\n'.join(report))
        f.close()
        
    def __encode(self, peaks):

        report = []

        for peak in peaks:

            number = str(peak.getNumber())

            name = str(peak.getSpectrum().getName())

            shift1 = peak.getProton1ChemicalShift()[0]
            if shift1 is None:
                shift1 = str(None)
            else:
                shift1 = '%.2f' % shift1

            shift2 = peak.getProton2ChemicalShift()[0]
            if shift2 is None:
                shift2 = str(None)
            else:
                shift2 = '%.2f' % shift2

            shift3 = peak.getHetero1ChemicalShift()[0]
            if shift3 is None:
                shift3 = str(None)
            else:
                shift3 = '%.2f' % shift3

            shift4 = peak.getHetero2ChemicalShift()[0]
            if shift4 is None:
                shift4 = str(None)
            else:
                shift4 = '%.2f' % shift4

            s = self.line % (number, name, shift1, shift2, shift3, shift4)
            report.append(s)

        return report
        
class PeakAssignerXMLPickler(XMLBasePickler):

    def _xml_state(self, x):
        e = XMLElement()

        return e

    def load_from_element(self, e):
        s = PeakAssignerSettings()
        
        return s

PeakAssignerSettings._xml_state = PeakAssignerXMLPickler()._xml_state
