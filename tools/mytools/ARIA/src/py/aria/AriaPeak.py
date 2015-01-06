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
import aria.CrossPeak as Xpk
from aria.xmlutils import XMLBasePickler

## constants

NA = 'N/A'

ASSIGNMENT_TYPE_DICT = {Xpk.CROSSPEAK_ASSIGNMENT_TYPE_MANUAL: 'M',
                        Xpk.CROSSPEAK_ASSIGNMENT_TYPE_SEMIAUTOMATIC: 'S',
                        Xpk.CROSSPEAK_ASSIGNMENT_TYPE_AUTOMATIC: 'A',
                        None: NA}

HEADER_PROJECT = \
"""
# Author:  %(author)s
# Date:    %(date)s
# Project: %(project)s
# Run:     %(run)s
# Path:    %(working_directory)s
"""[1:-1]

HEADER_ASSIGNMENT_TYPE = \
"""
# a_type:   Assignment type. Valid types are:
#           M (manual):         Restraint stems from existing (complete) cross-peak assignment.
#           S (semi-automatic): Reference cross-peak is only partially assigned. ARIA has
#                               completed the assignment.
#           A (automatic):      Reference cross-peak has not been assigned. ARIA has generated
#                               a list of probable assignments.
"""[1:-1]

HEADER_SEQUENCE_SEPARATION = \
"""
# sep:      sequence separation s: I:      s == 0  (intra-residual)
#                                  Q:      s == 1  (sequential)        
#                                  S: 2 <= s <= 3  (short)
#                                  M: 4 <= s <= 5  (medium)
#                                  L:      s >  5  (long)
#                                  i: inter-monomer
"""[1:-1]

HEADER_RESTRAINT_DEFINITION = \
"""
# ref_spec: Spectrum name
# ref_no:   Reference peak number (as in spectrum)
# id:       Internal ARIA restraint-id (unique over all spectra)
"""[1:-1]

HEADER_RESTRAINT_BOUNDS = \
"""                        
# d:        Restaint distance [A], after calibration.
# upper:    Upper bound [A]
# lower:    Lower bound [A]
"""[1:-1]

HEADER_RESTRAINT_ACTIVE = \
"""
# active:   yes(no): restraint has (not) been used in the calculation. A restraint is inactive
#           either if it was violated in the last iteration or if it has been merged with
#           another restraint.
"""[1:-1]

HEADER_DICT = {'project': HEADER_PROJECT,
               'assignment_type': HEADER_ASSIGNMENT_TYPE,
               'sequence_separation': HEADER_SEQUENCE_SEPARATION,
               'restraint_definition': HEADER_RESTRAINT_DEFINITION,
               'restraint_bounds': HEADER_RESTRAINT_BOUNDS,
               'restraint_active': HEADER_RESTRAINT_ACTIVE}

HEADER_ABBREVIATIONS = \
("""
#
# Abbreviations:
#
%(restraint_definition)s
%(restraint_active)s
%(restraint_bounds)s
#
# Statistical quantities which have been calculated with respect to the structures of
# this iteration:
#
# d_avg:    Average distance [A] found in structure ensemble.
# l_viol:   Average lower-bound violation [A].
# u_viol:   Average upper-bound violation [A].
# %%_viol:   Percentage of structures in which a restraint was violated.
# viol:     'yes' if the restraint has been violated. 'no' otherwise. Violated restraints
#           will not be used in the subsequent iteration.
#
%(assignment_type)s
#
""" % HEADER_DICT)[1:-1]

HEADER_AMBIG = \
"""
#
# List of ambiguous distance restraints.
#
# Created by Aria 2.0a, %(creation_date)s
#
%(project)s
#
# Restraints used during calculation: %(n_active)d
# Violated: %(n_violated)d
#
%(abbreviations)s
# n_c:      The number of contributions. (see noe_restraints.assignments for
#           explicit list of contributions).
#
# Statistical quantities are calculated with respect to a 'filtered' ensemble of structures
# (e.g. n lowest energy structures etc.) See <iteration> section in project-xml file.
#
"""[1:]

HEADER_UNAMBIG = \
"""
#
# List of unambiguous distance restraints.
#
# Created by Aria 2.0a, %(creation_date)s
#
%(project)s
#
# Restraints used during calculation: %(n_active)d
# Violated: %(n_violated)d
#
%(abbreviations)s
%(sequence_separation)s
# 
# All statistical quantities are calculated with respect to a 'filtered' ensemble of structures
# (e.g. n-lowest energy structures etc). See also <iteration> section in project-xml file.
#
"""[1:]

ABBREVIATIONS_CONTRIBUTIONS = \
"""
#
# List of assignments
#
# Created by Aria 2.0a, %(creation_date)s
#
%(project)s
#
# Abbreviations:
#
# Restraint:
#
%(restraint_definition)s
# d:        Average distance found in ensemble [A].
# u:        Upper bound [A]
# u_viol:   Average upper-bound violation [A]
# %%_viol:   Percentage of structures in which restraint was violated.
# viol:     yes/no. According to violation analysis, the restraint has (not) been violated.
%(assignment_type)s
#
# reliable: yes(no): Reference peak is (not) marked as 'reliable'. Reliable crosspeaks always
#           enter structure calculation.
#
# Contributions: residue1 atom1  -  residue2 atom2
#
# d:        Average distance found in the structure ensemble [A] +/- standard deviation.
#           If standard-deviation could not be calculated (e.g. if there is only 1 distance),
#           it is set to 'N/A'.
# weight:   The weight which is assigned to every contribution by the 'Partial Assignment'
#           step. The sum of weights over all contributions in a restraint is always greater
#           or equal to the parameter 'weight_threshold' which is defined in the section
#           <partial_assignment> in the project's xml file.
#

"""[1:]

CONTRIBUTION_SINGLE = \
"""
    %(spin_pair)s    d: %(d_avg)s +/- %(d_sd)s, weight: %(weight)s
"""[1:]

CONTRIBUTION_SPINPAIR = \
"""
    %(spin_pair)s
"""[1:]

CONTRIBUTIONS_HEADER = \
"""
ref_spec: %(spec_name)s, ref_peak: %(peak_no)s, id: %(id)s, d: %(d_avg)s, u: %(d_upper)s, u_viol: %(u_viol)s, %%_viol: %(p_viol)s, viol: %(violated)s, reliable: %(reliable)s, a_type: %(assignment_type)s
"""[1:]

HEADER_VIOLATIONS = \
"""
#
# Restraint list sorted by upper-bound violations.
#
# Created by Aria 2.0a, %(creation_date)s
#
%(project)s
#
# Abbreviations:
#
%(restraint_definition)s
%(restraint_active)s
%(restraint_bounds)s
#
# Statistical quantities which have been calculated with respect to the structures of
# this iteration:
#
# d_avg:    Average distance [A] found in structure ensemble.
# l_viol:   Average lower-bound violation [A].
# u_viol:   Average upper-bound violation [A].
# %%_viol:   Percentage of structures in which a restraint was violated.
# viol:     'yes' if the restraint has been violated. 'no' otherwise. Violated restraints
#           will not be used in the subsequent iteration.
#
%(assignment_type)s
#
# n_c:      For ambiguous restraints, n_c denotes the number of contributions
#           (see noe_restraints.assignments for explicit list of contributions).
#           For unambiguous restraints, n_c denotes the sequence separation, s,
#           of both contributions:
#                                  I:      s == 0  (intra-residual)
#                                  Q:      s == 1  (sequential)        
#                                  S: 2 <= s <= 3  (short)
#                                  M: 4 <= s <= 5  (medium)
#                                  L:      s >  5  (long)
# 
# All statistical quantities are calculated with respect to a 'filtered' ensemble of structures
# (e.g. n-lowest energy structures etc). See also <iteration> section in project-xml file.
#
"""[1:]

class TextPickler(AriaBaseClass):
    
    def convert(self, f, format_string, failed = 'n/a'):
        try:
            return str(format_string % f())
        except:
            return failed

    def indent(self, lines, n_times = 1, n_spaces = 4, first_line = None):

        if first_line is not None:
            n_spaces = len(first_line) + 1
            result = ['%s %s' % (first_line, lines[0])]
            lines = lines[1:]
        else:
            result = []
        
        prefix = ' ' * (n_times * n_spaces)
        
        return result + map(lambda s, p = prefix: p + s, lines)
    
    def format_output(self, table, header = None):

        import numpy

        if header is not None:
            table.insert(0, header)

        table = numpy.array(table, 'O')

        ## for each column: get largest entry

        lens = [max([len(x) for x in column]) for column in numpy.transpose(table)]

        ## set length of each item to maximum value plus some spacing

        formats = map(lambda l, s = str: '%' + s(l) + 's  ', lens)

        result = []

        for row in table:
            rr = [format % value for format, value in zip(formats, row)]
            result.append(''.join(rr))
            
##         ## TODO: check whether this needs too much memory

##         formats *= len(table)

##         result = map(lambda entry, format: format % entry,
##                      numpy.ravel(table), formats)

##         result = numpy.reshape(numpy.array(result, 'O'),
##                                  numpy.shape(table))

##         result = numpy.sum(result, 1)
        
        if header is not None:
            header = [result[0]]
            result = header + ['\n'] + list(result[1:])
        else:
            result = list(result)

        return result
    
class ContributionTextPickler(TextPickler):

    def encode_spinpair(self, sp):

        atom1, atom2 = sp.getAtoms()
                
        try:
            res1 = '%s' % atom1.getResidue().getName()
        except:
            res1 = NA
            
        try:
            res2 = '%s' % atom2.getResidue().getName()
        except:
            res2 = NA

        values = res1, atom1.getName(), atom1.getSegid(), res2, atom2.getName(), atom2.getSegid()

        return '%7s %4s %4s - %7s %4s %4s' % values

    def encode_contribs(self, contributions):

        text = ''

        for c in contributions:

            spin_pairs = c.getSpinPairs()

            d = {}

            ## statistics

            avg, sd = c.getAverageDistance()
            if avg is None:
                d['d_avg'] = NA
            else:
                d['d_avg'] = '%.2f' % avg
            if sd is None:
                d['d_sd'] = NA
            else:
                d['d_sd'] = '%.2f' % sd

            d['weight'] = '%.1f' % c.getWeight()
            d['spin_pair'] = self.encode_spinpair(spin_pairs[0])
            
            line = CONTRIBUTION_SINGLE % d
                
            if len(spin_pairs) > 1:
                for s in spin_pairs[1:]:
                    val = self.encode_spinpair(s)
                    line += CONTRIBUTION_SPINPAIR % {'spin_pair': val}

            line += '\n\n'
            text += line

        text = text.replace('\n\n', '\n')

        return text

    def encode(self, ap):

        d = {}

        ## reference peak

        ref_peak = ap.getReferencePeak()
        d['peak_no'] = ref_peak.getNumber()
        d['spec_name'] = ref_peak.getSpectrum().getName()
        d['id'] = ap.getId()

        ## restraint parameters

        d['d_upper'] = '%.2f' % ap.getUpperBound()

        AVG, SD = ap.analysis.getAverageDistance()
        
        if AVG is None:
            d['d_avg'] = NA
        else:
            d['d_avg'] = '%.2f' % AVG
            
        if SD is None:
            d['d_sd'] = NA
        else:
            d['d_sd'] = '%.2f' % SD

        
        ## statistics

        AVG, SD = ap.analysis.getLowerBoundViolation()

        if AVG is None:
            d['l_viol'] = NA
        else:
            d['l_viol'] = '%.2f' % AVG

        AVG, SD = ap.analysis.getUpperBoundViolation()

        if AVG is None:
            d['u_viol'] = NA
        else:
            d['u_viol'] = '%.2f' % AVG

        dov = ap.analysis.getDegreeOfViolation()

        if dov is None:
            d['p_viol'] = NA
        else:
            d['p_viol'] = '%.1f' % (dov * 100)

        violated = ap.analysis.isViolated()
        if violated:
            d['violated'] = YES
        elif violated == 0:
            d['violated'] = NO
        else:
            d['violated'] = NA

        ## misc.

        at = ref_peak.getAssignmentType()
        d['assignment_type'] = ASSIGNMENT_TYPE_DICT[at]

        reliable = ref_peak.isReliable()

        if reliable:
            d['reliable'] = YES
        else:
            d['reliable'] = NO
        
        ## contributions
        active = ap.getActiveContributions()

        contribs = self.encode_contribs(active)

        info = CONTRIBUTIONS_HEADER % d
        info += contribs

        ## clean-up

#        info = info.replace('\n\n', '\n')

        return info

    def dumps(self, ap):
        return self.encode(ap)

class AriaPeakTextPickler(TextPickler):

    def encode_common(self, ap):
        
        distance_format = '%.2f'

        number = '%d' % ap.getId()
        
        rp = ap.getReferencePeak()

        x = rp.getNumber()

        try:
            ref_peak_number = '%d' % x
        except:
            ref_peak_number = NA

        x = rp.getSpectrum().getName()

        try:
            ref_peak_spectrum = str(x)
        except:
            ref_peak_spectrum = NA

        x = ap.getDistance()

        try:
            distance = distance_format % x
        except:
            distance = NA

        x = ap.getLowerBound()

        try:
            lower = distance_format % x
        except:
            lower = NA

        x = ap.getUpperBound()

        try:
            upper = distance_format % x
        except:
            upper = NA

        x = ap.isActive()

        if x:
            active = YES
        else:
            active = NO

        at = rp.getAssignmentType()
        assignment_type = ASSIGNMENT_TYPE_DICT[at]
        
        ## analysis section

        ana = ap.analysis

        x = ana.getAverageDistance()[0]

        try:
            ana_avg_d = distance_format % x
        except:
            ana_avg_d = NA

        x = ana.getDegreeOfViolation() * 100

        try:
            ana_degree_of_viol = '%.1f' % x
        except:
            ana_degree_of_viol = NA

        x = ana.getLowerBoundViolation()[0]

        try:
            ana_l_viol = distance_format % x
        except:
            ana_l_viol = NA

        x = ana.getUpperBoundViolation()[0]

        try:
            ana_u_viol = distance_format % x
        except:
            ana_u_viol = NA

        x = ana.isViolated()
        if x:
            ana_violated = YES
        else:
            ana_violated = NO

        values = ref_peak_spectrum, ref_peak_number, number, \
                 active, distance, lower, upper, ana_avg_d, \
                 ana_l_viol, ana_u_viol, \
                 ana_degree_of_viol, ana_violated, assignment_type

        return list(values)

    def encode(self, ap):

        values = self.encode_common(ap)

        ## contributions

        contributions = ap.getContributions()

        ## take only active contributions

        contributions = ap.getActiveContributions()

        if len(contributions) == 1:

            ## get sequence separation
            ## in case of multuple spin-pairs,
            ## we just take the first one, since all are
            ## involve the same two residues

            atom1, atom2 = contributions[0].getSpinPairs()[0].getAtoms()
            seq_pos1 = atom1.getResidue().getNumber()
            seq_pos2 = atom2.getResidue().getNumber()

            # BARDIAUX 2.2
            s1, s2 = atom1.getSegid(), atom2.getSegid()
            if s1 <> s2:
                descr = 'i'
                values.append(descr)
                return values

            seq_sep = abs(seq_pos1 - seq_pos2)

            ## intra-residue

            if seq_sep == 0:
                descr = 'I'

            ## sequential
                
            elif seq_sep == 1:
                descr = 'Q'

            ## TODO: are these the correct values?

            ## short range

            elif seq_sep <= 3:
                descr = 'S'

            ## medium range

            elif seq_sep <= 5:
                descr = 'M'

            else:
                descr = 'L'

            values.append(descr)

        ## multiple contributions

        else:
            values.append(str(len(contributions)))

        return values

    def dumps(self, ap):
        return '\n'.join(self.encode(ap))

class AriaPeakListTextPicklerSettings(Settings):

    def create(self):

        from aria.Settings import ChoiceEntity
        
        d = {}

        choices = ('spectrum_name', 'distance_violation')
        descr = ('Specifies the criterion according to which' + \
                ' restraint-lists are sorted. Valid values are: %s.') \
                % str(choices)
        
        d['sort_criterion'] = ChoiceEntity(choices, description = descr)
        
        return d

    def create_default_values(self):
        d = {}
        d['sort_criterion'] = 'spectrum_name'

        return d

class AriaPeakListTextPickler(TextPickler):

    HEADER_COMMON = ['ref_spec', 'ref_no', 'id', 'active',
                     'd', 'lower', 'upper', 'd_avg', 'l_viol',
                     'u_viol', '%_viol', 'viol', 'a_type']

    COLUMNS = {'ambig': HEADER_COMMON + ['n_c'],
               'unambig': HEADER_COMMON + ['sep']}

    HEADER = {'ambig': HEADER_AMBIG,
              'unambig': HEADER_UNAMBIG}

    def __init__(self, settings):
        check_type(settings, 'AriaPeakListTextPicklerSettings')
        TextPickler.__init__(self, settings = settings)
    
    def get_column_header(self, _type):
        """
        _type is 'ambig' or 'unambig'
        """

        if not _type in ('ambig', 'unambig'):
            s = 'Header for peak-type "%s" not known.' % _type
            self.error(TypeError, s)

        return list(self.COLUMNS[_type])

    def encode(self, peak_list, header):

        pickler = AriaPeakTextPickler()
        all = map(pickler.encode, peak_list)
    
        ## add header

        if not len(all):
            return header

        if len(header) <> len(all[0]):
            s = 'Number of columns must match header-length.'
            self.error(Exception, s)

        header[0] = '# ' + header[0]
        
        ## show additional information

        active = [p for p in peak_list if p.isActive()]

        n_violated = len([p for p in active if p.analysis.isViolated()])

##         p_active = n_active * 100. / len(peak_list)

        d = self._compile_header_dict()
        
##         d['n_restraints'] = len(peak_list)
        d['n_violated'] = n_violated
##         d['p_violated'] = p_violated
        d['n_active'] = len(active)
##         d['p_active'] = p_active
        d['abbreviations'] = HEADER_ABBREVIATIONS

        text = self.format_output(all, header = header)

        ## add \n
        text = [line + '\n' for line in text]

        ## make string

        text = ''.join(text)

        return text, d
    
    def _write(self, s, filename, gzip = 0):
        import os
        
        if s is None:
            import aria.tools as tools
            
            tools.touch(filename)
            return

        if gzip:
            from aria.tools import gzip_open as open_func
        else:
            open_func = open

        filename = os.path.expanduser(filename)

        f = open_func(filename, 'w')
        f.write(s)
        f.close()

    def _compile_header_dict(self):
        from aria.Singleton import ProjectSingleton
        import time
        from copy import copy
        
        project = ProjectSingleton()
        project_settings = project.getSettings()

        infra = project.getInfrastructure()
        run_path = infra.get_run_path()

        d = {'date': project_settings['date'],
             'project': project_settings['name'],
             'run': project_settings['run'],
             'author': project_settings['author'],
             'working_directory': run_path}

        x = copy(HEADER_DICT)
        x['project'] %= d
        x['creation_date'] =time.ctime()

        return x

    def dump_assignments(self, peak_list, filename, gzip = 0):

        if peak_list:
        
            pickler = ContributionTextPickler()
            lines = [pickler.dumps(p) for p in peak_list]

            d = self._compile_header_dict()
            s = (ABBREVIATIONS_CONTRIBUTIONS % d)
            s += ''.join(lines)
#            s = s.replace('\n\n', '\n')
            
        else:
            s = None

        return self._write(s, filename, gzip)

    def dump_violations(self, peak_list, filename, gzip = 0):
        if peak_list:

            header = self.get_column_header('ambig')
            text, d = self.encode(peak_list, header)
            d.update(self._compile_header_dict())
            
            s = HEADER_VIOLATIONS % d
            s += text
#            s = s.replace('\n\n','\n') + text

        else:
            s = None

        return self._write(s, filename, gzip)

        return s

    def dump_ambiguous(self, peak_list, filename, gzip = 0):
        if peak_list:
            
            header = self.get_column_header('ambig')
            text, d = self.encode(peak_list, header)
            d.update(self._compile_header_dict())
            
            header = (self.HEADER['ambig'] % d)[1:]
            s = header + text
#            s = header.replace('\n\n','\n') + text

        else:
            s = None

        return self._write(s, filename, gzip)

    def dump_unambiguous(self, peak_list, filename, gzip = 0):
        if peak_list:

            header = self.get_column_header('unambig')
            text, d = self.encode(peak_list, header)
            d.update(self._compile_header_dict())
                     
            header = (self.HEADER['unambig'] % d)[1:]
            s = header + text
#            s = header.replace('\n\n','\n') + text

        else:
            s = None
        
        return self._write(s, filename, gzip)

class AriaPeakAnalysis(AriaBaseClass):

    """
    this container stores quantities which are calculated
    from structures generated based on values stored in
    the corresponding AriaPeak.
    """

    def __init__(self):
        AriaBaseClass.__init__(self)
        self.setDefaultParameters()

    def setDefaultParameters(self):

        from aria.Datum import Datum

        self.__average_distance = Datum(None)
        self.__degree_of_violation = None
        self.__lower_bound_violation = Datum(None)
        self.__upper_bound_violation = Datum(None)
        self.__is_violated = None
        self.__calculated_peak_size = Datum(None)
        self.__figure_of_merit = Datum(None)

    def setAverageDistance(self, d):
        check_type(d, 'Datum')
        self.__average_distance = d

    def getAverageDistance(self):
        return self.__average_distance

    def isViolated(self, v = None):
        """
        if argument is None, 1 if returned if
        peak is marked as violated, 0 otherwise.

        if argument is integer, the peak can be marked
        as violated.
        """

        check_type(v, INT, NONE)
        
        if v is None:
            return self.__is_violated

        else:
            self.__is_violated = v

    def setDegreeOfViolation(self, v):
        """
        fraction of structures in which restraint
        was violated.
        """
        check_type(v, FLOAT, NONE)
        self.__degree_of_violation = v

    def getDegreeOfViolation(self):
        return self.__degree_of_violation

    def setLowerBoundViolation(self, v):
        check_type(v, 'Datum')
        self.__lower_bound_violation = v

    def getLowerBoundViolation(self):
        return self.__lower_bound_violation

    def setUpperBoundViolation(self, v):
        check_type(v, 'Datum')
        self.__upper_bound_violation = v

    def getUpperBoundViolation(self):
        return self.__upper_bound_violation

    def setCalculatedPeaksize(self, v):
        check_type(v, 'Datum')
        self.__calculated_peak_size = v

    def getCalculatedPeaksize(self):
        return self.__calculated_peak_size

    def setFigureOfMerit(self, v):
        check_type(v, 'Datum')
        self.__figure_of_merit = v

    def getFigureOfMerit(self):
        return self.__figure_of_merit

    def __str__(self):
        
        s = \
"""AriaPeakAnalysis(average_distance=%s, model_peak_size=%s, lower_bound_violation=%s, upper_bound_violation=%s, degree_of_violation=%s, is_violated=%s, figure_of_merit=%s"""

        return s % (str(self.getAverageDistance()),
                    str(self.getCalculatedPeaksize()),
                    str(self.getLowerBoundViolation()),
                    str(self.getUpperBoundViolation()),
                    str(self.getDegreeOfViolation()),
                    str(self.isViolated()),
                    str(self.getFigureOfMerit()))

    __repr__ = __str__

class AriaPeakAnalysisXMLPickler(XMLBasePickler):

    order = ['degree_of_violation', 'violated', 'figure_of_merit',
             'average_distance', 'lower_bound_violation',
             'upper_bound_violation', 'model_peak_size']

    MSG_DATUM = 'aria_peak_analysis.%s: value/error type expected.'

    def _xml_state(self, x):

        from aria.xmlutils import XMLElement

        e = XMLElement(tag_order = self.order)

        e.average_distance = x.getAverageDistance()
        e.model_peak_size = x.getCalculatedPeaksize()
        e.lower_bound_violation = x.getLowerBoundViolation()
        e.upper_bound_violation = x.getUpperBoundViolation()
        e.degree_of_violation = x.getDegreeOfViolation()
        e.figure_of_merit = x.getFigureOfMerit()
        e.violated = x.isViolated()

        return e

    def load_from_element(self, e):

        from aria.Datum import Datum

        a = AriaPeakAnalysis()

        d = {e.average_distance: a.setAverageDistance,
             e.model_peak_size: a.setCalculatedPeaksize,
             e.lower_bound_violation: a.setLowerBoundViolation,
             e.upper_bound_violation: a.setUpperBoundViolation,
             e.figure_of_merit: a.setFigureOfMerit}

        for val, f in d.items():

            value, error =  val

            if value == '':
                value = None
            if error == '':
                error = None

            f(Datum(value, error))

        val = str(e.degree_of_violation)

        try:
            val = float(val)
        except:
            if val == '':
                val = None

        a.setDegreeOfViolation(val)
        
        val = str(e.violated)

        try:
            val = int(val)
        except:
            if val == '':
                val = None
                
        a.isViolated(val)

        return a
    
## BARDIAUX 2.2
class AbstractPeak(AriaBaseClass):
    """
    top class for holding AriaPeak and DistanceRestraints
    """

    def __init__(self):

        AriaBaseClass.__init__(self)

        self.__id = None
        self.__contributions = {}

        ## To indicate whether list of
        ## contributions had been modified
        
        self.__new_contributions = 0

        self.analysis = AriaPeakAnalysis()

        self.setDefaultParameters()

        self.restraint_weight = 1.
        
    def setDefaultParameters(self):

        self.__weight = 1.
        self.__active = 1
        
        self.__lower_bound = None
        self.__upper_bound = None
        self.__distance = None

        self.isMerged(0)
        self.setEquivalentPeak(None)

        self.analysis.setDefaultParameters()
        
        # BARDIAUX rMat
        self.__theoric_volume = 1.0
        self.__ispa = 1.0
        
        self._network = {}

    # BARDIAUX rMat
    def getTheoricVolume(self):
        return self.__theoric_volume

    def setTheoricVolume(self, v):
        check_float(v)
        self.__theoric_volume = v

    def getIspa(self):
        return self.__ispa

    def setIspa(self, v):
        check_float(v)
        self.__ispa = v
        
    def addContribution(self, value):

        check_type(value, 'Contribution')

        ## check if contribution already exists

        key = value.getId()

        if self.__contributions.has_key(key):
            s = 'Contribution "%d" does already exist.'
            raise KeyError, s % key

        self.__contributions[key] = value

        self.__new_contributions = 1

    def setContributions(self, c):

        check_type(c, LIST, TUPLE)
        check_elements(c, 'Contribution')

        self.__contributions = {}
        
        for cc in c:
            self.addContribution(cc)

        self.__new_contributions = 1
        
    def getContributions(self):

        ## Re-sort contributions with respect
        ## to their indices if contributions
        ## have changed

        contribs = self.__contributions.values()

        if self.__new_contributions:
            f = lambda a, b, c = cmp: c(a.getId(), b.getId())
            contribs.sort(f)
            self.__new_contributions = 0
        
        return contribs

    def getActiveContributions(self):
        active = [c for c in self.getContributions() if c.getWeight() > 0.]
        return active

    def isAmbiguous(self):
        """
        returns 1 if the number of active contributions 
        is greater than 1
        """

        active = self.getActiveContributions()
        return len(active) > 1

    def getId(self):
        return self.__id

    def setId(self, id):
        self.__id = id

    def setWeight(self, weight):
        check_type(weight, FLOAT, NONE)
        self.__weight = weight

    def getWeight(self):
        return self.__weight

    def isActive(self, a = None):
        if a is None:
            return self.__active
        else:
            check_int(a)
            self.__active = a

    ## distance, upper- and lower bound of an
    ## AriaPeak will be used to calculate structures

    def setDistance(self, d):
        check_float(d)
        self.__distance = d

    def getDistance(self):
        return self.__distance

    def setLowerBound(self, b):
        check_float(b)
        self.__lower_bound = b

    def getLowerBound(self):
        return self.__lower_bound

    def setUpperBound(self, b):
        check_float(b)
        self.__upper_bound = b

    def getUpperBound(self):
        return self.__upper_bound

    def isMerged(self, flag = None):

        if flag is None:
            return self.__merged
        else:
            check_int(flag)
            self.__merged = int(flag > 0)

    def setEquivalentPeak(self, peak):
        """
        Set the equivalent peak that will be used after merging.
        """
        check_type(peak, 'AriaPeak', NONE)

        self.__equivalent_peak = peak

    def getEquivalentPeak(self):
        return self.__equivalent_peak

    def __getitem__(self, key):

        if not self.__contributions.has_key(key):
            s = 'Contribution with the given index ("%d") does not exist.'
            raise KeyError, s % key
        
        return self.__contributions[key]

    def __len__(self):
        return len(self.__contributions)

# BARDIAUX 2.2
class AriaPeak(AbstractPeak):
    """
    An AriaPeak is organized as dictionary. It is given
    an id which should match the id of the corresponding CrossPeak.
    AriaPeak's elements are Contribution objects, indexed via its id.
    """

    ## TODO: call Id -> Number? or is it another thing?
    
    def __init__(self, id, reference_peak):

        ## TODO: reference_peak must not be None !

        AbstractPeak.__init__(self)

        check_int(id)

        self.setId(id)
        self._setReferencePeak(reference_peak)
        
    def getReferencePeak(self):
        return self.__reference_peak

    def _setReferencePeak(self, p):

        check_type(p, 'CrossPeak', NONE)
        self.__reference_peak = p
        

    def __copy__(self):
        """
        contributions are copied as well.
        """

        from copy import copy

        ## copy has *same* id!
        ## this is not to be confused with python's internal id!

        ## TODO: get rid of __dict__ stuff

        new = AriaPeak(self.getId(), self.getReferencePeak())
        new.__dict__ = copy(self.__dict__)

        ## set new analysis instance
        new.analysis = AriaPeakAnalysis()

        ## TODO: maybe there is a better solution than the one below

        ## copy contributions

        #attr_name = '_%s__contributions' % self.__class__.__name__
        #setattr(new, attr_name, {})
        attr_name = '_%s__contributions' % self.__class__.__bases__[0].__name__
        setattr(new, attr_name, {})

        [new.addContribution(copy(c)) for c in self.getContributions()]

        return new

    def __str__(self):
        
        s = 'AriaPeak(id=%d, weight=%s, distance=%s, ' + \
            'lower_bound=%s, upper_bound=%s, active=%s, merged=%s, ' + \
            'analysis=%s, contributions=%s, ref_peak=%s)'

        analysis_str = str(self.analysis)
        ref_peak = self.getReferencePeak()

        return s % (self.getId(),
                    str(self.getWeight()),
                    str(self.getDistance()),
                    str(self.getLowerBound()),
                    str(self.getUpperBound()),
                    str(self.isActive()),
                    str(self.isMerged()),
                    analysis_str,
                    str(self.getContributions()),
                    str(ref_peak))

    __repr__ = __str__

# BARDIAUX 2.2
class DistanceRestraint(AbstractPeak):
    """
    An AriaPeak is organized as dictionary. It is given
    an id which should match the id of the corresponding CrossPeak.
    AriaPeak's elements are Contribution objects, indexed via its id.
    """

    ## TODO: call Id -> Number? or is it another thing?

    counter = 0
    
    def __init__(self, id = None):

        AbstractPeak.__init__(self)

        #check_int(id)
        
        if id is None:
            self.setId(self.__class__.counter)
            self.__class__.counter += 1
        else:
            self.setId(id)

        self.__source = ''
        self.__number = None
        self.__combination = 0

        self.__reference_peak = None
    
    def getReferencePeak(self):
        return self.__reference_peak

    def setReferencePeak(self, p):

        check_type(p, 'CrossPeak', NONE)
        self.__reference_peak = p
        
    def setSource(self, s):
        check_type(s, STRING)
        self.__source = s

    def getSource(self):
        return self.__source

    def setNumber(self, n):
        check_type(n, INT)
        self.__number = n
        
    def getNumber(self):
        return self.__number

##     def getReferencePeak(self):
##         return self

##     def getSpectrum(self):
##         return self

##     def getName(self):
##         return self.getSource()

    def isCombination(self, c = None):
        if c is None:
            return self.__combination
        else:
            check_int(c)
            self.__combination = c       
        
    
    def __copy__(self):
        """
        contributions are copied as well.
        """

        from copy import copy

        ## copy has *same* id!
        ## this is not to be confused with python's internal id!

        ## TODO: get rid of __dict__ stuff

        new = DistanceRestraint(id = self.getId())
        new.__dict__ = copy(self.__dict__)

        ## set new analysis instance
        new.analysis = AriaPeakAnalysis()

        ## TODO: maybe there is a better solution than the one below

        ## copy contributions

        #attr_name = '_%s__contributions' % self.__class__.__name__
        attr_name = '_%s__contributions' % self.__class__.__bases__[0].__name__
        setattr(new, attr_name, {})

        [new.addContribution(copy(c)) for c in self.getContributions()]

        return new


    def __str__(self):
        
        s = 'DistanceRestraint(id=%d, weight=%s, distance=%s, ' + \
            'lower_bound=%s, upper_bound=%s, active=%s, merged=%s, ' + \
            'analysis=%s, contributions=%s)'

        analysis_str = str(self.analysis)

        return s % (self.getId(),
                    str(self.getWeight()),
                    str(self.getDistance()),
                    str(self.getLowerBound()),
                    str(self.getUpperBound()),
                    str(self.isActive()),
                    str(self.isMerged()),
                    analysis_str,
                    str(self.getContributions()))

    __repr__ = __str__

class AriaPeakXMLPickler:

    order = ['id', 'weight', 'distance', 'lower_bound',
             'upper_bound', 'active', 'merged', 'reference_peak',
             'contribution', 'analysis']

    def __init__(self):

        from aria.Datum import DatumXMLPickler
        from aria.Contribution import ContributionXMLPickler
        from aria.ShiftAssignment import SpinSystemXMLPickler
        from aria.Atom import AtomXMLPickler
        from aria.Datum import ChemicalShiftXMLPickler

        DP = DatumXMLPickler()

        self.analysis = AriaPeakAnalysisXMLPickler()
        self.peak = self
        self.contribution = ContributionXMLPickler()
##        self.spin_pair = SpinPairXMLPickler()
        self.spin_system = SpinSystemXMLPickler()
        self.atom = AtomXMLPickler()
        self.chemical_shift = ChemicalShiftXMLPickler()

        ## analysis section

        self.figure_of_merit = DP
        self.average_distance = DP
        self.model_peak_size = DP
        self.lower_bound_violation = DP
        self.upper_bound_violation = DP
        
    def _xml_state(self, x):
        from aria.xmlutils import XMLElement
        
        e = XMLElement(tag_order = self.order)
        
        e.id = x.getId()
        e.weight = x.getWeight()
        e.distance = x.getDistance()
        e.lower_bound = x.getLowerBound()
        e.upper_bound = x.getUpperBound()
        e.active = x.isActive()
        e.merged = x.isMerged()

        ref_peak = x.getReferencePeak()

        spec = ref_peak.getSpectrum()
        if spec is not None:
            spec_name = spec.getName()
        else:
            spec_name = ''
        
        e.reference_peak = XMLElement()
        e.reference_peak.spectrum = spec_name
        e.reference_peak.number = ref_peak.getNumber()
        
        e.contribution = tuple(x.getContributions())
        e.analysis = x.analysis

        return e

    def load_from_element(self, e):

        import aria.CrossPeak as CrossPeak
        import aria.NOESYSpectrum as NOESYSpectrum
        import aria.Datum as Datum
        from aria.tools import as_tuple

        ap = AriaPeak(int(e.id), None)
        
        ap.setWeight(float(e.weight))
        ap.setDistance(float(e.distance))
        ap.setLowerBound(float(e.lower_bound))
        ap.setUpperBound(float(e.upper_bound))
        ap.isMerged(int(e.merged))
        ap.isActive(int(e.active))

        ap.setContributions(as_tuple(e.contribution))
        ap.analysis = e.analysis

        ## reference peak

        ref = e.reference_peak

        ## spectrum

        ## TODO: hacked
        ## In order to keep track of all newly created spectra,
        ## we keep a list of them.

        try:
            self.__class__.spectra
        except:
            self.__class__.spectra = {}

        spec_name = str(ref.spectrum)

        if spec_name in self.__class__.spectra:
            spec = self.__class__.spectra[spec_name]
        else:
            spec = NOESYSpectrum.NOESYSpectrum(str(ref.spectrum), ())
            self.__class__.spectra[spec_name] = spec

        rp = [p for p in spec.getPeaks() if p.getNumber() == int(ref.number)]

        if len(rp) == 1:
            rp = rp[0]            
        elif len(rp) == 0:
            rp = CrossPeak.CrossPeak(int(ref.number),
                                     Datum.Datum(None, None), None)
            ## add cross-peak to spectrum.
            spec.addPeak(rp)
        else:
            self.error('multiple reference peaks found')
            
        ## set the cross-peaks's reference to the spectrum.

        ap._setReferencePeak(rp)

        return ap

AriaPeak._xml_state = AriaPeakXMLPickler()._xml_state
AriaPeakAnalysis._xml_state = AriaPeakAnalysisXMLPickler()._xml_state

