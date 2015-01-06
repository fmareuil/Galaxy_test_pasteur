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
from aria.AriaPeak import TextPickler
from numpy import *
from aria.mathutils import _average as average
from aria.mathutils import standardDeviation

NA = 'N/A'

REPORT = \
"""
## Report for iteration %(iteration_number)s
## 
## Created by ARIA %(version_string)s, %(creation_date)s
##
## Author:  %(author)s
## Date:    %(date)s
## Project: %(project)s
## Run:     %(run)s
## Path:    %(working_directory)s
##
## Restraints:
##
## Used during calculation: %(n)d
## Violated:                %(n_violated)d (%(p_violated)s)
## Merged:                  %(n_merged)d
##
## -----------------------------------------------------------------
##
%(spectrum)s
##
%(constraint_list)s
"""

UNAMBIG_REPORT = \
"""
## Unambiguous restraints: %(n)s
## Violated:               %(n_viol)s (%(p_viol)s%%)
"""

AMBIG_REPORT = \
"""
## Ambiguous restraints:   %(n)s
## Violated:               %(n_viol)s (%(p_viol)s%%)
## Contribution histogram: 2:    %(c_hist_2)s
##                         3:    %(c_hist_3)s
##                         4:    %(c_hist_4)s
##                         5:    %(c_hist_5)s
##                         6-10: %(c_hist_6-10)s
##                         >10:  %(c_hist_>10)s
"""

SPECTRUM_REPORT = \
"""
## Spectrum:                   %(name)s
## Calibration factor:         %(calibration_factor)s
##
## Restraints:
##
## Used during calculation: %(n)s
## Violated:                %(n_viol)s (%(p_viol)s%%)
##
%(unambig)s
##
%(ambig)s
## -----------------------------------------------------------------
##
"""

CONSTRAINT_LIST_REPORT = \
"""
## Constraint List:            %(name)s
## Calibration factor:         %(calibration_factor)s
##
## Restraints:
##
## Used during calculation: %(n)s
## Violated:                %(n_viol)s (%(p_viol)s%%)
##
%(unambig)s
##
%(ambig)s
## -----------------------------------------------------------------
##
"""

class IterationTextPickler(TextPickler):

    def encode_histogram(self, n_contribs):

        histogram = {}

        for i in range(2,6):
            key = 'c_hist_%d' % i
            histogram[key] = sum(equal(n_contribs, i))

        ## 6 .. 10

        m1 = greater(n_contribs, 5)
        m2 = less_equal(n_contribs, 10)

        histogram['c_hist_6-10'] = sum(logical_and(m1, m2))

        ## > 10

        histogram['c_hist_>10'] = sum(greater(n_contribs, 10))

        ## convert values into strings

        for k, v in histogram.items():
            histogram[k] = '%d' % v
        
        return histogram

    def encode_info(self, l):

        d = {}

        ## default values:

        if not l:
            d['n'] = '0'
            d['n_viol'] = '0'
            d['p_viol'] = '0.'
            d['n_active'] = '0'
            d['p_active'] = '0.'
##             d['avg_d_viol'] = NA
##             d['sd_d_viol'] = NA

            return d

        d['n'] = len(l)
        
        ## collect violated restraints
        violated = [r for r in l if r.analysis.isViolated()]
        active = [r for r in l if r.isActive()]
        
        d['n_viol'] = '%d' % len(violated)
        d['p_viol'] = '%.1f' % (len(violated) * 100. / len(l))

        return d

    def encode_unambig(self, restraint_list):

        unambig = [r for r in restraint_list if not r.isAmbiguous()]

        d = self.encode_info(unambig)
            
        return UNAMBIG_REPORT % d
    
    def encode_ambig(self, restraint_list):
        
        ambig = [r for r in restraint_list if r.isAmbiguous()]

        d = self.encode_info(ambig)

        if not ambig:
            keys = ('c_hist_2', 'c_hist_3', 'c_hist_4',
                    'c_hist_5', 'c_hist_6-10', 'c_hist_>10')

            for k in keys:
                d[k] = NA
            
        else:

            contributions = [r.getContributions() for r in restraint_list]
            
            n_active_contribs = []

            for contrib in contributions:
                active = [c for c in contrib if c.getWeight() > 0.]
                n_active_contribs.append(len(active))
            
            d.update(self.encode_histogram(n_active_contribs))

        report = AMBIG_REPORT % d

        return report
    
    def dumps(self, iteration):

        import time
        from aria.Singleton import ProjectSingleton
        from aria.ariabase import AriaBaseClass

        spectrum_reports = ''

        total_active = 0
        total_violated = 0
        total_merged = 0

        for spectrum, restraints in iteration.getPeaks().items():

            ## Forget all restraints which are inactive due
            ## to merging

            active = [r for r in restraints if r.isActive()]
            merged = [r for r in restraints if r.isMerged()]

            ## spectrum report

            d = {}
            d['name'] = spectrum.getName()
            
            try:
                factor = iteration.getCalibrationFactor(spectrum)
                d['calibration_factor'] = '%e' % factor
            except:
                d['calibration_factor'] = NA

            d['n'] = '%d' % len(active)

            violated = [r for r in active if r.analysis.isViolated()]

            total_active += len(active)
            total_violated += len(violated)
            total_merged += len(merged)
            
            d['n_viol'] = '%d' % len(violated)
            d['p_viol'] = '%.1f' % (len(violated) * 100. / len(active))
            d['unambig'] = self.encode_unambig(active)
            d['ambig'] = self.encode_ambig(active)

            spectrum_reports += SPECTRUM_REPORT % d

        constraint_list_reports = ''

        for spectrum, restraints in iteration.getDistanceRestraints().items():

            ## Forget all restraints which are inactive due
            ## to merging

            active = [r for r in restraints if r.isActive()]
            merged = [r for r in restraints if r.isMerged()]

            ## spectrum report

            d = {}
            d['name'] = spectrum.getName()
            
            try:
                factor = iteration.getCalibrationFactor(spectrum)
                d['calibration_factor'] = '%e' % factor
            except:
                d['calibration_factor'] = NA

            d['n'] = '%d' % len(active)

            violated = [r for r in active if r.analysis.isViolated()]

            total_active += len(active)
            total_violated += len(violated)
            total_merged += len(merged)
            
            d['n_viol'] = '%d' % len(violated)
            d['p_viol'] = '%.1f' % (len(violated) * 100. / len(active))
            d['unambig'] = self.encode_unambig(active)
            d['ambig'] = self.encode_ambig(active)

            constraint_list_reports += CONSTRAINT_LIST_REPORT % d

        ## generic info

        project = ProjectSingleton()
        project_settings = project.getSettings()

        infra = project.getInfrastructure()
        run_path = infra.get_run_path()

        d = {}
        
        d['iteration_number'] = '%d' % iteration.getNumber()
        d['creation_date'] = time.ctime()
        d['date'] = project_settings['date']
        d['project'] = project_settings['name']
        d['run'] = project_settings['run']
        d['author'] = project_settings['author']
        d['working_directory'] = run_path
        d['n'] = total_active
        d['n_violated'] = total_violated
        
        if total_active:
            d['p_violated'] = '%.1f' % (total_violated * 100. / total_active)

        else:
            d['p_violated'] = '-'
            
        d['n_merged'] = total_merged
        d['spectrum'] = spectrum_reports
        d['version_string'] = AriaBaseClass().get_version_string()

        d['constraint_list'] = constraint_list_reports

        report = REPORT % d

        ## remove blank lines

        report = report[1:].replace('\n\n', '\n')
        report = report.replace('\n\n', '\n')

        return report

    def dump(self, iteration, filename):
        import os

        filename = os.path.expanduser(filename)
        coded = self.dumps(iteration)

        f = open(filename, 'w')
        f.write(coded)
        f.close()

class Iteration(AriaBaseClass):

    ## TODO: shall Iteration behave like a
    ## dict? id's as keys, AriaPeaks as values?

    def __init__(self, number):

        check_int(number)

        AriaBaseClass.__init__(self)
        self.__number = number

        ## To indicate whether new peaks
        ## had been added.
        
        self.__new_peaks = 0

        self.setDefaultValues()

    def setDefaultValues(self):
        self.__peaks = {}
##         self.__peak_ids = {}
        self.__structure_ensemble = None
        self.__calibration_factors = {}

        ## BARDIAUX 2.2
        self.__distance_restraints = {}
        self.__combined_restraints = []
        
    def setStructureEnsemble(self, s):
        check_type(s, 'StructureEnsemble')
        self.__structure_ensemble = s

    def getStructureEnsemble(self):
        return self.__structure_ensemble

    def setCalibrationFactor(self, spectrum, c):
        check_type(spectrum, 'NOESYSpectrum')
        check_float(c)
        
        self.__calibration_factors[spectrum] = c

    def getCalibrationFactor(self, spectrum):
        check_type(spectrum, 'NOESYSpectrum')
        return self.__calibration_factors[spectrum]

    def getNumber(self):
        return self.__number

    def _setNumber(self, n):
        """
        private
        """
        check_int(n)
        
        self.__number = n

    def addPeak(self, a):

        check_type(a, 'AriaPeak')

        ref_peak = a.getReferencePeak()
        
        if ref_peak is None:
            self.error(ValueError, 'Peak has not reference-peak.')

        spectrum = ref_peak.getSpectrum()
        
        if spectrum is None:
            s = 'Reference peak is not known to any spectrum'
            self.error(ValueError, s)

        if not self.__peaks.has_key(spectrum):
            self.__peaks[spectrum] = []

        self.__peaks[spectrum].append(a)

        self.__new_peaks = 1

    def getPeaks(self, spectrum = None):
        """
        if spectrum is not None, its peaks are returned.
        otherwise, a dict with spectrum as key and the
        peak-lists as values is returned.
        """

        ## If new peaks have been added recently,
        ## re-sort peak-lists

        if self.__new_peaks:
            f = lambda a, b, c = cmp: c(a.getId(), b.getId())
            [l.sort(f) for l in self.__peaks.values()]
            self.__new_peaks = 0

        if spectrum is None:
            return self.__peaks
        else:
            return self.__peaks[spectrum]

    def getPeakList(self):
        """
        returns a list of all peaks (i.e. spectrum-wise
        peaks-lists are concatenated.
        """

        return reduce(lambda x,y: x+y, self.getPeaks().values(), [])

    ## BARDIAUX 2.2
    ## DistanceRestraints
    def addDistanceRestraint(self, a):

        check_type(a, 'DistanceRestraint')
        
        ref_peak = a.getReferencePeak()
        
        if ref_peak is None:
            self.error(ValueError, 'DistanceRestraint has not reference-peak.')

        constraint_list = ref_peak.getSpectrum()
        
        if constraint_list is None:
            s = 'Reference peak is not known to any constraint list'
            self.error(ValueError, s)
            
        
        if not self.__distance_restraints.has_key(constraint_list):
            self.__distance_restraints[constraint_list] = []

        self.__distance_restraints[constraint_list].append(a)

        #self.__new_peaks = 1

    def getDistanceRestraints(self, listname = None):

        ##return self.__distance_restraints
    
        if listname is None:
            return self.__distance_restraints
        else:
            return self.__distance_restraints[listname]

    def getDistanceRestraintsList(self):
        
        return reduce(lambda x,y: x+y, self.getDistanceRestraints().values(), [])
        
    def setDistanceRestraints(self, restraints):
        check_type(restraints, LIST)
        
        self.__distance_restraints = {}        
        self.__distance_restraints = restraints

    ## BARDIAUX 2.2
    ## CombinedRestraints from combination
        
    def setCombinedRestraints(self, restraints):
        check_type(restraints, LIST)
        
        self.__combined_restraints = []
        self.__combined_restraints = restraints        
    
    def getCombinedRestraints(self):
        return self.__combined_restraints
            


    def __copy__(self):
        """
        AriaPeaks are copied as well.
        """

        from copy import copy

        new = Iteration(self.getNumber())

        peaks = self.getPeakList()

        map(lambda peak, copy = copy, f = new.addPeak: \
            f(copy(peak)), peaks)

        ## BARDIAUX  2.2      
        restraints = self.getDistanceRestraintsList()

        map(lambda restraint, copy = copy, f = new.addDistanceRestraint: \
            f(copy(restraint)), restraints)

        return new

    def __str__(self):

        class_name = self.__class__.__name__

        n_peaks = len(self.getPeakList())
        n_constraints = len(self.getDistanceRestraintsListList())

        return '%s(number=%s, n_peaks=%s, n_constraints=%s)' % (class_name,
                                              str(self.getNumber()),
                                              str(n_peaks),
                                              str(n_constraints))
    __repr__ = __str__
    
