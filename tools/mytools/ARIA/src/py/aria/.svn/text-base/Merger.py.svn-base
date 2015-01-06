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
from aria.TypeChecking import check_type, check_elements, LIST, TUPLE
from aria.Settings import Settings
from aria.xmlutils import XMLBasePickler as _XMLBasePickler

REPORT_HEADER = \
"""
# List of merged peaks.
#
# Author:  %(author)s
# Date:    %(date)s
# Project: %(project)s
# Run:     %(run)s
# Path:    %(working_directory)s
#
# Total no. of merged peaks:    %(n)d
#
# Abbreviations:
#
# m_spec:  Name of spectrum.
# m_no:    Number of peak that has been merged.
# m_id:    Internal peak id (unique over all spectra).
#
# If two peaks M:=(m_spec, m_no) and P:=(spec, no), respectively, have been merged,
# peak M will not appear in the restraint list anymore. Peak P is used instead.
#
# spec:    Name of spectrum.
# no:      Number of peak that remains in restraint-list
# id:      Internal peak id (unique over all spectra).
#
"""

from aria.AriaPeak import TextPickler

class MergerTextPickler(TextPickler):

    HEADER = ['m_spec', 'm_no', 'm_id', 'spec', 'no', 'id']

    def create_header(self):
        from aria.Singleton import ProjectSingleton
        import time

        project = ProjectSingleton()
        s = project.getSettings()

        d = {}
        d['author'] = s['author']
        d['date'] = time.ctime()
        d['project'] = s['name']
        d['run'] = s['run']
        d['working_directory'] = s['working_directory']

        return d

    def dumps(self, restraints):
        
        l = []

        for r in restraints:
            
            m_id = r.getId()
            m_ref_peak = r.getReferencePeak()
            m_spec = m_ref_peak.getSpectrum().getName()
            m_no = m_ref_peak.getNumber()

            p = r.getEquivalentPeak()
            id = p.getId()
            ref_peak = p.getReferencePeak()
            spec = ref_peak.getSpectrum().getName()
            no = ref_peak.getNumber()
            
            x = [m_spec, m_no, m_id, spec, no, id]
            x = map(str, x)
            l.append(x)

        header = list(self.HEADER)
        header[0] = '# %s' % header[0]
        
        lines = self.format_output(l, header)
        lines = '\n'.join(lines)

        ## Insert global header
        d = self.create_header()
        d['n'] = len(restraints)
        header = REPORT_HEADER % d

        header += lines

        header = header[1:].replace('\n\n', '\n')

        return header

    def dump(self, restraints, filename, gzip = 0):
        coded = self.dumps(restraints)

        if gzip:
            from aria.tools import gzip_open as open_func
        else:
            open_func = open
        
        f = open_func(filename, 'w')
        f.write(coded)
        f.close()

class MergerSettings(Settings):

    def create(self):

        from aria.Settings import ChoiceEntity

        kw = {'method': ChoiceEntity(('standard', 'no_merging', 'combination'))}

        return kw

    def create_default_values(self):

        return {'method': 'standard'}

class Merger(AriaBaseClass):

    def __init__(self, settings):

        check_type(settings, 'MergerSettings')

        AriaBaseClass.__init__(self)

        self.setSettings(settings)

    def __call__(self, peaks):
        """
        Merges a set of aria-peak-list and a returns a list of
        nonredundant AriaPeaks.
        """

        check_type(peaks, LIST, TUPLE)
        check_elements(peaks, 'AriaPeak')

        from numpy import argmin

        import time

        t = time.clock()

        settings = self.getSettings()
            
        if not len(peaks):
            self.warning(ValueError, 'No data given.')
            return 0, []

        if settings['method'] == 'standard' or settings['method'] == 'combination':

            equiv_peaks = _find_equivalent_peaks(peaks)

            n = 0

            for pp in equiv_peaks:

                index = 0
                if len(pp) > 1:
                    values = map(lambda p: p.getDistance(), pp)
                    index = argmin(values)

                active_peak = pp[index]
                active_peak.isMerged(0)

                for peak in pp:

                    if peak == active_peak: continue

                    peak.isMerged(1)
                    peak.isActive(0)
                    peak.setEquivalentPeak(active_peak)
                    n += 1

            self.debug('Time: %ss' % str(time.clock() - t))

            #return n, []

            if settings['method'] == 'combination':

        
                long_range = _find_long_range_peaks(peaks)
                new_long = []

		from numpy.random import permutation

                ids = permutation(len(long_range))
           
                if len(ids) % 2 <> 0:
                    ids = ids[1:]

                for k in range(0, len(ids), 2):
                
                    i = ids[k]
                    j = ids[k+1]

                    to_combine = [long_range[x] for x in ids[k:k+2]]
                    p = self._combine(to_combine)

                    new_long.append(p)

                return n, new_long
        
            else:
                return n, []


        return -1

    def _combine(self, peaks):
        
        import numpy as N
        from aria.AriaPeak import DistanceRestraint

        d = DistanceRestraint()

        values = map(lambda p: p.getDistance(), peaks)
        
        t = N.power(N.sum(N.power(values, -6.)), -1./6.)
        t = max(values)
        
        u = t + 0.125 * t**2
        l = t - 0.125 * t**2

        d.setDistance(t)
        d.setLowerBound(l)
        d.setUpperBound(u)
        d.setWeight(1.)
        d.setNumber(peaks[0].getId())        
        d.setSource('combination')
        d.isCombination(1)

        for p in peaks:
            for c in p.getActiveContributions():
                d.addContribution(c)

        for p in peaks:
            p.isMerged(1)
            p.isActive(0)

        peaks[0].setEquivalentPeak(peaks[1])
        peaks[1].setEquivalentPeak(peaks[0])

        return d


def _find_equivalent_peaks(peaks):

    import numpy as num

    arrays = {}
    ids = {}
    peak_dict = {}
    
    for peak in peaks:

        sp_ids = [id(c.getSpinPairs()) for c in peak.getActiveContributions()]

        n = len(sp_ids)

        if not arrays.has_key(n):
            arrays[n] = []
            ids[n] = []

        arrays[n].append(sp_ids)
        ids[n].append(peak.getId())

        peak_dict[peak.getId()] = peak

    ## group peaks according to contribution content

    equivalent_peaks = []

    for n in arrays.keys():

        x = num.sort(num.array(arrays[n]), 1)
        i = num.array(ids[n])

        while len(x):

            mask = num.equal(num.sum(abs(x - x[0]), 1), 0)
            q = [peak_dict.get(ii) for ii in num.compress(mask, i)]
            equivalent_peaks.append(q)

            x = num.compress(num.logical_not(mask), x, 0)
            i = num.compress(num.logical_not(mask), i, 0)

    return equivalent_peaks

## BARDIAUX 2.2
## Constraint combination
def _get_range(c):
    
        atom1, atom2 = c.getSpinPairs()[0].getAtoms()
        seq_pos1 = atom1.getResidue().getNumber()
        seq_pos2 = atom2.getResidue().getNumber()

        return abs(seq_pos1 - seq_pos2)
    
def _find_long_range_peaks(peaks):

    import numpy as num
    
    longs = []
    
    for peak in peaks:

        f = _get_range
        
        contribs = peak.getActiveContributions()
        
        range = num.array([f(c) for c in contribs])
        
        mask_long = num.greater(range, 5)

        if num.sum(mask_long) == len(contribs):
            longs.append(peak)

    return longs

class MergerXMLPickler(_XMLBasePickler):

    def _xml_state(self, x):
        
        from aria.xmlutils import XMLElement

        e = XMLElement()        
        e.method = x['method']        

        return e

    def load_from_element(self, e):

        s = MergerSettings()        
        s['method'] = str(e.method)        

        return s

MergerSettings._xml_state = MergerXMLPickler()._xml_state
