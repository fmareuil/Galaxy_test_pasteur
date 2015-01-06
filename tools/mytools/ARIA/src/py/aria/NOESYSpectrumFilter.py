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
from aria.TypeChecking import *
from aria.Settings import *

## TODO: the whole class is more or less obsolete

class NOESYSpectrumFilterSettings(Settings):

    def create(self):

        from aria.Settings import NonNegativeFloat, PeakType, YesNoChoice

        s = 'Proton-shift error must be a non-negative float.'

        keywords = {'proton1_shift_err': NonNegativeFloat(error_message = s),
                    'proton2_shift_err': NonNegativeFloat(error_message = s),
                    'volume_or_intensity': PeakType(),
                    'filter_diagonal_peaks' : YesNoChoice(),
                    'filter_unassigned_peaks' : YesNoChoice()}

        return keywords

class NOESYSpectrumFilter(AriaBaseClass):

    def __init__(self, settings):

        check_type(settings, 'NOESYSpectrumFilterSettings')
        
        AriaBaseClass.__init__(self)

        self.setSettings(settings)

    def compile_statistic(self):
        """
        Returns a dictionary with counts for the different error
        types.
        """

        import aria.ChemicalShiftFilter as CSF
        import aria.CrossPeakFilter as CPF
        
        dict = {'proton1': 0, 'proton2': 0, 'size': 0, 'diagonal': 0, 'unassigned' : 0}

        for p, r in self.result.items():

            ## proton 1

            messages = r[p.getProton1ChemicalShift()]
            if CSF.INVALID_SHIFT_VALUE in messages or \
               CSF.NEGATIVE_SHIFT_ERROR in messages:
                dict['proton1'] += 1

            ## proton 2

            messages = r[p.getProton2ChemicalShift()]
            if CSF.INVALID_SHIFT_VALUE in messages or \
               CSF.NEGATIVE_SHIFT_ERROR in messages:
                dict['proton2'] += 1

            ## peak size

            if r['size'] in [CPF.NO_PEAK_SIZE, CPF.NEGATIVE_PEAK]:
                dict['size'] += 1

            ## diagonal peak

            if r['diagonal'] == CPF.DIAGONAL_PEAK:
                dict['diagonal'] += 1
                
            ## unassigned peak

            if r['unassigned'] == CPF.UNASSIGNED_PEAK:
                dict['unassigned'] += 1
                
        return dict
            
    def __call__(self, spectrum):

        check_type(spectrum, 'NOESYSpectrum')

        ## create cross-peak filter
        
        import aria.CrossPeakFilter as X
        
        settings = X.CrossPeakFilterSettings(default_settings = \
                                             self.getSettings())
        
        peak_filter = X.CrossPeakFilter(settings)

        result = {}

        for peak in spectrum.getPeaks():

            valid = peak_filter(peak)
            peak.is_valid(valid)

            if not valid: result.update(peak_filter.result)

        self.result = result

        return [p for p in spectrum.getPeaks() if p.is_valid()]


class NOESYSpectrumFilterTextPickler(AriaBaseClass):

    header = {'spectrum': 'spectrum',
              'peak': '## no.',
              'proton1': 'proton1',
              'error1' : 'error1',
              'proton2': 'proton2',
              'error2': 'error2',
              'size': 'size',
              'diagonal': 'diagonal',
              'unassigned' : 'unassigned'}

    line = '%6s%12s%10s%10s%10s%10s'
    line = '%(peak)6s%(spectrum)12s%(proton1)10s%(error1)10s' + \
           '%(proton2)10s%(error2)10s%(size)10s%(diagonal)10s%(unassigned)12s\n'

    def __init__(self):

        AriaBaseClass.__init__(self, name = 'NOESYSpectrumFilter.TP')

        self.summary = {}

    def dump(self, filter_results, filename):

        check_list(filter_results)
        check_string(filename)

        import os

        lines = [self.line % self.header, '\n']
        lines += [self.write_peak(*x) for x in filter_results]

        ## summary

        summary = []
        
        for name in self.summary.keys():
            
            t = self.summary[name]['total']
            f = self.summary[name]['filtered']
            p1 = self.summary[name]['proton1']
            p2 = self.summary[name]['proton2']
            s = self.summary[name]['size']
            d = self.summary[name]['diagonal']
            u = self.summary[name]['unassigned']
             
            summary.append('## Spectrum: %s, %d / %d (%.1f%%) filtered out\n'\
                           % (name, f, t, 100. * f / t))
            summary.append('## Details:\n')
            summary.append('## - diagonal peak         : %d / %d (%.1f%%)\n'\
                           % (d, t, 100. * d / t))
            summary.append('## - invalid proton 1 shift: %d / %d (%.1f%%)\n'\
                           % (p1, t, 100. * p1 / t))
            summary.append('## - invalid proton 2 shift: %d / %d (%.1f%%)\n'\
                           % (p2, t, 100. * p2 / t))
            summary.append('## - invalid peak size     : %d / %d (%.1f%%)\n'\
                           % (s, t, 100. * s / t))
            summary.append('## - unassigned peak       : %d / %d (%.1f%%)\n'\
                           % (u, t, 100. * u / t))
            summary.append('##\n')

        lines = summary + lines
                           
        filename = os.path.expanduser(filename)
        if os.path.exists(filename):
            self.warning('Path "%s" already exists and will be overwritten.'\
                         % filename)
        file = open(filename, 'w')
        file.writelines(lines)
        file.close()

        self.message('Spectrum filter report written to file "%s"' % filename)
    
    def write_peak(self, peak, report):

        import aria.ChemicalShiftFilter as CSF
        import aria.CrossPeakFilter as CPF
        from aria.ariabase import YES, NO
        
        dict = {}
        summary = self.summary
        
        spectrum = peak.getSpectrum()

        if spectrum is None:
            dict['spectrum'] = str(None)
        else:
            dict['spectrum'] = str(spectrum.getName())

        spectrum_name = dict['spectrum']
        if not summary.has_key(spectrum_name):
            summary[spectrum_name] = {'total': len(spectrum), 'filtered': 0,
                                      'proton1': 0, 'proton2': 0,
                                      'diagonal': 0, 'size': 0, 'unassigned' : 0}

        dict['peak'] = str(peak.getNumber())
        summary[spectrum_name]['filtered'] += 1
        
        ## proton 1

        messages = report[peak.getProton1ChemicalShift()]

        if CSF.INVALID_SHIFT_VALUE in messages:
            dict['proton1'] = 'invalid'            
        else:
            dict['proton1'] = 'valid'
            
        if CSF.NEGATIVE_SHIFT_ERROR in messages:
            dict['error1'] = 'negative'
        else:
            dict['error1'] = 'valid'

        summary[spectrum_name]['proton1'] += (dict['proton1'] <> 'valid') & \
                                             (dict['error1'] <> 'valid')

        ## proton 2

        messages = report[peak.getProton2ChemicalShift()]

        if CSF.INVALID_SHIFT_VALUE in messages:
            dict['proton2'] = 'invalid'
        else:
            dict['proton2'] = 'valid'

        if CSF.NEGATIVE_SHIFT_ERROR in messages:
            dict['error2'] = 'negative'
        else:
            dict['error2'] = 'valid'

        summary[spectrum_name]['proton2'] += (dict['proton2'] <> 'valid') & \
                                             (dict['error2'] <> 'valid')

        ## peak size

        if report['size'] == CPF.NO_PEAK_SIZE:
            dict['size'] = 'invalid'
            summary[spectrum_name]['size'] += 1
        elif report['size'] == CPF.NEGATIVE_PEAK:
            dict['size'] = 'negative'
            summary[spectrum_name]['size'] += 1            
        else:
            dict['size'] = 'valid'

        ## diagonal peak

        if report['diagonal'] == CPF.DIAGONAL_PEAK:
            dict['diagonal'] = YES
            summary[spectrum_name]['diagonal'] += 1
        else:
            dict['diagonal'] = NO
            
        ## unassigned peak

        if report['unassigned'] == CPF.UNASSIGNED_PEAK:
            dict['unassigned'] = YES
            summary[spectrum_name]['unassigned'] += 1
        else:
            dict['unassigned'] = NO

        return self.line % dict

