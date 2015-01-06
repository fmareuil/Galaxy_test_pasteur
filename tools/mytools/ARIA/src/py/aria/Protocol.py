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
from aria.xmlutils import XMLElement, XMLBasePickler
from numpy import *
from time import clock
from time import ctime

## BARDIAUX 280105
REPORT_UPDATED_SPECTRUM = '.assigned'

REPORT_SUMMARY = 'report'
REPORT_NOE_RESTRAINTS = 'noe_restraints'
MISSING_STRUCTURES = 'missing_structures'
KEYBOARD_INTERRUPT = 'keyboard_interrupt'


## local path where MOLMOL restraints files are stored
PATH_MOLMOL = 'molmol'

def files_exist(file_list):
    import os

    for fn in file_list:
        fn = os.path.expanduser(fn)
        if not os.path.exists(fn):
            return 0

    return 1

def dump_peaks_as_text(peak_list, filename, gzip, options = None):
    
    check_type(peak_list, LIST, TUPLE)
    check_string(filename)
    check_int(gzip)
    check_type(options, NONE, DICT)

    import aria.AriaPeak as AriaPeak

    ## write contributions as text file (both, ambiguous and
    ## unambiguous in one file)

    ## TODO: so far, text-pickler options are hard-coded.

    settings = AriaPeak.AriaPeakListTextPicklerSettings()

    ## set options

    if options is not None:
        settings.update(options)
    
    pickler = AriaPeak.AriaPeakListTextPickler(settings)

    ## sort peak list such, that active restraints come first,
    ## inactive last.

    active = [p for p in peak_list if p.isActive()]
    inactive = [p for p in peak_list if not p.isActive()]

    peak_list = active + inactive

    ## Dump assignments

    pickler.dump_assignments(peak_list, filename + '.assignments')

    ## decompose peak_list into ambiguous and unambiguous
    ## restraints

    ambig = [p for p in peak_list if p.isAmbiguous()]
    unambig = [p for p in peak_list if not p.isAmbiguous()]

    pickler.dump_ambiguous(ambig, filename + '.ambig', gzip)
    pickler.dump_unambiguous(unambig, filename + '.unambig', gzip)

    ## Dump restraint list, sorted by violations.
    ## Active restraints come first, inactive last
    
    violated = [p for p in peak_list if p.analysis.isViolated()]

    if violated:
        
        d_lower = []
        d_upper = []

        lower = []
        upper = []

        for r in violated:

            ## In case of a lower-bound violation,
            ## put restraint in list of lower-bound-violating
            ## restraints
            
            d = r.analysis.getLowerBoundViolation()[0]
            
            if d is not None and d > 0.:
                d_lower.append(d)
                lower.append(r)

            else:
                d = r.analysis.getUpperBoundViolation()[0]

                if d is not None and d > 0.:
                    d_upper.append(d)
                    upper.append(r)

        lower_indices = argsort(d_lower).tolist()
        lower_indices.reverse()

        upper_indices = argsort(d_upper).tolist()
        upper_indices.reverse()

        a = [lower[i] for i in lower_indices]
        b = [upper[i] for i in upper_indices]

        viol = a + b
            
    else:
        viol = []

    active = [p for p in viol if p.isActive()]
    inactive = [p for p in viol if not p.isActive()]
        
    viol = active + inactive
            
    pickler.dump_violations(viol, filename + '.violations', gzip)

class IterationSettings(Settings):

    def create(self):

        from aria.Settings import TypeEntity, ChoiceEntity, NonNegativeInt

        kk = {}

        kk['number'] = NonNegativeInt()
        kk['number_of_structures'] = NonNegativeInt()
        kk['number_of_best_structures'] = NonNegativeInt()
        kk['number_of_kept_structures'] = NonNegativeInt()

        choices = ['total_energy', 'restraint_energy','noe_violations', 'restraint_violations'] # BARDIAUX 2.2
        msg = 'Sort structures according to one of the following' + \
              'choices: %s' % ', '.join(choices)
        entity = ChoiceEntity(choices, error_message = msg, description = msg)
        kk['sort_criterion'] = entity

        msg = 'Argument has to be a %s instance.'

        d = {'calibrator_settings': 'CalibratorSettings',
             'peak_assigner_settings': 'PeakAssignerSettings',
             'violation_analyser_settings': 'ViolationAnalyserSettings',
             'contribution_assigner_settings': 'ContributionAssignerSettings',
             'merger_settings': 'MergerSettings',
             'network_anchoring_settings' : 'NetworkSettings'}

        for name, type in d.items():
            entity = TypeEntity(type, error_message = msg % type)
            kk[name] = entity

        return kk

    def create_default_values(self):

        import aria.Calibrator as Calibrator
        import aria.PeakAssigner as PeakAssigner
        import aria.ViolationAnalyser as VA
        import aria.ContributionAssigner as CA
        import aria.Merger as Merger
        ## BARDIAUX 2.2
        import aria.Network as Network

        d = {}

        d['calibrator_settings'] = Calibrator.CalibratorSettings()
        d['peak_assigner_settings'] = PeakAssigner.PeakAssignerSettings()
        d['violation_analyser_settings'] = VA.ViolationAnalyserSettings()
        d['contribution_assigner_settings'] = CA.ContributionAssignerSettings()
        d['merger_settings'] = Merger.MergerSettings()
        ## BARDIAUX 2.2
        d['network_anchoring_settings'] = Network.NetworkSettings()

        for value in d.values():
            value.reset()

        d['number'] = 0
        d['number_of_structures'] = 20
        d['number_of_best_structures'] = 7
        d['number_of_kept_structures'] = 10
        d['sort_criterion'] = 'total_energy'

        return d


class IterationSettingsXMLPickler(XMLBasePickler):

    order = ['number', 'n_structures', 'sort_criterion',
             'n_best_structures', 'n_kept_structures', 'assignment', 'merging',
             'calibration', 'violation_analysis',
             'partial_assignment', 'network_anchoring']
    
    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)

        e.number = x['number']
        e.n_structures = x['number_of_structures']
        e.sort_criterion = x['sort_criterion']
        e.n_best_structures = x['number_of_best_structures']

        e.assignment = x['peak_assigner_settings']
        e.merging = x['merger_settings']
        e.calibration = x['calibrator_settings']
        e.violation_analysis = x['violation_analyser_settings']
        e.partial_assignment = x['contribution_assigner_settings']
        
        # BARDIAUX 2.2
        e.network_anchoring = x['network_anchoring_settings']
        e.n_kept_structures = x['number_of_kept_structures']
        
        return e

    def load_from_element(self, e):

        s = IterationSettings()
        
        s['number'] = int(e.number)
        s['number_of_structures'] = int(e.n_structures)
        s['sort_criterion'] = str(e.sort_criterion)
        s['number_of_best_structures'] = int(e.n_best_structures)
        s['peak_assigner_settings'] = e.assignment
        s['merger_settings'] = e.merging
        s['calibrator_settings'] = e.calibration
        s['violation_analyser_settings'] = e.violation_analysis
        s['contribution_assigner_settings'] = e.partial_assignment
        
        ## BARDIAUX 2.2
        if hasattr(e, 'network_anchoring'):
            
            s['network_anchoring_settings'] = e.network_anchoring
        else:
            import aria.Network as Network
            nas = Network.NetworkSettings()
            nas.reset()
            s['network_anchoring_settings'] = nas

        ## BARDIAUX 2.2
        if hasattr(e, 'n_kept_structures'):
            s['number_of_kept_structures'] = int(e.n_kept_structures)
        else:
            s['number_of_kept_structures'] = 0

        return s

class ProtocolSettings(Settings):

    def create(self):

        from aria.Settings import TypeEntity, ChoiceEntity

        d = {}

        ## water-refinement parameters

        data_type = 'WaterRefinementParameters'
        msg = 'Argument has to be a %s instance.' % data_type
        entity = TypeEntity(data_type, error_message = msg)

        d['water_refinement'] = entity
        
        ## iteration settings

        entity = TypeEntity(DICT)
        entity.set({})

        d['iteration_settings'] = entity

        ## floating assignment flag

        choices = [YES, NO]
        msg = 'Floating assignment: %s?' % ' / '.join(choices) 
        entity = ChoiceEntity(choices, error_message = msg)

        d['floating_assignment'] = entity

        return d

    def create_default_values(self):
        return {'floating_assignment': YES}

    def addIterationSettings(self, settings):

        check_type(settings, 'IterationSettings')

        number = settings['number']
        it_settings = self['iteration_settings']

        if not it_settings.has_key(number):
            it_settings[number] = settings
        else:
            msg = 'IterationSettings instance with same number ' + \
                  '(%d) already stored in protocol settings.' % number
            self.error(KeyError, msg)

    def delIterationSettings(self, s):

        check_type(s, 'IterationSettings')
        
        number = s['number']
        it_settings = self['iteration_settings']

        if it_settings.has_key(number):
            del it_settings[number]
        else:
            msg = 'IterationSettings (number %d) d not exist.' % number
            self.error(KeyError, msg)

class Protocol(AriaBaseClass):

    def __init__(self, settings):

        check_type(settings, 'ProtocolSettings')

        AriaBaseClass.__init__(self, name = 'Protocol')
        self.setSettings(settings)

        self.__reports_written = 0
        self.constraintSetSerial = None

##        self.iterations = {}

    def getIterationSettings(self, number):

        check_int(number)
        return self.getSettings()['iteration_settings'][number]

    def getInfrastructure(self):
        from aria.Singleton import ProjectSingleton

        project = ProjectSingleton()
        return project.getInfrastructure()

    def initialize(self):
        """
        must be called after unpickling a protocol in
        order to finalize its initialization.
        create objects which implement the actual aria-method
        set parameters of first iteration.
        """

        from aria.NOEModel import NOEModel, ISPA
        import aria.Calibrator as Calibrator
        import aria.ViolationAnalyser as VA
        import aria.ContributionAssigner as CA
        import aria.Merger as Merger
        ## BARDIAUX 2.2
        import aria.Network as Network
        
        ## NOE model
        ## TODO: not nessecary to store model!

        #self.model = NOEModel()
        self.model = ISPA()

        ## calibrator

        cs = Calibrator.CalibratorSettings()
        self.calibrator = Calibrator.Calibrator(cs, self.model)

        ## violation-analyser

        vas = VA.ViolationAnalyserSettings()
        self.violation_analyser = VA.ViolationAnalyser(vas)

        ## contribution assigner

        cas = CA.ContributionAssignerSettings()
        self.contribution_assigner = CA.ContributionAssigner(cas)

        ## merger

        ms = Merger.MergerSettings()
        self.merger = Merger.Merger(ms)

        ## BARDIAUX 2.2
        ## network
        nas = Network.NetworkSettings()
        self.network_anchoring = Network.NetworkAnchoring(nas)        

    def _updateIterationSettings(self, iteration):
        """
        applies the parameters-settings in 'iteration' to
        sub-modules
        """

        check_type(iteration, 'Iteration')

        settings = self.getIterationSettings(iteration.getNumber())

        cs = settings['calibrator_settings']
        self.calibrator.getSettings().update(cs)

        vas = settings['violation_analyser_settings']
        self.violation_analyser.getSettings().update(vas)

        cas = settings['contribution_assigner_settings']
        self.contribution_assigner.getSettings().update(cas)

        ms = settings['merger_settings']
        self.merger.getSettings().update(ms)

        ## BARDIAUX 2.2
        nas = settings['network_anchoring_settings']
        self.network_anchoring.getSettings().update(nas)
        
        
    def _updateSpectrumSettings(self, spectrum):
        spectrum_data = spectrum.getDataSource()
        self.calibrator.getSettings().update(spectrum_data)
        self.violation_analyser.getSettings().update(spectrum_data)

    def mergePeakLists(self, iteration):
        
        ## BARDIAUX 2.2
        # remove previous combined peaks
        iteration.setCombinedRestraints([])
        
        if self.merger.getSettings()['method'] == 'no_merging':
            return

        peaks = [p for p in iteration.getPeakList() if p.isActive()]

        ## merge peaks. peaks which are merged (away), i.e.
        ## are marked via peak.isMerged(1)

        n, combined = self.merger(peaks)

        # add combined restraints
        iteration.setCombinedRestraints(combined)

        if n:
            msg = 'Peaks merged: %d / %d peaks (%.1f %%) (see report for details).'
            self.message(msg % (n, len(peaks), 100. * n / len(peaks)))

        if combined:
            msg = 'Peaks combined: %d / %d peaks (%.1f %%).'
            self.message(msg % (len(combined), len(peaks), 100. * len(combined) / len(peaks)))

##             self.message('old %d new %d.' % (len(peaks), len([p for p in iteration.getPeakList() if p.isActive()])))   
            
        
    def compileFileList(self, number, **kw):
        """
        Compiles a tuple of filenames which must exist after
        iteration 'number' has been successfully completed.
        If the method is called with keywords, the returned tuple
        contains only a sub-set of all filenames. The sub-set depends
        on the particular keyword(s)
        """

        check_int(number)

        STRUCTURES = 'structures'
        FLOAT_FILES = 'float_files'

        from os.path import join

        keywords = (STRUCTURES, FLOAT_FILES)

        unknown = [key for key in kw if key not in keywords]
            
        if unknown:
            s = 'Unknown sub-list identifier ("%s"). Known identifier ' + \
                'are: %s'
            self.error(s % (str(unknown), str(keywords)))

        if kw:
           keywords = kw.keys() 
        
        infra = self.getInfrastructure()

        ## structures (pdb-files)

        it_settings = self.getIterationSettings(number)
        n_structures = it_settings['number_of_structures']
            
        path = infra.get_iteration_path(number)

        l = []

        if STRUCTURES in keywords:
            name_template = infra.get_file_root() + '_%d.pdb'
            l += [join(path, name_template % j) \
                  for j in range(1, n_structures + 1)]

        if FLOAT_FILES in keywords:
            name_template = infra.get_file_root() + '_%d.float'
            l += [join(path, name_template % j) \
                  for j in range(1, n_structures + 1)]

        return tuple(l)

    def findFirstIteration(self):
        """
        For every iteration, this method checks if the required
        PDB-files exist. If so, it seeks to the next iteration until
        it ends up at an iteration with missing files. If files are
        present for all iterations, None is returned. Otherwise the
        number of the iteration with missing files is returned.
        """

        infra = self.getInfrastructure()
        n_iterations = infra.getSettings()['n_iterations']
        
        for i in range(n_iterations):

            ## Get filenames of all structures for iteration i

            pdb_files = self.compileFileList(i, structures=1)
            if not files_exist(pdb_files):
                return i

        return None

    def _writePDBFileList(self, iteration):
        import os

        infra = self.getInfrastructure()
        
        iteration_path = infra.get_iteration_path(iteration.getNumber())
        molmol_path = os.path.join(iteration_path, PATH_MOLMOL)
        
        ok = 1
        
        if not os.path.exists(molmol_path):
                
            try:
                os.makedirs(molmol_path)

            except Exception, msg:
                import aria.tools as tools
                self.warning(tools.last_traceback())
                msg = 'Could not create directory for MOLMOL output.'
                self.warning(msg)
                ok = 0
                
            if not ok:
                return None

        ## file.nam

        ## files are sorted according to 'sort_criterion'
        ## which is part of the <structure_generation>
        ## section of project xml.

        ensemble = iteration.getStructureEnsemble()
        file_list = ensemble.getFiles()
            
        if file_list is None:
            s = 'Could not write MolMol file.nam: No PDB-files' + \
                'in structure ensemble.'
            self.warning(s)
            return None

        filename = os.path.join(molmol_path, 'file.nam')
        
        try:
            f = open(filename, 'w')
        except:
            self.warning('Could not create %s.' % filename)
            return None
        
        s = '\n'.join(file_list)
        
        try:
            f.write(s)
            f.close()
            
        except:
            self.warning('Could not write PDB-file list (%s).' % filename)
            return None

        return molmol_path

    def structure_calc_done(self, engine, molecule, iteration):

        ## Check if there are missing structures

        if engine.missingStructures():
            self.__done = MISSING_STRUCTURES
            return
        
        self.message('Structure calculation done.')

        it_settings = self.getIterationSettings(iteration.getNumber())
        ensemble = engine.getEnsemble(it_settings, molecule)

        ## store structure-ensemble in iteration

        iteration.setStructureEnsemble(ensemble)

        self.__done = 1

    def startStructureCalculation(self, iteration, molecule, kept_structures = []):

        from aria.Singleton import ProjectSingleton

        ## get list of all restraints
        
        peaks = iteration.getPeakList()


        ## BARDIAUX 2.2 Add user-distance CCPN restraints
        restraints = iteration.getDistanceRestraintsList()
        peaks += restraints
        # add also combined restraints
        combined = iteration.getCombinedRestraints()
        peaks += combined
        
        project = ProjectSingleton()
        engine = project.getStructureEngine()

        ## callback that will be called when pdb-files have been generated.

        f = lambda engine = engine, molecule = molecule, it = iteration, \
            c = self.structure_calc_done: c(engine, molecule, it)
        
        engine.set_callback(f)
        ## < Mareuil
        MODE = self.getMODE()
        engine.go(peaks, self.getSettings(), iteration.getNumber(), kept_structures = kept_structures, MODE=MODE)
        ## Mareuil >

    def doSolventRefinement(self, iteration, molecule):

        from aria.Singleton import ProjectSingleton
        import aria.StructureEnsemble as SE
            
        ## If the given iteration lacks a structure ensemble,
        ## load PDB-files from last iteration.
        ## < Mareuil
        MODE = self.getMODE()
        ## Mareuil >
        ensemble = iteration.getStructureEnsemble()

        if ensemble is None:
        
            ## create structure-ensemble from given 'filenames'
            number = iteration.getNumber()
            it_settings = self.getIterationSettings(number)
            se_settings = SE.StructureEnsembleSettings()

            ## Set iteration-specific settings
            
            se_settings.update(it_settings)
            
            ensemble = SE.StructureEnsemble(se_settings)
            
            filenames = self.compileFileList(number, structures=1)
            naming_convention = 'cns'
            ensemble.read(filenames, molecule, naming_convention)
            
            iteration.setStructureEnsemble(ensemble)
        
        project = ProjectSingleton()
        engine = project.getStructureEngine()
        ## < Mareuil
        file_list = engine.solvent_refine(self.getSettings(), ensemble, MODE=MODE)
        ## Mareuil >
        ## cleanup (for solvent refinement and last iteration)

        infra = self.getInfrastructure()
        last_it = infra.getSettings()['n_iterations'] - 1
        
        engine.cleanup(infra.get_refinement_path())
        engine.cleanup(infra.get_iteration_path(last_it))

    def doViolationAnalysis(self, restraints, ensemble, store_analysis = 0):
        """
        'restraints': list of AriaPeaks. The bounds of every
        restraint in that list is checked against distances found
        in the 'ensemble'.
        
        'targets': list of AriaPeaks. The violationAnalyser will
        store all intermediate results in their analysis-section.
        Note: we assume, that peaks[i] corresponds to results[i]
        for all i !. If a restraint has been violated, the
        corresponding 'target'_restraint will be marked as violated.
        """
        
        f = self.violation_analyser.analysePeak

        violated = []
        non_violated = []

        ## get theshold for current iteration
        settings = self.violation_analyser.getSettings()
        threshold = settings['violation_threshold']
	#Mareuil
        if settings['sigma_mode'] == 'auto':
            ecars = []
            for restraint in restraints:
                temp_ecars = self.violation_analyser.tolerance(ensemble, restraint, store_analysis)
                for ecar in temp_ecars:
                    ecars.append(ecar)
	    
            ecar_avg = sum(ecars)/len(ecars)
            tol = 0
            for n_ecar in ecars:
                tol = tol + power(n_ecar - ecar_avg,2)
            print "TOLERANCE ", power(tol/len(ecars),1./2)
            tol =  power(tol/len(ecars),1./2) * settings['violation_tolerance']
            print "AVG RESTRAINT ", tol
	  
        else:
            tol = None
	#Mareuil    
        for restraint in restraints:
            R_viol = f(ensemble, restraint, tol, store_analysis, sig_mode = settings['sigma_mode'])
            
            ##
            ## If a restraint has been violated in too many structures
            ## (according to 'threshold'), mark is a violated.
            ##

            if R_viol > threshold:
                restraint.analysis.isViolated(1)
                violated.append(restraint)
                
            else:
                restraint.analysis.isViolated(0)
                non_violated.append(restraint)

        ## For violated restraints: if bound-correction is enabled,
        ## repeat violation-analysis with modified bounds.

        if settings['lower_bound_correction']['enabled'] == YES:
            new_lower = settings['lower_bound_correction']['value']
        else:
            new_lower = None

        if settings['upper_bound_correction']['enabled'] == YES:
            new_upper = settings['upper_bound_correction']['value']
        else:
            new_upper = None

        if new_lower is not None or new_upper is not None:

            ## We forget 'store_analysis' here, since it has already
            ## been stored (if set).
            
            R_viol = [f(ensemble, r, lower_correction = new_lower,
                        upper_correction = new_upper) for r in violated]

            ## List of restraint-indices which are no longer
            ## violated after bound modification.

            indices = flatnonzero(less(R_viol, threshold))
            new_non_violated = [violated[i] for i in indices]
            
            [r.analysis.isViolated(0) for r in new_non_violated]

        else:
            new_non_violated = None

        return violated, non_violated, new_non_violated

    def done(self):
        self.message('finished.')

    def _dump_peak_list(self, peaks, path):

        check_list(peaks)
        check_string(path)
        
        import os
        from aria.Singleton import ProjectSingleton
        import aria.Merger as Merger

        project = ProjectSingleton()
        s = project.getReporter()['noe_restraint_list']

        ## Dump only non-merged restraitns

        not_merged = [r for r in peaks if not r.isMerged()]

        if s['xml_output'] in (YES, GZIP):
            
            from aria.AriaXML import AriaXMLPickler as pickler
            
            filename = os.path.join(path, REPORT_NOE_RESTRAINTS + '.xml')
            t = clock()
            p = pickler()

            # BARDIAUX 2.2: It's not possible to dump
            # CCPN restraints as XML
            not_merged = [r for r in not_merged if not is_type(r, 'DistanceRestraint')]
            
            gzip = {YES: 0, GZIP: 1}[s['xml_output']]
            p.dump(not_merged, filename, gzip = gzip)
            self.message('NOE-restraint list (xml) written (%s).' % filename)
            self.debug('Time: %ss' % str(clock()-t))

        if s['text_output'] in (YES, GZIP):
            t = clock()
            filename = os.path.join(path, REPORT_NOE_RESTRAINTS)

            gzip = {YES: 0, GZIP: 1}[s['text_output']]
            dump_peaks_as_text(not_merged, filename, gzip)
            self.message('NOE-restaint list (text) written (%s).' % filename)
            self.debug('Time: %ss' % str(clock()-t))
            
        if s['pickle_output'] in (YES, GZIP):

            from aria.tools import Dump
            
            filename = os.path.join(path, REPORT_NOE_RESTRAINTS + '.pickle')
            t = clock()

            gzip = {YES: 0, GZIP: 1}[s['pickle_output']]
            Dump(not_merged, filename, gzip = gzip)
            self.message('NOE-restraint list (python pickle) ' + \
                         'written (%s).' % filename)
            self.debug('Time: %ss' % str(clock()-t))

        ## Write Merger report

        t = clock()

        merged = [r for r in peaks if r.isMerged()]
        filename = os.path.join(path, REPORT_NOE_RESTRAINTS + '.merged')
        
        p = Merger.MergerTextPickler()
        p.dump(merged, filename)
        self.message('Report for Merging-step written (%s).' % filename)
        self.debug('Time: %ss' % str(clock() - t))


    def _dump_spectra(self, peaks, spectrum, filename):

        check_list(peaks)
        check_type(spectrum, 'NOESYSpectrum')
        check_string(filename)
        
        import os
        from aria.Singleton import ProjectSingleton
        from aria.NOESYSpectrum import NOESYSpectrum
        import aria.Assignment as Assignment
        from copy import deepcopy
        from aria.AriaXML import AriaXMLPickler as pickler

        project = ProjectSingleton()
        
        s = project.getReporter()['spectra']
            

        if s['write_unambiguous_only'] == YES:
            peaks = [pk for pk in peaks if not pk.isAmbiguous()]
        
        if peaks and s['write_assigned'] == YES:

            spectrum = deepcopy(spectrum)
                
            for pk in peaks:

                ref_pk = spectrum.findPeak(pk.getReferencePeak().getNumber())

                if ref_pk is None:
                    self.error('Inconsistency: Could not find reference peak %d in copy of spectrum "%s"' %
                               (ref_pk.getNumber(), spectrum.getName()))
                
                active_contribs = [c for c in pk.getContributions() \
                                   if c.getWeight() > 0.]
                
                if not ref_pk.isAssigned() or \
                       (ref_pk.isAssigned() and s['write_assigned_force'] == YES) and\
                       active_contribs:
                    
                    h = active_contribs[0]

                    ss = h.getSpinSystems()

                    for sp_sys in ss:

                        dim = ref_pk.getDimension(sp_sys)

                        AUTO = Assignment.ASSIGNMENT_TYPE_AUTOMATIC

                        A = Assignment.Assignment(sp_sys.getAtoms(), AUTO)

                        if dim == 1:                       
                            ref_pk.setProton1Assignments((A,))
                        elif dim == 2:
                            ref_pk.setProton2Assignments((A,))

                    for h in active_contribs[1:]:
                    
                        ss = h.getSpinSystems()

                        for sp_sys in ss:

                            dim = ref_pk.getDimension(sp_sys)
                            
                            AUTO = Assignment.ASSIGNMENT_TYPE_AUTOMATIC
                            
                            A = Assignment.Assignment(sp_sys.getAtoms(), AUTO)

                            if dim == 1:
                                atoms_assi = [a.getAtoms() for a in ref_pk.getProton1Assignments()]
                                if not A.getAtoms() in atoms_assi:
                                    ref_pk.addProton1Assignment(A)
                            elif dim == 2:
                                atoms_assi = [a.getAtoms() for a in ref_pk.getProton2Assignments()]
                                if not A.getAtoms() in atoms_assi:
                                    ref_pk.addProton2Assignment(A)
                                
                
            t = clock()
            p = pickler()

            p.dump(spectrum, filename, gzip = 0)
            self.debug('Time: %ss' % str(clock()-t))

            self.message('Assigned spectrum "%s" written to file %s"' \
                         % (spectrum.getName(), filename))

    def _dump_ccpn(self, iteration, path, is_water_refinement=False):

        from aria.Singleton import ProjectSingleton

        project = ProjectSingleton()

        if project.ccpn_project_instance is None:
            return

##         # BARDIAUX 2.3 : when solvent molec present, abort CCPN PDB export
##         write_solvent = self.getSettings()['water_refinement']['write_solvent_molecules']
##         if is_water_refinement and write_solvent == YES:
##             msg = 'CCPN export: Solvent refined structures contain solvent molecules.' + \
##                   'PDB files will not be exported.'
##             self.warning(msg)
##             return

        # BARDIAUX 2.3
        try:
            import aria.exportToCcpn as ccpnExport
        except:
            self.warning('CCPNmr Analysis not found. CCPN export aborted.')
            return
        
        from aria.importFromCcpn import getKeysFromString

        project_settings = project.getSettings()
        run_name = project_settings['run'].strip().replace(' ', '_')
        file_root = project_settings['file_root'].strip().replace(' ', '_')

        ccpn_project = project.ccpn_project_instance
        nmrProject = ccpn_project.currentNmrProject or ccpn_project.findFirstNmrProject()
        
        ## TBD: proper handling of NMR projects, e.g. the case in which
        ##      multiple NMR projects exist.

        if not nmrProject:
            self.message('CCPN export: No NMR project found, creating new one.')

            name = 'ARIA2_run_%s' % run_name

            nmrProject = ccpn_project.newNmrProject(name=name)
            
        s = project.getReporter()['ccpn']

        export = 0

        # TJS: Includes inactive ones as these are
        # now passed back to CCPN (and marked) 
        restraints = iteration.getPeakList()

        # BARDIAUX: Update for multiple chains
        aria_chain = project.getMolecule().get_chains()[0]
        chains = ccpnExport.getChainsFromAria2(restraints, ccpn_project,
                                               aria_chain=aria_chain)

        # TJS: Get Run object for CCPN to group objects
        # If data comes from CCPN, this would have been made
        # in Extend-NMR GUI or at load time
        if chains:                                           
            ccpnMolSystem = chains[0].molSystem                                     
            ccpnAriaRun = ccpnExport.getAriaRun(ccpnMolSystem)

	    ask_water = self.getSettings()['water_refinement']['enabled'] == YES
            is_last_iteration = iteration.getNumber() == project.getSettings()['n_iterations']-1
            
            if ccpnAriaRun:
                # TJS: Completed runs will not be used again after this
                if (is_water_refinement and ask_water) or \
                       (is_last_iteration and not ask_water):
                    ccpnAriaRun.status = 'completed'
	    
        else:
            ccpnAriaRun = None

        # TJS: Export structures before restraints & peaks thus  
        # we have to link NmrConstraintStore to the structure generation
        # at the very end
        
        struct_gen = None
        if s['export_structures'] == YES and \
               (is_last_iteration or is_water_refinement):

            export = 1
            
            # TJS: Moved above beacause chains needed for CCPn Run Object
            # chains = ccpnExport.getChainsFromAria2(restraints, ccpn_project,
            #                                        aria_chain=aria_chain)

            if not chains:
                self.warning('CCPN export: No molecular system found.')
                return

            structures = ccpnExport.getStructuresFromAria2Dir(path, chains)

            if not structures:
                self.warning('CCPN export: Unable to load any structures from iteration directory %s' % path)

            else:
                self.message('CCPN export: PDB files exported.')
                

                if not is_water_refinement:
                    name = 'ARIA2_run_%s_it%d' % (run_name, iteration.getNumber())
                    details = 'Structures created by ARIA2, %s run "%s", iteration %d.' % \
                                                (file_root, run_name, iteration.getNumber())

                else:
                    name = 'ARIA2_run_%s_water_refine' % (run_name,)
                    details = 'Structures created by ARIA2, %s run "%s", water refinement.' % \
                                                                           (file_root, run_name)
                
                
                # TJS: Structure generation is the old way of storing results
                struct_gen = ccpnExport.makeStructureGeneration(structures, None)
                struct_gen.name = name
                struct_gen.details = details
                    
                # TJS: Store ARIA generated data as output of CCPN NmrCalc run
                if ccpnAriaRun:
                    dataObj = ccpnAriaRun.newStructureEnsembleData(ensembleId=structures.ensembleId,
                                                                   molSystemCode=chains[0].molSystem.code,
                                                                   ioRole='output', name=name,
                                                                   details=details)
                                                             
        else:
            structures = None


        if not is_water_refinement and (s['export_noe_restraint_list'] == 'all' or \
                                        (s['export_noe_restraint_list'] == 'last' and \
                                         is_last_iteration)):
        
            export = 1


            if self.constraintSetSerial is not None:
                constraintSet = nmrProject.findFirstNmrConstraintStore(serial=self.constraintSetSerial)
                self.message('CCPN export: Using existing constraint set.')

            else:
                self.message('CCPN export: Creating new constraint set.')
                constraintSet = ccpnExport.makeNmrConstraintStore(nmrProject)
                self.constraintSetSerial = constraintSet.serial
            
            # TJS: Link back to any structure generation object
            if struct_gen:
              struct_gen.nmrConstraintStore = constraintSet


            # BARDIAUX: Update for multiple chains
            chains = ccpnExport.getChainsFromAria2(restraints, ccpn_project,
                                                   aria_chain=aria_chain)

            if restraints:
              # TJS expert rejected (not active) restraints back to CCPN for analysis
              constrList, rejectList, violList = ccpnExport.getConstraintsFromAria2(restraints, chains, constraintSet,
                                                          structures)

              constrList.name = 'ARIA2_run%s_NOEs_it%d' % (run_name, iteration.getNumber())
              constrList.details = 'ARIA2 NOEs, project "%s", run "%s", iteration %d.' % \
                                               (file_root, run_name, iteration.getNumber())

              # TJS: Store ARIA generated data as output of CCPN NmrCalc run
              if ccpnAriaRun:
                  cSet = constraintSet.serial
                  dataObj = ccpnAriaRun.newConstraintStoreData(constraintStoreSerial=cSet,
                                            constraintListSerials=[constrList.serial],
                                            ioRole='output', name=constrList.name,
                                            details=constrList.details)
                  
              # TJS constraint list contaning the inactive restraints
              if rejectList is not None:
                  rejectList.name = 'ARIA2_REJECT_run%s_NOEs_it%d' % (run_name, iteration.getNumber())
                  rejectList.details = 'ARIA2 *INACTIVE* NOEs, project "%s", run "%s", iteration %d.' % \
                                                              (file_root, run_name, iteration.getNumber())
 
                  # TJS: Store ARIA generated data as output of CCPN NmrCalc run
                  if ccpnAriaRun:
                      cSet = constraintSet.serial
                      dataObj = ccpnAriaRun.newConstraintStoreData(constraintStoreSerial=cSet,
                                                    constraintListSerials=[rejectList.serial],
                                                    ioRole='output', name=rejectList.name,
                                                    details=rejectList.details)
 
              self.message('CCPN export: NOE restraint list exported.')

            # BARDIAUX export Constraints from CCPN
            for orig_list, distance_constraints in iteration.getDistanceRestraints().items():
                

                # TJS adjust; we want the inactive too. ;-)
                # distance_constraints = [r for r in distance_constraints if r.isActive()]
                constrList, rejectList, violList = ccpnExport.getConstraintsFromAria2(distance_constraints, chains,
                                                                                      constraintSet, structures)
                                                                                      
                 # TJS adjust
                ccpnKeys = getKeysFromString(orig_list.getDataSource().get('ccpn_id'))
                list_id  = '%s|%s' % (ccpnKeys[-2],ccpnKeys[-1])# Project id is too verbose and redundant inside same project
                
                constrList.name = 'ARIA2_run%s_Dists_%s_it%d' % (run_name, list_id, iteration.getNumber())
                constrList.details = 'ARIA2 Distances %s, project "%s", run "%s", iteration %d.' % \
                                        (orig_list.getName(), file_root, run_name, iteration.getNumber())

                # TJS: Store ARIA generated data as output of CCPN NmrCalc run
                if ccpnAriaRun:
                    cSet = constraintSet.serial
                    dataObj = ccpnAriaRun.newConstraintStoreData(constraintStoreSerial=cSet,
                                                  constraintListSerials=[constrList.serial],
                                                  ioRole='output', name=constrList.name,
                                                  details=constrList.details)
                  
                # TJS constraint list contaning the inactive restraints
                if rejectList is not None:
                    rejectList.name = 'ARIA2_REJECT_run%s_Dists_%s_it%d' % (run_name, list_id, iteration.getNumber())
                    rejectList.details = 'ARIA2 *INACTIVE* Distances %s, project "%s", run "%s", iteration %d.' % \
                                                   (orig_list.getName(), file_root, run_name, iteration.getNumber())
 
                    # TJS: Store ARIA generated data as output of CCPN NmrCalc run
                    if ccpnAriaRun:
                        cSet = constraintSet.serial
                        dataObj = ccpnAriaRun.newConstraintStoreData(constraintStoreSerial=cSet,
                                                      constraintListSerials=[rejectList.serial],
                                                      ioRole='output', name=rejectList.name,
                                                      details=rejectList.details)

            self.message('CCPN export: DistanceConstraints list exported.')
            
        if not is_water_refinement and (s['export_assignments'] == YES  and is_last_iteration):

            import aria.DataContainer as DC
            
            export = 1

            ## compile dict that maps peaklist keys to name of new peaklist

            ccpn_spectra = [x for x in project.ccpn_data_sources.get(DC.DATA_SPECTRUM, ())]

            names_dict = {}

            template = 'ARIA2_NOE_Peaks_run%s_it%d_%s'

            for x in ccpn_spectra:

                ccpn_id = x['peaks']['ccpn_id']

                args = run_name, iteration.getNumber(), ccpn_id
                names_dict[ccpn_id] = template % args

            peak_list_keys = ccpnExport.getPeakAssignmentsFromAria2(ccpn_project, restraints, namesDict=names_dict,
                                                                    aria_chain=aria_chain)

            # TJS: Store ARIA generated data as output from run
            if ccpnAriaRun and peak_list_keys:
                from aria.importFromCcpn import getCcpnPeakList
                
                for peak_list_key in peak_list_keys:
                    ccpnPeakList = getCcpnPeakList(ccpn_project, peak_list_key)
                   
                    if ccpnPeakList:
                       spectrum = ccpnPeakList.dataSource
                       eSerial = spectrum.experiment.serial
                       sSerial = spectrum.serial
                       pSerial = ccpnPeakList.serial
                       details = ccpnPeakList.details
                       plId = '%d|%d|%d' % (eSerial, sSerial, pSerial)
                       name = template % (run_name, iteration.getNumber(), plId)
                       
                       ccpnAriaRun.newPeakListData(experimentSerial=eSerial,
                                                   dataSourceSerial=sSerial,
                                                   peakListSerial=pSerial,
                                                   ioRole='output',
                                                   details=details,
                                                   name=name[:80])
                      
                   

            self.message('CCPN export: NOE assignments exported.')

 

        if export:

            try:
                ccpn_project.saveModified()

            except Exception, msg:

                self.warning('CCPN export: Could not save project. Error message was: %s' % msg) #by AWSS, msg is an object not a string


            self.message('CCPN project saved.')

    def _dump_iteration(self, iteration, path):

        import os
        import aria.Iteration as Iteration

        filename = os.path.join(path, REPORT_SUMMARY)
        pickler = Iteration.IterationTextPickler()
        pickler.dump(iteration, filename)

    ## BARDIAUX 2.2
    def _dump_rms(self, peaks, iteration_number):
        
        import aria.RmsReport as RmsReport
        infra = self.getInfrastructure()
        text_path = infra.get_iteration_path(iteration_number)
        graphics_path = infra.get_iteration_graphics_path(iteration_number)
        
        rp = RmsReport.RmsReport(peaks, iteration_number, text_path, graphics_path)
        rp.go()
        

    def dumpIteration(self, iteration):

        import os
        from aria.Singleton import ProjectSingleton
        
        project = ProjectSingleton() 
        check_type(iteration, 'Iteration')

        infra = self.getInfrastructure()
        path = infra.get_iteration_path(iteration.getNumber())
        
        # BARDIAUX 2.2
        # dump also CCPN DistanceRestraints
        # all_peaks = iteration.getPeaks()
        all_peaks = {}
        all_peaks.update(iteration.getPeaks())
        all_peaks.update(iteration.getDistanceRestraints())
        
        order = all_peaks.keys()
        order.sort()

        peak_list = []

        for spectrum in order:

            peaks = all_peaks[spectrum]
            
            ## sort peak_list wrt ref_peak number

            peaks.sort(lambda a, b, c = cmp: \
                       c(a.getReferencePeak().getNumber(), \
                         b.getReferencePeak().getNumber()))

            
            peak_list += peaks
            
            ## BARDIAUX 25/04/05
            # Dump updated spectrum

            if is_type(spectrum, 'ConstraintList'):
                # Not for CCPN ConstraintList
                continue
            
            spec_name = spectrum.getName()
            spec_name = spec_name.replace(' ','_')
            
            filename = os.path.join(path, spec_name + \
                                    REPORT_UPDATED_SPECTRUM +'.xml')

       
            s = project.getReporter()['spectra']
            
            if s['iteration'] == 'all'  or \
                   (s['iteration'] == 'last' and  \
                    iteration.getNumber() == project.getSettings()['n_iterations']-1):
                
                self._dump_spectra(peaks, spectrum, filename)


        self._dump_peak_list(peak_list, path)
        self._dump_iteration(iteration, path)
        self._dump_ccpn(iteration, path)
        
        ## BARDIAUX 2.2
        ## RMS report
        
        if peak_list and iteration.getNumber() > 0:
            
            try:
                self._dump_rms(peak_list, iteration.getNumber())
            #except:
            #    self.message('Error during RMS analysis/graphics generation.')
            except Exception, msg:
                import aria.tools as tools
                self.warning(tools.last_traceback())
                msg = 'Error during RMS analysis/graphics generation.'
                self.warning(msg)
        
        self.__reports_written = 1

    def reportsWritten(self):
        """
        returns 1 if any report files have been written.
        """

        return self.__reports_written

    def _get_peak_sizes(self, peaks):
        """
        depending on the calibrator-setting 'volume_or_intensity',
        the method returns volumes or intensities
        """
            
        ref_peaks = [p.getReferencePeak() for p in peaks]

        if self.calibrator.getSettings()['volume_or_intensity'] == 'volume':
            peak_sizes = [p.getVolume()[0] for p in ref_peaks]
        else:
            peak_sizes = [p.getIntensity()[0] for p in ref_peaks]

        return peak_sizes

    # Bardiaux rMat
    def _get_peak_theoric_volumes(self, peaks):

        """
        get the theoric volume from the spin diffusin correction.
        """
        return array([p.getTheoricVolume() for p in peaks])

    # BARDIAUX 2.2
    # get settings of DistanceRestraint sourec list
    def __getListSource(self, p):
        return p.getReferencePeak().getSpectrum().getListSource()
    
    def setModelIntensities(self, restraints, ensemble, calibration_factor):

        from aria.Datum import Datum
        
        f = self.model.calculatePeaksize

        model_peak_sizes = array([f(r, ensemble) for r in restraints])
        calculated_peak_sizes = model_peak_sizes * calibration_factor

        for i in range(len(restraints)):
            d = Datum(calculated_peak_sizes[i], None)
            restraints[i].analysis.setCalculatedPeaksize(d)

    def calculateBounds(self, factor, peaks, bound_corrected = None, ensemble = None):
        """
        calculate lower- and upper bounds for every peak using
        the calibration 'factor'. values are stored.
        'bound_corrected': list of restraints which are classified
        as correct after bound-modification.
        """
        # BARDIAUX 2.2
        # ConstraintList can bypass the calibration
        if is_type(peaks[0], 'DistanceRestraint'):
            calibrate = self.__getListSource(peaks[0])['calibrate']
            if calibrate == NO or \
                   (self.findFirstIteration() == 0 and calibrate == 'all_iterations_except_first'):                
                return
            
        factor = power(factor, 1./6)

        peak_sizes = self._get_peak_sizes(peaks)
        
        # Malliavin/Bardiaux rMat
        cs = self.calibrator.getSettings()
        if cs['relaxation_matrix'] == YES and ensemble is not None: 

##             from NOEModel import ISPA
##             ispa = ISPA()
##             f = ispa.calculatePeaksize

##             ispa_peak_sizes = array([f(r, ensemble) for r in peaks])
            
            ispa_peak_sizes =   array([p.getIspa() for p in peaks])
            peak_theoric_vol = self._get_peak_theoric_volumes(peaks)

            ratio = ispa_peak_sizes / peak_theoric_vol
            distances = factor * power(peak_sizes * ratio, -1. / 6)


        else:

            distances = factor * power(peak_sizes, -1. / 6)
        
        ## TODO: hard-coded 0.125

        if cs['error_estimator'] == 'intensity':
            errors = 0.125 *  power((factor * power(peak_sizes, -1. / 6)), 2.)

        else:
            errors = 0.125 * power(distances, 2.)

##         # BARDIAUX : To be implemented
##         if cs['error_estimator'] == 'no_errors':
##             errors = [0.] * len(distances)
            
#        errors = [0.] * len(distances)
            
##         distances = factor * power(peak_sizes, -1. / 6)
##         ## TODO: hard-coded 0.125
##         errors = 0.125 * power(distances, 2.)
        
        ## lower bounds are >= 0.

        lower_bounds = clip(distances - errors, 0., 1.e10)
        upper_bounds = distances + errors

        for i in range(len(peaks)):

            peak = peaks[i]
                
            peak.setDistance(distances[i])
            peak.setLowerBound(lower_bounds[i])
            peak.setUpperBound(upper_bounds[i])

        ## Set new (fixed) bounds for bound-corrected restraints

        if bound_corrected:
            va_settings = self.violation_analyser.getSettings()
            
            if va_settings['lower_bound_correction']['enabled'] == YES:
                new_bound = va_settings['lower_bound_correction']['value']
                [r.setLowerBound(new_bound) for r in bound_corrected]

            if va_settings['upper_bound_correction']['enabled'] == YES:
                new_bound = va_settings['upper_bound_correction']['value']
                [r.setUpperBound(new_bound) for r in bound_corrected]

    def doDumboCalibration(self, peaks):
        
        # BARDIAUX 2.2
        # ConstraintList can bypass calibration
        if is_type(peaks[0], 'DistanceRestraint'):
            do_cal = self.__getListSource(peaks[0])['calibrate'] == 'all_iterations_except_first'
            if do_cal:
                return 1.
            
        peak_sizes = self._get_peak_sizes(peaks)

        ## Assume that average distance of atoms
        ## causing an NOE is 3.0A

        d_calib = 3.0
        
        sum_noe_calc = len(peaks) * (d_calib ** -6)
        
        factor = sum(peak_sizes) / sum_noe_calc
        
        return factor

    def doCalibration(self, restraints, ensemble, store_analysis = 0):

        do_cal = 1
        
        # BARDIAUX 2.2
        # ConstraintList can bypass calibration

        if is_type(restraints[0], 'DistanceRestraint'):
            do_cal = self.__getListSource(restraints[0])['calibrate'] == NO
            if do_cal:
                return 1.
        
        if ensemble is None:            
            factor = self.doDumboCalibration(restraints)
        else:
            f = self.calibrator.calculateEstimator
            factor = f(restraints, ensemble, store_analysis = store_analysis)

            ## if factor is None, i.e. if no NOEs stronger than a
            ## certain cutoff were found, set the cutoff to zero
            ## and calibrate again.

            if factor is None:

                s = self.calibrator.getSettings()
                d_cutoff = s['distance_cutoff']

                s = 'Could not perform 1st calibration, since ' + \
                    'no distances less than %.1f A were found in the ' + \
                    'ensemble. Omitting distance-cutoff and ' + \
                    'calibrating again...'

                self.message(s % d_cutoff)

                factor = f(restraints, ensemble, \
                           store_analysis = store_analysis, use_cutoff = 0)

        return factor

    def finalize_engine(self, molecule):
        from aria.Singleton import ProjectSingleton
    
        project = ProjectSingleton()
        engine = project.getStructureEngine()
        engine.prepare(self.getSettings(), molecule)

    def doAnalysis(self, restraints, ensemble):

        ##
        ## 1st calibration
        ##
        ## Store results.
        ##

        factor = self.doCalibration(restraints, ensemble, store_analysis = 1)

        ##
        ## Violation analysis
        ##
        ## Store results:
        ##
        ## - average violation distance
        ## - whether a restraint is violated
        ## - degree of violation etc.
        ##

        violated, non_violated, new_non_violated = \
                  self.doViolationAnalysis(restraints, ensemble,
                                           store_analysis = 1)

        ## If possible, enlarge set of non-violated restraints

        if new_non_violated:
            non_violated += new_non_violated

        ##
        ## 2nd calibration
        ##

        if non_violated:
            factor = self.doCalibration(non_violated, ensemble)

        ## Calculate model peak-sizes

        self.setModelIntensities(restraints, ensemble, factor)

    def export_water_refined_structures(self, iteration):

        if self.getSettings()['water_refinement']['enabled'] == NO:
            return

        infra = self.getInfrastructure()
        path = infra.get_refinement_path()

        self._dump_ccpn(iteration, path, is_water_refinement=True)
        
    def writeReports(self, iteration):

        t = clock()

        s = 'Performing analysis on calculated structures...'
        self.message(s)

        ensemble = iteration.getStructureEnsemble()

        # BARDIAUX 2.2
        # add CCPN DistanceRestraints
        all_restraints = {}
        all_restraints.update(iteration.getPeaks())
        all_restraints.update(iteration.getDistanceRestraints())

        for spectrum, restraints in all_restraints.items():
        #for spectrum, restraints in iteration.getPeaks().items():
            
            ## Set spectrum-specific settings
            self._updateSpectrumSettings(spectrum)
            self.doAnalysis(restraints, ensemble)

        self.message('Analysis done.')
        self.debug('Time: %ss' % str(clock() - t))

        self.dumpIteration(iteration)

    def readPDBFiles(self, filenames, molecule, float_files = None,
                     format = 'cns'):
        """
        'filenames' is a dict or a tuple. In case of a dict,
        keys are filenames, values are format-strings. In case of
        a tuple, the argument 'format' is used as format-string.
        Reads PDB-files 'filenames' and the respective
        'float_files'. Returns a StructureEnsemble.
        """

        check_type(filenames, TUPLE, DICT)
        check_type(molecule, 'Molecule')
        check_type(float_files, TUPLE, NONE)
        check_string(format)

        import aria.StructureEnsemble as SE

        if float_files == ():
            float_files = None

        if float_files is not None:
            
            from aria.FloatFile import FloatFile

            parser = FloatFile()
            swapped_atoms = [parser.parse(f) for f in float_files]
            
        else:
            swapped_atoms = None

        ## create structure-ensemble from given 'filenames'

        se_settings = SE.StructureEnsembleSettings()
        ensemble = SE.StructureEnsemble(se_settings)
        ensemble.read(filenames, molecule, format, swapped_atoms)

        return ensemble

    def start(self, iteration, molecule):
        """
        setup the first iteration. calculate distance-estimates
        using either a template structure or by doing a dumbo-
        calibration.
        """

        check_type(molecule, 'Molecule')
        check_type(iteration, 'Iteration')
        
        ## set-up first iteration

        import os
        from aria.Singleton import ProjectSingleton
        from aria.DataContainer import DATA_TEMPLATE_STRUCTURE
        
        t_iteration = clock()

        ## create seed-assignment

        ## if a template-structure has been specified,
        ## use that structure to create a better seed-assignment.
        ## otherwise we use an extended chain (in case of CNS,
        ## it is created by using generate_template.inp

        infra = self.getInfrastructure()
        project = ProjectSingleton()

        templates = project.getData(DATA_TEMPLATE_STRUCTURE)
        templates = [t for t in templates if t['enabled'] == YES]
        
        ## If template structures have been specified, we
        ## create an initial structure ensemble which will
        ## be used in the following steps to calculate
        ## an initial calibration factor etc.

        if templates:

            filenames = [t['filename'] for t in templates]
            formats = [t['format'] for t in templates]
            
            s = 'Template structure(s) specified (%s). Using template(s)' + \
                ' for calibration / seed-assignment generation.'

            self.message(s % ', '.join(map(os.path.basename, filenames)))

            ## Compile local filenames (template PDB-files are
            ## read from local data-directory)

            filenames = [infra.get_local_filename(f, DATA_TEMPLATE_STRUCTURE) \
                         for f in filenames]

            ## Build dict with filenames as keys and formats as
            ## values.

            d = {}

            for i in range(len(formats)):
                d[filenames[i]] = formats[i]

            ensemble = self.readPDBFiles(d, molecule)
            ensemble.getSettings()['number_of_best_structures'] = 'all'
            iteration.setStructureEnsemble(ensemble)

        self.message('Initial iteration created.', verbose_level = VL_LOW)
        self.debug('Time: %ss' % str(clock() -t_iteration))

        ## register initial iteration and start protocol
        
##         self.iterations[iteration.getNumber()] = iteration

        ## Return last iteration
        return self.run_protocol(molecule, iteration)
    ## < Mareuil
    def go(self, iteration, molecule, MODE = None):
        self.setMODE(MODE)
        ## Mareuil >
        check_type(molecule, 'Molecule')
        check_type(iteration, 'Iteration')

        ## Check which iterations have already been
        ## calculated, so that we can skip them.

        first_iteration = self.findFirstIteration()

        ## If no iteration has been calculated yet, run the
        ## full protocol.
        
        if first_iteration == 0:
            iteration = self.start(iteration, molecule)

        ## If some iterations already exist, calculate the remaing ones.
            
        elif first_iteration is not None:

            t = clock()
            
            s = 'Existing iterations found. Resuming with iteration ' + \
                '%d.'
            
            self.message(s % first_iteration)

            iteration._setNumber(first_iteration - 1)

            ## Get filenames of structures which ought exist
            ## after parent-iteration has been successfully completed.

            CFL = self.compileFileList
            
            filenames = CFL(first_iteration - 1, structures = 1)
            float_files = CFL(first_iteration - 1, float_files = 1)

            it_settings = self.getIterationSettings(iteration.getNumber())

            ## Default PDB-format is 'cns'
            
            ensemble = self.readPDBFiles(filenames, molecule,
                                         float_files)


            ## Apply structure-ensemble settings of current
            ## iteration.

            se_settings = ensemble.getSettings()
            se_settings.update(it_settings)

            ensemble.settingsChanged()
            
            iteration.setStructureEnsemble(ensemble)

            s = 'PDB-files / FLOAT-files for iteration %d read: %ss'
            self.debug(s % (first_iteration, str(clock() - t)))

            iteration = self.run_protocol(molecule, iteration)
            
        else:

            from aria.Iteration import Iteration
            
            self.message('Iterations are already complete.')

            ## TODO: this is a bit fiddled!
            ## Since the solvent-refinement step does not refer to any
            ## restraint, we simply create an empty iteration
            ## and set its number to the last iteration's.

            infra = self.getInfrastructure()
            last_iteration = infra.getSettings()['n_iterations'] - 1
            iteration = Iteration(last_iteration)
            
        ## Start water refinement
        self.doSolventRefinement(iteration, molecule)

        ## Export structures to CCPN

        self.export_water_refined_structures(iteration)

        ## End of ARIA protocol.
        self.done()

        ## return last iteration

        return iteration

    ## BARDIAUX 2.2
    def _checkSpectrumData(self, spectrum):

        k = spectrum.getExperimentData().values()
        return len([v for v in k if v <> 0.0]) == len(k)
    
    def run_iteration(self, molecule, last_iteration, iteration,
                      is_first_iteration = 0):
        """
        last_iteration is introduced, to provide time-shifted
        pickling of restraint lists (i.e. while another
        iteration is on its way)
        """

        from time import sleep
        
        check_type(molecule, 'Molecule')
        check_type(last_iteration, 'Iteration')

        ## if all iterations shall be stored in memory ...

##         self.iterations[iteration.getNumber()] = iteration

        ## apply settings of target iteration to all sub-modules
        self._updateIterationSettings(iteration)

        ## get structure ensemble calculated in last iteration
        ensemble = last_iteration.getStructureEnsemble()

        ## We now analyse the structure ensemble of the last
        ## iteration to generate a new set of restraints which
        ## will then be used to calculate a new ensemble of
        ## structures for the target iteration.

        spectra = last_iteration.getPeaks().keys()
        spectra.sort()

        # BARDIAUX 2.2
        # add CCPN DistanceRestraints
        constraint_lists = last_iteration.getDistanceRestraints().keys()
        constraint_lists.sort()

        # Malliavin/Bardiaux rMat
        cs = self.calibrator.getSettings()

        for s in spectra:
            if cs['relaxation_matrix'] <> NO and not self._checkSpectrumData(s):
                self.warning("Missing experimental data for spectrum \"%s\". "  % s.getName())
                self.warning("Spin-diffusion correction will be disabled.\n")

                cs['relaxation_matrix'] = NO
                #continue
            
        if cs['relaxation_matrix'] == YES and ensemble is not None:
            
            from aria.NOEModel import SpinDiffusionCorrection
            
            self.model = SpinDiffusionCorrection()
            self.calibrator.setModel(self.model)
            self.model.prepare(molecule, ensemble)

            # initialize matrix
            for spectrum in spectra:
                self.model.setIntensityMatrix(spectrum)

        spectra += constraint_lists

        for spectrum in spectra:

            self.message('Calibrating spectrum "%s"...' % spectrum.getName())

            ## get peaks stored in current iteration
            ## BARDIAUX 2.2 ConstraintList 
            if is_type(spectrum, 'ConstraintList'):
                peaks = iteration.getDistanceRestraints(spectrum)
            else:
                peaks = iteration.getPeaks(spectrum)

            ##
            ## Calibration
            ##
            ## If we do not have a seed structure ensemble, we
            ## perform a very simple calibration
            ##

            ## set spectrum-specific parametes for Calibrator
            
            
            self._updateSpectrumSettings(spectrum)

            if ensemble is None:
                factor = self.doCalibration(peaks, None)
                new_non_violated = None

            else:

                ## Calculate initial calibraton factor.
            
                factor = self.doCalibration(peaks, ensemble)

                ## Calculate upper/lower bounds for restraints of current
                ## iteration. 

                t = clock()

                self.calculateBounds(factor, peaks, ensemble = ensemble)

                s = '1st calibration and calculation of new ' + \
                    'distance-bounds done (calibration factor: %e)'
                self.message(s % factor, verbose_level = VL_LOW)
                self.debug('Time: %ss' % str(clock() - t))
                
                ##
                ## Violation Analysis
                ##
                ## Assess every restraint regarding its degree of violation.
                ##
                ## Violated restraints will be disabled for the
                ## current iteration and thus will not be used during
                ## structure calculation.
                ##

                t = clock()

                violated, non_violated, new_non_violated = \
                          self.doViolationAnalysis(peaks, ensemble)

                ## Augment set of non-violated restraints

                if new_non_violated:
                    non_violated += new_non_violated

                n = len(peaks)
                n_viol = len(violated)
                p_viol = n_viol * 100. / n

                s = 'Violation analysis done: %d / %d restraints ' + \
                    '(%.1f %%) violated.'

                self.message(s % (n_viol, n, p_viol))

                if new_non_violated:
                    s = 'Number of valid restraints has been increased ' + \
                        'by %d (%.1f%%) after applying a bound-correction.'

                    p_new = len(new_non_violated) * 100. / n

                    self.message(s % (len(new_non_violated), p_new))

                self.debug('Time: %ss' % str(clock()-t))

                ##
                ## 2nd calibration - wrt to non-violated restraints.
                ## If no restraints have been violated, we use the
                ## 1st calibration factor.
                ## Again, we do not store results.
                ##

                if non_violated:
                        factor = self.doCalibration(non_violated, ensemble)

                ##
                ## Activate restraints explicitly.
                ## We consider a restraint as active, if it has
                ## not been violated or if its reference cross-peak is
                ## 'reliable'.
                ##

                for r in peaks:
                    if not r.analysis.isViolated() or \
                           r.getReferencePeak().isReliable():
                        r.isActive(1)
                    else:
                        r.isActive(0)

            ## Store final calibration factor for current iteration.
            iteration.setCalibrationFactor(spectrum, factor)

            ## Calculate upper/lower bounds for restraint-list
            ## used in the current iteration. I.e. these bounds
            ## will be used to calculated the structures.

            t = clock()

            self.calculateBounds(factor, peaks, new_non_violated, ensemble = ensemble)

            s = 'Final calibration and calculation of new distance-bounds' + \
                ' done (calibration factor: %e).' % factor
            self.message(s, verbose_level = VL_LOW)
            self.debug('Time: %ss' % str(clock() - t))

            ##
            ## Partial assignment for restraint-list used in
            ## current iteration, i.e. for all restraints:
            ## Calculate weight for every contribution
            ## and (depends possibly on partial analyser
            ## settings) throw away 'unlikely' contributions.
            ##
            ## If we do not have an ensemble, all contributions
            ## are activated.
            ##

            t = clock()
            
            # BARDIAUX 2.2
            # if we are using the NA, we do not filter the contributions yet
            # if the spectrum is a ConstraintList, check if needs to be filtered
            NA = self.network_anchoring
            ns = NA.getSettings()
            do_filter = ns['enabled'] == NO
            
            if is_type(spectrum, 'ConstraintList'):
                do_filter = do_filter or \
                            spectrum.getListSource()['run_network_anchoring'] == NO
                
                do_filter = do_filter and \
                            spectrum.getListSource()['filter_contributions'] == YES
                        
            self.contribution_assigner.assign(peaks, ensemble, filter_contributions = do_filter)
            
            self.message('Partial assignment done.', verbose_level = VL_LOW)
            self.debug('Time: %ss' % str(clock() - t))
            
        # BARDIAUX 2.2 NA
        NA = self.network_anchoring
        ns = NA.getSettings()
        
        if ns['enabled'] == YES:
            NA.run(iteration)
            
            # filter contributions
            for spectrum in spectra:

                # CCPN ConstraintList may not be filtered
                if is_type(spectrum, 'ConstraintList'):
                    if spectrum.getListSource()['filter_contributions'] == NO:
                        continue
                    peaks = iteration.getDistanceRestraints(spectrum)
                else:
                    peaks = iteration.getPeaks(spectrum)
                    
                self.contribution_assigner.filter_contributions(peaks)

                
        ##
        ## Merge spectra
        ##
                
        self.mergePeakLists(iteration)
        
        ## BARDIAUX
        # test maxn remove peak
        maxn = self.contribution_assigner.getSettings()['max_contributions']
        for p in iteration.getPeakList():
            if len(p.getActiveContributions()) > maxn and \
               not p.getReferencePeak().isReliable():
                p.isActive(0) 
        
        ## to known when the structure generation has been completed.
        self.__done = 0

        ## Display some cache statistics 
        if AriaBaseClass.cache and ensemble is not None:
            cache = ensemble._StructureEnsemble__cache
            hit_rate = float(cache['hit']) * 100. / cache['total']
            self.debug('StructureEnsemble cache hit-rate: %.1f%%' % hit_rate)
        
        ##
        ## Run structure calculation
        ##

        # test BARDIAUX 2.2 kept_structures
        kept_structures = []
        if ensemble is not None:
            kept_structures = ensemble.getFiles()
        self.startStructureCalculation(iteration, molecule, kept_structures = kept_structures)

        ##
        ## Perform analysis wrt to last iteration
        ## Dump report files
        ##

        if not is_first_iteration:
            
            ## apply settings of last iteration to sub-modules
            self._updateIterationSettings(last_iteration)

            self.writeReports(last_iteration)
            
            ## Restore old settings
            self._updateIterationSettings(iteration)

        ## Wait for completion of structure generation.
        self.message('Waiting for completion of structure calculation...')

        try:
        
            while not self.__done:
                sleep(1.)

        except KeyboardInterrupt:
            self.__done = KEYBOARD_INTERRUPT
            
        if self.__done in (MISSING_STRUCTURES, KEYBOARD_INTERRUPT):

            ## shut down job-manager

            from aria.Singleton import ProjectSingleton

            project = ProjectSingleton()
            engine = project.getStructureEngine()
            engine.getJobScheduler().shutdown()

            if self.__done == MISSING_STRUCTURES:

                ## BARDIAUX 2.2
                # get list of missing structures
                missing  =  engine.missingStructures()
                msg = 'Structure calculation failed for structure %d.\n' #Please check your setup and the CNS output files for errors.'
                err_msg = ''
                for m in missing:
                    err_msg += msg % m

                err_msg += 'Please check your setup and the CNS output files for errors.'

##                 err_msg = 'Structure calculation failed. ' + \
##                           'Some structures are missing. ' + \
##                           'Please check your setup and/or CNS output files.'

            else:
                err_msg = 'Interrupted by user.'
                
            self.error(StandardError, err_msg)
        
        ## Return current, most recently calculated iteration
        return iteration

    def run_protocol(self, molecule, iteration):

        check_type(molecule, 'Molecule')
        check_type(iteration, 'Iteration')

        from copy import copy
        from aria.Singleton import ProjectSingleton
        
        infra = self.getInfrastructure()
        project = ProjectSingleton()
        engine = project.getStructureEngine()

        msg = '---------------------- Iteration %d -----------------------\n'

        is_first_it = 1

        n_iterations = len(self.getSettings()['iteration_settings'])

        t_main_protocol = clock()

        n_first = iteration.getNumber()

        for it_number in range(n_first, n_iterations - 1):

            ## create new iteration

            t = clock()

            target = copy(iteration)
            target._setNumber(iteration.getNumber() + 1)

            self.debug('New iteration created: %ss' % str(clock() - t))

            self.message(msg % (target.getNumber()))
            
            ## Run calculation

            t = clock()
            
            iteration = self.run_iteration(molecule, iteration, target,
                                           is_first_it)

            ## cleanup

            if it_number < n_iterations - 1:
                it_path = infra.get_iteration_path(it_number)
                engine.cleanup(it_path)
            
            is_first_it = 0

            self.message('Iteration %s done.' % target.getNumber())
               
            self.debug('Time: %ss' % str(clock() - t))

        self.debug('Total time for all iterations: %ss' % \
                   str(clock() - t_main_protocol))

        ## Write reports for final iteration
        self.writeReports(iteration)
        self.writeFinalReports(iteration)

        ## Return last iteration
        return iteration

    def writeFinalReports(self, iteration):

        check_type(iteration, 'Iteration')
        
        import os
        import aria.MolMol as MolMol
        from aria.Singleton import ProjectSingleton

        ## write MOLMOL file-list "file.nam"

        molmol_path = self._writePDBFileList(iteration)

        if molmol_path is None:
            msg = 'Could not write MOLMOL file-list.'
            self.warning(msg)
            return

        ## if that worked, attempt to write restraint-lists

        project = ProjectSingleton()
        s = project.getReporter()['molmol']

        if s['enabled'] in (YES, GZIP):
            
            ## restraint files

            restraints = iteration.getPeakList()
            not_merged = [r for r in restraints if not r.isMerged()]
           
            gzip = {YES: 0, GZIP: 1}[s['enabled']]

            MolMol.write_noe_restraints(not_merged, molmol_path,
                                        gzip = gzip)
            self.message('MOLMOL lower and upper bound (.lol, .upl) ' + \
                         'files written.')

class ProtocolXMLPickler(XMLBasePickler):

    order = ['floating_assignment', 'iteration',
             'water_refinement']
    
    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)

        s = x.getSettings()

        e.iteration = tuple(s['iteration_settings'].values())
        e.water_refinement = s['water_refinement']
        e.floating_assignment = s['floating_assignment']

        return e

    def load_from_element(self, e):

        from aria.tools import as_tuple

        s = ProtocolSettings()

        it_settings = as_tuple(e.iteration)

        [s.addIterationSettings(i) for i in it_settings]
        s['water_refinement'] = e.water_refinement
        s['floating_assignment'] = str(e.floating_assignment)

        d = s['iteration_settings']

        it_numbers = d.keys()

        if min(it_numbers) <> 0:

            self.warning('Iteration numbering does not start at 0, renumbering ...')

            _min = min(it_numbers)

            for number, it_settings in d.items():

                new_number = number - _min

                it_settings['number'] = new_number
                d[new_number] = it_settings

                del d[number]

        p = Protocol(s)

        return p

IterationSettings._xml_state = IterationSettingsXMLPickler()._xml_state
Protocol._xml_state = ProtocolXMLPickler()._xml_state
