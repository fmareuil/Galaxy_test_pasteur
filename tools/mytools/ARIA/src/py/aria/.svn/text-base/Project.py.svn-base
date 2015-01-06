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
from aria.Singleton import Singleton
from threading import Thread
from aria.xmlutils import XMLElement, XMLBasePickler

QUALITY_CHECKS_HEADER = \
"""

Summary for quality checks:
---------------------------
Created: %(time)s

"""

def checksum(filename):

    import md5

    f = open(filename)
    l = f.read()
    f.close()

    m = md5.new()
    m.update(l)

    return m.digest()

class ProjectSettings(Settings):

    def create(self):

        from aria.Settings import Path, String, PositiveFloat, AbsolutePath
        from aria.Settings import PositiveInteger, YesNoChoice
        import aria.Settings as Settings

        d = {'topology_xml': Path(),
             'name': String(),
             'version': PositiveFloat(),
             'author': String(),
             'date': String(),
             'description': String(),
             'comment': String(),
             'references': String()}

        descr = 'Defines the root directory of your project. Every run is then stored in a separate subdirectory PROJECT_PATH/RUNxxx.'
             
        d['working_directory'] = AbsolutePath(description = descr, exists = 0)
        d['path_name_iterations'] = Path(exists = 0)

        descr = 'For every run, ARIA creates a directory [TEMP_PATH]/aria_temp.@xxxxx to store all temporary files (e.g. CNS output). If ommited, the current directory will be used. When calculating on multiple machines, the temporary directory *must* be accessible from all machines.'
             
        d['temp_root'] = AbsolutePath(description = descr, exists = 0)

        d['aria_temp_path'] = AbsolutePath(exists = 0)
        
        d['n_iterations'] = PositiveInteger()

        descr = "Specifies the current RUN of a project. Every RUN is stored in a separate directory with the project directory as its root, i.e. PROJECT_PATH/RUNxxxx"
        
        d['run'] = String(description = descr)

        descr = 'The file_root is used to build filenames. E.g., if the file_root is set to "my_structure", the structure PDB files are called my_structure_1.pdb etc.'
        
        d['file_root'] = Settings.NonEmptyString(description = descr,
                                                 mandatory = 1)
        
        descr = "If enabled, ARIA caches all XML data files (i.e. molecule definition and spectra). If a RUN is executed more than once, ARIA just accesses the cache instead of reading the complete XML files; this considerably speeds-up the reading process. If the original data (i.e. those that are stored in the run's local data directory PROJECT_PATH/RUNxxx/data) are modified, the cache is automatically invalidated so that the data are read from their original XML files."
             
        d['cache'] =YesNoChoice(description = descr)

        descr = 'Enable this option if you want ARIA to remove the temporary directory TEMP_PATH/aria.xxxxx. Of course, the main temporary directory, TEMP_PATH, is not removed.'
        
        d['cleanup'] = YesNoChoice(description = descr)

        return d

    def create_default_values(self):

        from aria.ariabase import YES, NO

        d = {}

        d['version'] = 1.0
        d['cache'] = NO
        d['run'] = str(1)
        d['cleanup'] = YES

        return d

class Project(AriaBaseClass):

    def __init__(self, settings):

        check_type(settings, 'ProjectSettings')

        AriaBaseClass.__init__(self, name = 'Project')
        
        self.setSettings(settings)
        self.setDefaultValues()

    def setDefaultValues(self):
        self.__infra = None
        self.__data = []
        self.__protocol = None
        self.__analysis_parameters = None
        self.__engine = None
        self.__reporter = None
        self.__analyser = None
        
        self.chemical_shift_list_filter = None
        self.noesy_spectrum_filter = None
        self.peak_assigner = None
        self.molecule = None
        self.ccpn_model = None
        self.ccpn_project_instance = None
        ## BARDIAUX 2.2
        self.__ref_segid = "    "

        self.distance_restraints = []

    def initialize(self):
        """
        must be called after unpickling a project to finish
        the initialization.
        """

        ## infrastructure settings

        infra = self.getInfrastructure()
        infra_settings = infra.getSettings()

        protocol = self.getProtocol()

        n_iterations = len(protocol.getSettings()['iteration_settings'])
        infra_settings['n_iterations'] = n_iterations

        ## TODO: hard-coded
        
        infra_settings['path_name_iterations'] = 'structures/it'

        self.createObjects()

    def setInfrastructure(self, infra):
        check_type(infra, 'Infrastructure')
        self.__infra = infra

    def getInfrastructure(self):
        return self.__infra

    def setStructureEngine(self, engine):
        check_type(engine, 'StructureEngine')
        self.__engine = engine

    def getStructureEngine(self):
        return self.__engine

    def setReporter(self, r):
        check_type(r, 'ReportSettings')
        self.__reporter = r

    def getReporter(self):
        return self.__reporter

    def addData(self, data):
        check_type(data, 'DataContainer')
        self.__data.append(data)

    def delData(self, data):
        check_type(data, 'DataContainer')
        
        if data not in self.__data:
            s = 'Data "%s" does not exist.' % str(data)
            self.error(ValueError, s)

        self.__data.remove(data)

    def getData(self, type):

        from aria.DataContainer import DATA_TYPES
        
        if not type in DATA_TYPES:
            self.error(TypeError, 'Data type "%s" not known.' % str(type))
            
        return filter(lambda d, t = type: d.getType() == t, self.__data)

    def setProtocol(self, p):
        check_type(p, 'Protocol')
        self.__protocol = p

    def getProtocol(self):
        return self.__protocol

    def setAnalyser(self, a):
        check_type(a, 'Analyser')
        self.__analyser = a

    def getAnalyser(self):
        return self.__analyser
    
    def getChemicalShiftListFilter(self):
        return self.chemical_shift_list_filter

    def getNOESYSpectrumFilter(self):
        return self.noesy_spectrum_filter

    def getPeakAssigner(self):
        return self.peak_assigner

    def getMolecule(self):
        return self.molecule

    def getExperiments(self):
        return self.experiments

    def createObjects(self):

        import aria.ShiftAssignmentFilter as SA
        import aria.NOESYSpectrumFilter as NOESY
        import aria.PeakAssigner as PA
        
        ## TODO: parameters should come from project.xml

        filter = SA.ChemicalShiftListFilter(\
            SA.ChemicalShiftListFilterSettings())
        filter.getSettings()['ppm_type'] = (FLOAT, NONE)
        filter.getSettings()['max_no_shifts'] = 2

        self.chemical_shift_list_filter = filter
        
        assigner_settings = PA.PeakAssignerSettings()
        self.peak_assigner = PA.PeakAssigner(assigner_settings)

        noesy_settings = NOESY.NOESYSpectrumFilterSettings(default_settings = \
                                                           assigner_settings)
        
        self.noesy_spectrum_filter = NOESY.NOESYSpectrumFilter(noesy_settings)

    def setup(self, force = 0):
        """
        Setup infrastructure (data directories, copy data locally), setup
        CNS environment. If force is nonzero, the project will be
        re-setup.
        """

        self.getInfrastructure().setup(force)

    def check(self, probe_commands = 1):
        """
        has to be called prior to project.go(). 
        """

        import os, sys
        
        self.getInfrastructure().finalize()

        ## Check whether run-path exists
        ## as a simple test whether project
        ## has been setup.

        run_path = self.getInfrastructure().get_run_path()

        if not os.path.exists(run_path):
            run = self.getSettings()['run']
            s = 'The directory for run "%s" (%s) does not exist. ' + \
                'Has the project been set up?'
            self.message(s % (run, run_path))
            sys.exit(1)

        ## Finalize job-manager

        engine = self.getStructureEngine()
        job_manager = engine.getJobScheduler()
        jm_settings = job_manager.getSettings()

        ## In case of an empty host-list, use localhost

        if not jm_settings['host_list']:

            from aria.JobManager import HostSettings
            
            try:
                pipe = os.popen('hostname')
                localhost = pipe.read().strip()
                pipe.close()
            except:
                s = """List of available machines is empty and querying the name of the local host has failed. Please specify a host list explicitly (<structure_generation> section in project-xml)."""
                self.error(StandardError, s)

            hs = HostSettings()
            hs.reset()

            ## use job-managers default command
            ## for job-dispatch
            
            default_cmd = jm_settings['default_command'].strip()

            ## just in case
            
            if default_cmd.strip() == '':
                s = 'Job manager\'s <default_command> has not been set. Using "csh -f".'
                self.message(s)
                default_cmd = 'csh -f'
                
            elif default_cmd.lower() in ('ssh', 'rsh'):
                s = 'Job manager\'s <default_command> is set to "%s" which makes little sense without additional options. Using "csh -f" instead.'
                self.message(s % default_cmd)
                default_cmd = 'csh -f'
            ## < Mareuil
            hs['submit_command'] = default_cmd
            ## Mareuil >
            ## Get default binary from structure engine

            engine = self.getStructureEngine()
            executable = engine.getSettings()['local_executable']

            ## Check whether local executable exists

            if not os.path.exists(executable):
                s = 'CNS executable not found: %s'
                self.error(IOError, s % executable)
            
            hs['executable'] = executable

            jm_settings['host_list'] = [hs]

            s = 'List of available machines is empty. Running ' + \
                'calculations on localhost "%s", assuming 1 CPU. ' + \
                'Using %s as structure engine.'
            
            self.message(s % (localhost, executable))

        ## check whether specified commands are valid
        ## < Mareuil
        if probe_commands and jm_settings['default_command'] != "GRID":
        ## Mareuil >
            self.message('Checking host list ...')

            passed = job_manager.probeHosts(self.getInfrastructure())

            if not passed:
                msg = 'At least one command could not be executed. Please check your host list setup.'
                self.warning(msg)
                self.halt()

        else:
            ## < Mareuil
            if jm_settings['default_command'] == "GRID":
                self.message('No host list check in GRID mode')
            else:
                self.message('Host list check has been disabled.')
            ## Mareuil >
	    
    def cache_data(self, molecule, experiments):

        from aria.tools import Dump
        from aria.Singleton import AtomFactory

        ## set-up check-sum for sequence source file

        filename = molecule.getDataSource()['filename']

        try:
            cs = checksum(filename)
        except:
            s = 'Could not generate checksum for file %s.'            
            self.error(IOError, s % filename)
        
        cache = {}

        cache['atom_factory'] = AtomFactory()
        cache['molecule'] = {molecule.getDataSource()['filename']: \
                             {'object': molecule, 'checksum': cs}}

        cache['experiments'] = {}

        for e in experiments:

            shifts_data = e.getDataSource()['shifts']
            peaks_data = e.getDataSource()['peaks']

            ## set-up checksums for shift- and spectrum source files

            shifts_cs = checksum(shifts_data['filename'])
            peaks_cs = checksum(peaks_data['filename'])

            d = {shifts_data['filename']: {'object': e.getShiftList(),
                                           'checksum': shifts_cs},
                 peaks_data['filename']: {'object': e.getSpectrum(),
                                          'checksum': peaks_cs}}

            cache['experiments'].update(d)

        cache_filename = self.getInfrastructure().get_cache_filename()

        try:
            Dump(cache, cache_filename)
        except:
            s = 'Could not save cached data files ("%s").'
            self.error(IOError, s % cache_filename)

    def check_cached_molecule(self, cache):
        from aria.DataContainer import DATA_SEQUENCE

        cache = cache['molecule']

        ## TODO: [0] ...

        molecule_data = self.getData(DATA_SEQUENCE)[0]
        filename, format = molecule_data.getLocation()

        ## Check whether cache knows about sequence

        if not filename in cache:
            return 0

        check_sum = checksum(filename)

        return check_sum == cache[filename]['checksum']

    def check_cached_spectra(self, cache):
        from aria.DataContainer import DATA_SPECTRUM

        cache = cache['experiments']
        spectrum_data = self.getData(DATA_SPECTRUM)

        shift_files = [sd['shifts'].getLocation()[0] for \
                        sd in spectrum_data]

        spectrum_files = [sd['peaks'].getLocation()[0] for \
                          sd in spectrum_data]

        local_filenames = shift_files + spectrum_files

        local_filenames.sort()

        ## If list of experiment files does not match list
        ## of cached files, reload the whole thing.

        names_cache = cache.keys()
        names_cache.sort()

        ok = local_filenames == names_cache

        ## if both lists are equal, check checksums

        if ok:
            for fn in local_filenames:
                c = checksum(fn)
                if cache[fn]['checksum'] <> c:
                    ok = 0
                    break
                
        return ok

    def read_molecule(self, cache = None):
        """
        Molecule data are read from its local directory.
        """

        from aria.DataContainer import DATA_SEQUENCE

        if DATA_SEQUENCE in self.ccpn_data_sources:
            return

        molecule_data = self.getData(DATA_SEQUENCE)[0]
        filename, format = molecule_data.getLocation()

        s = 'Reading molecule definition %s.'
        self.message(s % filename)

        if cache is not None:
            molecule = cache['molecule'][filename]['object']

        if cache is None:
            
            import aria.AriaXML as AriaXML
            
            pickler = AriaXML.AriaXMLPickler()
            molecule = pickler.load(filename)

        molecule.setDataSource(molecule_data)
        
        self.molecule = molecule

    def REREAD(self, x):

        from aria.AriaXML import AriaXMLPickler
        from time import clock
        import tempfile, os

        tmp = tempfile.mktemp()

        filename = os.path.split(tmp)[1]

        path = self.getInfrastructure().get_temp_path()

        tmp = os.path.join(path, filename)

        pickler = AriaXMLPickler()

        t0 = clock()

        s = pickler.dump(x, tmp)

        obj = pickler.load(tmp)

        self.debug('REREAD: %ss' % str(clock() - t0))

        return obj

    def read_spectra(self, cache = None):
        """
        Spectrum data are read from their local directories.
        """
        from aria.Experiment import Experiment
        from aria.DataContainer import DATA_SPECTRUM

        spectrum_data = self.getData(DATA_SPECTRUM)

        experiments = []
        
        if cache is None:
            import aria.AriaXML as AriaXML
            pickler = AriaXML.AriaXMLPickler()

            spec_id = 0
            spectrum_names = []

        ccpn_spectra = [x for x in self.ccpn_data_sources.get(DATA_SPECTRUM, ())]

        for sd in spectrum_data:

            if sd in ccpn_spectra:
                continue
            
            ## Forget the format, must be XML anyway

            shifts, format = sd['shifts'].getLocation()
            spec, format = sd['peaks'].getLocation()

            if sd['enabled'] <> YES:
                m = 'Spectrum "%s" disabled.'
                self.message(m % spec)
                continue

            s = 'Reading spectrum %s.'
            self.message(s % spec)
            s = 'Reading chemical shift list %s.'
            self.message(s % shifts)

            if cache is not None:
                assignment = cache['experiments'][shifts]['object']
                spectrum = cache['experiments'][spec]['object']

            if cache is None:

                assignment = pickler.load(shifts)
                spectrum = pickler.load(spec)

                if spectrum.getName() is None:
                    spectrum._setName('spec_%d' % spec_id)
                    s = 'Spectrum %s has no name. Name set to "%s".'
                    self.warning(s % (spec, spectrum.getName()))
                    spec_id += 1

                ## Ensure that spectra have unique names

                spec_name = spectrum.getName()
                
                if spec_name in spectrum_names:
                    s = 'Spectra must have unique names.'
                    self.error(StandardError, s)
                else:
                    spectrum_names.append(spec_name)
                    
            E = Experiment(spectrum, assignment)
            E.setDataSource(sd)
            experiments.append(E)

        self.experiments += experiments

    def read_ccpn_data(self):
        
        from memops.general.Io import loadProject
        from aria.Experiment import Experiment
        import aria.DataContainer as DC
        import aria.importFromCcpn as ccpnImport

        msg = '-' * 12 + ' Reading data from CCPN data model ' + '-' * 12
        self.message(msg)

        ccpn_fn = self.ccpn_model['filename']

        self.message('Opening CCPN data model "%s"...' % ccpn_fn)
        
        ccpn_project = loadProject(ccpn_fn)

        self.ccpn_project_instance = ccpn_project
 
        self.message('CCPN data model opened.')

        X = self.ccpn_data_sources

        if DC.DATA_SEQUENCE in X:
            
            seq_data = X[DC.DATA_SEQUENCE][0]

            molsystem_id = tuple(ccpnImport.getKeysFromString(seq_data['ccpn_id']))

            codes = ",".join(molsystem_id[1:])
            self.message('Retrieving molecular system "%s", chain code(s) "%s"...' % (molsystem_id[0], codes))
            
        else:
            molsystem_id = None
            
        ## BARDIAUX DistanceConstraint
        restraints_names = {}
        restraints_types = [DC.DATA_AMBIGUOUS,DC.DATA_UNAMBIGUOUS, DC.DATA_HBONDS, DC.DATA_SSBONDS,
                            DC.DATA_DIHEDRALS, DC.DATA_KARPLUS, DC.DATA_RDCS]
    
        for d in restraints_types:
            if d in X:
                restraints_names[d] = []
                for c in X[d]:
                    if c['enabled'] == YES:
                        restraints_names[d].append(c['ccpn_id'])
                        self.message('Retrieving Constraint List "%s" of type "%s" ...' % (c['ccpn_id'], d))
        ##
                        
        ## BARDIAUX Models
        initial_structure_names = []
        if DC.DATA_INITIAL_STRUCTURE in X:
          for init_struc_data in X[DC.DATA_INITIAL_STRUCTURE]:
            if init_struc_data['enabled'] == YES:
              initial_structure_names.append(init_struc_data['ccpn_id'])
          if initial_structure_names:
              self.message('Retrieving initial structure(s) "%s" ...' % ", ".join(initial_structure_names))

        template_structure_names = []
        if DC.DATA_TEMPLATE_STRUCTURE in X:
          for template_struc_data in X[DC.DATA_TEMPLATE_STRUCTURE]:
            if template_struc_data['enabled'] == YES:
              template_structure_names.append(template_struc_data['ccpn_id'])
          if template_structure_names:
              self.message('Retrieving template structure(s) "%s" ...' % ", ".join(template_structure_names))
        ##
              
        exp_names = []

        if DC.DATA_SPECTRUM in X:

            for spec_data in X[DC.DATA_SPECTRUM]:

                if spec_data['enabled'] <> YES:
                    #m = 'Spectrum "%s" disabled.'
                    #self.message(m % spec)
                    continue
            
                peakListString  = spec_data['peaks']['ccpn_id']
                shiftListString = spec_data['shifts']['ccpn_id']
            
                shiftListKeys = ccpnImport.getKeysFromString(shiftListString)
                peakListKeys  = ccpnImport.getKeysFromString(peakListString)
                
                exp_names.append((peakListKeys,shiftListKeys))

                pText = '(project=%s, experiment=%d, data_source=%d, serial=%d)' %  tuple(peakListKeys)

                self.message('Retrieving peak list "%s", shift list "%d" ...' % \
                            (pText, shiftListKeys[-1]))

        if exp_names or molsystem_id:
            
            ccpChains         = ccpnImport.getCcpnChains(ccpn_project, molsystem_id)#[0]

            constraint_lists = ccpnImport.getCcpnConstraintLists(ccpn_project, restraints_names)

            ccpn_initial_models = ccpnImport.getCcpnModels(ccpn_project, initial_structure_names)

            ccpn_template_models = ccpnImport.getCcpnModels(ccpn_project, template_structure_names)
            
            from aria.exportToCcpn import getAriaRun
            # TJS: Get Run object for CCPN to group objects
            # Run object may come from Extend-NMR GUI or a
            # previously intialised project
            # An old run is only used if it has no output
            if constraint_lists and constraint_lists.values()[0]:
              constraintSet = constraint_lists.values()[0][0].nmrConstraintStore
            else:
              constraintSet = None
            
            # TJS: This will be None if the CCPN version is pre 2.1  
            ccpnAriaRun = getAriaRun(ccpChains[0].molSystem)
             
            peak_shift_lists = ccpnImport.getCcpnPeakAndShiftLists(ccpn_project, ccpChains[0].molSystem, exp_names)

            if molsystem_id:
                self.molecule = ccpnImport.makeAriaMolecule(ccpChains[0].molSystem, ccpChains)

            from aria.Singleton import AtomFactory

            atom_factory = AtomFactory(__new_instance__ = 1)

            ## conversion ... sorry for that
            
            self.molecule = self.REREAD(self.molecule)

            ## BARDIAUX 2.2
            self.map_chains()
            
            atom_factory.freeze()

            ## set up peakList/shiftList pairs

            experiments = []


            for peakList, shiftList in peak_shift_lists:
                
                # TJS: Store CCPN objects as input to run if absent
                # This is a record of what was actually used.
                if ccpnAriaRun:
                  dataObj = ccpnAriaRun.findFirstData(className='PeakListData',
                                                      peakList=peakList,
                                                      ioRole='input')
                
                  if not dataObj:
                    spectrum = peakList.dataSource
                    eSerial = spectrum.experiment.serial
                    sSerial = spectrum.serial
                    pSerial = peakList.serial
 
                    ccpnAriaRun.newPeakListData(experimentSerial=eSerial,
                                                dataSourceSerial=sSerial,
                                                peakListSerial=pSerial,
                                                ioRole='input')

                aria_shift_list = ccpnImport.makeAriaShiftList(shiftList, ccpChains, self.molecule)
                aria_spectrum   = ccpnImport.makeAriaSpectrum(peakList, self.molecule)
                
                ## test dimer 2.2
                sym_settings = self.getData(DC.DATA_SYMMETRY)[0]
        
                if sym_settings['symmetry_type'] in ["C2", "C3", "C5"]:
                   
                   for p in aria_spectrum.getPeaks():

                       p1 = p.getProton1Assignments()
                       h1 = p.getHetero1Assignments()

                       newp1 = [a for a in p1 if a.getAtoms()[0].getSegid() == self.__ref_segid]
                       newh1 = [a for a in h1 if a.getAtoms()[0].getSegid() == self.__ref_segid]

                       p.setProton1Assignments(newp1)
                       p.setHetero1Assignments(newh1)

                
                spectrum = self.REREAD(aria_spectrum)
                shifts   = self.REREAD(aria_shift_list)

                E = Experiment(spectrum, shifts)

                spec_data = None

                for x in X[DC.DATA_SPECTRUM]:
                   if ccpnImport.getObjectKeyString(peakList) == x['peaks']['ccpn_id']:
                      if ccpnImport.getObjectKeyString(shiftList) == x['shifts']['ccpn_id']:
                         spec_data = x
                         break

                if spec_data is None:
                    self.error(KeyError, 'Could not retrieve CCPN objects.')

                E.setDataSource(spec_data)

                experiments.append(E)

            self.experiments += experiments
            
            ## BARDIAUX CCPN Restraints
            self.distance_restraints = []

            for cType, constraintLists in constraint_lists.items():

                data_dir = self.getInfrastructure().get_data_directories()[cType]

                for i, constraintList in enumerate(constraintLists):


                    if X[cType][i]['enabled'] == YES:
                    
                        # TJS: Store CCPN objects as input to run if absent
                        # This is a record of what was actually used.
                        if ccpnAriaRun:
                          cSet = constraintList.nmrConstraintStore.serial
                          dataObj = ccpnAriaRun.findFirstData(className='ConstraintStoreData',
                                                              constraintStoreSerial=cSet,
                                                              constraintListSerials=[constraintList.serial],
                                                              ioRole='input')
 
                          if not dataObj:
                            ccpnAriaRun.newConstraintStoreData(constraintStoreSerial=cSet,
                                                               constraintListSerials=[constraintList.serial],
                                                               ioRole='input')
 
                        convert_to_aria = cType in [DC.DATA_AMBIGUOUS, DC.DATA_UNAMBIGUOUS] and \
                                          (X[cType][i]['add_to_network'] == YES or \
                                           X[cType][i]['calibrate'] != NO or \
                                           X[cType][i]['filter_contributions'] == YES)
                        
                        if convert_to_aria:
                            ## CCPN Restraint added to distance restraints list
                            ## TODO: DC.DATA_HBONDS ?

                            constraints, list_of_constraints = ccpnImport.getAriaDistanceRestraintsList(constraintList,
                                                                                           cType, self.molecule)
                            list_of_constraints.setListSource(X[cType][i])

                            # set all CCPN constraints as reliable if we don't filter them
                            if X[cType][i]['filter_contributions'] == NO:
                                [r.getReferencePeak().isReliable(1) for r in constraints]
                            
                            self.distance_restraints += constraints
                          
                            #self.distance_restraints += ccpnImport.getAriaDistanceRestraintsList(constraintList,
                            #                                                                    type, self.molecule)
                            
                            m = 'Contraints list %s added for analysis.' 
                            self.message(m % restraints_names[cType][i])
                                
                        else:
                            
                            m = 'CCPN Contraints list %s written to %s.'
                            
                            # TJS: Added molSytem specification - FormatConverter was asking
                            # if more than one was avail, even it though must be consistent
                            file_name = ccpnImport.dumpRestraintsList(constraintList, ccpChains[0].molSystem,
                                                                      data_dir, cType)
                            X[cType][i]['filename'] = file_name

                            self.message(m %(restraints_names[cType][i], file_name))


            ## Initial structures
            models_dict = {}
            molSystemCode = ccpChains[0].molSystem.code
            
            data_dir = self.getInfrastructure().get_data_directories()[DC.DATA_INITIAL_STRUCTURE]
            
            for i, structure in enumerate(ccpn_initial_models):
              
              file_name = ccpnImport.dumpModel(structure, ccpChains[0].molSystem,
                                               data_dir)
              
              structureEnsembleId = structure.structureEnsemble.ensembleId
              
              models_dict.setdefault(structureEnsembleId, [])
              if structure.serial not in models_dict[structureEnsembleId]:
                  models_dict[structureEnsembleId].append(structure.serial)

              X[DC.DATA_INITIAL_STRUCTURE][i]['filename'] = file_name
              X[DC.DATA_INITIAL_STRUCTURE][i]['format'] = 'cns'

            ## template structures
            data_dir = self.getInfrastructure().get_data_directories()[DC.DATA_TEMPLATE_STRUCTURE]
            
            for i, structure in enumerate(ccpn_template_models):
              
              file_name = ccpnImport.dumpModel(structure, ccpChains[0].molSystem,
                                               data_dir)
              
              structureEnsembleId = structure.structureEnsemble.ensembleId
              
              models_dict.setdefault(structureEnsembleId, [])
              if structure.serial not in models_dict[structureEnsembleId]:
                  models_dict[structureEnsembleId].append(structure.serial)

              X[DC.DATA_TEMPLATE_STRUCTURE][i]['filename'] = file_name
              X[DC.DATA_TEMPLATE_STRUCTURE][i]['format'] = 'cns'
            

            if ccpnAriaRun:
              for structureEnsembleId, structureSerials in models_dict.items():
                
                dataObj = ccpnAriaRun.findFirstData(className='StructureEnsembleData',
                                                    ensembleId=structureEnsembleId,
                                                    molSystemCode = molSystemCode,
                                                    modelSerials=structureSerials,
                                                    ioRole='input')
 
                if not dataObj:
                  ccpnAriaRun.newStructureEnsembleData(ensembleId=structureEnsembleId,
                                                       molSystemCode = molSystemCode,
                                                       modelSerials=structureSerials,
                                                       ioRole='input')

            
            self.message('CCPN import done.')

    def get_ccpn_data_sources(self):

        import aria.DataContainer as DC

        d = {}
        
        ## BARDIAUX 2.2
        ## some data_types are not related to CCPN for the moment
        not_ccpn_data = [DC.DATA_SYMMETRY, DC.DATA_SSBRIDGE, DC.DATA_HISPATCH, DC.DATA_CISPROPATCH]
        
        for DT in DC.DATA_TYPES:

            x = self.getData(DT)

            for xx in x:

                val = None

                if DT == DC.DATA_SPECTRUM:

                    ## this assumes, that both peak and shift list
                    ## have a ccpn_id; in principle one could also
                    ## read one list from a ccpn data model, the
                    ## other one from an ARIA xml file. However,
                    ## this is not supported for now.

                    if xx['peaks']['format'] == 'ccpn' or \
                       xx['shifts']['format'] == 'ccpn':

                        val = xx

                ## BARDIAUX 2.2
                elif DT in not_ccpn_data:
                    continue
                
                elif xx.get('format') == 'ccpn':
                    val = xx

                if not val:
                    continue

                if not DT in d:
                    d[DT] = []

                d[DT].append(val)

        return d

    def map_chains(self):

        ## BARDIAUX 2.2
        import aria.DataContainer as DC
        sym_settings = self.getData(DC.DATA_SYMMETRY)[0]
        ##print sym_settings
        
        #if sym_settings['symmetry_type'] <> "C2":
        #   return

        if sym_settings['enabled'] <> YES:
           return
        
        segids = [c.getSegid() for c in self.molecule.get_chains()]
        chains = [] 
        self.__ref_segid = segids[0]
            
        for s in segids:
            chain = [a for c in self.molecule.get_chains() for r in c.getResidues() for a in r.getAtoms() if c.getSegid() == s]
            chain.sort(lambda a, b: cmp(a.getId(), b.getId()))
            chains.append(chain)

        if len(segids) == 1:
            # monomer:
            pass
        elif len(segids) > 1:
            if sym_settings['symmetry_type'] in ['C2','C3','C5']:
                from aria.tools import circular_permutation
                for atoms in zip(*chains):
                    cp = circular_permutation(atoms)
                    for s in cp:
                        for homo in s[1:]:
                            s[0].addHomologous(homo)
                            
            elif sym_settings['symmetry_type'] == "D2":
                from aria.tools import circular_permutation
                for atoms in zip(*chains):
                    cp = circular_permutation(atoms)
                    for i in range(0,4):
                        if (i % 2) == 0:
                            atoms[i].setHomologous(cp[i][1:])
                        else:
                            x = cp[i][1:]
                            x.reverse()
                            atoms[i].setHomologous(x)

        del chain, chains

    def read_data(self):

        from time import clock
        from aria.tools import Load
        import os
        import aria.DataContainer as DC

        ## there is only one CCPN data model ...

        ccpn = self.ccpn_model['filename'] <> ''
        infra = self.getInfrastructure()
        cache_on = infra.getSettings()['cache'] == YES and not ccpn

        ## check which data sources are stored in a ccpn model

        if ccpn:
            self.ccpn_data_sources = self.get_ccpn_data_sources()
        else:
            self.ccpn_data_sources = {}

        if not self.ccpn_data_sources:
            
            msg = '-' * 20 + ' Reading data ' + '-' * 20
            self.message(msg)

        ## Set the directory of all data-sources to the local
        ## data-directory. So, local data-files are loaded.

        f = infra.get_local_filename

        ## Spectra

        ## replace filename with local filenames

        for sd in self.getData(DC.DATA_SPECTRUM):

            if sd in self.ccpn_data_sources.get(DC.DATA_SPECTRUM, ()):
                continue

            fn =  sd['shifts']['filename']
            sd['shifts']['filename'] = f(fn, DC.DATA_SPECTRUM)
            fn = sd['peaks']['filename']
            sd['peaks']['filename'] = f(fn, DC.DATA_SPECTRUM)

        ## Molecule

        if not DC.DATA_SEQUENCE in self.ccpn_data_sources.keys():

            ## TODO: [0]...
            
            molecule_data = self.getData(DC.DATA_SEQUENCE)[0]
            fn = molecule_data['filename']
            molecule_data['filename'] = f(fn, DC.DATA_SEQUENCE)

        ## Read data

        t = clock()

        d = {1: 'Cache is enabled.',
             0: 'Cache is disabled.'}

        self.message(d[cache_on])

        ## Attempt to load pickles

        if cache_on:

            name = infra.get_cache_filename()

            if os.path.exists(name):
                try:
                    cache = Load(name)

                except:
                    self.message('Could not load cache file')
                    cache = None
            else:
                s = 'Cache file does not exist. Creating new file.'
                self.message(s)
                cache = None
        else:
            cache = None

        ## If cache has been loaded successfully, first check 
        ## data integrity

        if cache is not None:

            data_ok = self.check_cached_molecule(cache)

            if not data_ok:
                self.message('Molecule definition has changed.')

            if data_ok:
                
                data_ok = self.check_cached_spectra(cache)

                if not data_ok:
                    self.message('Peak/shift listse have changed.')
                    
        else:
            data_ok = 0

        if cache_on and cache is not None:

            d = {1: 'Loading cached data files...',
                 0: 'Cache is not up to date. Reloading...'}

            self.message(d[data_ok])

        ## If none of the data have been modified, use pickles.

        ## Molecule

        if not DC.DATA_SEQUENCE in self.ccpn_data_sources.keys():

            if data_ok:

                ## Instead of instantiating new AtomFactory,
                ## get pickled version. Factory is already frozen.

                atom_factory = cache['atom_factory']

            else:
                ## Otherwise, create a new AtomFactory

                from aria.Singleton import AtomFactory

                atom_factory = AtomFactory(__new_instance__ = 1)
                cache = None

            ## Load and set molecule definition and freeze factory

            self.read_molecule(cache)

            ## BARDIAUX 2.2
            self.map_chains()
            
            atom_factory.freeze()

        self.experiments = []

        ## Load and set experiments
        
        self.read_spectra(cache)

        if not self.ccpn_data_sources:

            self.message('Data files read.')
            self.debug('Time: %ss' % str(clock() - t))
        
        ## Re-cache data source-files is necessary.

        if cache_on and not data_ok:
            
            t = clock()
            self.cache_data(self.molecule, self.experiments)
            self.message('Data files cached.')
            self.debug('Time: %ss' % str(clock()-t))
        
        ## Load remaining data from CCPN data model

        if ccpn and self.ccpn_data_sources:
            self.read_ccpn_data()
        
    def apply_filters(self, experiment):

        from time import clock

        msg = 'Applying filters to spectrum/shiftlist "%s":'
        self.message(msg % experiment.getName())
            
        ## filter shift-assignments
        
        t = clock()
        
        shifts = experiment.getShiftList()

        filtered_shifts = self.chemical_shift_list_filter(shifts)

        m = 'Chemical shift list filtered: %d / %d shifts (%.2f %%) ' + \
            'removed.'

        a = len(filtered_shifts)
        b = len(shifts)

        self.message(m % (b - a, b, 100. - 100. * a / b))
        self.debug('Time: %ss' % str(clock() - t))

        ## filter cross-peaks

        t = clock()

        spectrum_data = experiment.getDataSource()

        spectrum_filter = self.noesy_spectrum_filter
        spectrum_filter.getSettings().update(spectrum_data['peaks'])
        # BARDIAUX 2.2
        # diagonal peaks filter
        spectrum_filter.getSettings().update(spectrum_data)
            
        spectrum = experiment.getSpectrum()
        filtered_peaks= spectrum_filter(spectrum)
        statistic = spectrum_filter.compile_statistic()

        m = 'NOESY spectrum "%s" filtered: %d / %d peaks (%.2f %%) ' + \
            'removed. Details:'

        a = len(filtered_peaks)
        b = len(spectrum)

        self.message(m % (experiment.getName(), b - a, b,
                          100. - 100. * a / b))
        self.message('- no. of invalid proton 1 shifts: %d' \
                     % statistic['proton1'])
        self.message('- no. of invalid proton 2 shifts: %d' \
                     % statistic['proton2'])
        self.message('- no. of invalid peak sizes:      %d' \
                     % statistic['size'])
        self.message('- no. of diagonal peaks:          %d' \
                     % statistic['diagonal'])
        self.message('- no. of unassigned peaks:        %d' \
                     % statistic['unassigned'])        

        ## TODO: fraction hard coded, support fraction via setting
        if float(a) / float(b) < 0.1:
            m = '%.0f %% of the peaks of spectrum "%s" have been removed. ' + \
                'This might be due to a wrong choice of the peak size type '+ \
                '("volume" or "intensity") that will be used for ' + \
                'calibration and the definition of distance bounds.'
            self.error(ValueError, m % (100.*(1.- float(a)/float(b)),
                                        experiment.getName()))
            
        self.debug('Time: %ss' % str(clock() - t))

        return spectrum_filter.result

    def assign_spectra(self, experiments):

        from aria.AriaPeak import AriaPeak
        from aria.PeakAssigner import PeakAssignerTextPickler
        from time import clock
        from aria.DataContainer import DATA_SPECTRUM, DATA_SYMMETRY
        import os

        msg = '-' * 20 + ' Assigning spectra ' + '-' * 20
        self.message(msg)
        
        t = clock()
        id = 0
        aS = self.peak_assigner.getSettings()
        aria_peaks = []

        unassigned_peaks = {}

        for experiment in experiments:

            spectrum_data = experiment.getDataSource()

            use_assignments = spectrum_data['use_assignments'] == YES

            m = 'Creating seed assignment for spectrum "%s" ...'

            if not use_assignments:
                m += ' Manual assignments will not be used.'

            self.message(m % experiment.getName())

            spectrum = experiment.getSpectrum()

            unassigned_peaks[spectrum] = []

            shifts = experiment.getFilteredShiftAssignments()
            cross_peaks = experiment.getFilteredPeaks()

            ## set frequency-windows

            aS.update(spectrum_data['peaks'])

            ## set spectrum-specific assigner settings
            
            shift_data = spectrum_data['shifts']
            aS['default_shift_error'] = shift_data['default_shift_error']
            aS['use_assignments'] = spectrum_data['use_assignments']
            ## BARDIAUX 2.2
            aS['ref_segid'] = self.__ref_segid
            aS['spec_type'] = spectrum.getExperimentData()['ambiguity_type']
            aS['sym_type'] = self.getData(DATA_SYMMETRY)[0]['symmetry_type']
            aS['structural_rules_enabled'] = spectrum_data['structural_rules_enabled']
            
            #self.message('Using Settings:\n'+str(aS))

            self.peak_assigner.setShiftAssignments(shifts)

            ## If crosspeak is already assigned, classify
            ## it as 'reliable' if requested.

            if spectrum_data['trust_assigned_peaks'] == YES:

                msg = 'Trusting fully assigned peaks.'
                self.message(msg, verbose_level = VL_SETTINGS)
                
                [p.isReliable(p.isAssigned()) for p in cross_peaks]

            msg = 'Using the following settings:\n%s' % str(aS)
            self.message(msg, verbose_level = VL_SETTINGS)

            ## switch if spectrum intra or inter/?

            for peak in cross_peaks:

                contributions = self.peak_assigner.assign(peak)

                if not contributions:
                    unassigned_peaks[spectrum].append(peak)
                    continue

                ap = AriaPeak(id, peak)

                for contribution in contributions:
                    ap.addContribution(contribution)

                if self.use_restraint_weights:

                    import numpy

                    weights = [c.spectral_weight for c in contributions]

                    ap.restraint_weight = numpy.sum(weights) / \
                                          float(len(weights))
                    
                aria_peaks.append(ap)

                id += 1

            m = 'Done. %d / %d peaks (%.2f %%) could not be ' + \
                'assigned (s. peak assigner report file for details).'

            a = len(unassigned_peaks[spectrum])
            b = len(cross_peaks)

            ## emergency brake: if more than 90% of all peaks could
            ## not be assigned, we stop. this could be due to a wrong
            ## nuclei - freq. dimension mapping, erroneous window sizes etc.

            if (float(a) / b > 0.9):
                msg = """More than 90% of the cross peaks could not be assigned. This might be due to interchanged frequency dimensions or undersized frequency windows. Please check your setup."""
                self.warning(msg)
                self.halt()
            
            self.message(m % (a, b, 100. * a / b))
            self.debug('Time: %ss' % str(clock() - t))

        ## dump unassigned peaks
        ## TODO: destination hard-coded

        dest = self.getInfrastructure().get_data_directories()[DATA_SPECTRUM]
        dest = os.path.join(dest, 'peak_list.unassigned')

        pickler = PeakAssignerTextPickler()

        ## generate generic information

        info = []

        for experiment in experiments:
            spectrum = experiment.getSpectrum()
            spec_name = spectrum.getName()
            n_peaks = len(experiment.getFilteredPeaks())
            n_unassigned = len(unassigned_peaks[spectrum])
            
            s = 'Spectrum: %s, %d / %d (%.1f%%) not assigned'
            
            info.append(s % (spec_name, n_unassigned, n_peaks,
                             n_unassigned * 100 / float(n_peaks)))
        
        pickler.set_info(info)
        
        for peaks in unassigned_peaks.values():
            peaks.sort(lambda a, b, cmp = cmp:
                       cmp(a.getNumber(), b.getNumber()))
            pickler.add_peaks(peaks)

        pickler.dump(dest)

        return aria_peaks

    def createFirstIteration(self):

        from aria.Iteration import Iteration

        aria_peaks = self.assign_spectra(self.experiments)

        iteration = Iteration(-1)
        [iteration.addPeak(p) for p in aria_peaks]

        ## BARDIAUX 2.2 Distance Restraints
        [iteration.addDistanceRestraint(r) for r in self.distance_restraints]
        
        return iteration

    def load_and_preprocess_data(self):
        
        self.read_data()
        self.preprocess_data()

    def preprocess_data(self):

        from aria.NOESYSpectrumFilter import NOESYSpectrumFilterTextPickler
        from aria.DataContainer import DATA_SPECTRUM

        import os
        
        ## apply filters to all spectra

        msg = '-' * 19 + ' Filtering input data ' + '-' * 19
        self.message(msg)

        filter_results = []

        for experiment in self.experiments:
            filter_results.append(self.apply_filters(experiment))

        ## sort entries

        sorted_results = []
        for d in filter_results:

            k = d.keys()
            k.sort(lambda a, b: cmp(a.getNumber(), b.getNumber()))

            sorted_results += [(kk, d[kk]) for kk in k]

        ## dump information on filtered peaks
        ## TODO: destination still hard-coded

        dest = self.getInfrastructure().get_data_directories()[DATA_SPECTRUM]
        dest = os.path.join(dest, 'peak_list.filtered')

        report = NOESYSpectrumFilterTextPickler()
        report.dump(sorted_results, dest)

    def do_analysis(self, last_iteration):
        """
        performs a list of tests/checks for structures
        generated in the last iteration and the solvent-
        refinement.
        """

        import os, time

        analyser = self.getAnalyser()
        infra = self.getInfrastructure()
        protocol = self.getProtocol()
        protocol_settings = protocol.getSettings()
        
        ## Run CNS analysis scripts
        analyser.cns_analyses(last_iteration)

        # BARDIAUX 2.2
        ## Run CNS analysis scripts water refinement
        if protocol_settings['water_refinement']['enabled'] == YES:
            analyser.cns_analyses(last_iteration, is_water=1)

        ## Run quality checks
        last_it = last_iteration.getNumber()
        path_last_it = infra.get_iteration_path(last_it)

        ## TODO: better solution

        s = self.getProtocol().getSettings()['iteration_settings']
        n_best_structures = s[last_iteration.getNumber()]['number_of_best_structures']
        pdb_directories = [(path_last_it, n_best_structures)]

        if protocol_settings['water_refinement']['enabled'] == YES:
            n_structures =protocol_settings['water_refinement']['n_structures']
            pdb_directories.append((infra.get_refinement_path(),
                                    n_structures))

        [analyser.quality_checks(*s) for s in pdb_directories]

        ## Append the quality-checks summary to our summary report.

        if protocol.reportsWritten():

            ## name of summary file

            from legacy.QualityChecks.QualityChecks import FILENAME_REPORT
            src = os.path.join(path_last_it, FILENAME_REPORT)

            ## summary report
            from aria.Protocol import REPORT_SUMMARY
            dst = os.path.join(path_last_it, REPORT_SUMMARY)

            if os.path.exists(src) and os.path.exists(dst) and \
                   analyser.isEnabled():

                d = {'time': time.ctime()}

                f = open(src)
                lines = f.readlines()
                f.close()

                ## exchange header

                lines = lines[1:]
                lines = QUALITY_CHECKS_HEADER % d + ''.join(lines).strip()

                f = open(dst, 'a')
                f.write(lines)
                f.close()


    def doExport(self, iteration):
        """
        Write ARIA results in other data formats
        """
        
        reporter = self.getReporter()


    def finalize(self):

        from time import ctime
        
        if self.getMolecule() is None:
            self.error(StandardError, 'No molecule set. Have the data been loaded?')
        

        msg = '-' * 16 + ' Preparing structure engine ' + 16 * '-'
        self.message(msg)

        self.getStructureEngine().setup_data()

        protocol = self.getProtocol()
        protocol.finalize_engine(self.getMolecule())

        ## start ARIA

        s = 'Starting ARIA main protocol on %s'
        self.message(s % str(ctime()))

        return self.createFirstIteration()

    def go(self, use_condor=False):

        ## the condor stuff should go into the project file at some point

        self.getStructureEngine().set_condor(use_condor)
        
        first_iteration = self.finalize()
        self.run_protocol(first_iteration)
   
    ## < Mareuil
    def setCOMMAND(self):
        engine = self.getStructureEngine()
        job_manager = engine.getJobScheduler()
        jm_settings = job_manager.getSettings()
        return [jm_settings['default_command'], jm_settings['job_management']]
       
    def run_protocol(self, iteration):

        from time import ctime
        
        MODE = self.setCOMMAND()
        last_iteration = self.getProtocol().go(iteration, self.getMolecule(), MODE)
        ## Mareuil >
        ## run analysis on last iteration and solvent-refined
        ## structures (if enabled)

        self.do_analysis(last_iteration)

        ## Run exporters
        self.doExport(last_iteration)

        ## remove temporary directories

        self.getInfrastructure().cleanup()

        ## STOP here

        s = 'ARIA run completed at %s'
        self.message(s % str(ctime()))

        if AriaBaseClass.log_file is not None:
            AriaBaseClass.log_file.close()


class ProjectSingleton(Project, Singleton):

    def __init__(self, *args, **kw):
        Singleton.__init__(self)
        Project.__init__(self, *args, **kw)

    def __getstate__(self):
        return self.__dict__, Singleton.__getstate__(self)

    def __setstate__(self, s):
        self.__dict__, singleton_state = s
        Singleton.__setstate__(self, singleton_state)
        
        self.initialize()

class ProjectThread(Thread):

    def __init__(self, project):
        Thread.__init__(self)
        self.p = project

    def run(self):
        self.p.go()

class ProjectXMLPickler(XMLBasePickler):

    order = ['name', 'version', 'author', 'date', 'description',
             'comment', 'references', 'working_directory', 'temp_root',
             'run', 'file_root', 'cache', 'cleanup',
             'data', 'structure_generation',
             'protocol', 'analysis', 'report']

    def __init__(self):

        from aria.Analyser import AnalyserXMLPickler
        
        self.project = self
        self.analysis = AnalyserXMLPickler()

        self.setReportPickler()
        self.setStructGenPickler()
        self.setCNSPickler()
        self.setProtocolPickler()
        self.setDataPickler()
        self.setIterationPickler()

    def setReportPickler(self):
        import aria.Report as R
        
        self.report = R.ReportXMLPickler()
        self.noe_restraint_list = R.NOEListSettingsXMLPickler()
        self.ccpn = R.CCPNSettingsXMLPickler()
        self.molmol = R.MolMolSettingsXMLPickler()
        self.spectra = R.UpSpecSettingsXMLPickler()
        
    def setStructGenPickler(self):

        from aria.cns import CNSXMLPickler
        import aria.JobManager as JobManager
        
        self.cns = CNSXMLPickler()
        self.job_manager = JobManager.JobManagerXMLPickler()
        self.host = JobManager.HostSettingsXMLPickler()

    def setCNSPickler(self):

        import aria.DataContainer as DC

        self.annealing_parameters = DC.AnnealingParametersXMLPickler()
        self.ambiguous_restraints = DC.AmbiguousParametersXMLPickler()
        self.unambiguous_restraints = DC.UnambiguousParametersXMLPickler()
        self.dihedral_restraints = DC.DihedralParametersXMLPickler()
        self.karplus_restraints = DC.KarplusParametersXMLPickler()
        self.rdc_restraints = DC.RDCParametersXMLPickler()
        self.flat_bottom_harmonic_wall = DC.FBHWParametersXMLPickler()
        self.md_parameters = DC.MDParametersXMLPickler()
        self.hbond_restraints = DC.HBondParametersXMLPickler()
        ## BARDIAUX 2.2
        self.symmetry_restraints = DC.SymmetryParametersXMLPickler()
        ## BERNARD 2.3
        self.logharmonic_potential = DC.LogHarmonicParametersXMLPickler()
        
    def setProtocolPickler(self):
        
        import aria.Protocol as Protocol
        import aria.DataContainer as DC

        self.protocol = Protocol.ProtocolXMLPickler()
        self.iteration = Protocol.IterationSettingsXMLPickler()
        self.water_refinement = DC.WaterRefinementXMLPickler()

    def setDataPickler(self):

        import aria.DataContainer as DC
        #from importFromCcpn import CCPNDataXMLPickler

        self.molecule = DC.SequenceDataXMLPickler()
        self.shifts = DC.ShiftDataXMLPickler()
        self.peaks = DC.PeakDataXMLPickler()
        ## BARDIAUX rMat
        self.experiment_data = DC.ExperimentDataXMLPickler()

        self.ccpn_model = DC.CCPNDataXMLPickler()

        f = DC.BoundCorrectionXMLPickler
        
        self.lower_bound_correction = f(DC.LowerBoundCorrection)
        self.upper_bound_correction = f(DC.UpperBoundCorrection)
        self.template_structure = DC.TemplateDataXMLPickler()
        self.initial_structure = DC.InitialStructureDataXMLPickler()
        self.hbonds = DC.HBondDataXMLPickler()
        self.dihedrals = DC.DihedralDataXMLPickler()
        self.jcouplings = DC.KarplusDataXMLPickler()
        self.rdcs = DC.RDCDataXMLPickler()

        self.cyspatch = DC.CysPatchXMLPickler()
                
        self.ssbonds = DC.SSBondDataXMLPickler()

        
        self.ssbridge = DC.SSBridgeXMLPickler()
        self.hispatch = DC.HisPatchXMLPickler()
        ## BARDIAUX 2.2
        self.symmetry =  DC.SymmetryXMLPickler()
        self.cispropatch = DC.CisProPatchXMLPickler()        
        self.znpatch = DC.ZnPatchXMLPickler()
        self.other_data = DC.OtherDataXMLPickler()
        
        self.spectrum = DC.SpectrumDataXMLPickler()
        self.ambiguous_distance_restraints = \
                                       DC.AmbiguousDistanceDataXMLPickler()
        self.unambiguous_distance_restraints = \
                                       DC.UnambiguousDistanceDataXMLPickler()

    def setIterationPickler(self):

        import aria.PeakAssigner as PA
        import aria.Merger as Merger
        import aria.Calibrator as Calibrator
        import aria.ViolationAnalyser as VA
        import aria.ContributionAssigner as CA
        ## BARDIAUX 2.2
        import aria.Network as Network

        self.assignment = PA.PeakAssignerXMLPickler()
        self.merging = Merger.MergerXMLPickler()
        self.calibration = Calibrator.CalibratorXMLPickler()
        self.violation_analysis = VA.ViolationAnalyserXMLPickler()
        self.partial_assignment = CA.ContributionAssignerXMLPickler()
        ## BARDIAUX 2.2
        self.network_anchoring = Network.NetworkXMLPickler()

    def _xml_data_state(self, x):

        import aria.DataContainer as DC

        order = ['ccpn_model', 'molecule', 'spectrum', 'unambiguous_distance_restraints',
                 'ambiguous_distance_restraints', 'jcouplings', 'rdcs',
                 'hbonds', 'dihedrals', 'ssbonds', 'ssbridge', 'hispatch', 'cispropatch', 'znpatch', 'symmetry', 'other_data', # BARDIAUX 2.2
                 'initial_structure', 'template_structure']

        e = XMLElement()

        f = x.getData

        e.ccpn_model = x.ccpn_model

        e.molecule = f(DC.DATA_SEQUENCE)
        e.spectrum = f(DC.DATA_SPECTRUM)
        e.jcouplings = f(DC.DATA_KARPLUS)
        e.rdcs = f(DC.DATA_RDCS)
        e.hbonds = f(DC.DATA_HBONDS)
        e.dihedrals = f(DC.DATA_DIHEDRALS)
        e.ssbonds = f(DC.DATA_SSBONDS)
        e.ssbridge = f(DC.DATA_SSBRIDGE)
        e.hispatch = f(DC.DATA_HISPATCH)
        ## BARDIAUX 2.2

        e.cispropatch = f(DC.DATA_CISPROPATCH)
        e.symmetry = f(DC.DATA_SYMMETRY)
        e.znpatch = f(DC.DATA_ZNPATCH)
        e.other_data = f(DC.DATA_OTHER)        

        
        e.template_structure = f(DC.DATA_TEMPLATE_STRUCTURE)
        e.initial_structure = f(DC.DATA_INITIAL_STRUCTURE)
        e.ambiguous_distance_restraints = f(DC.DATA_AMBIGUOUS)
        e.unambiguous_distance_restraints = f(DC.DATA_UNAMBIGUOUS)

        optional = ('jcouplings', 'rdcs', 'hbonds', \
                    'dihedrals', 'ssbonds', 'ssbridge', 'hispatch', 'cispropatch', 'znpatch', 'other_data',# BARDIAUX 2.2
                    'template_structure', 'ambiguous_distance_restraints',
                    'unambiguous_distance_restraints')

        for tag in optional:

            if not getattr(e, tag):
                delattr(e, tag)
                order.remove(tag)

        e.set_tag_order(tuple(order))

        return e

    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)

        ## data

        e.data = self._xml_data_state(x)

        ## structure generation

        engine = x.getStructureEngine()

        order = ('engine', 'cns', 'job_manager')
        
        struct_gen = XMLElement(tag_order = order)
        struct_gen.engine = engine.getName()
        struct_gen.cns = engine
        struct_gen.job_manager = engine.getJobScheduler()
        
        e.structure_generation = struct_gen

        ## protocol

        e.protocol = x.getProtocol()

        ## analysis

        e.analysis = x.getAnalyser()

        ## report
        
        e.report = x.getReporter()

        ## attributes

        s = x.getSettings()

        e.name = s['name']
        e.version = s['version']
        e.author = s['author']
        e.date = s['date']
        e.description = s['description']
        e.comment = s['comment']
        e.references = s['references']
        e.working_directory = s['working_directory']
        e.temp_root = s['temp_root']
        e.run = s['run']
        e.file_root = s['file_root']

        e.cache = s['cache']
        e.cleanup = s['cleanup']


        return e

    def get_ccpn_data(self, e):

        #from importFromCcpn import CCPNData
        from aria.DataContainer import CCPNData
        
        ## try to create CCPN data structure
        ## compatibility code; ccpn stuff has been added to
        ## version 2.1

        from aria.tools import as_tuple
        
        if hasattr(e, 'ccpn_model'):

            ccpn = as_tuple(e.ccpn_model)

            if len(ccpn) <> 1:

                self.error(StandardError, 'Only one element "ccpn_model" expected. %d given.' %len(ccpn))

            else:
                ccpn = ccpn[0]

        else:

            ## either ARIA XML format error, or version < 2.1

            self.message('XML element "ccpn_model" is missing and has been added to your project file. This XML element is mandatory from ARIA version 2.1 on.')

            ccpn = CCPNData()

            ccpn.reset()

        return ccpn

    def add_data(self, project, e):

        from aria.tools import as_tuple
        import aria.DataContainer as DC

        p = project

        ccpn = self.get_ccpn_data(e)

        p.ccpn_model = ccpn

        p.addData(e.molecule)

        initial_structures = as_tuple(e.initial_structure)
        
        if len(initial_structures) > 1:
            s = 'Only 1 initial structure is supported, %d given.'
            self.error(ValueError, s % len(initial_structures))
        
        p.addData(initial_structures[0])

        optional = ('spectrum', 'jcouplings', 'rdcs', 'hbonds', \
                    'dihedrals', 'ssbonds', 'ssbridge', 'hispatch', 'cispropatch','znpatch', 'other_data', # BARDIAUX 2.2
                    'template_structure', 'ambiguous_distance_restraints',
                    'unambiguous_distance_restraints')
        
        # BARDIAUX 2.2
        if hasattr(e, 'symmetry'):
            p.addData(e.symmetry)

        else:
            from aria.DataContainer import Symmetry
            z = Symmetry()
            z.reset()

            self.message("No symmetry data specified. Symmetry will be disabled.")
            p.addData(z)


        for tag in optional:
            try:
                value = as_tuple(getattr(e, tag))
            except:
                continue

            [p.addData(d) for d in value]

    def load_from_element(self, e):

        from aria.Singleton import ProjectSingleton
        from aria.Infrastructure import Infrastructure

        s = ProjectSettings()
        s['name'] = str(e.name)
        s['version'] = float(e.version)
        s['author'] = str(e.author)
        s['date'] = str(e.date)
        s['description'] = str(e.description)
        s['comment'] = str(e.comment)
        s['references'] = str(e.references)
        s['working_directory'] = str(e.working_directory)
        s['temp_root'] = str(e.temp_root)
        s['run'] = str(e.run)

        if self.relaxed:

            entity = s.getEntity('file_root')
            
            m = entity.is_mandatory()
            entity.mandatory(0)

        s['file_root'] = str(e.file_root)

        if self.relaxed:
            entity.mandatory(m)

        s['cache'] = str(e.cache)
        s['cleanup'] = str(e.cleanup)

        p = ProjectSingleton(s, __new_instance__ = 1)

        ## TODO: remove ...

        infra = Infrastructure(s)
        p.setInfrastructure(infra)

        ## structure generation

        struct_gen = e.structure_generation
        engine = struct_gen.cns
        engine._setName(struct_gen.engine)
        engine.setJobScheduler(struct_gen.job_manager)
        
        p.setStructureEngine(engine)

        ## protocol
        
        p.setProtocol(e.protocol)

        ## add data

        self.add_data(p, e.data)
        
        ## analysis

        p.setAnalyser(e.analysis)

        ## reporter

        p.setReporter(e.report)

        ## do final initialization

        p.initialize()
        p.getProtocol().initialize()

        return p

Project._xml_state = ProjectXMLPickler()._xml_state
