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
from exceptions import Exception
from aria.xmlutils import XMLElement, XMLBasePickler

## TODO: how shall we handle this?

CACHE_FILENAME = 'cache'

class AriaDirectoryCreationError(Exception):

    def __str__(self):
        return 'Cannot create directory "%s"' % Exception.__str__(self)

class Infrastructure(AriaBaseClass):
    def __init__(self, settings = None):

        ## TODO: merge Infrastructure and Project

        check_type(settings, 'ProjectSettings')

        AriaBaseClass.__init__(self)
        self.setSettings(settings)
        self.checked_directories = []

        ## TODO: Infrastructure will become Project anyway
        self._set_name('Project')

    def get_aria_path(self):
        return AriaBaseClass.install_path
        
    def get_conversion_table(self):

        if not hasattr(self, 'conversion_table'):
            from aria.ConversionTable import ConversionTable
            self.conversion_table = ConversionTable()

        return self.conversion_table

    def get_project_path(self):
        return self.getSettings()['working_directory']

    def get_run_path(self):
        import os
        
        return os.path.join(self.get_project_path(),
                            'run' + self.getSettings()['run'])

    def get_local_filename(self, filename, t):
        """
        Depending on the data-type, t, this method returns
        the full local path of 'filename'. It is not checked
        whether the local actually exists.
        """

        import os

        local_path = self.get_data_directory(t)
        base = os.path.basename(filename)
        local_name = os.path.join(local_path, base)

        return local_name

    def get_data_directory(self, t):
        d = self.get_data_directories()
        if not t in d:
            s = 'Data-type "%s" not known. Valid types are: %s'
            self.error(KeyError, s % (str(s), str(d.keys())))

        return d[t]

    def get_data_directories(self):
        """
        Directory names, where Aria stores all data.
        """
        import aria.DataContainer as DC
        import os

        d = {DC.DATA_SEQUENCE: 'data/sequence',
             DC.DATA_SPECTRUM: 'data/spectra',
             DC.DATA_HBONDS: 'data/hbonds',
             DC.DATA_DIHEDRALS: 'data/dihedrals',
             DC.DATA_KARPLUS: 'data/jcouplings',
             DC.DATA_RDCS: 'data/rdcs',
             DC.DATA_SSBONDS: 'data/ssbonds',
             DC.DATA_AMBIGUOUS: 'data/distances',
             DC.DATA_UNAMBIGUOUS: 'data/distances',
             DC.DATA_TEMPLATE_STRUCTURE: 'data/templates',
             DC.DATA_INITIAL_STRUCTURE: 'data/begin',
             DC.DATA_OTHER : 'data/other'}

        base = self.get_run_path()

        for key, value in d.items():
            d[key] = os.path.join(base, value)

        return d

    def get_temp_root(self):
        return self.getSettings()['temp_root']

    def get_temp_path(self):
        return self.getSettings()['aria_temp_path']

    def get_file_root(self):
        return self.getSettings()['file_root']

    def get_psf_file(self):
        return self.get_file_root() + '.psf'

    def get_iteration_path(self, iteration_number):
        """
        returns the absolute path for iteration with number
        'iteration_number'
        """
        import os
        
        iteration_root = self.getSettings()['path_name_iterations']
        return os.path.join(self.get_run_path(),
                            '%s%d' % (iteration_root, iteration_number))  

    def get_refinement_path(self):
        import os

        d = self.getSettings()['path_name_iterations']
        return os.path.join(self.get_run_path(), os.path.split(d)[0], 'refine')

    ## BARDIAUX 2.2
    def get_iteration_graphics_path(self, iteration_number):
        import os
        
        it_path = self.get_iteration_path(iteration_number)
        return os.path.join(it_path, 'graphics')
        

    def get_cache_filename(self):
        import os

        return os.path.join(self.get_run_path(), CACHE_FILENAME)

    ## CNS stuff

    def get_cns_path(self):
        import os
        return os.path.join(self.get_run_path(), 'cns')

    def get_cns_templates_path(self):
        """
        The path where templates for CNS scripts are stored.
        """
        return AriaBaseClass.data_path

    def get_cns_protocols_path(self):
        import os
        return os.path.join(self.get_cns_path(), 'protocols')

    def get_cns_analysis_path(self):
        import os
        return os.path.join(self.get_cns_protocols_path(), 'analysis')

    def get_cns_topology_path(self):
        import os
        return os.path.join(self.get_cns_path(), 'toppar')

    def get_cns_data_directory(self, t):
        return self.get_cns_data_directories(t)

    def get_cns_data_directories(self, data_type = None):
        """
        CNS data directories.
        """

        import aria.DataContainer as DC
        import os

        if data_type is not None:
            check_string(data_type)
        
        d = {DC.DATA_SEQUENCE: 'data/sequence',
             DC.DATA_HBONDS: 'data/hbonds',
             DC.DATA_DIHEDRALS: 'data/dihedrals',
             DC.DATA_KARPLUS: 'data/jcouplings',
             DC.DATA_RDCS: 'data/rdcs',
             DC.DATA_SSBONDS: 'data/ssbonds',
             DC.DATA_AMBIGUOUS: 'data/distances',
             DC.DATA_UNAMBIGUOUS: 'data/distances',
             DC.DATA_INITIAL_STRUCTURE: 'begin',
             DC.DATA_OTHER: 'data/other'}

        base = self.get_cns_path()
        
        for key, value in d.items():
            d[key] = os.path.join(base, value)

        if data_type is not None:
            if data_type not in d:
                s = 'Data type %s not known.'
                self.error(KeyError, s % data_type)

            d = d[data_type]
            
        return d

    def cleanup(self):
        """
        Remove temporary working directory.
        """

        if self.getSettings()['cleanup'] == 'no':
            return

        import os
        
        os.system('rm -rf %s' % self.get_temp_path())

        self.message('Temporary directory "%s" removed.' \
                     % self.get_temp_path())
        
    def create_directory(self, name):
        import os

        try:
            os.makedirs(name)
        except Exception, message:
            self.warning(message)
            raise AriaDirectoryCreationError, name

    def __check_directory(self, name, msg):
        
        import tempfile, os, time

        tempfile.tempdir = name
        filename = tempfile.mktemp()
        filename += str(int(time.time()))

        try:
            f = open(filename, 'w')
            f.write('test')
            f.close()
        except Exception, msg:
            s = 'No write permission or write-error in directory %s' % name
            self.error(IOError, s)

        ## check whether we have read permissions

        try:
            f = open(filename, 'r')
            s = f.read()
            f.close()
        except Exception, msg:
            s = 'No read permission or read-error in directory %s' % name
            self.error(IOError, s)

        ## try to remove file
            
        try:
            os.unlink(filename)
        except Exception, msg:
            s = 'Could not remove test-file in directory %s' % filename
            self.error(IOError, s)

    def create_and_check_directory(self, name, msg = None):
        import os

        if name in self.checked_directories:
            return

        if msg is None:
            msg = 'Directory "%s" does already exist.' % name

        if os.path.exists(name):
            self.warning(msg)
        else:
            self.create_directory(name)
            self.checked_directories.append(name)

        ## check whether a file can be created,
        ## i.e. whether we have write permission

        self.__check_directory(name, msg)

    def finalize(self):

        import os, tempfile, time
        
        ## create temp directory

        tempfile.tempdir = ''

        local_temp = tempfile.mktemp()
        local_temp += str(int(time.time()))
        
        path = os.path.join(self.get_temp_root(),
                            'aria_temp.%s' % local_temp)
                            
        self.create_and_check_directory(path)

        settings = self.getSettings()
        settings['aria_temp_path'] = path

        s = 'Temporary directory has been set to %s'
        self.message(s % path)

        self.message('Finalized.', verbose_level = VL_LOW)

    def setup(self, force = 0):
        """
        Creates all directories (except the aria_temp path).
        All data files are copied to their respective destinations.
        If force is non-zero, overwriting of files is enforced.
        """

        import os

        settings = self.getSettings()

        ## create project-path

        path = self.get_project_path()
        warning = 'Project path "%s" does already exist.' % path
        self.create_and_check_directory(path, warning)

        ## create path for current run

        path = self.get_run_path()
        warning = 'Path of current run "%s" does already exist.' % path
        self.create_and_check_directory(path, warning)

        ## create paths for data used in the current run

        for path in self.get_data_directories().values():
            warning = 'Data path "%s" for current run does already exist.'
            self.create_and_check_directory(path, warning % path)
            
        ## create iteration sub-directories

        n_iterations = settings['n_iterations']
        for i in range(n_iterations):

            path = self.get_iteration_path(i)            
            self.create_and_check_directory(path)

        ## create 'analysis' directory for last iteration

        analysis_path = os.path.join(self.get_iteration_path(n_iterations - 1),
                                     'analysis')
        self.create_and_check_directory(analysis_path)
        
        ## create directory for solvent refinement

        self.create_and_check_directory(self.get_refinement_path())

        ## BARDIAUX 2.2
        ## create 'analysis' directory for solvent refinement
        analysis_path = os.path.join(self.get_refinement_path(),
                                     'analysis')
        
        self.create_and_check_directory(analysis_path)
        
        ## BARDIAUX 2.2
        ## create graphics directory in each iteration sub-dir
        
        n_iterations = settings['n_iterations']
        for i in range(n_iterations):
            
            path = self.get_iteration_path(i)
            graphics_path = os.path.join(path, 'graphics')
            self.create_and_check_directory(graphics_path)
        graphics_path = os.path.join(self.get_refinement_path(), 'graphics')    
        self.create_and_check_directory(graphics_path)
        
        ## create cns directory

        path = self.get_cns_path()
        warning = 'Path for cns calculations "%s" does already exist.' % path
        self.create_and_check_directory(path, warning)

        ## create 'begin' directory; needed by cns protocol

        path = os.path.join(self.get_cns_path(), 'begin')
        self.create_and_check_directory(path)

        ## create 'data' sub-tree

        for path in self.get_cns_data_directories().values():
            self.create_and_check_directory(path)

        ## copy cns protocols

        self.copy_protocols(force)
        
        ## create toppar directory and copy files

        path = os.path.join(self.get_cns_path(), 'toppar')
        self.create_and_check_directory(path)

        source_path = AriaBaseClass.cns_directories['toppar']
                            
        os.system('cp %s/*.* %s' % (source_path, path))

        self.message('Directory tree created.')

        ## copy data specified in the project to the current run dir

        self.message('Copying data files into local data-directory...')
        self.copy_data(force)

        ## In forced setup: remove existing psf-file.

        if force:
            import aria.DataContainer as DC
            import os
            
            path = self.get_cns_data_directories(DC.DATA_INITIAL_STRUCTURE)
            psf_file = os.path.join(path, self.get_psf_file())

            if os.path.exists(psf_file):
                try:
                    os.unlink(psf_file)
                except:
                    s = 'Could not remove existing psf-file'
                    self.message(s)

    def copy_protocols(self, force = 0):
        """
        Copies CNS core and analysis protocols.
        If force is non-zero, existing protocols will be overwritten
        """

        import os
        import aria.cns as cns
        from aria.tools import copy_file

        cns_dirs = AriaBaseClass.cns_directories
        source_path = cns_dirs['protocols']

        ## create path for cns protocols and copy them

        dest_path = self.get_cns_protocols_path()
        warning = 'CNS protocol path "%s" does already exist.' % dest_path
        self.create_and_check_directory(dest_path, warning)
        
        analysis_dest = self.get_cns_analysis_path()
        warning = ('CNS analysis protocols path "%s" does ' + \
                  'already exist.') % analysis_dest
        self.create_and_check_directory(analysis_dest, warning)

        ## set correct path for CNS analysis protocols
        analysis_protocols = [os.path.join('analysis', p) \
                              for p in cns.PROTOCOLS_ANALYSIS]

        cns_protocols = cns.PROTOCOLS_CORE + \
                        tuple(analysis_protocols)

        protocols_copied = 0

        for filename in cns_protocols:

            src = os.path.join(source_path, filename)
            dst = os.path.join(dest_path, filename)

            if force or not os.path.exists(dst):
                try:
                    copy_file(src, dst)
                except Exception, msg:
                    self.error(Exception, msg)
                    
                protocols_copied = 1

        ## Copy run.cns template to local path

        src = os.path.join(self.get_cns_templates_path(), cns.RUN_CNS)
        dst = os.path.join(dest_path, cns.RUN_CNS)

        if force or not os.path.exists(dst):
            try:
                copy_file(src, dst)
            except Exception, msg:
                msg += 'Could not copy run.cns template to local ' + \
                       'path "%s"' % dst
                self.error(Exception, msg)

            protocols_copied = protocols_copied or 1

        if protocols_copied:
            if force:
                s = 'Existing CNS protocols have been replaced ' + \
                    'with default protocols.'
                self.message(s)
                
            self.message('Protocols copied.')
            
    def copy_data(self, force = 0):
        """
        Copies all data specified in the project to the current run
        directory to have everything in one place.
        if force is non-zero, data-files are copied in any case.
        """
        from aria.Singleton import ProjectSingleton
        import aria.DataContainer as DC
        from aria.tools import copy_file
        from os.path import basename, exists, join

        project = ProjectSingleton()
        used_bases = {}

        ## TODO: take care that only one sequence file has been set

        for data_type, dst in self.get_data_directories().items():

            used_bases[data_type] = {}

            for data in project.getData(data_type):

                if data_type == DC.DATA_SPECTRUM:
                    files = [data['shifts'].getLocation()[0],
                             data['peaks'].getLocation()[0]]

                elif data_type == DC.DATA_INITIAL_STRUCTURE and \
                    not data['enabled'] == YES:
                    continue

                else:
                    files = [data.getLocation()[0]]


                for file in files:

                    if file == '':
                        continue
                    
                    base = basename(file)

                    if used_bases[data_type].has_key(base) and \
                       file not in used_bases[data_type][base]:
                        m = 'Basenames for %s data files have ' % data_type + \
                            'to be unique because they are copied ' + \
                            'to the same data directory.'
                        self.error(ValueError, m)

                    else:
                        used_bases[data_type][base] = []
                        
                    used_bases[data_type][base].append(file)
                    
                    if exists(join(dst, base)) and not force:
                    
                        m = 'File "%s" for data of type ' % base + \
                            '"%s" already exists in local ' % data_type + \
                            'data directory; local copy will be used.'
                        self.warning(m)
                        continue
                    
                    try:
                        copy_file(file, dst)
                    except IOError, message:
                        m = 'Copying file "%s" of data-type "%s" to ' + \
                            'destination "%s" failed: %s'
                        self.error(IOError, m % (file, data_type,dst,message))

                ## Copy misc. files

                if data_type == DC.DATA_SEQUENCE:
                    self.copy_linkage_files(data, force)
                    self.copy_topology_files(data, force)
                    self.copy_parameter_files(data, force)

    def copy_linkage_files(self, data, force = 0):

        from aria.DataContainer import SequenceData as SD
        from os.path import join, basename, exists
        from aria.tools import copy_file

        name = data['linkage_name']
        path = self.get_cns_topology_path()

        if name <> SD.USER_DEFINED:
            return

        src = data['linkage_filename']
        dst = join(path, basename(src))

        if not force and exists(dst):

            msg = 'User-defined Linkage definition ' + \
                  ' file "%s" already exists locally. ' + \
                  'Using local copy.'

            self.message(msg % src)

        else:
            try:
                copy_file(src, dst)
                msg = 'User-defined Linkage definition ' + \
                      ' file "%s" copied.'
                self.message(msg % src)

            except IOError, message:
                m = 'Could not copy user-defined ' + \
                    'Linkage file from "%s" to "%s".'
                self.error(IOError, m % (src, dst))

    def copy_topology_files(self, data, force = 0):

        from aria.DataContainer import SequenceData as SD
        from os.path import join, basename, exists
        from aria.tools import copy_file

        name = data['topology_name']
        path = self.get_cns_topology_path()        

        if name <> SD.USER_DEFINED:
            return

        src = data['topology_filename']
        dst = join(path, basename(src))

        if not force and exists(dst):

            msg = 'User-defined Topology definition ' + \
                  ' file "%s" already exists locally. ' + \
                  'Using local copy.'

            self.message(msg % src)

        else:
            try:
                copy_file(src, dst)
                msg = 'User-defined Topology definition ' + \
                      ' file "%s" copied.'
                self.message(msg % src)

            except IOError, message:
                m = 'Could not copy user-defined ' + \
                    'Topology file from "%s" to "%s".'

                self.error(IOError, m % (src, dst))

    def copy_parameter_files(self, data, force = 0):

        from aria.DataContainer import SequenceData as SD
        from os.path import join, basename, exists
        from aria.tools import copy_file

        name = data['parameter_name']
        path = self.get_cns_topology_path()
        
        if name <> SD.USER_DEFINED:
            return

        src = data['parameter_filename']
        dst = join(path, basename(src))

        if not force and exists(dst):

            msg = 'User-defined Parameter definition ' + \
                  ' file "%s" already exists locally. ' + \
                  'Using local copy.'

            self.message(msg % src)

        else:
            try:
                copy_file(src, dst)
                msg = 'User-defined Parameter definition ' + \
                      ' file "%s" copied.'
                self.message(msg % src)

            except IOError, message:
                m = 'Could not copy user-defined ' + \
                    'Parameter file from "%s" to "%s".'

                self.error(IOError, m % (src, dst))
