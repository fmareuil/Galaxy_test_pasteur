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




## set correct python-path

ARIA_ENV = 'ARIA2'
PATH_MODULES = 'src/py/aria2'

## program-wide constants

YES = 'yes'
NO = 'no'
GZIP = 'gzip'
ALWAYS = 'always'
PROJECT_TEMPLATE = 'project_template.xml'
CCPN_EXTENSION = '.ccpn'

## verbose level constants

VL_STANDARD = 0
VL_SETTINGS = 2
VL_LOW = 1

def get_path():
    
    print 'Aria environment variable (%s) missing. ' % ARIA_ENV + \
          'Please specify root path:'

    return raw_input()

def get_aria_root():

    import os

    if not os.environ.has_key(ARIA_ENV):
        missing = 1
    else:
        aria_path = os.environ[ARIA_ENV]
        missing  = not os.path.exists(aria_path)

    while missing:
        aria_path = get_path()
        missing = not os.path.exists(aria_path)
        if not missing:
            os.environ[ARIA_ENV] = aria_path

    return aria_path
        
from aria.TypeChecking import *

class AriaBaseClass:

    use_restraint_weights = 0

    VERSION = 2.3
    VERSION_RELEASE = 1
    VERSION_FORMAT = '%.1f.%d'
    VERSION_STRING = VERSION_FORMAT % (VERSION, VERSION_RELEASE)

    display_warnings = 1
    display_messages = 1
    display_deprecated = 1
    display_debug = 0
    
    warnings_as_errors = 0
    
    wrap_lines = 1
    line_length = 80
    description_length = 100
    cache = 1
    save_memory = 1
    log_file = None
    # BB extendNmr
    log_gui = None
    log_stdout = 1
    # BB
    check_type.active = 0
    verbose_level = 0

    ## if Aria root path has not been set yet, do it now

    try:
        install_path
        data_path
        
        has_root = 1
    except:
        has_root = 0

    if not has_root:

        import sys, os
        
        install_path = get_aria_root()

        ## misc. files (templates, etc) are stored here
        data_path = os.path.join(install_path, 'src/py/data')

        ## CNS specific files
        toppar_path = os.path.join(install_path, 'cns/toppar')
        protocols_path = os.path.join(install_path, 'cns/protocols')
        analysis_path = os.path.join(protocols_path, 'analysis')
        
        cns_directories = {'toppar': toppar_path,
                           'protocols': protocols_path,
                           'analysis': analysis_path}
        del os
        del sys
    ## < Mareuil    
    def __init__(self, settings = None, name = None, MODE = None):
    ## Mareuil >
        """
        'name': every class which is inherited from AriaBaseClass can
        be assigned its own name. it will be used when displaying
        warnings, errors, messages etc.
        """
        
        if settings is not None:
            self.setSettings(settings)
        else:
            self.__settings = None

        if name is not None:
            self._set_name(name)

    def get_version_string(self):
        return self.VERSION_FORMAT % (self.VERSION, self.VERSION_RELEASE)
    
    ## < Mareuil
    def setMODE(self,v):
        check_string(v)
        self.default_MODE = v
    
    def getMODE(self):
        return self.default_MODE
    ## Mareuil >

    def _set_name(self, name):
        check_string(name)
        self._name = name

    def deprecated(self, msg):
        """
        can be used to display information about code-changes etc.
        """

        if self.display_deprecated:
            prefix = self.__compile_name('DEPRECATED')
            print self.__format(prefix, msg)

    def debug(self, msg):
        """
        checks the class variable display_debug. if non-zero
        'msg' is displayed.
        """

        if self.__class__.display_debug:
            self.message(msg, prefix = 'DEBUG')

    def error(self, exception = None, error = '', msg = None):

        import inspect

        if exception is None:
            exception = Exception

        ## shut-down job manager

        self.shutdown()

        try:
            frame = inspect.currentframe().f_back.f_back
            have_frame = 1
        except:
            have_frame = 0

        if have_frame:
            code = frame.f_code
            func_name = code.co_name
            filename = code.co_filename
            lineno = frame.f_lineno

            descr = 'File "%s", line %d in %s\n%s'
            msg = descr % (filename, lineno, func_name, error)

        if msg is None:
            msg = ''

        msg = 'USER ERROR <%s> ' % str(self.__class__) + msg

        self.__log(msg)

        raise exception, msg

    def __format(self, tag, msg):

        check_string(tag)
        check_string(msg)
        
        import aria.tools as tools

        if self.wrap_lines:
            lines = tools.make_block(msg, self.line_length - len(tag))
        else:
            lines = [msg]
            
        lines = tools.indent(lines, tag)

        return lines

    def __compile_name(self, prefix):
        if hasattr(self, '_name'):
            name = self._name
        else:
            name = self.__class__.__name__
            
        return '%s [%s]: ' % (prefix, name)

    def __print(self, prefix, msg, verbose_level):

        if verbose_level <= self.verbose_level or self.display_debug:
            lines = self.__format(prefix, msg)
            if self.log_gui:
                self.log_gui.write(lines + '\n')
            if self.log_stdout:
                print lines
            self.__log(lines)

    def __log(self, s):
        if self.log_file is not None:
            self.log_file.write(s + '\n')
            self.log_file.flush()

    def message(self, msg, prefix = 'MESSAGE', verbose_level = VL_STANDARD):
        if self.display_messages:
            prefix = self.__compile_name(prefix)
            msg = str(msg)
            self.__print(prefix, msg, verbose_level)

    def warning(self, msg, verbose_level = VL_STANDARD):
        
        if self.display_warnings:

            msg = str(msg)
            
            if self.warnings_as_errors:
                self.error(None, msg)
            else:
                prefix = self.__compile_name('WARNING')
                self.__print(prefix, msg, verbose_level)

    def halt(self):
        """
        aborts ARIA
        """

        import time, sys

        msg = '\nARIA halted at %s.\n' % str(time.ctime())

        ## close log-file

        if self.log_file is not None:
            self.log_file.write(msg)
            self.log_file.close()

        print msg

        sys.exit(0)

    def shutdown(self):

        return
        
        from aria.Singleton import ProjectSingleton

        try:
            project = ProjectSingleton()
            job_manager = project.getStructureEngine().getJobScheduler()
            job_manager.shutdown()

        except:
            import aria.tools as tools
            print tools.last_traceback()
            print 'Could not shutdown job-manager.'
            s = 'If this message occurs more than once, the method' + \
                ' AriaBaseClass.error has not been used properly!' + \
                ' E.g. within a try-except statement.'
            print s

    def setSettings(self, s):

        check_type(s, 'Settings')
        self.__settings = s

    def getSettings(self):

        if self.__settings is None:
            s = '%s: Settings are None.' % self.__class__.__name__
            self.error(s)

        return self.__settings
