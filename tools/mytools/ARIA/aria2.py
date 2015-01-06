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
$Date: 2010/03/23 15:27:16 $
"""


import sys, os

ARIA_ENV = 'ARIA2'

CCPN_ENV = 'CCPNMR_TOP_DIR'

PATH_MODULES = ['src/py']

MODULES = ('aria.Analyser', 'aria.ariabase', 'aria.AriaPeak', 'aria.AriaXML', \
           'aria.Atom', 'aria.Calibrator', 'aria.Chain', 'aria.ChemicalShiftFilter', \
           'aria.ChemicalShiftList', 'aria.Contribution', 'aria.ConversionTable', \
           'aria.PDBReader', 'aria.Report', 'aria.mathutils', 'aria.ContributionAssigner', \
           'aria.CrossPeak', 'aria.CrossPeakFilter', 'aria.DataContainer', 'aria.Datum', \
           'aria.Experiment', 'aria.Factory', 'aria.Infrastructure', 'aria.Iteration', \
           'aria.JobManager', 'aria.Merger', 'aria.NOEModel', 'aria.NOESYSpectrum', \
           'aria.NOESYSpectrumFilter', 'aria.OrderedDict', 'aria.PeakAssigner', \
           'aria.Project', 'aria.Protocol', 'aria.Residue', 'aria.Settings', \
           'aria.ShiftAssignment', 'aria.ShiftAssignmentFilter', 'aria.Singleton', \
           'aria.SpinPair', 'aria.StructureEnsemble', 'aria.Topology', 'aria.TypeChecking', \
           'aria.ViolationAnalyser', 'aria.cns', 'aria.conversion', 'aria.mathutils', \
           'aria.tools', 'aria.xmlutils', 'aria.Molecule', 'aria.FloatFile', \
#           'Scientific.IO.PDB', 'aria.MolMol', 'aria.WhatifProfile', 'aria.RmsReport', \
           'aria.scientific.PDB', 'aria.MolMol', 'aria.WhatifProfile', 'aria.RmsReport', \
           'aria.Network', 'aria.CovalentDistances', 'aria.Relaxation')

USAGE = \
"""Invalid arguments. Use --help for help

Usage: aria2 [options] [project XML file]
"""

HELP = \
'''     
Options:

   --convert [-t] conversion.xml
   
      Reads in a conversion XML file, say "conversion.xml", and
      converts all files specified into various output-formats.
      Valid formats are: ansig, nmrview, xeasy. 

      Options:
      
        -t  Creates a template (i.e. pre-formatted) conversion XML file and
            stores it in "conversion.xml". In order to run the conversion,
            the template-file must be edited/completed by the user.

   -g, --gui
   
      Start the graphical project editor.

   --help
   
      This help screen.

   --no-test
   
      Disables the initial dry run of the commands (specified in the host
      list) used to launch a structure calculation.
      
   --output=filename
   
      Stores ARIAs output in the file "filename". 

   --overview
   
      Gives a brief overview of the given project: Date, author, description,
      data sources, etc.

   --project_template
   
      Use this option if you want ARIA to create a template XML file for a
      new project. NOTE: setting-up a new project is very much simpified by
      using ARIAs GUI.

   -s, --setup [-f]
   
      To setup a project, the project XML file is read and checked for errors.
      After reading the project file, the following actions are performed:

      1. Create directory tree
      2. Copy all data files from its source locations to the
         project directory.
      3. Copy CNS protocols

      Switch: -f
 
      If the project already exists, neither data files in the project
      directory nor CNS protocols will be overwritten. To enforce
      overwriting use the switch -f.

   --verbose-level=n
   
      Set verbose level to n. Default: n=0. If n=1, additional information
      will be shown.

How to create, setup and run a typical ARIA project:

1. Data conversion
------------------

  ARIA 2.1 uses a unified XML description for all data sources. Unless your
  data are already stored in ARIAs XML format, they first need to be converted
  into ARIAs XML format. Version 2.1 supports XML document definitions for the
  following data sources:

     1. Molecule definition
     2. Peak lists
     3. Chemical shift lists

  Currently, ARIA supports the file formats used by the programs Ansig,
  NMRView, and XEasy (as well as several other formats, e.g. Pronto and Sparky,
  if you use the CCPN format converter; see below). J couplings, H bonds or
  RDCs data can directly be incorporated by providing them in CNS ".tbl"
  format. In order to perform the conversion, some information on your data
  needs to be specified in a simple XML file, "conversion.xml", say. You may
  create an empty conversion template XML file by using the command

      aria2 --convert -t conversion.xml               

  Here, "-t" means "template". This command creates a new (text) file,
  "conversion.xml", formatted in XML. It contains the essential information
  needed to perform a data conversion. The next step is to specify those
  informations for your data, for instance the format or your spectrum files,
  proton / hetero dimensions etc. This is done by
  filling in some fields in the newly created "conversion.xml" file.

  Once your completed "conversion.xml" file has been stored, run ARIA again
  to perform the actual conversion:

      aria2 --convert conversion.xml

  This command first loads and parses the conversion file. It then reads in
  your data files (i.e spectra, shifts-lists, sequence-file), converts them
  into XML and stores them. The new data XML files can then be edited using a
  XML or text-editor.

2. The project file 
-------------------

  The complete definition of a project, i.e. data sources, parameters for the
  minimization protocol, parameters for the other sub-modules of ARIA etc., is
  encapsulated in the project\'s XML file.
   
  If a project file has been specified in the conversion XML file (cf.
  <project> in conversion XML), a project XML file has been generated
  during data conversion. In that case, the project file already references the
  converted data.

  If a project filename has been omitted, you can create an (empty) project
  template XML file by invoking the command

      aria --project_template [filename]

  To complete the project template file, you then have to fill in the fields
  referencing your data (represented as XML, cf. (1)) by hand.

  In order to setup an ARIA project (cf. (3)), it is necessary to provide
  some additional information:

      - working directory               The root directory of your project
      
      - file_root                       The \'nickname\' of your project.
                                        For instance, the file_root serves as
                                        basename for all PDB-files: e.g. if
                                        file_root is set to "bpti", PDB-files
                                        will be called bpti_1.pdb ... .
                                        
      - temporary path
      
      - list of available machines      In particular if you want to run the
                                        structure calculation on several
                                        machines simultaniously, you need
                                        to specify a list of those machines

                                        If the list is empty, ARIA uses the
                                        local host.
      - Path of CNS binary                                   


2.1 Editing a project
---------------------

  In general you can use ARIAs graphical user interface (GUI), an XML editor or
  a text editor to display and edit your project. However, ARIA\'s GUI is
  intended to streamline the project setup and further provides brief
  descriptions/help for most of the parameter settings. For starting the GUI
  with the project "your_project.xml" use the command

      aria2 --gui your_project.xml    or
      
      aria2 -g your_project.xml

  If the project file is omitted, the GUI starts without loading a project.


3. Project setup
----------------

  Once the project XML file has been completed, the project is ready for
  setup. To setup a project "new_project.xml", say, run ARIA again:

      aria2 -s new_project.xml

  This command reads in the XML file "new_project.xml" and checks it for
  errors. If the project has been loaded successfully, the following actions
  are performed:

  1. Directory-tree creation

     It creates the full directory tree needed by ARIA and CNS using
     {WORKING_DIRECTORY}/{RUNxxx} as its root (WORKING_DIRECTORY and RUNxxx
     are specified in the project file).
             
  2. Data setup

     Copies all (XML) data-files from its source locations (as specified in the
     project file) into the local "data" directory of the directory tree.

     When running an ARIA project, all data-files are read from the *local*
     directory. In other words, only modifications of the data stored in the
     local "data" directory are considered when re-running a project. The
     source-locations specified in the project file are only used for project
     setup.
             
  3. Copy CNS-specific files

     All CNS-specific files (protocols, topology/parameters files etc) are
     copied from their source location (ARIAs installation path) into the local
     directory "cns/...".
          
  If the project has already been set up, re-setup skips all existing files.
  To enforce overwriting of existing files, use the option -f:

         aria2 -sf new_project.xml

  Forced setup overwrites/updates the following files:

     a) Data files in the project\'s local data directory.
     
     b) CNS-specific files.


4. Running ARIA
---------------

  If the project "new_project.xml" has been successfully set up, start ARIA by
  invoking the command:

      aria2 new_project.xml

  Use the GUI to modify protocol- or minimization-parameters. As long as no 
  data-files or CNS protocols have been changed, it is not necessary to setup
  the project again. 
   
      
Example
--------

Data-conversion, setup and running ARIA for a typical project:

1) aria2 --convert -t conversion.xml    Creates an empty conversion XML
                                        template file which has to be filled
                                        in on order to convert your data.
                                        
2) aria2 --convert conversion.xml       Convert the data specified in
                                        "conversion.xml". We assume, that the
                                        conversion XML does not contain a
                                        project-filename. Otherwise step (3)
                                        can be skipped.

3) aria2 --project_template project.xml Creates an empty project XML template.

4) aria2 -g project.xml                 Starts ARIAs GUI to edit/complete the
                                        project "project.xml".

5) aria2 -s project.xml                 Setup "project.xml".

6) aria2 project.xml                    Run ARIA



Further examples
----------------

Running ARIA:

   aria2 a_new_project.xml              Runs ARIA on a_new_project.xml. The
                                        project must be already setup
                                        (see below.)
   aria2 --output=aria.out project.xml  Runs ARIA on project.xml and stores
                                        the output in aria.out

Starting ARIAs GUI:

   aria2 -g                             start GUI
   aria2 --gui your_project.xml         start GUI and edit your_project.xml

Project setup etc.
   
   aria2 -s another_project.xml         setup another_project.xml
   aria2 -sf existing_project.xml       re-setup existing_project.xml
   aria2 --overview test_project.xml    display project info

Conversion:

   aria2 --convert my_conversion.xml    Reads the XML file 'my_conversion.xml'
                                        and converts all files specified.

   aria2 --convert -t template.xml      Creates a raw conversion XML file which
                                        has to be completed in order to run the
                                        actual conversion.
'''

def print_sequence(data):
    print 'Filename:', data['filename']
    print 'Format:  ', data['format']

def print_spectrum(data):
    print 'Chemical shifts:'
    print '  Filename:', data['shifts']['filename']
    print '  Format:  ', data['shifts']['format']
    print 'Peaks:'
    print '  Filename:', data['peaks']['filename']
    print '  Format:  ', data['peaks']['format']
    print
    print 'Enabled:', data['enabled']
    print 'Proton1 freq. window: %.2f ppm' % data['peaks']['proton1_shift_err']
    print 'Proton2 freq. window: %.2f ppm' % data['peaks']['proton2_shift_err']
    print 'Hetero1 freq. window: %.2f ppm '% data['peaks']['hetero1_shift_err']
    print 'Hetero2 freq. window: %.2f ppm' % data['peaks']['hetero2_shift_err']

def print_template_structure(data):
    print 'Filename:', data['filename']
    print 'Format:  ', data['format']
    print 'Enabled: ', data['enabled']

def default_printer(data):
    print data

def show_help():
    from aria.tools import make_block

    block = make_block(HELP )
    block = '\n'.join(block)
                       
    print block

def get_path():
    
    print "ARIA's environment variable (%s) missing. " % ARIA_ENV + \
          'Please specify ARIA root path:'

    return raw_input()

def get_aria_root():

    if not os.environ.has_key(ARIA_ENV):
        missing = 1
    else:
        aria_root = os.environ[ARIA_ENV]
        missing = not os.path.exists(aria_root)

    while missing:
        aria_root = get_path()
        missing = not os.path.exists(aria_root)
        if not missing:
            os.environ[ARIA_ENV] = aria_root

    return aria_root

def welcome():

    import aria.ariabase as ariabase

    message = \
"""
ARIA Version %s. Authors: Benjamin Bardiaux, Michael Habeck, Jens Linge,
Therese Malliavin, Sean O'Donoghue, Wolfgang Rieping, and Michael Nilges.

If you use this software, please quote the following reference(s):

Rieping W., Habeck M., Bardiaux B., Bernard A., Malliavin T.E.,
Nilges M.(2007) ARIA2: automated NOE assignment and data integration in NMR
structure calculation. Bioinformatics 23:381-382
"""
    print message % ariabase.AriaBaseClass().get_version_string()

def check_python():
    import sys

    version = float(sys.version[:3])

    if version < 2.4:
        print 'Python version 2.4 or higher required.'
        sys.exit(1)
        
def check_numeric_slice():

    # BARDIAUX
    msg_NUMERIC = "\nThe version of Numeric (%s) is known to be incompatible with ARIA.\nConsider reverting to a more stable version (like 23.8).\n"
    msg_NUMPY   = "\nThe version of numpy (%s) is known to be incompatible with ARIA.\nConsider reverting to a more stable version.\n"
    
    msg = {'numpy'   : msg_NUMPY,
           'numeric' : msg_NUMERIC}
    
    numerix = 'numpy'
    
    try:
        from numpy import ones
        from numpy import __version__ as NUMERIC_VERSION        
    except:
        print "No numpy found."
        sys.exit(1)
    
    if not len(ones(10)[2:]):
        print msg[numerix] % NUMERIC_VERSION
        sys.exit(1)
        
def invalid_args():

    import sys

    print USAGE
    sys.exit(1)

def get_project_name(s):

    import os

    path, filename = os.path.split()
    filename, ext = os.path.splitext()

    if ext == '':
        ext = '.xml'

    return '%s/%s%s' % (path, filename, ext)
    
def parse_args(args):

    import getopt, sys

    options = 'gs:f:t'
    long_options = ('gui', 'setup=', 'convert', 'overview', \
                    'project_template', 'output=', 'help', 'debug',
                    'no-test', 'verbose-level=', 'condor',
                    'convert_ccpn')
    
    try:
        opt_list, args = getopt.getopt(args, options, long_options)
    except:
        invalid_args()

    d = {}

    for key, value in opt_list:
        if key in d:
            invalid_args()

        if value == '':
            value = None
            
        d[key] = value

    if not d and not args:
        invalid_args()

    not_none = ('-s', '--setup', '--output', '--verbose-level')

    has_args = ('--convert', '--convert_ccpn', '--project_template',
                '--debug', '--no-test', '--verbose-level')

#    has_args = ('--convert', '--project_template',
#                '--debug', '--no-test', '--verbose-level')
    
    for option in not_none:
        if option in d:
            if d[option] is None:
                invalid_args()

    for option in has_args:
        if option in d:
            if not args:
                invalid_args()

    if '-s' in d:
        if d['-s'] in ('-f', 'f') or d['-s'] is None:
            if not args:
                invalid_args()

    elif '--setup' in d:
        if d['--setup'] in ('-f', None) or d['--setup'] is None:
            if not args:
                invalid_args()

    elif '-g' in d:
        if args:
            d['-g'] = args[0]
            
    elif '--gui' in d:
        if args:
            d['--gui'] = args[0]

    elif '-t' in d:
        if not args or not '--convert' in d:
            invalid_args()

    elif '--overview' in d:
        if args:
            if d['--overview'] is not None:
                invalid_args()
        else:
            invalid_args()

    elif '--help' in d:
        if args:
            invalid_args()

    if len(args):
        args = args[0]

    return d, args
    
def run_gui(project_file):
    from aria.gui import gui
            
    gui.go(project_file)

def load_project(filename):
    import aria.AriaXML as AriaXML
    import os
    import sys
    from xml.parsers import expat
    from aria.TypeChecking import is_type

    filename = os.path.expanduser(filename)

    if not os.path.exists(filename):
        print 'Could not load project', filename
        sys.exit(0)

    pickler = AriaXML.AriaXMLPickler()

    if AriaBaseClass.log_stdout:
        print '\nLoading project "%s"...' % os.path.basename(filename)

    try:
        project = pickler.load(filename)
    except expat.error, msg:
        print 'XML format error.'
        raise expat.error, msg

    if not is_type(project, 'Project'):
        s = 'Specified XML file %s is not an ARIA project file.'
        print s % filename
        sys.exit(1)

    ## debug: otherwise, project instance would not
    ## be accessible.

    global __project
    __project = project

    return project

def convert_data(filename, create_template=0, use_ccpn=0):
    from aria.conversion import Converter
    import os, sys

    if use_ccpn and not CCPN_ENV in os.environ:

        print 'CCPN environment variable (CCPNMR_TOP_DIR) not set.'

        import sys
        sys.exit(1)

    filename = os.path.expanduser(filename)

    if not create_template and not os.path.exists(filename):
        print 'Could not open', filename
        sys.exit(0)

    converter = Converter()

    if create_template:

        ok = not os.path.exists(filename)

        if not ok:

            a = raw_input('\nFile %s exists. Overwrite (y/n)? ' \
                          % filename)

            ok = a.lower() == 'y'

        if ok:
            converter.write_template(filename)
        
    else:

        if use_ccpn:
            print '\nStarting data conversion using CCPNMR ...\n'
            converter.convert_ccpn(filename)

        else:
            print '\nStarting data conversion ...\n'
            converter.convert(filename)

def write_project_template(filename):

    try:
        dst = open(filename, 'w')
    except:
        print 'Could not create %s.' % filename
        sys.exit(0)

    from aria.ariabase import PROJECT_TEMPLATE
    from os.path import join

    fn = join(AriaBaseClass.data_path, PROJECT_TEMPLATE)
    src = open(fn)
    lines = src.readlines()
    src.close()

    dst.writelines(lines)
    dst.close()

    if AriaBaseClass.log_stdout:
        s = 'Project xml-template %s created.'
        print s % filename
    
def setup_aria(opt, args):

    import os

    if opt in ('-f', 'f'):
        force = 1
        project_file = args
    else:
        force = 0
        project_file = opt

    project = load_project(project_file)

    infra = project.getInfrastructure()
    run_path = infra.get_run_path()
    path, filename = os.path.split(project_file)

    run = project.getSettings()['run']

    if AriaBaseClass.log_stdout:
        print 'Setting-up %s (run: %s)' % (project_file, run)

    project.setup(force = force)

    if AriaBaseClass.log_stdout:
        print '\nDone.\n'
        print 'For running ARIA: aria2 %s' % project_file
        print

def project_info(args):

    ## TODO: could (should) be implemented as __str__ method
    ## of ProjectSettings.
    
    import aria.tools as tools
    import os
    import aria.DataContainer as DC
    from aria.OrderedDict import OrderedDict

    project = load_project(args)

    settings = project.getSettings()

    keys = ['name', 'version', 'author', 'description',
            'comment', 'references']

    max_key_len = max([len(k) for k in keys])

    LINE_LENGTH = 80

    print '\nInformation on project %s:\n' % os.path.basename(args)

    ## generic 

    for k in keys:
        new_key = k + ': ' + ' ' * (max_key_len - len(k))
        
        value = settings[k]

        prefix = new_key.upper()

        value = tools.wrap_string(str(value), LINE_LENGTH - len(prefix))
        value = tools.indent([value], prefix)
        
        print value

    print

    ## data

    keys = OrderedDict()

    keys[DC.DATA_SEQUENCE] = 'Sequence'
    keys[DC.DATA_SPECTRUM] = 'Spectra'
    keys[DC.DATA_HBONDS] = 'H-Bonds'
    keys[DC.DATA_DIHEDRALS] = 'Dihedrals'
    keys[DC.DATA_KARPLUS] = 'Karplus'
    keys[DC.DATA_RDCS] = 'Residual dipolar couplings'
    keys[DC.DATA_SSBONDS] = 'Disulfide bonds'
    keys[DC.DATA_SSBRIDGE] = 'Covalent disulfide bridges'
    keys[DC.DATA_UNAMBIGUOUS] = 'Unambiguous distance restraints'
    keys[DC.DATA_AMBIGUOUS] = 'Ambiguous distance restraints'
    keys[DC.DATA_TEMPLATE_STRUCTURE] = 'Template structures'

    printer = {DC.DATA_SEQUENCE: print_sequence,
               DC.DATA_SPECTRUM: print_spectrum,
               DC.DATA_TEMPLATE_STRUCTURE: print_template_structure}

    for key, label in keys.items():
        
        data = project.getData(key)
        
        if not data:
            print '%s: None\n' % label
            continue

        print ('%s:\n'+ '-' * (len(label) + 1)) % label

        p = printer.get(key, default_printer)

        for d in data:
            p(d)
            print
            
def run_aria(project_file, log_file = None, debug = 0, test_commands = 1,
             verbose_level = 0, use_condor=False):

    if log_file is not None:

        import os

        try:
            log_file = os.path.expanduser(log_file)
            f = open(log_file, 'w')
        except:
            print 'Could not create log-file: %s', log_file
            f = None

        AriaBaseClass.log_file = f

    if debug:
        from aria.TypeChecking import check_type
        AriaBaseClass.display_debug = 1
        check_type.active = 1

        print 'Running in --debug mode. Type checking has been turned on.'

    AriaBaseClass.verbose_level = verbose_level

    import os
    import aria.Settings as Settings

    ## Turn off path-checking. This is done for the following
    ## reasons:
    ##
    ## 1. All files/directories are already checked when the
    ##    project is set up.
    ## 2. Once a project has been set up, all files (data etc.)
    ##    have been copied to their local locations in the
    ##    project directory tree. If the user deletes any of
    ##    them, tough luck. If the user instead deletes source
    ##    files (i.e. the files which are specified in project-xml),
    ##    ARIA should not complain. Therefore, path-checking
    ##    is turned of once a project has been set up.
    
    Settings.Path.global_mandatory = 0
    
    project = load_project(project_file)

    project.check(test_commands)
    project.load_and_preprocess_data()

    project.go(use_condor)

def convert(source, dest, format):
    s = 'Converting source file %s (%s format) into xml-output %s'
    print s % (source, dest, format)
    
## check python version

check_python()

## BARDIAUX: check numeric
check_numeric_slice()

## add ARIA modules to python path

aria_root = os.path.join(get_aria_root(), 'src/py')

sys.path.insert(0, aria_root)

## modules = [os.path.join(aria_root, p) for p in PATH_MODULES \
##            if not p in sys.path]

## sys.path = modules + sys.path

from aria.ariabase import AriaBaseClass

if __name__ == '__main__':

    welcome()

    ## check command-line args

    a, args = parse_args(sys.argv[1:])

    if '--debug' in a:

        ## check whether all modules are importable .. except
        ## for tools.py

        import aria.tools as tools

        ## remove doubles ... just in case

        d = {}

        for m in MODULES:
            d[m] = 0

        MODULES = tuple(d.keys())

        tools.check_modules(MODULES)

    from aria.ariabase import AriaBaseClass
    AriaBaseClass.log_stdout = 1
    

    run_options = ('--output', '--debug', '--no-test', '--verbose-level',
                   '--condor')

    is_run_option = len([1 for o in run_options if o in a])

    if '-g' in a:
        run_gui(a['-g'])
    elif '--gui' in a:
        run_gui(a['--gui'])
    elif '-s' in a:
        setup_aria(a['-s'], args)

    elif '--setup' in a:
        setup_aria(a['--setup'], args)

    elif not a or is_run_option:

    ##    import profile
    ##    profile.run('run_aria(args)')

        debug = '--debug' in a
        test_commands = '--no-test' not in a
        condor = '--condor' in a #False

        verbose_level = a.get('--verbose-level', 0)

        try:
            verbose_level = int(verbose_level)
        except:
            import sys
            print 'Invalid arguments: verbose-level must be number (0-1).'
            sys.exit(0)

        run_aria(args, a.get('--output', None), debug, test_commands,
                 verbose_level, use_condor=condor)

    elif '--overview' in a:
        project_info(args)

    elif '--convert' in a:
        if '-t' in a:
            convert_data(args, create_template = 1)
        else:
            convert_data(args)

    elif '--convert_ccpn' in a:
        convert_data(args, use_ccpn=1)

    elif '--project_template' in a:
        write_project_template(args)

    elif '--help' in a:
        show_help()

    else:
        print 'Internal error. Options not parsed:', a, args
        invalid_args()

