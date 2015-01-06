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

## CNS core protocols (SA etc.)
PROTOCOLS_CORE = ('boxtyp20.pdb', 'define_methyls_all.cns',
                  'define_methyls_ini.cns', 'dmso.pdb',
                  'flags_new.cns', 'generate.inp',
                  'generate_dmso.cns', 'generate_template.inp',
                  'generate_water.cns', 'print_coorheader.cns',
                  'read_data.cns', 'refine.inp',
                  'refine_dmso.inp', 'refine_water.inp',
                  'sa_l_cool1.cns', 'sa_l_hightemp.cns',
                  'sa_l_randomchain.cns', 'sa_l_reduced.cns',
                  'sa_ls_cool2.cns', 'sa_ltad_cool1.cns',
                  'sa_ltad_hightemp4.cns', 'setup_swap_init.cns',
                  'swap.cns', 'torsiontop.cns', 'xplortodiana3.inp',
                  'write_symnoe.cns', 'write_symnoe2.cns', 'newsymmetry.cns', 'sa_rc_hightemp.cns', # BARDIAUX
                  'freezn.cns', 'logn_auto_weight.cns', 'print_coorheader_auto.cns') # BARDIAUX
		   

## CNS analysis scripts. cf. run_analysis method
PROTOCOLS_ANALYSIS = ('wellordered.inp', 'rotares.cns', 'minimize.inp',
                      'rmsave.inp', 'print_noes.inp', 'ensemble_rmsd.inp',
                      'rmsd.inp', 'print_geom.inp', 'print_dih.inp',
                      'print_coup.inp', 'print_sani.inp', 'noe_violations.inp',
                      'cop.inp', 'energy.inp')

RUN_CNS = 'run.cns'

SECONDARY_STRUCTURE = 'secondarystructure.cns'

DIHEDRALS = 'dihedrals.tbl'
DIHEDRALS_CSI = 'dihedrals_csi.tbl'
DIHEDRALS_TALOS = 'dihedrals_talos.tbl'

HBONDS = 'hbonds.tbl'
HBONDS_CSI = 'hbonds_csi.tbl'

PLAN = 'plan.tbl'

CNS_OUTPUT_PATH_NAME = 'cns'

## < Mareuil
CSH_SCRIPT_REFINE = \
'''
# PBS facility
#PBS -N %(sge_job_name)s
#PBS -S %(sge_job_shell)s

# SGE facility
#$ -N %(sge_job_name)s
#$ -S %(sge_job_shell)s

## results will be stored here
setenv NEWIT %(iteration_path)s

## pdb path
setenv PATHPDB %(initial_path_pdb)s

## project path
setenv RUN %(cns_resource_path)s

## individual run.cns is stored here
setenv RUN_CNS %(cns_working_dir)s

## CNS working directory
cd %(cns_working_dir)s

## solves some NFS sync problems

cat %(cns_input_file)s > /dev/null

## command line
%(cns_executable)s < %(cns_input_file)s >! %(cns_output_file)s

touch done
'''

CSH_SCRIPT_SGE = \
'''
# PBS facility
#PBS -N %(sge_job_name)s
#PBS -S %(sge_job_shell)s

# SGE facility
#$ -N %(sge_job_name)s
#$ -S %(sge_job_shell)s
#$ -t 1-%(sge_cpu_num)s

## results will be stored here
setenv NEWIT %(iteration_path)s

## pdb path
setenv PATHPDB %(initial_path_pdb)s

## project path
setenv RUN %(cns_resource_path)s

## individual run.cns is stored here
setenv RUN_CNS %(cns_working_dir_sge)s$SGE_TASK_ID

## CNS working directory
cd %(cns_working_dir_sge)s$SGE_TASK_ID

## solves some NFS sync problems

cat %(cns_input_file)s > /dev/null

## command line
%(cns_executable)s < %(cns_input_file)s >! %(cns_output_file)s

touch done
'''

CONDOR_SCRIPT = \
'''
universe = standard

environment =  NEWIT=%(iteration_path)s; RUN=%(cns_resource_path)s; RUN_CNS=%(cns_working_dir)s; ANALYSIS=%(analysis_output_path)s; PATHPDB=%(initial_path_pdb)s

executable = %(cns_executable)s

+Joblength = "Short"
log = condor.log
error = condor.error
notification = NEVER
Requirements = Arch == "X86_64" && OpSys == "LINUX" && Memory > 200
'''

CSH_SCRIPT_REFINE_CONDOR = \
'''
## results will be stored here
setenv NEWIT %(iteration_path)s

## pdb path
setenv PATHPDB %(initial_path_pdb)s

## project path
setenv RUN %(cns_resource_path)s

## individual run.cns is stored here
setenv RUN_CNS %(cns_working_dir)s

## CNS working directory
cd %(cns_working_dir)s

## solves some NFS sync problems

cat %(cns_input_file)s > /dev/null

## command line
condor_submit condor.job
condor_wait condor.log
touch done
'''

CSH_SCRIPT_ANALYSIS = \
"""
## results will be stored here
setenv ANALYSIS %(analysis_output_path)s

## results from last iteration
setenv NEWIT %(iteration_path)s

## project path
setenv RUN %(cns_resource_path)s

## individual run.cns and output-files are stored here
setenv RUN_CNS %(cns_working_dir)s

## working directory
cd %(cns_working_dir)s

## command line
"""

FILE_CNS_TEMPLATE = \
"""
module(
filenames;
iterations;
)

evaluate (&iterations.Structures=%(n_structures)s)
evaluate (&iterations.Assignstruct=%(n_best_structures)s)
"""

GRID_JDL_SCRIPT = \
"""
Executable = "%(csh_script_grid)s";
Requirements = (other.GlueHostArchitecturePlatformType == "x86_64");
Rank = other.GlueCEStateEstimatedResponseTime;
RetryCount = 0;
ShallowRetryCount = 5;
InputSandbox = {"%(csh_working_dir)s", "%(input_tar_dir)s"};
OutputSandbox = {"%(input_tar_grid)s"}
"""

GRID_CSH_SCRIPT = \
"""#!/bin/csh -f
## base directory
setenv BASE `pwd`
setenv BASE_CNS /tmp/`mktemp tmp.XXXXX`
ln -s $BASE $BASE_CNS

## decompression and removing of archive
tar -xzf $BASE/%(input_tar_grid)s
rm $BASE/%(input_tar_grid)s

## results will be stored here
setenv NEWIT $BASE_CNS%(iteration_path)s

## pdb path
setenv PATHPDB $BASE_CNS%(initial_path_pdb)s

## project path
setenv RUN $BASE_CNS%(cns_resource_path)s

## individual run.cns is stored here
setenv RUN_CNS $BASE_CNS%(cns_working_dir)s

## CNS working directory
cd $BASE%(cns_working_dir)s

## command line
mv $BASE%(cns_working_dir)s%(cns_executable_grid)s $BASE%(cns_executable_grid)s
chmod 700 $BASE%(cns_executable_grid)s
$BASE%(cns_executable_grid)s < $BASE%(cns_input_file)s >! %(cns_output_file)s
touch done
cd $BASE
tar -czf %(input_tar_grid)s .%(cns_working_dir)s .%(run_resource_path)s
rm -rf $BASE%(cns_executable_grid)s $BASE_CNS
"""
## Mareuil >

## constants

NAME_NOE_RESTRAINTS_AMBIG = 'ambig.tbl'
NAME_NOE_RESTRAINTS_UNAMBIG = 'unambig.tbl'

CNS_TRUE = 'true'
CNS_FALSE = 'false'

class CNSSettings(Settings):

    def create(self):

        from aria.Settings import Path, GZipChoice, ChoiceEntity

        d = {'local_executable': \
             Path(description = 'Local CNS executable')}

        descr = 'Specifies whether CNS output-files shall be copied to the project directory [PROJECT_PATH]/[RUN]/structures/itx/cns. This affects CNS output files for the structure calculation, the solvent-refinement and the CNS analyses.'

        d['keep_output'] = GZipChoice(description = descr)

        descr = \
"""If set to '%s', ARIA runs a script that automatically creates a proper PSF file using the user-defined XML sequence file (specified in the project-xml file).

If set to 'always', any existing PSF-file is replaced with a new one. In particular if the molecule definition has been changed, this option *must* be set to 'always' in order to create a new PSF-file - otherwise the existing one will be used.

If set to '%s', the user must provide a proper PSF-file which has to be stored/named in [PROJECT_PATH]/[RUN]/cns/begin/[FILEROOT].psf (e.g. .../my_project/run1/cns/begin/bpti.psf).""" % (str(YES), str(NO))

        choices = (YES, NO, ALWAYS)
        d['create_psf_file'] = ChoiceEntity(choices, description = descr)

        descr = \
"""If set to '%s', ARIA will run a script that automatically creates a PDB-file [PROJECT_PATH]/[RUN]/cns/begin/[FILEROOT]_template.pdb (e.g. .../my_project/run1/cns/begin/bpti_template.pdb). The structure will serve as initial structure used by the minimization protocol in all structure calculation runs.

If set to 'always', the initial structure is created (and stored as specified above) anyway; the existing one is replaced. In particular if the molecule definition has been changed, this option *must* be set to 'always' in order to create a new initial structure - otherwise the existing one will be used.

If set to '%s', the user must provide an initial PDB-file. It can be defined in the <Protocol> section.""" % (str(YES), str(NO))

        d['generate_template'] = ChoiceEntity(choices, description = descr)

        nb_names = ('PROLSQ', 'PARMALLH6', 'PARALLHDG', 'OPLSX')
        descr = \
"""ARIA supports 4 parametrizations for non-bonded interactions."""

        d['nonbonded_parameters'] = ChoiceEntity(nb_names,
                                                 description = descr)

        descr = 'ARIA writes temporary CNS restraint files, ambig.tbl / unambig.tbl. They define all distance restraints that are used during structure calculation. If this option is set to "no", they are deleted after an iteration has been compelted. If set to "gzip", they are kept an gzipped.'

        d['keep_restraint_files'] = GZipChoice(description = descr)

        return d

    def create_default_values(self):

        d = {}

        d['keep_output'] = YES
        d['create_psf_file'] = YES
        d['generate_template'] = YES
        d['nonbonded_parameters'] = 'PROLSQ'
        d['keep_restraint_files'] = YES

        return d

class StructureEngine(AriaBaseClass):
    pass

## TODO: 'refine' script-name is hard-coded

class CNS(StructureEngine):

    _atom_selection_template = '(segid "%s" and resid %4d and name %4s)'
    _or_template = '    or %s'

    def __init__(self, settings, name = None):


        """
        'name': must be the same as project-xml's
                structure_generation:engine attribute.
                if the CNS structure-engine is unpickled automatically
                by ARIA's xml-pickler, name is set to 'cns'
        """

        check_type(settings, 'CNSSettings')
        check_type(name, STRING, NONE)

        self.setSettings(settings)
        self.__name = name
        self.set_callback(None)
        self.setDefaultValues()

        self.__cache = {}
        self.__missing_structures = []
        ## BARDIAUX 2.2
        self.__nrestraints = 0

        self.use_condor = False

    def set_condor(self, x):
        self.use_condor = x

    def setDefaultValues(self):
        self.__protocol = None
        self.__scheduler = None
        self.__annealing_parameters = None
        self.__md_parameters = None

    def _setName(self, name):
        check_string(name)
        self.__name = name

    def getName(self):
        return self.__name

    def setAnnealingParameters(self, p):
        check_type(p, 'AnnealingParameters')
        self.__annealing_parameters = p

    def getAnnealingParameters(self):
        return self.__annealing_parameters
    
    ## < Mareuil
    def revpath(self,path,start):
        """
        return the relative path of an absolute path
        """
        import os, string
        real = start
        pathtemp = string.split(path,"/")
        realpath = ""
        realpathtemp = pathtemp[int(pathtemp.index(os.path.basename(real)))+1:]
        for i in realpathtemp:
            realpath = realpath +"/"+ i
        return realpath
    ## Mareuil >

    def setMDParameters(self, p):
        check_type(p, 'MDParameters')
        self.__md_parameters = p

    def getMDParameters(self):
        return self.__md_parameters

    def setJobScheduler(self, js):
        check_type(js, 'JobScheduler')
        self.__scheduler = js

        ## subsitute 'default' for 'executable'-setting
        ## to engine's executable

        try:
            hosts = js.getSettings()['host_list']
        except:
            self.warning('Job manager has empty host-list.')
            return

        executable = self.getSettings()['local_executable']

        for host in hosts:
            if host['executable'] == 'default':
                host['executable'] = executable

    def getJobScheduler(self):
        return self.__scheduler

    def getInfrastructure(self):

        from aria.Singleton import ProjectSingleton

        project = ProjectSingleton()

        return project.getInfrastructure()

    def get_local_filename(self, fn, t):
        """
        Depending on the data-type, t, this method returns
        the full local path of 'filename' with respect to the
        CNS sub-branch of the project directory. It is not checked
        whether the local actually exists.
        """

        import os

        infra = self.getInfrastructure()
        local_path = infra.get_cns_data_directory(t)
        base = os.path.basename(fn)
        local_name = os.path.join(local_path, base)

        return local_name

    def getConversionTable(self):

        return self.getInfrastructure().get_conversion_table()

    def getPDBReader(self):

        from aria.PDBReader import PDBReader

        table = self.getConversionTable()

        return PDBReader(table)

    def _create_atom_selection(self, atom):

        if AriaBaseClass.cache:
            key = atom.getId()
            if key in self.__cache:
                return self.__cache[atom.getId()]

        from aria.PDBReader import BASE_TYPES

        residue = atom.getResidue()

        base_type = BASE_TYPES[residue.getChain().getType()]

        segid = atom.getSegid()

        args = (residue.getType(), atom.getName(), 'iupac', 'cns', base_type)

        atom_name = self.getConversionTable().convert_atom(*args)

        selection = self._atom_selection_template \
                    % (segid, residue.getNumber(), atom_name)

        if AriaBaseClass.cache:
            self.__cache[key] = selection

        return selection

    def _create_peak_selections(self, peak):

        """
        the current implementation of an ADR in CNS,
        averages over all atom-pairs (i.e. the information
        from which contribution a set of atom-pairs
        is coming is neglected.). we therefore loop over all
        contributions, then over all spinpairs and build
        the list of atom-selections.
        """

        atom_selections = []

        for c in peak.getContributions():

            ## TODO: introduce isActive or so?

            if c.getWeight() == 0.:
                continue

            for spin_pair in c.getSpinPairs():

                atom1, atom2 = spin_pair

                sel1 = self._create_atom_selection(atom1)
                sel2 = self._create_atom_selection(atom2)

                line = '%s %s' % (sel1, sel2)

                atom_selections.append(line)

        if not atom_selections:
            return None, None

        ## connect atom pair selections via 'or' statement

        first = atom_selections[0]
        others = [self._or_template % s for s in atom_selections[1:]]

        return first, others

    ## BARDIAUX 2.2
    # dump homolgous pekas for dimers  even multimers)
    def _create_homo_peak_selections(self, peak, n):

        """
        the current implementation of an ADR in CNS,
        averages over all atom-pairs (i.e. the information
        from which contribution a set of atom-pairs
        is coming is neglected.). we therefore loop over all
        contributions, then over all spinpairs and build
        the list of atom-selections.

        create the homoguous restraints
        n is the number of the homologuous atom
        """

        atom_selections = []

        for c in peak.getContributions():

            ## TODO: introduce isActive or so?

            if c.getWeight() == 0.:
                continue

            for spin_pair in c.getSpinPairs():

                atom1, atom2 = spin_pair

                atom1 = atom1.getHomologous()[n]
                atom2 = atom2.getHomologous()[n]

                sel1 = self._create_atom_selection(atom1)
                sel2 = self._create_atom_selection(atom2)

                line = '%s %s' % (sel1, sel2)

                atom_selections.append(line)

        if not atom_selections:
            return None, None

        ## connect atom pair selections via 'or' statement

        first = atom_selections[0]
        others = [self._or_template % s for s in atom_selections[1:]]

        return first, others

    def setup_data(self):

        from aria.DataContainer import DATA_HBONDS, DATA_DIHEDRALS, DATA_KARPLUS, \
             DATA_RDCS, DATA_SSBONDS, DATA_INITIAL_STRUCTURE, \
             DATA_AMBIGUOUS, DATA_UNAMBIGUOUS, DATA_OTHER

        from aria.Singleton import ProjectSingleton
        from aria.Chain import TYPE_PROTEIN
        from aria.tools import cat_files

        import os

        project = ProjectSingleton()
        infra = self.getInfrastructure()
        cns_dirs = infra.get_cns_data_directories()
        data_dirs = infra.get_data_directories()

        reader = self.getPDBReader()

        ## helper dict

        CA = 'category'
        FN = 'filenames'

        STD = 'standard'
        CSI = 'csi'
        CLASS = 'class'
        TYPE = 'type'

        dict = {DATA_HBONDS: {CA: TYPE,
                              FN: {STD: HBONDS,
                                   CSI: HBONDS_CSI}},

                DATA_AMBIGUOUS: {CA: None,
                                 FN: {None: NAME_NOE_RESTRAINTS_AMBIG}},

                DATA_UNAMBIGUOUS: {CA: None,
                                   FN: {None: NAME_NOE_RESTRAINTS_UNAMBIG}},

                DATA_DIHEDRALS: {CA: TYPE,
                                 FN: {STD: DIHEDRALS,
                                      'talos': DIHEDRALS_TALOS,
                                      CSI: DIHEDRALS_CSI}},

                DATA_KARPLUS: {CA: CLASS,
                               FN: {1: 'c1.tbl', 2: 'c2.tbl',
                                    3: 'c3.tbl', 4: 'c4.tbl',
                                    5: 'c5.tbl'}},

                DATA_RDCS: {CA: CLASS,
                            FN: {1: 'rdc1.tbl', 2: 'rdc2.tbl',
                                 3: 'rdc3.tbl', 4: 'rdc4.tbl',
                                 5: 'rdc5.tbl'}},

                DATA_SSBONDS: {CA: None,
                               FN: {None: 'ssbonds.tbl'}},

                DATA_OTHER: {CA: TYPE,
                             FN: {'planarity': PLAN}}}


        for data_type in dict.keys():

            files = {}

            ## remove existing data

            for file in dict[data_type][FN].values():

                dst = os.path.join(cns_dirs[data_type], file)
                if os.path.exists(dst):
                    os.unlink(dst)
                    m = 'Existing copy "%s" has been removed ' % dst + \
                        'and will be replaced by a new copy.'
                    self.warning(m)

            ## categorize data resources

            for data in project.getData(data_type):

                ##if data['enabled'] == NO: continue
                if data['enabled'] == NO or data['filename'] == '': continue # BARDIAUX 2.2

                category_key = dict[data_type][CA]
                if category_key is not None:
                    category = data[category_key]
                else:
                    category = None
                if not files.has_key(category): files[category] = []
                files[category].append(data)

            ## copy / cat resource files

            for category, data in files.items():

                dst = os.path.join(cns_dirs[data_type],
                                   dict[data_type][FN][category])

                src_dir = data_dirs[data_type]

                src = map(lambda d: d.getLocation()[0], data)
                src = map(os.path.basename, src)
                src = map(lambda f, d = src_dir, os = os:
                          os.path.join(d, f), src)

                try:
                    cat_files(src, dst)
                except IOError, message:
                    print src, dst
                    self.error(IOError, message)

        ## If user has specified an initial pdb-file (i.e. the
        ## pdb-file which will be used as initial structures for
        ## the SA protocol), we convert it into IUPAC format and it
        ## them in CNS' 'begin' directory.

        initial_structures = project.getData(DATA_INITIAL_STRUCTURE)
        src_dir = data_dirs[DATA_INITIAL_STRUCTURE]
        dst_dir = cns_dirs[DATA_INITIAL_STRUCTURE]

        for structure in initial_structures:

            if structure['enabled'] == NO: continue

            base = os.path.basename(structure['filename'])
            name = os.path.join(src_dir, base)

            ## TODO: chain type hard coded

            chain = reader.read(name, format = structure['format'])

            dst = os.path.join(dst_dir, base)

            if os.path.exists(dst):

                m = 'Initial structure PDB-file "%s" does already ' % dst + \
                    'exist and will be overwritten.'
                self.warning(m)

            reader.write(chain, dst, 'iupac', 'cns')

    def getEnsemble(self, it, molecule):
        """
        returns a 'StructureEnsemble' instance.
        """

        ## TODO: maybe we should place this method somewhere
        ## else. on the one hand, the structure-ensemble could
        ## be rather specific to cns-generated pdb-files
        ## (e.g. parser for energies), on the other hand,
        ## iteration-specific parameters have to be set.

        check_type(it, 'IterationSettings')
        check_type(molecule, 'Molecule')

        from time import clock
        import aria.StructureEnsemble as SE
        from aria.FloatFile import FloatFile

        import os

        infra = self.getInfrastructure()

        pdb_path = infra.get_iteration_path(it['number'])
        file_root = infra.get_file_root()
        number = it['number_of_structures']

        pdbfiles = [os.path.join(pdb_path, file_root + '_%d.pdb' % i)
                    for i in range(1, number + 1)]

        floatfiles = [os.path.join(pdb_path, file_root + '_%d.float' % i)
                      for i in range(1, number + 1)]

        parser = FloatFile()

        if len(floatfiles):
            swapped_atoms = [parser.parse(file) for file in floatfiles]
        else:
            swapped_atoms = None

        ## read structure ensemble

        t = clock()

        ## create structure-ensemble
        ## settings will be set later, when ensemble is read

        es = SE.StructureEnsembleSettings(default_settings = it)
        ensemble = SE.StructureEnsemble(es)

        ensemble.read(pdbfiles, molecule, float_files = swapped_atoms)

        if len(floatfiles):

            n_swapped = [len(x) for x in swapped_atoms]
            self.message('%d PDB-files and float-files read.' \
                         % (len(pdbfiles)))

            [self.message('%s: %d atom pairs swapped' % (f, n))
             for f, n in zip(floatfiles, n_swapped)]

        else:
            self.message('%d PDB-files read.' % len(pdbfiles))

        self.debug('Time: %ss' % str(clock() - t))

        return ensemble

    def dumps_ariapeak_list(self, object):

        check_type(object, LIST, TUPLE)
        #check_elements(object, 'AriaPeak')
        check_elements(object, 'AbstractPeak')

        import numpy

        ambiguous = {}
        unambiguous = {}

        assign_template = 'assign %s %.3f %.3f %.3f weight %.3f'
        comment_template = '%s ! spec=%s, no=%d, id=%d, vol=%.6e'

        ## BARDIAUX 2.2
        ## multiple homlogous restraints written on the fly
        from aria.Singleton import ProjectSingleton
        project = ProjectSingleton()
        import aria.DataContainer as DC
        symmetry = project.getData(DC.DATA_SYMMETRY)[0]
        self.__nrestraints = 0

        if self.use_restraint_weights:

            weight_normaliser = [x.restraint_weight for x in object \
                                 if x.isActive()]

            weight_normaliser = numpy.sum(weight_normaliser) / \
                                len(weight_normaliser)

        else:
            weight_normaliser = 1.

        for peak in object:

            if not peak.isActive():
                continue

            ## encode contributions, i.e. create atom-selection
            ## statements

            first_selection, other_selections = \
                             self._create_peak_selections(peak)

            ## in case of empty atom selection, skip restraints
            if first_selection is None:
                continue

            ## get distance information

            d = peak.getDistance()

            if d is None:
                s = 'AriaPeak %d: distance is None.'
                self.error(ValueError, s % peak.getId())

            lower = peak.getLowerBound()

            if lower is None:
                s = 'AriaPeak %d: lower bound is None.'
                self.error(ValueError, s % peak.getId())
            else:
                lower = d - lower

            upper = peak.getUpperBound()

            if upper is None:
                s = 'AriaPeak %d: upper bound is None.'
                self.error(ValueError, s % peak.getId())
            else:
                upper -= d

            ## add distance information

            weight = peak.restraint_weight / weight_normaliser

            first_selection = assign_template % (first_selection,
                                                 d, lower, upper, weight)

            ## add comment

            peak_id = peak.getId()
            ref_xpk = peak.getReferencePeak()

            spec = ref_xpk.getSpectrum().getName()
            if len(spec) > 11:
                spec = spec[:4] + '...' + spec[-4:]

            peak_no = ref_xpk.getNumber()

            volume = ref_xpk.getVolume().getValue() or \
                     ref_xpk.getIntensity().getValue()

            first_selection = comment_template % \
                              (first_selection, spec, peak_no, peak_id, volume)

            assign_statement = '\n'.join([first_selection] + other_selections)
            self.__nrestraints += len([first_selection] + other_selections)


            if symmetry['enabled'] == YES:
                for n in range(0, symmetry['n_monomers'] -1):

                    first_homo_selection, other_homo_selections = \
                                          self._create_homo_peak_selections(peak, n)

                    if first_homo_selection is not None:

                        first_homo_selection = assign_template % (first_homo_selection,
                                                              d, lower, upper, weight)


                        first_homo_selection = comment_template % \
                                               (first_homo_selection, spec, peak_no, peak_id, volume)

                        assign_homo_statement = '\n'.join([first_homo_selection] + other_homo_selections)

                        assign_statement = assign_statement + '\n' + assign_homo_statement
                        self.__nrestraints += len([first_homo_selection] + other_homo_selections)


            ## depending on no. of active contributions,
            ## add assign-statement to ambig or unambig set

            key = '%s %6d' % (spec, peak_no)

            if peak.isAmbiguous():
                ambiguous[key] = assign_statement
            else:
                unambiguous[key] = assign_statement

        ## sort wrt spec and ref. peak number

        keys = unambiguous.keys()
        keys.sort()

        unambiguous = [unambiguous[key] for key in keys]

        keys = ambiguous.keys()
        keys.sort()

        ambiguous = [ambiguous[key] for key in keys]

        return '\n'.join(unambiguous) + '\n', \
               '\n'.join(ambiguous) + '\n'

    # OBSOLETE
    def dumps_restraint_list(self, object):

        check_type(object, LIST, TUPLE)
        #check_elements(object, 'AriaPeak')
        check_elements(object, 'DistanceRestraint') # BARDIAUX 2.2

        import numpy

        ambiguous = {}
        unambiguous = {}

        assign_template = 'assign %s %.3f %.3f %.3f weight %.3f'
        comment_template_restraint = '%s ! list=%s, no=%d, id=%d vol=%.6e' # BARDIAUX 2.2

        if self.use_restraint_weights:

            weight_normaliser = [x.restraint_weight for x in object \
                                 if x.isActive()]

            weight_normaliser = numpy.sum(weight_normaliser) / \
                                len(weight_normaliser)

        else:
            weight_normaliser = 1.

        for restraint in object:

            if not restraint.isActive():
                continue

            ## encode contributions, i.e. create atom-selection
            ## statements

            first_selection, other_selections = \
                             self._create_peak_selections(restraint)


            ## in case of empty atom selection, skip restraints
            if first_selection is None:
                continue

            ## get distance information

            d = restraint.getDistance()

            if d is None:
                s = 'DistanceRestraint %d: distance is None.'
                self.error(ValueError, s % restraint.getId())

            lower = restraint.getLowerBound()

            if lower is None:
                s = 'DistanceRestraint %d: lower bound is None.'
                self.error(ValueError, s % restraint.getId())
            else:
                lower = d - lower

            upper = restraint.getUpperBound()

            if upper is None:
                s = 'DistanceRestraint %d: upper bound is None.'
                self.error(ValueError, s % restraint.getId())
            else:
                upper -= d

            ## add distance information
            weight = restraint.restraint_weight / weight_normaliser

            first_selection = assign_template % (first_selection,
                                                 d, lower, upper, weight)


            ## BARDIAUX 2.2
            peak_id = restraint.getId()
            ref_xpk = restraint.getReferencePeak()


            list_name = ref_xpk.getSpectrum().getName()

            if len(list_name) > 11:
                 list_name= list_name[:4] + '...' + list_name[-4:]

            peak_no = ref_xpk.getNumber()

            volume = ref_xpk.getVolume().getValue() or \
                     ref_xpk.getIntensity().getValue()


            first_selection = comment_template_restraint % \
                              (first_selection, source, peak_no, peak_id, volume)


            assign_statement = '\n'.join([first_selection] + other_selections)
            self.__nrestraints += len([first_selection] + other_selections)


            ## depending on no. of active contributions,
            ## add assign-statement to ambig or unambig set

            key = '%s %6d' % (list_name, peak_no)

            if restraint.isAmbiguous():
                ambiguous[key] = assign_statement
            else:
                unambiguous[key] = assign_statement

        ## sort wrt spec and ref. peak number

        keys = unambiguous.keys()
        keys.sort()

        unambiguous = [unambiguous[key] for key in keys]


        keys = ambiguous.keys()
        keys.sort()

        ambiguous = [ambiguous[key] for key in keys]


        return '\n'.join(unambiguous) + '\n', \
               '\n'.join(ambiguous) + '\n'


    # BARDIAUX 2.2
    def write_peaklists(self, object, name_unambiguous, name_ambiguous):

        check_string(name_unambiguous)
        check_string(name_ambiguous)

        from aria.tools import cat_files
        from aria.DataContainer import DATA_AMBIGUOUS, DATA_UNAMBIGUOUS

        import os

        name_unambiguous = os.path.expanduser(name_unambiguous)
        name_ambiguous = os.path.expanduser(name_ambiguous)

        if os.path.exists(name_unambiguous):
            s = 'File %s does already exist and will be overwritten'
            self.warning(s % name_unambiguous)

        if os.path.exists(name_ambiguous):
            s = 'File %s does already exist and will be overwritten'
            self.warning(s % name_ambiguous)

        unambig, ambig = self.dumps_ariapeak_list(object)

        f = open(name_unambiguous, 'w')
        f.write(unambig)
        f.close()

        f = open(name_ambiguous, 'w')
        f.write(ambig)
        f.close()

        ## cat distance restraints

        dirs = self.getInfrastructure().get_cns_data_directories()

        distances = os.path.join(dirs[DATA_UNAMBIGUOUS], NAME_NOE_RESTRAINTS_UNAMBIG)

        if os.path.exists(distances):
            cat_files([name_unambiguous, distances], name_unambiguous)

        distances = os.path.join(dirs[DATA_AMBIGUOUS], NAME_NOE_RESTRAINTS_AMBIG)

        if os.path.exists(distances):
            cat_files([name_ambiguous, distances], name_ambiguous)

    def write_molecule(self, molecule, filename):

        check_type(molecule, 'Molecule')
        check_string(filename)

        import os
        from aria.OrderedDict import OrderedDict

        reader = self.getPDBReader()

        filename = os.path.expanduser(filename)

        if os.path.exists(filename):
            s = 'File "%s" does already exist and will be overwritten.'
            self.warning(s % filename)

        #molecule_dict = {}
        # BARDIAUX 2.2 keep chain order
        molecule_dict = OrderedDict()

        for chain in molecule.get_chains():

            molecule_dict[chain.getSegid()] = {}

            chain_dict = molecule_dict[chain.getSegid()]
            chain_dict['chain_type'] = chain.getSettings()['type']

            for residue in chain.getResidues():

                residue_dict = {}

                residue_dict['residue_type'] = residue.getType()

                for atom in residue.getAtoms():

                    residue_dict[atom.getName()] = [0., 0., 0.]

                chain_dict[residue.getNumber()] = residue_dict

        reader.write(molecule_dict, filename, 'iupac', 'cns')

    def _params_as_cns_dict(self, params):
        """
        Helper function
        """
        from aria.DataContainer import DATA_ANNEALING_UNAMBIG, \
             DATA_ANNEALING_AMBIG, DATA_ANNEALING_HBOND, \
             DATA_ANNEALING_DIHEDRAL, DATA_ANNEALING_KARPLUS, \
             DATA_ANNEALING_RDC, DATA_ANNEALING_FBHW, \
             DATA_DYNAMICS, DATA_ANNEALING_SYM, DATA_ANNEALING_LOGHARMONIC

        params_type = params.getType()
        dict = {}

        if params_type == DATA_ANNEALING_UNAMBIG:
            for key, value in params.items():
                dict['unambig_%s' % key] = value

        elif params_type == DATA_ANNEALING_AMBIG:
            for key, value in params.items():
                dict['ambig_%s' % key] = value

        elif params_type == DATA_ANNEALING_HBOND:
            for key, value in params.items():
                dict['hbond_%s' % key] = value

        elif params_type == DATA_ANNEALING_DIHEDRAL:
            for key, value in params.items():
                dict['dihedral_%s' % key] = value

        elif params_type == DATA_ANNEALING_FBHW:
            for key, value in params.items():
                dict['fbhw_%s' % key] = value

        elif params_type == DATA_DYNAMICS:
            for key, value in params.items():
                dict[key] = value

        elif params_type == DATA_ANNEALING_KARPLUS:
            c = params['class']
            for key, value in params.items():
                dict['c%d_%s' % (c, key)] = value

        elif params_type == DATA_ANNEALING_RDC:
            c = params['class']
            for key, value in params.items():
                dict['rdc%d_%s' % (c, key)] = value

        # BARDIAUX 2.2
        elif params_type == DATA_ANNEALING_SYM:
            for key, value in params.items():
                dict[key] = value

        # BERNARD 2.3
        elif params_type == DATA_ANNEALING_LOGHARMONIC:
            d = {YES: CNS_TRUE, NO: CNS_FALSE}
            for key, value in params.items():
                if key in ['enabled', 'use_auto_weight']:
                    value = d[value]
                dict["logharmonic_%s" % key] = value
                
        return dict

    def _update_run_cns_dict(self, dict):
        """
        Adds data-specific parameters to 'dict'.
        """

        import aria.DataContainer as DC
        import os

        params = self.getAnnealingParameters().items()
        params.append((None, self.getMDParameters()))

        ## update parameter settings in cns-dict

        for p_type, p in params:

            if is_type(p, 'DataContainer'):
                dict.update(self._params_as_cns_dict(p))

            elif is_type(p, DICT):

                ## TODO: hard-coded: 5 classes

                for i in range(1, 6):

                    if p.has_key(i):
                        pp = p[i]

                    else:
                        if p_type == DC.DATA_ANNEALING_KARPLUS:
                            pp = DC.KarplusParameters()
                        elif p_type == DC.DATA_ANNEALING_RDC:
                            pp = DC.RDCParameters()
                        else:
                            self.error('Something went wrong.')

                        pp.reset()
                        pp['class'] = i

                    dict.update(self._params_as_cns_dict(pp))

##         ## BARDIAUX 2.2
##         sym_enabled = ['packing_enabled','ncs_enabled']
##         for i in dimers_enabled:
##             if dict[i] == YES :
##                 dict[i] = CNS_TRUE
##             else:
##                 dict[i] = CNS_FALSE

        ## set enabled flags according to wether the local CNS copy
        ## of a data file exists or not

        dirs = self.getInfrastructure().get_cns_data_directories()

        files = {DC.DATA_HBONDS:
                 (('hbonds.tbl', 'hbond_enabled'),
                  ('hbonds_csi.tbl', 'hbond_csi_enabled')),

                 DC.DATA_SSBONDS:
                 (('ssbonds.tbl', 'ssbonds_enabled'),),

                 DC.DATA_DIHEDRALS:
                 (('dihedrals.tbl', 'dihedral_enabled'),
                  ('dihedrals_csi.tbl', 'dihedral_csi_enabled'),
                  ('dihedrals_talos.tbl', 'dihedral_talos_enabled')),

                 DC.DATA_KARPLUS:
                 (('c1.tbl', 'c1_enabled'), ('c2.tbl', 'c2_enabled'),
                  ('c3.tbl', 'c3_enabled'), ('c4.tbl', 'c4_enabled'),
                  ('c5.tbl', 'c5_enabled')),

                 DC.DATA_RDCS:
                 (('rdc1.tbl', 'rdc1_enabled'), ('rdc2.tbl', 'rdc2_enabled'),
                  ('rdc3.tbl', 'rdc3_enabled'), ('rdc4.tbl', 'rdc4_enabled'),
                  ('rdc5.tbl', 'rdc5_enabled'))}

        for data_type in files.keys():

            dir = dirs[data_type]

            for file, flag in files[data_type]:
                if os.path.exists(os.path.join(dir, file)):
                    dict[flag] = CNS_TRUE
                else:
                    dict[flag] = CNS_FALSE

        return dict

    def _ssbridge_to_cns(self, ssbridge, number):
        """
        Helper function that converts an SSBridge-instance into
        the corresponding CNS commands for defining this bridge.
        """

        check_type(ssbridge, 'SSBridge')
        check_int(number)

        from aria.tools import string_to_segid

        cns_command = 'evaluate (&toppar.ss_i_resid_%(n)d=%(residue1)d)\n' + \
                      'evaluate (&toppar.ss_i_segid_%(n)d="%(segid1)s")\n' + \
                      'evaluate (&toppar.ss_j_resid_%(n)d=%(residue2)d)\n' + \
                      'evaluate (&toppar.ss_j_segid_%(n)d="%(segid2)s")\n'

        d = {'n': number + 1,
             'segid1': string_to_segid(ssbridge['segid1']),
             'segid2': string_to_segid(ssbridge['segid2']),
             'residue1': ssbridge['residue1'],
             'residue2': ssbridge['residue2']}

        return cns_command % d

    # BARDIAUX 2.2
    def _cyspatch_to_cns(self, cyspatch1, cyspatch2, number):
        """
        Helper function that converts an cyspatch-instance into
        the corresponding CNS commands for patching cysteins with
        the DISN patch (HG is removed and charges change)
        """

        check_type(cyspatch1, 'CysPatch')
        check_type(cyspatch2, 'CysPatch')
        check_int(number)

        from aria.tools import string_to_segid

        cns_command = 'evaluate (&toppar.cyspatch_i_resid_%(n)d=%(residue1)d)\n' + \
                      'evaluate (&toppar.cyspatch_i_segid_%(n)d="%(segid1)s")\n' +\
                      'evaluate (&toppar.cyspatch_j_resid_%(n)d=%(residue2)d)\n' + \
                      'evaluate (&toppar.cyspatch_j_segid_%(n)d="%(segid2)s")\n'

        d = {'n': number + 1,
             'segid1': string_to_segid(cyspatch1['segid']),
             'residue1': cyspatch1['residue'],
             'segid2': string_to_segid(cyspatch2['segid']),
             'residue2': cyspatch2['residue']}

        return cns_command % d

    def _hispatch_to_cns(self, hispatch, number):

        from aria.tools import string_to_segid

        check_type(hispatch, 'HisPatch')
        check_int(number)

        cns_command = 'evaluate (&toppar.%(proton)s_resid_%(n)d=' + \
                      '%(residue)d)\n' + \
                      'evaluate (&toppar.%(proton)s_segid_%(n)d=' + \
                      '"%(segid)s")\n'

        d = {'n': number + 1, 'residue': hispatch['residue'],
             'segid': string_to_segid(hispatch['segid']),
             'proton': hispatch['proton'].lower()}

        return cns_command % d

    # BARDIAUX 2.2
    def _cispropatch_to_cns(self, cispropatch, number):

        from aria.tools import string_to_segid

        check_type(cispropatch, 'CisProPatch')
        check_int(number)

        cns_command = 'evaluate (&toppar.cispro_resid_%(n)d=' + \
                      '%(residue)d)\n' + \
                      'evaluate (&toppar.cispro_segid_%(n)d=' + \
                      '"%(segid)s")\n'

        d = {'n': number + 1, 'residue': cispropatch['residue'],
             'segid': string_to_segid(cispropatch['segid'])}

        return cns_command % d

    # BARDIAUX 2.3
    def _znpatch_to_cns(self, znpatch, number):

        from aria.tools import string_to_segid

        check_type(znpatch, 'ZnPatch')
        check_int(number)

        cns_command = 'evaluate (&toppar.%(type)s_1_resid_%(n)d=%(resid1)d)\n' + \
                      'evaluate (&toppar.%(type)s_1_segid_%(n)d="%(segid1)s")\n'+ \
                      'evaluate (&toppar.%(type)s_2_resid_%(n)d=%(resid2)d)\n' + \
                      'evaluate (&toppar.%(type)s_2_segid_%(n)d="%(segid2)s")\n'+ \
                      'evaluate (&toppar.%(type)s_3_resid_%(n)d=%(resid3)d)\n' + \
                      'evaluate (&toppar.%(type)s_3_segid_%(n)d="%(segid3)s")\n'+ \
                      'evaluate (&toppar.%(type)s_4_resid_%(n)d=%(resid4)d)\n' + \
                      'evaluate (&toppar.%(type)s_4_segid_%(n)d="%(segid4)s")\n'+ \
                      'evaluate (&toppar.%(type)s_5_resid_%(n)d=%(resid_zn)d)\n' + \
                      'evaluate (&toppar.%(type)s_5_segid_%(n)d="%(segid_zn)s")\n'

        d = {'n': number + 1, 'type' : znpatch['type'].lower(),
             'resid1': znpatch['residue1'],
             'segid1': string_to_segid(znpatch['segid1']),
             'resid2': znpatch['residue2'],
             'segid2': string_to_segid(znpatch['segid2']),
             'resid3': znpatch['residue3'],
             'segid3': string_to_segid(znpatch['segid3']),
             'resid4': znpatch['residue4'],
             'segid4': string_to_segid(znpatch['segid4']),
             'resid_zn': znpatch['residue_zn'],
             'segid_zn': string_to_segid(znpatch['segid_zn'])}

        return cns_command % d

    # BARDIAUX 2.3
    def _dnapatch_to_cns(self, dnapatch, number):

        from aria.tools import string_to_segid

        check_int(number)

        cns_command = 'evaluate (&toppar.dna_segid_%(n)d="%(segid)s")'

        d = {'n': number + 1, 'segid' : string_to_segid(dnapatch.getSegid())}

        return cns_command % d

    def _update_chain_data(self, chain_data):
        """
        Augments run.cns dict with chain-specific parameters
        """

        import aria.DataContainer as DC
        import aria.Chain as Chain
        import os
        from aria.Singleton import ProjectSingleton

        ## TODO: clumsy. Up to now, we do not
        ## have a rule how/where to store the
        ## actual data (not the DataSources -
        ## they can be accessed vio project.getData)

        project = ProjectSingleton()
        molecule = project.getMolecule()

        ## Default Linkage/Topology files

        d_linkage = {Chain.TYPE_PROTEIN: 'topallhdg5.3.pep',
                     Chain.TYPE_RNA: 'dna-rna.link',
                     Chain.TYPE_DNA: 'dna-rna.link',
                     Chain.TYPE_NONPOLYMER: 'ions.link'}

        d_topology = {Chain.TYPE_PROTEIN: 'topallhdg5.3.pro',
                      Chain.TYPE_RNA: 'dna-rna-allatom.top',
                      Chain.TYPE_DNA: 'dna-rna-allatom.top',
                      Chain.TYPE_NONPOLYMER: None}

        d_parameters = {Chain.TYPE_PROTEIN: 'parallhdg5.3.pro',
                        Chain.TYPE_RNA: 'dna-rna-allatom.param',
                        Chain.TYPE_DNA: 'dna-rna-allatom.param',
                        Chain.TYPE_NONPOLYMER: None}

        ## Linkage, topology, parameters definition file

        dicts = {'linkage': d_linkage,
                 'topology': d_topology,
                 'parameter': d_parameters}

        run_cns_keys = {}

        for base, dict in dicts.items():

            ## Compile entity names and
            ## key that is used in run.cns

            name = chain_data[base + '_name']
            filename = base + '_filename'
            key = base + '_file'
            run_cns_keys[key + '_1'] = ""
            run_cns_keys[key + '_2'] = ""

            if name == DC.SequenceData.USER_DEFINED:
                filename = chain_data[filename]

            ## If automatic: linkage-filename depends on chain type
            ## BARDIAUX 2.2 extend
            ## upgraded for hetero-complex
            elif name == DC.SequenceData.AUTOMATIC:
                filenames = []
                for a in [dict[a] for a in molecule.getTypes()]:
                    if a and a not in filenames:
                        filenames.append(a)
                # apparently dna_rna.link should be first...
                filenames.sort()
                if len(filenames) > 1:
                    run_cns_keys[key + '_2'] = os.path.basename(filenames[1])
                filename  = filenames[0]

            else:
                filename = name

            run_cns_keys[key + '_1'] = os.path.basename(filename)


        ## TODO: hack => BARDIAUX 2.2 : Done
        ##run_cns_keys['segid'] = molecule.get_chains()[0].getSegid()
        ## BARDIAUX 2.2 , segid value for run.cns

        segid_key = ['segid_1','segid_2', 'segid_3' ,\
                     'segid_4', 'segid_5']
        chains_key = [ch.getSegid() for ch in molecule.get_chains()]
        empty = len(segid_key) - len(chains_key)

        for u in segid_key:
            run_cns_keys[u] = ""

        for i in xrange(len(chains_key)):
            run_cns_keys[segid_key[i]] = chains_key[i]

        # BARDIAUX mononum 2.2
        # then, allow pack restraint for hetero dimers
        mononum = len(molecule.get_chains())
        run_cns_keys["mononum"] = mononum

        # BARDIAUX DNA chains
        dna = [c for c in molecule.get_chains() if c.getType() == Chain.TYPE_DNA]
        dna_chains = [self._dnapatch_to_cns(dna[i], i) for i in range(len(dna))]
        run_cns_keys['dna_chains'] = '\n'.join(dna_chains)
        run_cns_keys['number_of_dna_chains'] = len(dna)

        return run_cns_keys

    def _write_run_cns(self, cns_working_dir, cns_out_dir,
                       iteration_number = None, counter = 1,
                       initial_structure = None, protocol_settings = None):
        """
        Private method for preparation of the run.cns script from
        its template.
        """
	
        check_string(cns_working_dir)
        check_string(cns_out_dir)
        check_int(counter)
        check_type(initial_structure, NONE, STRING)
        check_type(protocol_settings, 'ProtocolSettings')

        from aria.Singleton import ProjectSingleton
        import aria.DataContainer as DC

        import os

        project = ProjectSingleton()
        infra = self.getInfrastructure()
        cns_data_dirs = infra.get_cns_data_directories()

        ## create dict to patch run.cns template

        if iteration_number is None:
            iteration_number = 0

        initial_pdbs = project.getData(DC.DATA_INITIAL_STRUCTURE)
        initial_pdbs = [s for s in initial_pdbs if s['enabled'] == YES]

        ## If user has not specified an initial structure, we use
        ## the standard template file which has been generated
        ## by the 'generate.inp' script.

        if not len(initial_pdbs):
            initial_pdb = '%s_template.pdb' % infra.get_file_root()
        else:
            ## Should be only 1 structure anyway.
            initial_pdb = os.path.basename(initial_pdbs[0]['filename'])

        filename = os.path.join(cns_data_dirs[DC.DATA_INITIAL_STRUCTURE],
                                initial_pdb)

        ## TODO: to specify the initial PDB file explictely might
        ## not be the most elegant way but is necessary for water
        ## refinement.

        ## Main chain dihedral angles are randomized only if
        ## we start from an extended state.

        if initial_structure is None:
            initial_structure = filename
            randomize_mainchain = CNS_TRUE
            use_template = CNS_TRUE

        elif not os.path.exists(initial_structure):
            m = 'Initial PDB-file "%s" does not exist.' % initial_structure
            self.error(ValueError, m)

        else:
            randomize_mainchain = CNS_FALSE
            use_template = CNS_FALSE
            s = 'Using %s as initial structure.'
            self.message(s % initial_structure)

        if protocol_settings:
            x = protocol_settings['floating_assignment']
            floating_assignment = {YES: CNS_TRUE, NO: CNS_FALSE}[x]
        else:
            floating_assignment = CNS_FALSE



        ## BARDIAUX 2.2
        ## With the symmetry (multimers) number of restraints
        ## in multiply by the number of conformers.
        ## so noe nrestraints is allocated on the fly
        ## we add 20.000 for the distances restraints added by user.
        n_rest = self.__nrestraints + 20000

        cns_settings = self.getSettings()
        ## < Mareuil
        run_cns_dict = {'run_path': infra.get_run_path(),
                        'psf_file': infra.get_psf_file(),
                        'sequence_pdb': '%s.pdb' % infra.get_file_root(),
                        'file_root': infra.get_file_root(),
                        'initial_pdb': os.path.basename(initial_structure),
                        'out_dir': cns_out_dir,
                        'iteration': iteration_number,
                        'randomize_mainchain': randomize_mainchain,
                        'counter': counter,
                        'floating_assignment': floating_assignment,
                        'nonbonded_parameters': \
                        cns_settings['nonbonded_parameters'],
                        'nrestraints' : n_rest,
                        'use_template' : use_template}
        ## Mareuil >
        ## for all data types: add specific parameters to dict
        self._update_run_cns_dict(run_cns_dict)

        ## handle ss-bridges
        data = project.getData(DC.DATA_SSBRIDGE)

        run_cns_dict['number_of_ssbridges'] = len(data)
        ssbridges = map(lambda d, n, f = self._ssbridge_to_cns: f(d, n),
                        data, range(len(data)))
        run_cns_dict['ssbridges'] = '\n'.join(ssbridges)

        ## handle ambiguous ss-bonds/cyspatch
        ## BARDIAUX 2.2
        run_cns_dict['ssbonds_ambig'] = "unambiguous"

        data = []
        for d in project.getData(DC.DATA_SSBONDS):
            for c in d['cyspatch']:
                data.append(c)

        cyspatches = []
        if len(data) % 2 <> 0:
            cyspatches = [self._cyspatch_to_cns(data[0], data[-1], 0)]
            data = data[1:]

        for k in range(0, len(data), 2):
            cyspatches.append(self._cyspatch_to_cns(data[k], data[k+1], len(cyspatches)))

        run_cns_dict['number_of_cyspatches'] = len(cyspatches)
        run_cns_dict['cyspatches'] = '\n'.join(cyspatches)

        if len(cyspatches):
            run_cns_dict['ssbonds_ambig'] = "ambiguous"

        ## handle his-patches
        data = project.getData(DC.DATA_HISPATCH)

        hisd_patches = [d for d in data if d['proton'] == 'HISD']
        run_cns_dict['number_of_hisd_patches'] = len(hisd_patches)
        hisd_patches = [self._hispatch_to_cns(hisd_patches[i], i) for \
                        i in range(len(hisd_patches))]
        run_cns_dict['hisd_patches'] = '\n'.join(hisd_patches)

        hise_patches = [d for d in data if d['proton'] == 'HISE']
        run_cns_dict['number_of_hise_patches'] = len(hise_patches)
        hise_patches = [self._hispatch_to_cns(hise_patches[i], i) for \
                        i in range(len(hise_patches))]
        run_cns_dict['hise_patches'] = '\n'.join(hise_patches)

        ## handle cispro-patches
        ## BARDIAUX 2.2
        data = project.getData(DC.DATA_CISPROPATCH)

        run_cns_dict['number_of_cispro_patches'] = len(data)
        cispro_patches = [self._cispropatch_to_cns(data[i], i) for \
                        i in range(len(data))]
        run_cns_dict['cispro_patches'] = '\n'.join(cispro_patches)

        ## handle zn-patches
        ## BARDIAUX 2.3
        data = project.getData(DC.DATA_ZNPATCH)
        ssss_patches = [d for d in data if d['type'] == 'SSSS']
        ssse_patches = [d for d in data if d['type'] == 'SSSE']
        sssd_patches = [d for d in data if d['type'] == 'SSSD']

        run_cns_dict['number_of_ssss_patches'] = len(ssss_patches)
        run_cns_dict['number_of_ssse_patches'] = len(ssse_patches)
        run_cns_dict['number_of_sssd_patches'] = len(sssd_patches)

        zn_patches = [self._znpatch_to_cns(ssss_patches[i], i) for \
                        i in range(len(ssss_patches))]
        zn_patches += [self._znpatch_to_cns(ssse_patches[i], i) for \
                        i in range(len(ssse_patches))]
        zn_patches += [self._znpatch_to_cns(sssd_patches[i], i) for \
                        i in range(len(sssd_patches))]
        run_cns_dict['zn_patches'] = '\n'.join(zn_patches)

        ## BARDIAUX : symmetry parameters statement 2.2
        d = {"None": 1, "C2" : 2, "C3" : 3, "D2" : 222, "C5" : 5}
        symmetry = project.getData(DC.DATA_SYMMETRY)[0]

        run_cns_dict['symmetry'] = d[symmetry['symmetry_type']]

        d = {YES: CNS_TRUE, NO: CNS_FALSE}
        ncs = symmetry['ncs_enabled']
        run_cns_dict['ncs_enabled'] = d[ncs]
        pack = symmetry['packing_enabled']
        run_cns_dict['packing_enabled'] = d[pack]

        ## solvent refinement

        solvent_params = protocol_settings['water_refinement']

        d = {YES: CNS_TRUE, NO: CNS_FALSE}

        keep_mol = solvent_params['write_solvent_molecules']
        run_cns_dict['write_solvent_molecules'] = d[keep_mol]

        ## Molecule-specific parameters

        chain_data = project.getData(DC.DATA_SEQUENCE)[0]
        run_cns_dict.update(self._update_chain_data(chain_data))

        ## patch and store local run.cns
        if not os.path.exists(cns_working_dir):
            try:
                os.makedirs(cns_working_dir)
            except:

                from aria.Infrastructure import AriaDirectoryCreationError

                m = 'Could not create working directory ' + \
                    '"%s" for current CNS job.' % cns_working_dir
                self.error(AriaDirectoryCreationError, m)

        template = os.path.join(infra.get_cns_protocols_path(), RUN_CNS)
        file = open(template)
        script = file.read() % run_cns_dict
        file.close()

        run_cns = os.path.join(cns_working_dir, RUN_CNS)
        file = open(run_cns, 'w')
        file.write(script)
        file.close()
        ## < Mareuil
        return initial_structure
        ## Mareuil >

    def _write_condor_script(self, parameters, condor_append = None):

        """
        creates script that is used by condor to launch a CNS calculation
        """

        check_dict(parameters)

        import os

        script = CONDOR_SCRIPT % parameters
        script += condor_append % parameters

        ## csh-script will be stored in same directory as run.cns

        cns_working_dir = parameters['cns_working_dir']

        script_name = os.path.join(cns_working_dir, 'condor.job')

        f = open(script_name, 'w')
        f.write(script)
        f.close()
    ## < Mareuil
    
    def _write_csh_script_async(self, cns_working_dir, iteration_path, cns_script, job_name,
                              nb_structs, n_cpu, use_default_executable, initial_structure = None, Mode = None):
        
    ## Mareuil >
        """
        if use_default_executable is false, the executable will be
        patched into the
        csh-script later (see JobManager). if it is true, we
        use cns' default_executable instead. it has to be handled
        in that way, since cns-jobs are executed on different
        machines, each having its individual cns binary.
        """

        check_string(cns_working_dir)
        check_string(iteration_path)
        check_string(cns_script)
        import os
	## < Mareuil 0.1
        initial_path_structure = os.path.dirname(initial_structure)
        ## Mareuil 0.1 >

        infra = self.getInfrastructure()
        cns_input = os.path.join(infra.get_cns_protocols_path(),
                                 '%s.inp' % cns_script)

        ## name of local output file
        cns_output = '%s.out' % os.path.basename(cns_script)

        if use_default_executable or self.use_condor:
            cns_executable = self.getSettings()['local_executable']
            cns_grid = '/' + os.path.basename(cns_executable)
        else:
            cns_executable = '%(executable)s'
            cns_grid = '/%(executable)s'
        ## < Mareuil
        script_name = os.path.basename(cns_script)
        script_name = os.path.join(cns_working_dir, '%s.csh' % script_name)
        ## create different csh scripts depending on whether we shall
	
	if Mode == None:
            name_dir = job_name + "_"
	    cns_working_dir_sge = os.path.join(infra.get_temp_path(), name_dir)
            script_sge = os.path.join(infra.get_temp_path(), '%s.csh' % os.path.basename(cns_script))
            d = {'cns_input_file': cns_input,
                 'cns_output_file': cns_output,
                 'cns_executable': cns_executable,
                 'cns_working_dir': cns_working_dir,
		 'cns_working_dir_sge': cns_working_dir_sge,
                 'initial_path_pdb': initial_path_structure,
                 'cns_resource_path': infra.get_cns_path(),
                 'iteration_path': iteration_path,
                 'analysis_output_path': ''}
            # some statements for SGE
            # job_name for queue (not mandatory but usefull)
            d['sge_job_name'] = job_name
            d['sge_job_shell'] = '/bin/csh'
            d['sge_job_num'] = nb_structs
            SGE = CSH_SCRIPT_SGE % d
            f = open(script_sge, 'w')
            f.write(SGE)
            f.close()
	    
        elif Mode[0] == "GRID":
            script_jdl = os.path.basename(cns_script)
            script_jdl = os.path.join(cns_working_dir, '%s.jdl' % script_jdl)
            input_tar_grid = "%s_%s.tar.gz" % (os.path.basename(infra.get_run_path()),os.path.basename(cns_working_dir))
            run_path = os.path.split(os.path.split(cns_working_dir)[0])[0]
            input_tar = os.path.join(run_path, input_tar_grid)
            d = {'input_tar_dir': input_tar,
                 'input_tar_grid': input_tar_grid,
                 'csh_script_grid': os.path.basename(script_name),
                 'csh_working_dir': script_name,
                 'cns_input_file': self.revpath(cns_input, run_path),
                 'initial_path_pdb': self.revpath(initial_path_structure, run_path),
                 'cns_output_file': cns_output,
                 'cns_executable': cns_executable,
                 'cns_executable_grid': cns_grid,
                 'cns_working_dir': self.revpath(cns_working_dir, run_path),
                 'cns_resource_path': self.revpath(infra.get_cns_path(), run_path),
                 'run_resource_path': self.revpath(infra.get_run_path(), run_path),
                 'iteration_path': self.revpath(iteration_path, run_path),
                 'analysis_output_path': ''}
            JDL = GRID_JDL_SCRIPT % d
            script = GRID_CSH_SCRIPT % d
            g = open(script_jdl, 'w')
            g.write(JDL)
            g.close()
            f = open(script_name, 'w')
            f.write(script)
            f.close()  
            
        elif Mode[0] == "LOCAL" or Mode[0] == "CLUSTER":
            name_dir = job_name + "_"
	    cns_working_dir_sge = os.path.join(infra.get_temp_path(), name_dir)
            script_sge = os.path.join(infra.get_temp_path(), '%s.csh' % os.path.basename(cns_script))
            d = {'cns_input_file': cns_input,
                 'cns_output_file': cns_output,
                 'cns_executable': cns_executable,
                 'cns_working_dir': cns_working_dir,
		 'cns_working_dir_sge': cns_working_dir_sge,
                 'initial_path_pdb': initial_path_structure,
                 'cns_resource_path': infra.get_cns_path(),
                 'iteration_path': iteration_path,
                 'analysis_output_path': ''}
            # some statements for SGE
            # job_name for queue (not mandatory but usefull)
            d['sge_job_name'] = job_name
            d['sge_job_shell'] = '/bin/csh'
            d['sge_job_num'] = nb_structs
	    d['sge_cpu_num'] = n_cpu
            SGE = CSH_SCRIPT_SGE % d
            f = open(script_sge, 'w')
            f.write(SGE)
            f.close()
        ## csh-script will be stored in same directory as run.cns
        return script_sge    
    
    def _write_csh_script_sync(self, cns_working_dir, iteration_path, cns_script,
                          use_default_executable, initial_structure = None, Mode = None):
    ## Mareuil >
        """
        if use_default_executable is false, the executable will be
        patched into the
        csh-script later (see JobManager). if it is true, we
        use cns' default_executable instead. it has to be handled
        in that way, since cns-jobs are executed on different
        machines, each having its individual cns binary.
        """
        
        check_string(cns_working_dir)
        check_string(iteration_path)
        check_string(cns_script)

        import os
	## < Mareuil 0.1
        initial_path_structure = os.path.dirname(initial_structure)
        ## Mareuil 0.1 >

        infra = self.getInfrastructure()

        cns_input = os.path.join(infra.get_cns_protocols_path(),
                                 '%s.inp' % cns_script)

        ## name of local output file
        cns_output = '%s.out' % os.path.basename(cns_script)

        if use_default_executable or self.use_condor:
            cns_executable = self.getSettings()['local_executable']
            cns_grid = '/' + os.path.basename(cns_executable)
        else:
            cns_executable = '%(executable)s'
            cns_grid = '/%(executable)s'
        ## < Mareuil
        script_name = os.path.basename(cns_script)
        script_name = os.path.join(cns_working_dir, '%s.csh' % script_name)
        ## create different csh scripts depending on whether we shall
        ## use condor or ssh or glite to dispatch jobs
        if self.use_condor:
            d = {'cns_input_file': cns_input,
                 'cns_output_file': cns_output,
                 'cns_executable': cns_executable,
                 'cns_working_dir': cns_working_dir,
                 'initial_path_pdb': initial_path_structure,
                 'cns_resource_path': infra.get_cns_path(),
                 'iteration_path': iteration_path,
                 'analysis_output_path': ''}
            condor_queue = "input = %(cns_input_file)s\noutput = %(cns_output_file)s\nqueue\n"
            self._write_condor_script(d, condor_append = condor_queue)
            script = CSH_SCRIPT_REFINE_CONDOR % d
            f = open(script_name, 'w')
            f.write(script)
            f.close()
	    
        if Mode == None:
            d = {'cns_input_file': cns_input,
                'cns_output_file': cns_output,
                'cns_executable': cns_executable,
                'cns_working_dir': cns_working_dir,
                'initial_path_pdb': initial_path_structure,
                'cns_resource_path': infra.get_cns_path(),
                'iteration_path': iteration_path,
                'analysis_output_path': ''}
            # some statements for SGE
            # job_name for queue (not mandatory but usefull)
            sge_job_name = os.path.basename(cns_working_dir)
            d['sge_job_name'] = sge_job_name
            d['sge_job_shell'] = '/bin/csh'
	    
            script = CSH_SCRIPT_REFINE % d
            f = open(script_name, 'w')
            f.write(script)
            f.close()
	    
        elif Mode[0] == "GRID":
            script_jdl = os.path.basename(cns_script)
            script_jdl = os.path.join(cns_working_dir, '%s.jdl' % script_jdl)
            input_tar_grid = "%s_%s.tar.gz" % (os.path.basename(infra.get_run_path()),os.path.basename(cns_working_dir))
            run_path = os.path.split(os.path.split(cns_working_dir)[0])[0]
            input_tar = os.path.join(run_path, input_tar_grid)
            d = {'input_tar_dir': input_tar,
                 'input_tar_grid': input_tar_grid,
                 'csh_script_grid': os.path.basename(script_name),
                 'csh_working_dir': script_name,
                 'cns_input_file': self.revpath(cns_input, run_path),
                 'initial_path_pdb': self.revpath(initial_path_structure, run_path),
                 'cns_output_file': cns_output,
                 'cns_executable': cns_executable,
                 'cns_executable_grid': cns_grid,
                 'cns_working_dir': self.revpath(cns_working_dir, run_path),
                 'cns_resource_path': self.revpath(infra.get_cns_path(), run_path),
                 'run_resource_path': self.revpath(infra.get_run_path(), run_path),
                 'iteration_path': self.revpath(iteration_path, run_path),
                 'analysis_output_path': ''}
            JDL = GRID_JDL_SCRIPT % d
            script = GRID_CSH_SCRIPT % d
            g = open(script_jdl, 'w')
            g.write(JDL)
            g.close()
            f = open(script_name, 'w')
            f.write(script)
            f.close()
            
        elif Mode[0] == "LOCAL" or Mode[0] == "CLUSTER":
            d = {'cns_input_file': cns_input,
                 'cns_output_file': cns_output,
                 'cns_executable': cns_executable,
                 'cns_working_dir': cns_working_dir,
                 'initial_path_pdb': initial_path_structure,
                 'cns_resource_path': infra.get_cns_path(),
                 'iteration_path': iteration_path,
                 'analysis_output_path': ''}
            # some statements for SGE
            # job_name for queue (not mandatory but usefull)
            sge_job_name = os.path.basename(cns_working_dir)
            d['sge_job_name'] = sge_job_name
            d['sge_job_shell'] = '/bin/csh'
	    
            script = CSH_SCRIPT_REFINE % d
            f = open(script_name, 'w')
            f.write(script)
            f.close()
        ## csh-script will be stored in same directory as run.cns
        return script_name

    def header_cns(self, cns_wdir, iteration_path, iteration_number = None,
                   counter = 1, initial_structure = None, protocol_settings = None):
		   
        check_type(protocol_settings, NONE, 'ProtocolSettings')
        check_string(cns_wdir)
        check_string(iteration_path)
        check_type(iteration_number, INT, NONE)
        check_int(counter)
	
        return self._write_run_cns(cns_wdir, iteration_path, iteration_number,
                            counter, initial_structure, protocol_settings)	   
		   
    def setup_job_csh(self, cns_wdir, iteration_path, cns_protocol, job_name,
                  nb_structs = 1, n_cpu = 1, use_default_executable = 1, pdb_path = None, Mode = None):

        check_string(cns_wdir)
        check_string(iteration_path)
        check_string(cns_protocol)
        check_int(nb_structs)
			    
        return self._write_csh_script_async(cns_wdir, iteration_path, cns_protocol, job_name,
                                      nb_structs, n_cpu, use_default_executable, pdb_path, Mode)				      
				      
    def setup_job(self, cns_wdir, iteration_path, cns_protocol,
                  iteration_number = None, counter = 1,
                  initial_structure = None, use_default_executable = 1,
                  protocol_settings = None, Mode = None):

        check_type(protocol_settings, NONE, 'ProtocolSettings')
        check_string(cns_wdir)
        check_string(iteration_path)
        check_string(cns_protocol)
        check_type(iteration_number, INT, NONE)
        check_int(counter)
	
        pdb_path = self._write_run_cns(cns_wdir, iteration_path, iteration_number,
                            counter, initial_structure, protocol_settings)
			    
        return self._write_csh_script_sync(cns_wdir, iteration_path, cns_protocol,
                                      use_default_executable, pdb_path, Mode)
				      
    ## Mareuil >
    def create_psf_file(self, protocol_settings, replace):

        check_type(protocol_settings, 'ProtocolSettings')

        import os
        import aria.DataContainer as DC

        infra = self.getInfrastructure()

        if not replace and self.psf_file_exists():
            s = 'PSF-file already exists. Using existing file: %s'
            self.message(s % infra.get_psf_file())

            return

        ## check dependencies
        ## check sequence PDB file

        base = infra.get_file_root() + '.pdb'
        sequence_pdb_name = self.get_local_filename(base, DC.DATA_SEQUENCE)

        if not os.path.exists(sequence_pdb_name):
            s = 'Could not create PSF file: ' +\
                'Template PDB file (%s) does not exist.'
            self.error(IOError, s % sequence_pdb_name)

        cns_protocol = 'generate'
        cns_wdir = os.path.join(infra.get_temp_path(), cns_protocol)
        it_path = os.path.join('.', '')
        csh_script = self.setup_job(cns_wdir, it_path, cns_protocol,
                                    counter = 0,
                                    protocol_settings = protocol_settings)

        ## run cns script

        os.system('csh -f %s' % csh_script)

        ## check whether psf-file has been created

        ## TODO: assumed: initial structure and psf-file are in
        ## same directory.

        psf_file = infra.get_psf_file()
        psf_file = self.get_local_filename(psf_file, \
                                           DC.DATA_INITIAL_STRUCTURE)

        if self.getSettings()['keep_output'] == YES:
            copy_cns_output = 1
        elif self.getSettings()['keep_output'] == GZIP:
            copy_cns_output = 2
        else:
            copy_cns_output = 0

        if copy_cns_output:

            ## copy cns-output file into same directory as
            ## the psf and name it 'generate.out'

            src = os.path.join(cns_wdir, cns_protocol + '.out')
            dest = self.get_local_filename(cns_protocol + '.out', \
                                           DC.DATA_INITIAL_STRUCTURE)

            os.system('cp %s %s' % (src, dest))

            if copy_cns_output == 2: os.system('gzip -f %s' % dest)

            s = 'CNS output file has been copied to "%s".'
            self.message(s % dest)

        if not os.path.exists(psf_file):

            ## TODO: would be nice to know a bit more about the
            ## error. parsing of cns output file?

            s = 'An error occurered during creation of the PSF-file: ' + \
                '"%s" could not be found.'
            self.error(IOError, s % psf_file)

        self.message('PSF-file has been created.')

    def generate_template_pdb(self, protocol_settings, replace):

        check_type(protocol_settings, 'ProtocolSettings')

        import aria.DataContainer as DC
        import os

        infra = self.getInfrastructure()

        ## check whether initial PDB file already exists
        name = '%s_template.pdb' % infra.get_file_root()
        name = self.get_local_filename(name, DC.DATA_INITIAL_STRUCTURE)

        if not replace and os.path.exists(name):
            s = 'Template PDB file (%s) does already exist.' + \
                ' Using existing file.'
            self.message(s % name)
            return

        ## check dependencies

        if not self.psf_file_exists():
            s = 'Cannot generate template PDB file: ' + \
                'PSF file does not exist.'
            self.error(IOError, s)

        cns_protocol = 'generate_template'
        cns_wdir = os.path.join(infra.get_temp_path(), cns_protocol)
        it_path = os.path.join('.', '')
        csh_script = self.setup_job(cns_wdir, it_path, cns_protocol,
                                    counter = 0,
                                    protocol_settings = protocol_settings)

        ## run cns script

        os.system('csh -f %s' % csh_script)

        ## copy cns output

        if self.getSettings()['keep_output'] == YES:
            copy_cns_output = 1
        elif self.getSettings()['keep_output'] == GZIP:
            copy_cns_output = 2
        else:
            copy_cns_output = 0

        if copy_cns_output:

            ## copy cns-output file into same directory as
            ## the psf and name it 'generate.out'

            src = os.path.join(cns_wdir, cns_protocol + '.out')
            dest = self.get_local_filename(cns_protocol + '.out', \
                                           DC.DATA_INITIAL_STRUCTURE)

            os.system('cp %s %s' % (src, dest))

            if copy_cns_output == 2: os.system('gzip -f %s' % dest)

            s = 'CNS output file has been copied to %s.'
            self.message(s % dest)

        ## check whether template pdb-file has been created

        if not os.path.exists(name):
            s = 'Creation of template PDB-file was not successful.'
            self.error(IOError, s)
        else:
            self.message('Template PDB-file has been created.')

    def create_sequence_pdb(self, molecule):

        import os
        import aria.DataContainer as DC

        infra = self.getInfrastructure()

        ## write molecule as pdb-file

        dst = infra.get_cns_data_directories()[DC.DATA_SEQUENCE]
        sequence_pdb_name = infra.get_file_root() + '.pdb'
        dst = os.path.join(dst, sequence_pdb_name)

        self.write_molecule(molecule, dst)

        ## check whether PDB file has been created

        if not os.path.exists(dst):
            s = 'PDB sequence file has not been created.'
            self.error(IOError, s)

        self.message('Sequence PDB-file written.')

    def prepare(self, protocol_settings, molecule):

        check_type(molecule, 'Molecule')
        check_type(protocol_settings, 'ProtocolSettings')

        ## create sequence PDB file

        self.create_sequence_pdb(molecule)

        settings = self.getSettings()

        ## create psf-file

        if settings['create_psf_file'] in (YES, ALWAYS):
            replace = settings['create_psf_file'] == ALWAYS
            self.create_psf_file(protocol_settings, replace)

        ## generate template file (extended structure)

        if settings['generate_template'] in (YES, ALWAYS):
            replace = settings['generate_template'] == ALWAYS
            self.generate_template_pdb(protocol_settings, replace)

    def psf_file_exists(self):
        """
        returns 1 if psf-file exists. 0 otherwise
        """

        from aria.DataContainer import DATA_INITIAL_STRUCTURE
        import os

        infra = self.getInfrastructure()
        psf_file = infra.get_psf_file()
        psf_file = self.get_local_filename(psf_file, DATA_INITIAL_STRUCTURE)

        return os.path.exists(psf_file)

    def check_environment(self):
        """
        checks the cns environment.

        - psf-file

        returns 0 or raises an error if some file etc. could
        not be found.
        """

        if not self.psf_file_exists():
            self.error(IOError, 'PSF-file could not be found.')

        return 1

    def analysis(self, protocol_settings, ensemble, is_water=0):
        """
        Runs a list of analysis scripts on the last
        iteration.
        """

        import os
        from time import clock
        import aria.DataContainer as DC

        cmd_template = '%(cns_executable)s < %(cns_input_file)s ' + \
                       '>! %(cns_output_file)s\n'

        cmd_template_condor = "condor_submit condor.job; condor_wait condor.log\n"

        condor_queue = "input = %(cns_input_file)s\noutput = %(cns_output_file)s\nqueue\n"

        ## Add additional analysis script here

        protocols = ('wellordered.inp', 'minimize.inp', 'rmsave.inp',
                     'print_noes.inp', 'ensemble_rmsd.inp', 'rmsd.inp',
                     'print_geom.inp', 'print_dih.inp', 'print_coup.inp',
                     'print_sani.inp', 'noe_violations.inp', 'cop.inp',
                     'energy.inp')

        infra = self.getInfrastructure()

        ## check whether cns-environment is ok.
        ## (i.e. psf-file etc.)

        self.check_environment()

        ## Create file.cns file

        ## Get all PDB-files that have been
        ## calculated in last iteraton. The list is
        ## sorted with respect to total energy

        se_settings = ensemble.getSettings()
        se_settings['sort_criterion'] = 'total_energy'
        se_settings['number_of_best_structures'] = 'all'

        structures = ensemble.getFiles()
        n_structures = len(structures)

        ## Create temporary 'file.cns' which is used
        ## many scripts.

        template = 'evaluate (&filenames.bestfile_%d="NEWIT:%s")\n'

        lines = []

        for i in range(n_structures):
            base = os.path.basename(structures[i])
            lines.append(template % (i + 1, base))

        file_cns = ''.join(lines)

        ## Patch settings for last iteration into
        ## file.cns
        last_it = infra.getSettings()['n_iterations'] - 1

        if is_water:
            s = protocol_settings['water_refinement']
            n_best_structures = s['n_structures']

        else:
            s = protocol_settings['iteration_settings']
            n_best_structures = s[last_it]['number_of_best_structures']


        if n_best_structures > n_structures:
            n_best_structures = n_structures

        d = {'n_structures': n_structures,
             'n_best_structures': n_best_structures}

        header = FILE_CNS_TEMPLATE % d

        ## Store file.cns in temporary path (same
        ## location where run.cns is stored)

        temp_path = infra.get_temp_path()
        working_dir = os.path.join(temp_path, 'analysis')
        filename = os.path.join(working_dir, 'file.cns')

        ## Create working directory

        if not os.path.exists(working_dir):
            try:
                os.makedirs(working_dir)
            except Exception, msg:
                s = 'Could not create working directory for ' + \
                    'CNS analysis scripts ("%s").'
                self.error(s % working_dir)

        f = open(filename, 'w')
        f.write(header + file_cns)
        f.close()

        ## Create run.cns
        if is_water:
            iteration_path = infra.get_refinement_path()
        else:
            iteration_path = infra.get_iteration_path(last_it)

        self._write_run_cns(working_dir, iteration_path,
                            last_it, 1, None,
                            protocol_settings)

        ## create csh script

        cns_executable = self.getSettings()['local_executable']

        ## Analyses are store here.
        analysis_output = os.path.join(iteration_path, 'analysis')

        csh_dict = {'iteration_path': iteration_path,
                    'cns_resource_path': infra.get_cns_path(),
                    'analysis_output_path': analysis_output,
                    'cns_working_dir': working_dir,
                    'cns_executable': cns_executable}

        csh_script = CSH_SCRIPT_ANALYSIS % csh_dict

        condor_append = '\n'

        ## CNS analysis-script source path.

        protocol_path = infra.get_cns_analysis_path()

        msg = 'Preparing analysis script: %s'

        cns_output_filenames = {}

        for protocol in protocols:

            self.message(msg % protocol)

            protocol_name = os.path.join(protocol_path, protocol)

            output_name = os.path.splitext(protocol)[0] + '.out'
            output_name = os.path.join(working_dir, output_name)

            cns_output_filenames[protocol] = output_name

            csh_dict['cns_input_file'] = protocol_name
            csh_dict['cns_output_file'] = output_name

            if self.use_condor:

                condor_append += condor_queue % csh_dict

            else:

                csh_script += cmd_template % csh_dict

        if self.use_condor:

            csh_script += cmd_template_condor

            self._write_condor_script(csh_dict, condor_append = condor_append)

        script_name = os.path.join(working_dir, 'analysis.csh')

        ## Create some misc. files needed by
        ## some of the analysis scripts

        misc = {SECONDARY_STRUCTURE: DC.DATA_SEQUENCE,
                DIHEDRALS: DC.DATA_DIHEDRALS,
                DIHEDRALS_CSI: DC.DATA_DIHEDRALS,
                DIHEDRALS_TALOS: DC.DATA_DIHEDRALS,
                HBONDS: DC.DATA_HBONDS}

        for filename, data_type in misc.items():

            local_name = self.get_local_filename(filename, data_type)

            if not os.path.exists(local_name):
                from aria.tools import touch
                touch(local_name)

        ## Store csh-script in working directory

        f = open(script_name, 'w')
        f.write(csh_script)
        f.close()

        self.message('Starting analyses...')

        os.system('csh ' + script_name)

        ## Copy CNS output-files

        if self.getSettings()['keep_output'] == YES:
            copy_cns_output = 1
        elif self.getSettings()['keep_output'] == GZIP:
            copy_cns_output = 2
        else:
            copy_cns_output = 0

        if copy_cns_output:

            ## Create cns-output directory in iteration
            ## path if necessary.

            cns_output_dst = os.path.join(analysis_output,
                                          CNS_OUTPUT_PATH_NAME)

            if not os.path.exists(cns_output_dst):
                os.makedirs(cns_output_dst)

            for protocol in protocols:
                src = cns_output_filenames[protocol]
                path, name = os.path.split(src)
                dst = os.path.join(cns_output_dst, name)

                ## TODO: use shutil instead.

                os.system('cp %s %s' % (src, dst))

                if copy_cns_output == 2:
                    os.system('gzip %s' % dst)

            self.message('CNS Analysis output-files copied.')

        s = 'Analyses done. Results are in %s'
        self.message(s % analysis_output)
    ## < Mareuil
    def solvent_refine(self, protocol_settings, ensemble, MODE = None):
    ## Mareuil >
        """
        water-refinement
        """

        check_type(protocol_settings, 'ProtocolSettings')
        check_type(ensemble, 'StructureEnsemble')

        ## check whether cns-environment is ok.
        ## (i.e. psf-file etc.)

        self.check_environment()

        from aria.PDBReader import PDBReader
        from numpy import argsort

        import aria.JobManager as JM
        import time, os

        params = protocol_settings['water_refinement']

        if params['enabled'] == NO:
            return

        infra = self.getInfrastructure()
        ## get settings

        n_structures = params['n_structures']
        solvent_type = params['solvent']

        cns_protocol_name = 'refine_' + solvent_type

        ## get file names of sorted structures
        files = ensemble.getFiles()[:n_structures]
        ## compile name-template for solvent-refined PDB files

        pdb_name_template = os.path.join(infra.get_refinement_path(),
                                         infra.get_file_root() + '_%d' + \
                                         '_%s.pdb' % solvent_type)

        ## run refinment script of each structure.
        ## if PDB-file already exists, skip it.

        jobs = []
        counters = []
        structures = []

        ## < Mareuil
        if MODE == None:
            for i in range(len(files)):

                src = files[i]
                counter = os.path.basename(src)
                counter = counter.replace(infra.get_file_root() + '_', '')
                counter = counter.replace('.pdb', '')
                counter = int(counter)
                if os.path.exists(pdb_name_template % counter):
        	    structures.append(pdb_name_template % counter)
        	    continue

                counters.append(counter)

                local_path = 'refine_%d' % counter

                cns_wdir = os.path.join(infra.get_temp_path(), local_path)

                iteration_number = self.getInfrastructure().getSettings()['n_iterations'] - 1
                csh_script = self.setup_job(cns_wdir, infra.get_refinement_path(),
        	 			    cns_protocol_name, iteration_number,
        		 		    counter, src,
        			 	    use_default_executable=0,
        				    protocol_settings = protocol_settings, Mode = MODE)

                job_settings = JM.JobSettings()
                job_settings['script'] = csh_script
                job_settings['working_directory'] = cns_wdir
                job_settings['default_command'] = MODE
                job_settings['iteration_path'] = infra.get_refinement_path()
                if src:
         	    job_settings['pdb'] = src
                job_settings['run'] = infra.get_run_path()
                job_settings['run_cns'] = infra.get_cns_path()
                job = JM.Job(job_settings)

                jobs.append(job)
		
        elif MODE[0] == "CLUSTER" and MODE[1] == "Asynchro":
            
            for i in range(len(files)):
                src = files[i]
                counter = os.path.basename(src)
                counter = counter.replace(infra.get_file_root() + '_', '')
                counter = counter.replace('.pdb', '')
                counter = int(counter)
                if os.path.exists(pdb_name_template % counter):
        	    structures.append(pdb_name_template % counter)
        	    continue

                counters.append(counter)

                local_path = 'refine_%d' % int(i+1)

                cns_wdir = os.path.join(infra.get_temp_path(), local_path)

                iteration_number = self.getInfrastructure().getSettings()['n_iterations'] - 1
                job_name = "refine"
		PDBPATH = self.header_cns(cns_wdir, infra.get_refinement_path(), iteration_number,
                                          counter, src, protocol_settings = protocol_settings)
					  
            if len(structures) < len(files):
                csh_script = self.setup_job_csh(cns_wdir, infra.get_refinement_path(), cns_protocol_name, job_name,
            				len(files), len(files), use_default_executable=0, pdb_path = PDBPATH, Mode = MODE)

                job_settings = JM.JobSettings()
                job_settings['script'] = csh_script
                job_settings['working_directory'] = cns_wdir
                job_settings['default_command'] = MODE[0]
                job_settings['job_management'] = MODE[1]
                job_settings['iteration_path'] = infra.get_refinement_path()
                #if src:
                #	job_settings['pdb'] = src
                job_settings['run'] = infra.get_run_path()
                job_settings['run_cns'] = infra.get_cns_path()
                job_settings['number_of_structures'] = str(n_structures)
                job_settings['n_cpu'] = str(len(files))
	        job_settings['job_name'] = job_name
                job_settings['pdb_name'] = infra.get_file_root()
                job = JM.Job(job_settings)

                jobs.append(job)
		
        else:
            for i in range(len(files)):

                src = files[i]
                counter = os.path.basename(src)
                counter = counter.replace(infra.get_file_root() + '_', '')
                counter = counter.replace('.pdb', '')
                counter = int(counter)
                if os.path.exists(pdb_name_template % counter):
        	    structures.append(pdb_name_template % counter)
        	    continue

                counters.append(counter)

                local_path = 'refine_%d' % counter

                cns_wdir = os.path.join(infra.get_temp_path(), local_path)

                iteration_number = self.getInfrastructure().getSettings()['n_iterations'] - 1
                csh_script = self.setup_job(cns_wdir, infra.get_refinement_path(),
        	 			    cns_protocol_name, iteration_number,
        		 		    counter, src,
        			 	    use_default_executable=0,
        				    protocol_settings = protocol_settings, Mode = MODE)

                job_settings = JM.JobSettings()
                job_settings['script'] = csh_script
                job_settings['working_directory'] = cns_wdir
                job_settings['default_command'] = MODE[0]
		job_settings['job_management'] = MODE[1]
                job_settings['iteration_path'] = infra.get_refinement_path()
                if src:
         	    job_settings['pdb'] = src
                job_settings['run'] = infra.get_run_path()
                job_settings['run_cns'] = infra.get_cns_path()
                job = JM.Job(job_settings)

                jobs.append(job)
        ## Mareuil >
        ## start scheduler

        copy_cns_output = 0

        if jobs:
            if params['write_solvent_molecules'] == YES:
                s = 'PDB-files will include solvent molecules.'
                self.message(s)

            ## Copy noe-restraint files of last iteration into
            ## solvent-refinement directory.

            ## Assumed: solvent-refinement with respect to
            ## last iteration !

            it_number = infra.getSettings()['n_iterations'] - 1
            it_path = infra.get_iteration_path(it_number)

            tbl_names = (NAME_NOE_RESTRAINTS_AMBIG,
                         NAME_NOE_RESTRAINTS_UNAMBIG)

            src_files = [os.path.join(it_path, fn) for fn in tbl_names]

            dst = infra.get_refinement_path()

            for fn in src_files:
                if os.path.exists(fn):
                    os.system('cp %s %s' % (fn ,dst))
                else:
                    s = 'Could not find (un)ambig.tbl of last iteration. ' + \
                        'Solvent refinement aborted.'
                    self.message(s)
                    return

            scheduler = self.getJobScheduler()
            scheduler.go(jobs)

            ## wait until all jobs are done.

            while not scheduler.is_done():
                time.sleep(1.)

            ## check wether pdb-files exist

            msg = 'Refinement failed for structure %d. File %s does not exist. Please check you setup and the CNS output files for errors.'
            if self.getSettings()['keep_output'] == YES:
                copy_cns_output = 1
            elif self.getSettings()['keep_output'] == GZIP:
                copy_cns_output = 2
            else:
                copy_cns_output = 0
            ## Create directory for cns output files
            ## if necessary

            if copy_cns_output:
                refinement_path = infra.get_refinement_path()
                cns_output_dst = os.path.join(refinement_path,
                                              CNS_OUTPUT_PATH_NAME)

                if not os.path.exists(cns_output_dst):
                    os.makedirs(cns_output_dst)

            for i in range(len(counters)):
                j = counters[i]
                if not os.path.exists(pdb_name_template % j):
                    self.error(IOError, msg % (j, pdb_name_template % j))

                structures.append(pdb_name_template % j)
		
                if MODE == None:
                    j = counters[i]
                   # dest = os.path.join(cns_output_dst,
                   #                     cns_protocol_name + '_%d.out' % j)
                elif MODE[0] == 'CLUSTER' and MODE[1] == "Asynchro":
                    j = i+1
                    #dest = os.path.join(cns_output_dst,
                    #                    cns_protocol_name + '_%d.out' % j)
                else:
                    j = counters[i]
                    #dest = os.path.join(cns_output_dst,
                    #                    cns_protocol_name + '_%d.out' % j)

                if copy_cns_output:

                    ## path where the i-th structure has been stored

                    src = os.path.join(infra.get_temp_path(),
                                       'refine_%d' % j,
                                       cns_protocol_name + '.out')
                    dest = os.path.join(cns_output_dst,
                                        cns_protocol_name + '_%d.out' % j)
                    ## copy cns-output file into same directory as
                    ## pdb-files and name it 'refine_<pdb_number>.out'

                    os.system('cp %s %s' % (src, dest))

                    if copy_cns_output == 2: os.system('gzip %s' % dest)

            if copy_cns_output:
                s = 'CNS output files have been copied to %s.'
                self.message(s % cns_output_dst)

        else:
            self.message('Solvent-refined PDB-files do already exist.')

        if copy_cns_output:
            s = 'Output of solvent-refinement protocols copied.'
            self.message(s)

        self.message('Solvent refinement done.')

        ## sort files according to energy criterion

        reader = PDBReader()
        criterion = ensemble.getSettings()['sort_criterion']
        energies = [reader.get_energy(file, criterion) for file in structures]

        return [structures[index] for index in argsort(energies)]

    def cleanup(self, path):
        """
        attemps to delete ambig.tbl / unambig.tbl. however, some care must
        be taken, not to delete any of the files in the last iteration when
        solvent-refinement is turned on since it relies on both files.
        therefore, we do not delete the files in the last iteration.
        """

        check_string(path)

        keep_restraint_files = self.getSettings()['keep_restraint_files']

        if keep_restraint_files == YES:
            return

        gzip = keep_restraint_files == GZIP

        import os

        tbl_files =  (NAME_NOE_RESTRAINTS_AMBIG,
                      NAME_NOE_RESTRAINTS_UNAMBIG)

        action_performed = 0

        for tbl_file in tbl_files:

            fn = os.path.join(path, tbl_file)

            if not os.path.exists(fn):
                continue

            action_performed = 1

            if gzip:
                os.system('gzip -f ' + fn)
            else:
                try:
                    os.unlink(fn)
                except:
                    pass

        if action_performed:

            msg = {1: '(un)ambig.tbl files in %s have been gzipped',
                   0: '(un)ambig.tbl files in %s have been deleted'}

            self.message(msg[gzip] % path)

    def refine_done(self, d):
        """
        is called when all structure calculation jobs have terminated.
        """

        import os

        ## depending on 'keep_output' settings, cns-output
        ## files are copied as well.
        if self.getSettings()['keep_output'] == YES:
            copy_cns_output = 1
        elif self.getSettings()['keep_output'] == GZIP:
            copy_cns_output = 2
        else:
            copy_cns_output = 0

        infra = self.getInfrastructure()
        cns_protocol = d['cns_protocol_name']

        cns_output_dst = os.path.join(d['it_path'], CNS_OUTPUT_PATH_NAME)
        if not os.path.exists(cns_output_dst):
            os.makedirs(cns_output_dst)

        ## check whether new pdb-files exist.
        msg = 'Structure calculation failed for structure %d. Please check your setup and the CNS output files for errors.'

        for i in d['counters']:

            if copy_cns_output:

                ## path where the i-th structure has been stored

                src = os.path.join(infra.get_temp_path(),
                                   'run_cns_%d' % i, cns_protocol + '.out')

                ## copy cns-output file into same directory as
                ## pdb-files and name it 'refine_<pdb_number>.out'

                dest = os.path.join(cns_output_dst, \
                                    cns_protocol + '_%d.out' % i)

                os.system('cp %s %s' % (src, dest))

                if copy_cns_output == 2: os.system('gzip -f %s' % dest)

        if copy_cns_output:
            s = 'CNS output files have been copied to %s'
            self.message(s % cns_output_dst)

        ## Check whether all structures have been created.

        for i in d['counters']:
            if not os.path.exists(d['pdb_file_template'] % i):

                ## some structure is missing
                self.__missing_structures.append(i)

        self.__done = 1

        if self.__callback is not None:
            self.__callback()

    def done(self):
        return self.__done

    def missingStructures(self):
        return self.__missing_structures
    
    ## < Mareuil
    def go(self, peaks, protocol_settings, iteration_number, restraints = None, kept_structures = [], MODE = None):
        """
        runs cns. nessecary settings must be in the ProtocolSettings
        instance.
        """
        import os
        import time
        import aria.JobManager as JM

        self.setMODE(MODE)
        ## Mareuil >

        check_type(peaks, LIST, TUPLE)
        # BARDIAUX 2.2
        #check_elements(peaks, 'AriaPeak')
        check_elements(peaks, 'AbstractPeak')
        check_type(protocol_settings, 'ProtocolSettings')
        check_int(iteration_number)

        ## check whether cns-environment is ok

        self.check_environment()

        it_settings = protocol_settings['iteration_settings'][iteration_number]

        infra = self.getInfrastructure()
        it_path = infra.get_iteration_path(iteration_number)
        pdb_file = os.path.join(it_path, infra.get_file_root() + '_%d.pdb')

        ## check whether pdb-files already exist
        counters = range(1, it_settings['number_of_structures'] + 1)
        counters = filter(lambda i, e = os.path.exists, f = pdb_file:
                          not e(f % i), counters)

        if not len(counters):
            self.message('PDB-files do already exist for iteration %d.' \
                         % iteration_number)

            ## signal that we are done.

            if self.__callback is not None:
                self.__callback()

            return

        ## write peak-lits

        t = time.clock()

        ambig_name= os.path.join(it_path, NAME_NOE_RESTRAINTS_AMBIG)
        unambig_name = os.path.join(it_path, NAME_NOE_RESTRAINTS_UNAMBIG)

        self.write_peaklists(peaks, unambig_name, ambig_name)

        self.message('Restraint files written.', verbose_level = VL_LOW)
        self.debug('Time: %ss' % str(time.clock() - t))

        ## run structure generation jobs

        jobs = []

        cns_protocol = 'refine'

        ## TODO: it's a bit weired to be forced to do this
        ## as _write_run_cns is already accessing the initial
        ## pdb-file

        from aria.Singleton import ProjectSingleton

        import aria.DataContainer as DC

        project = ProjectSingleton()
        cns_dirs = infra.get_cns_data_directories()

        initial_structures = project.getData(DC.DATA_INITIAL_STRUCTURE)
        initial_structures = [s for s in initial_structures
                              if s['enabled'] == YES]

        if len(initial_structures):
            initial_structure = initial_structures[0]['filename']
            initial_structure = os.path.join(\
                cns_dirs[DC.DATA_INITIAL_STRUCTURE], \
                os.path.basename(initial_structure))

        else:
            initial_structure = None

        # BARDIAUX 2.3
        default_initial_structure = initial_structure

        if iteration_number > 0:
            nb_kept_struc = it_settings['number_of_kept_structures']
        else:
            nb_kept_struc = 0

        ## < Mareuil
        if MODE == None:
	    for i in counters:

                #BARDIAUX 2.3
                if i <= nb_kept_struc:
                    initial_structure = kept_structures[i-1]
                else:
                    initial_structure = default_initial_structure
		    
                cns_wdir = os.path.join(infra.get_temp_path(), 'run_cns_%d' % i)
            ## create csh-script and store script-name

                csh_script = self.setup_job(cns_wdir, it_path, cns_protocol,
                                        iteration_number, i,
                                        use_default_executable = 0,
                                        initial_structure = initial_structure,
                                        protocol_settings = protocol_settings, Mode = MODE)

                job_settings = JM.JobSettings()
                job_settings['script'] = csh_script
                job_settings['working_directory'] = cns_wdir
                job_settings['default_command'] = MODE[0]
                job_settings['job_management'] = MODE[1]
                job_settings['run'] = infra.get_run_path()
                job_settings['run_cns'] = infra.get_cns_path()
                if initial_structure:
                    job_settings['pdb'] = initial_structure
                job_settings['iteration_path'] = it_path
                job = JM.Job(job_settings)

                jobs.append(job)
            ## start scheduler

            scheduler = self.getJobScheduler()

            ## a bit too pragmatic: store some things which
            ## have to be accessed later (c.f. refine_done)

            vars = {'counters': counters,
                'pdb_file_template': pdb_file,
                'cns_protocol_name': cns_protocol,
                'it_path': it_path,
                'current_iteration': iteration_number}

            callback = lambda v = vars, f = self.refine_done: f(v)
            scheduler.set_callback(callback)

            self.__done = 0
            self.__missing_structures = []

            scheduler.go(jobs)

        elif MODE[0] == "CLUSTER" and MODE[1] == "Asynchro":
            #BARDIAUX 2.3
            scheduler = self.getJobScheduler()
            n_cpu = scheduler.getSettings()['host_list'][0]['n_cpu']
	    n_cpu_counters = range(1,n_cpu+1)
	    n_cpu_counters = filter(lambda i, e = os.path.exists, f = pdb_file:
                          not e(f % i), n_cpu_counters)

            for i in n_cpu_counters:
                if i <= nb_kept_struc:
                    initial_structure = kept_structures[i-1]
                else:
                    initial_structure = default_initial_structure

                cns_wdir = os.path.join(infra.get_temp_path(), 'run_cns_%d' % int(n_cpu_counters.index(i)+1))
		PDBPATH = self.header_cns(cns_wdir, it_path, iteration_number, 
                                          i, initial_structure = initial_structure, protocol_settings = protocol_settings)

        ## create csh-script and store script-name
            if len(counters) != 0:
                job_name = "run_cns"
                csh_script = self.setup_job_csh(cns_wdir, it_path, cns_protocol, job_name, 
            			        len(counters), len(n_cpu_counters), use_default_executable = 0, pdb_path = PDBPATH, Mode = MODE)

                job_settings = JM.JobSettings()
                job_settings['script'] = csh_script
                job_settings['working_directory'] = cns_wdir
                job_settings['default_command'] = MODE[0]
                job_settings['job_management'] = MODE[1]
                job_settings['run'] = infra.get_run_path()
                job_settings['run_cns'] = infra.get_cns_path()
	        job_settings['number_of_structures'] = str(it_settings['number_of_structures'])
		job_settings['n_cpu'] = str(len(n_cpu_counters))
	        job_settings['job_name'] = job_name
		job_settings['pdb_name'] = infra.get_file_root()
                #if initial_structure:
                #    job_settings['pdb'] = initial_structure
                job_settings['iteration_path'] = it_path
                job = JM.Job(job_settings)

                jobs.append(job)

            ## a bit too pragmatic: store some things which
            ## have to be accessed later (c.f. refine_done)
            vars = {'counters': range(1, it_settings['number_of_structures'] + 1),
                'pdb_file_template': pdb_file,
                'cns_protocol_name': cns_protocol,
                'it_path': it_path,
                'current_iteration': iteration_number}

            callback = lambda v = vars, f = self.refine_done: f(v)
            scheduler.set_callback(callback)

            self.__done = 0
            self.__missing_structures = []

            error_list = scheduler.go(jobs)
            print error_list    
	    
	else:
            for i in counters:

                #BARDIAUX 2.3
                if i <= nb_kept_struc:
                    initial_structure = kept_structures[i-1]
                else:
                    initial_structure = default_initial_structure
		    
                cns_wdir = os.path.join(infra.get_temp_path(), 'run_cns_%d' % i)
            ## create csh-script and store script-name
                csh_script = self.setup_job(cns_wdir, it_path, cns_protocol,
                                        iteration_number, i,
                                        use_default_executable = 0,
                                        initial_structure = initial_structure,
                                        protocol_settings = protocol_settings, Mode = MODE)

                job_settings = JM.JobSettings()
                job_settings['script'] = csh_script
                job_settings['working_directory'] = cns_wdir
                job_settings['default_command'] = MODE[0]
                job_settings['job_management'] = MODE[1]
                job_settings['run'] = infra.get_run_path()
                job_settings['run_cns'] = infra.get_cns_path()
                if initial_structure:
                    job_settings['pdb'] = initial_structure
                job_settings['iteration_path'] = it_path
                job = JM.Job(job_settings)

                jobs.append(job)

            ## start scheduler

            scheduler = self.getJobScheduler()

            ## a bit too pragmatic: store some things which
            ## have to be accessed later (c.f. refine_done)

            vars = {'counters': counters,
                'pdb_file_template': pdb_file,
                'cns_protocol_name': cns_protocol,
                'it_path': it_path,
                'current_iteration': iteration_number}

            callback = lambda v = vars, f = self.refine_done: f(v)
            scheduler.set_callback(callback)

            self.__done = 0
            self.__missing_structures = []

            scheduler.go(jobs)
            ## Mareuil >
        return

    def set_callback(self, f):
        """
        the callback will be called after the job-manager
        has executed all jobs
        """

        self.__callback = f

    def get_callback(self):
        return self.__callback

class CNSXMLPickler(XMLBasePickler):

    order = ('annealing_parameters', 'md_parameters',
             'local_executable', 'keep_output',
             'keep_restraint_files', 'create_psf_file',
             'generate_template', 'nonbonded_parameters')

    def _xml_state(self, x):

        e = XMLElement(tag_order = self.order)

        ## annealing / md parameters

        e.annealing_parameters = x.getAnnealingParameters()
        e.md_parameters = x.getMDParameters()

        s = x.getSettings()

        e.local_executable = s['local_executable']
        e.keep_output = s['keep_output']
        e.create_psf_file = s['create_psf_file']
        e.generate_template = s['generate_template']
        e.nonbonded_parameters = s['nonbonded_parameters']
        e.keep_restraint_files = s['keep_restraint_files']

        return e

    def load_from_element(self, e):

        s = CNSSettings()

        ## relaxed mode? executable need not exist

        if self.relaxed:
            entity = s.getEntity('local_executable')
            entity.mandatory(0)

        s['local_executable'] = str(e.local_executable)

        if self.relaxed:
            entity.mandatory()

        s['keep_output'] = str(e.keep_output)
        s['create_psf_file'] = str(e.create_psf_file)
        s['generate_template'] = str(e.generate_template)
        s['nonbonded_parameters'] = str(e.nonbonded_parameters)

        ## TODO: remove in release version

        if hasattr(e, 'keep_restraint_files'):
            value = str(e.keep_restraint_files)
        else:
            value = YES

        s['keep_restraint_files'] = value

        cns = CNS(s)
        cns.setMDParameters(e.md_parameters)
        cns.setAnnealingParameters(e.annealing_parameters)

        return cns

CNS._xml_state = CNSXMLPickler()._xml_state

