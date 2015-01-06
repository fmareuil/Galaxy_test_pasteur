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
from aria.xmlutils import XMLBasePickler
from aria.xmlutils import XMLElement

import aria.Settings as Settings

class AnalysisParameters(Settings.Settings):

    def create(self):

        ## TODO: delete hack for paths; enable check for existence

        d = {'procheck_executable': Settings.Path(exists = 0),
             'procheck_enabled': Settings.YesNoChoice(),
             'whatif_executable': Settings.Path(exists = 0),
             'whatif_enabled': Settings.YesNoChoice(),
             'prosa_executable': Settings.Path(exists = 0),
             'prosa_enabled': Settings.YesNoChoice(),
             'clashlist_executable': Settings.Path(exists = 0),
             'clashlist_enabled': Settings.YesNoChoice()}

        descr = 'If enabled, ARIA runs a number of scripts to analyse ' + \
                'the structures in last iteration and solvent refinement.'

        d['cns_analyses'] = Settings.YesNoChoice(description = descr)

        return d

    def create_default_values(self):

        defaults = {'procheck_enabled': YES,
                    'whatif_enabled': YES,
                    'prosa_enabled': YES,
                    'clashlist_enabled': YES,
                    'procheck_executable': 'procheck',
                    'prosa_executable': 'prosaII',
                    'whatif_executable': 'whatif',
                    'clashlist_executable': 'clashlist',
                    'cns_analyses': YES}

        return defaults

class Analyser(AriaBaseClass):

    def __init__(self, settings):
        AriaBaseClass.__init__(self, settings = settings, name = 'Analyser')

    def isEnabled(self):
        """
        convenience function. returns 1 if any of the analyses
        is enabled.
        """

        entities = ('procheck_enabled', 'whatif_enabled',
                    'prosa_enabled')

        s = self.getSettings()

        for e in entities:
            if s[e] == YES:
                return 1

        return 0

    def quality_checks(self, pdb_path, n_structures):
        from aria.Singleton import ProjectSingleton
        from legacy.QualityChecks import QualityChecks as checker
        from os.path import exists 

        # BARDIAUX 2.2
        from aria.WhatifProfile import WhatifProfile
        
        check_string(pdb_path)

        project = ProjectSingleton()
        infra = project.getInfrastructure()

        temp_path = infra.get_temp_path()

        s = self.getSettings()

        procheck_enabled = s['procheck_enabled'] == YES
        prosa_enabled = s['prosa_enabled'] == YES
        whatif_enabled = s['whatif_enabled'] == YES
        clashlist_enabled = s['clashlist_enabled'] == YES
                

        ## check whether executables exist

        msg = 'Could not find executable for %s ("%s"). Quality-check ' + \
            'has been disabled.'

        if procheck_enabled and not exists(s['procheck_executable']):
            self.message(msg % ('Procheck', s['procheck_executable']))
            procheck_enabled = 0

        if prosa_enabled and not exists(s['prosa_executable']):
            self.message(msg % ('Prosa II', s['prosa_executable']))
            prosa_enabled = 0

        if whatif_enabled and not exists(s['whatif_executable']):
            self.message(msg % ('WhatIf', s['whatif_executable']))
            whatif_enabled = 0

        if clashlist_enabled and not exists(s['clashlist_executable']):
            self.message(msg % ('Clashlist', s['clashlist_executable']))
            clashlist_enabled = 0
            
        if not(procheck_enabled or prosa_enabled or whatif_enabled or clashlist_enabled):
            return

        msg = 'Running quality checks on structures in %s'
        self.message(msg % pdb_path)

        try:

            checker.runChecks(workingDirectory = pdb_path,
                              trashDirectory = temp_path,
                              procheckOnOff = procheck_enabled,
                              prosaOnOff = prosa_enabled,
                              whatifOnOff = whatif_enabled,
                              clashlistOnOff = clashlist_enabled,
                              procheckExe = s['procheck_executable'],
                              prosaExe = s['prosa_executable'],
                              whatIfExe = s['whatif_executable'],
                              clashlistExe = s['clashlist_executable'],
                              howManyPdb = n_structures)
            
            ## BARDIAUX 2.2
            if whatif_enabled:
                w_prof = WhatifProfile(pdb_path, whatIfExe = s['whatif_executable'])
                w_prof.makeProfiles()
                
                
        except Exception, msg:
            from aria.tools import last_traceback
            self.message(last_traceback())
            s = 'Could not run quality-checks: %s'
            self.warning(s % msg)

        self.message('Quality checks done.')
        
    # BARDIAUX 2.2
    def _getSolventEnsemble(self, iteration, settings):
        """
        return the StructureEnsemble of the refine structures
        """

        from aria.Singleton import ProjectSingleton
        import os
        project = ProjectSingleton()
        infra = project.getInfrastructure()
        
        n_structures = settings['n_structures']
        solvent_type = settings['solvent']
        
        ensemble = iteration.getStructureEnsemble()
        files = ensemble.getFiles()[:n_structures]
        
        pdb_name_template = os.path.join(infra.get_refinement_path(),
                                         infra.get_file_root() + '_%d' + \
                                         '_%s.pdb' % solvent_type)

        float_name_template = os.path.join(infra.get_refinement_path(),
                                         infra.get_file_root() + '_%d' + \
                                         '_%s.float' % solvent_type)  
        counters = []
        structures = []
        float_files = []
        
        for i in range(len(files)):

            src = files[i]
            
            counter = os.path.basename(src)
            counter = counter.replace(infra.get_file_root() + '_', '')
            counter = counter.replace('.pdb', '')
            counter = int(counter)
            
            if os.path.exists(pdb_name_template % counter):
                structures.append(pdb_name_template % counter)
                continue

            if os.path.exists(float_name_template % counter):
                float_files.append(float_name_template % counter)
                continue

            counters.append(counter)

        if float_files == []:
            float_files = None

        if float_files is not None:
            
            from aria.FloatFile import FloatFile

            parser = FloatFile()
            swapped_atoms = [parser.parse(f) for f in float_files]
            
        else:
            swapped_atoms = None

        ## create structure-ensemble from given 'structures'
        import aria.StructureEnsemble as SE
        
        molecule = project.getMolecule()
        
        se_settings = SE.StructureEnsembleSettings()
        ensemble = SE.StructureEnsemble(se_settings)
        ensemble.read(structures, molecule, 'cns', swapped_atoms)

        return ensemble        
            
    def cns_analyses(self, iteration, is_water = 0): # BARDIAUX 2.2
        
        check_type(iteration, 'Iteration')

        if self.getSettings()['cns_analyses'] <> YES:
            return

        if is_water:
            msg = 'water refinement'
        else:
            msg = 'iteration %d' %  iteration.getNumber()
        s = 'Running CNS analyses for %s.'
        self.message(s % msg)
        
        from aria.Singleton import ProjectSingleton

        project = ProjectSingleton()
        
        protocol_settings = project.getProtocol().getSettings()
        engine = project.getStructureEngine()
        
        if is_water:
            ensemble = self._getSolventEnsemble(iteration, \
                                                protocol_settings['water_refinement'])

    
        else:
            ensemble = iteration.getStructureEnsemble()
        
        engine.analysis(protocol_settings, ensemble, is_water)
        
##     def cns_analyses(self, iteration):
##         check_type(iteration, 'Iteration')

##         if self.getSettings()['cns_analyses'] <> YES:
##             return

##         s = 'Running CNS analyses for iteration %d.'
##         self.message(s % iteration.getNumber())
        
##         from Singleton import ProjectSingleton

##         project = ProjectSingleton()
        
##         protocol_settings = project.getProtocol().getSettings()
##         engine = project.getStructureEngine()
##         ensemble = iteration.getStructureEnsemble()
        
##         engine.analysis(protocol_settings, ensemble)
        
class AnalyserXMLPickler(XMLBasePickler):

    order = ('structures_analysis', 'procheck', 'prosa', 'whatif', 'clashlist')

    def _xml_state(self, x):

        x = x.getSettings()

        procheck = XMLElement()
        procheck.enabled = x['procheck_enabled']
        procheck.executable = x['procheck_executable']

        whatif = XMLElement()
        whatif.enabled = x['whatif_enabled']
        whatif.executable = x['whatif_executable']

        prosa = XMLElement()
        prosa.enabled = x['prosa_enabled']
        prosa.executable = x['prosa_executable']

        clashlist = XMLElement()
        clashlist.enabled = x['clashlist_enabled']
        clashlist.executable = x['clashlist_executable']        

        cns = XMLElement()
        cns.enabled = x['cns_analyses']

        e = XMLElement(tag_order = self.order)

        e.structures_analysis = cns
        e.procheck = procheck
        e.whatif = whatif
        e.prosa = prosa
        e.clashlist = clashlist

        return e

    def load_from_element(self, e):

        s = AnalysisParameters()
        
        s['procheck_executable'] = str(e.procheck.executable)
        s['procheck_enabled'] = str(e.procheck.enabled)
        s['whatif_executable'] = str(e.whatif.executable)
        s['whatif_enabled'] = str(e.whatif.enabled)
        s['prosa_executable'] = str(e.prosa.executable)
        s['prosa_enabled'] = str(e.prosa.enabled)
        s['cns_analyses'] = str(e.structures_analysis.enabled)

        if hasattr(e, 'clashlist'):
            s['clashlist_executable'] = str(e.clashlist.executable)
            s['clashlist_enabled'] = str(e.clashlist.enabled)
        else:
            s['clashlist_executable'] = ''
            s['clashlist_enabled'] = NO

        return Analyser(s)

Analyser._xml_state = AnalyserXMLPickler()._xml_state
