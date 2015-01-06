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


import Tix
from Tkconstants import *
import sys, os
from aria.ariabase import *
import panels
from aria.tools import Dump, Load, last_traceback
from aria.JobManager import JobSchedulerSettings as JSSettings

TCL_DONT_WAIT		= 1<<1
TCL_WINDOW_EVENTS	= 1<<2
TCL_FILE_EVENTS		= 1<<3
TCL_TIMER_EVENTS	= 1<<4
TCL_IDLE_EVENTS		= 1<<5
TCL_ALL_EVENTS		= 0

ARIA_ENV = 'ARIA2'
PATH_GUI_MODULES = 'src/py/aria/gui'
GUI_VERSION = '1.4'

GUI_MODULES = ('decorators.py', 'panels.py', 'widgets.py')

TK_OPTIONS_TEMPLATE = \
"""
*font:          Arial %(font_size)d
*Label*font:    Arial %(font_size)d 
*Entry*background: white
*Text*background: white
"""

GUI_SETTINGS_FILENAME = '~/.aria_gui'

CANCEL = 'cancel'

class GUISettings(JSSettings):

    def __init__(self):

        JSSettings.__init__(self)

    def create(self):
        from aria.Settings import PositiveInteger, TypeEntity, String
        from aria.Settings import Integer, Float, YesNoChoice
             
        font_size_descr = 'The font size in points. A new font size will take effect when restarting the application.'

        cns_exe_descr = 'The default location of the CNS executable. For a new project, it is also used as default executable for every machine in the host-list.'

        ccpn_names_descr = 'Show detailed names for CCPN Data (this can make the display of Data panels slower)'

        d = {'n_recent_files': PositiveInteger(),
             'recent_files': TypeEntity(LIST),
             'geometry': String(),
             'output_panel_size': Integer(),
             'tree_panel_size': Integer(),
             'font_size': PositiveInteger(description = font_size_descr),
             'text_editor': String(),
             'pdb_viewer': String(),
             'cns_executable': \
             String(description=cns_exe_descr),
             'prosa_executable': \
             String(description='Default Prosa II executable'),
             'procheck_executable': \
             String(description='Default Procheck executable'),
             'whatif_executable': \
             String(description='Default WhatIf executable'),
             'host_list': TypeEntity(LIST),
             ## < Mareuil
             'default_command': String(),
             'job_management': String(),
             ## Mareuil >
             'author': String(),
             'use_tooltips': String(),
             'no_main_window_tooltips': String(),
             'show_ccpn_names' : YesNoChoice(description = ccpn_names_descr)}
        return d

    def create_default_values(self):

        d = {'n_recent_files': 5,
             'recent_files': [],
             'geometry': '1000x800+11+33',
             'output_panel_size': 200,
             'tree_panel_size': 300,
             'font_size': 10,
             'text_editor': 'emacs',
             'pdb_viewer': 'rasmol',
             'cns_executable': '',
             'prosa_executable': '',
             'procheck_executable': '',
             'whatif_executable': '',
             'host_list': [],
             ## < Mareuil
             'default_command': '',
             'job_management': '',
             ## Mareuil >
             'author': '',
             'use_tooltips': '1',
             'no_main_window_tooltips': '0',
             'show_ccpn_names' : YES}
        return d
        
def popup(p):

    w = Tix.PopupMenu(p, title = 'popup')
    ## < Mareuil
    w.menu.add_command(label = 'submit_command')
    ## Mareuil >

    return w

def DiscardDialog(type = 'YesNo', title = 'Warning'):

    import widgets

    text = \
"""
Some entry fields are missing or have not been
saved yet. Do you want to continue anyway?
"""

    if type == 'YesNo':
        dialog = widgets.YesNoDialog(title = title, text = text, modal = 1)
        return dialog.get()
    else:
        dialog = widgets.OKDialog(title = title, text = text, modal = 1)

def EntriesMissingDialog(title = 'Warning'):

    import widgets

    text = \
"""
Some entry fields are missing or
have not been saved yet.
"""

    dialog = widgets.OKDialog(title = title, text = text, modal = 1)

def NotSavedDialog():

    import widgets

    text = 'Project has not been saved yet.\n' + \
           'Save now?'

    dialog = widgets.YesNoCancelDialog(title = 'Warning', text = text,
                                       modal = 1)

    return dialog.get()

def AboutBox():

    import widgets, Tkinter
    from aria.ariabase import AriaBaseClass

    imageFile = os.path.join(AriaBaseClass.install_path, PATH_GUI_MODULES, 'logo.gif')
    logo = None
    if os.path.exists(imageFile):
        logo = Tkinter.PhotoImage(file=imageFile)    

    text = """
ARIA Version %s, GUI version %s.

Copyright Benjamin Bardiaux, Michael Habeck, Jens P. Linge, Therese E. Malliavin,
Wolfgang Rieping, and Michael Nilges.

All rights reserved.

If you use this software, please quote the following reference(s):

Rieping W., Habeck M., Bardiaux B., Bernard A., Malliavin T.E., Nilges M. (2007)
ARIA2: automated NOE assignment and data integration in NMR structure calculation
Bioinformatics 23:381-382
""" % (AriaBaseClass().get_version_string(), GUI_VERSION)

    dialog = widgets.OKDialog(title = 'About', text = text, modal = 1, image = logo)

class GUI(AriaBaseClass):

    messages = \
             {'log_window': 'Log window. Python output and other information is stored here.',
              'tree_panel': 'This is the project tree. To select a node, <click> it; to expand/collaps a node, <click> the small "plus"/"minus" symbol at its left. All entry fields are self-validating. In other words, an entry field would complain, for instance, if you try to enter a string when a number is expected. If an input-mask (panel) could not be validated (due to erroneous or missing input-data), the corresponding folder-icon changes to a red exclamation mark. A project cannot be saved unless all input-masks are valid. Switching from an invalid or incomplete panel to another input-mask, renders the former panel invalid (i.e. its icon becomes an exclamation mark).\n\nThis message can be turned off. Check the menu "Help".'}

    tkoptions_filename = '~/.aria_gui_tkoptions'
    
    def __init__(self, top):

        AriaBaseClass.__init__(self)

        ## turn-off line-wrapping for aria errors/warnings/messges

        AriaBaseClass.wrap_lines = 0
        
        self.root = top
        self.__exit = -1

        self.balloon = None			# balloon widget
        self.useBalloons = Tix.StringVar()
        self.no_tooltips_main_window = Tix.StringVar()
        self.statusbar = None			# status bar widget

        self.project = {'file': None, 'data': None}

    def load_settings(self):
        
        import os

        ## load settings

        self.settings_file = GUI_SETTINGS_FILENAME

        try:
            settings = Load(self.settings_file)
            ok = 1
            
        except Exception, message:
            self.message(str(message))
            self.message('Could not load settings (%s).' % self.settings_file)
            self.message('Using default settings.')
            settings = GUISettings()

            ## set current version
            
            self.version = GUI_VERSION
            
            settings.reset()

            ok = 0

        if ok:
            
            ## check version. In new versions of the GUI
            ## we always dump the tuple (version, gui_settings)

            if type(settings) == type(()):
                self.version, settings = settings

            ## TODO: remove for release version.
            ## old version

            else:
                for host in settings['host_list']:
                    ## < Mareuil
                    host['submit_command'] = ''
                    ## Mareuil >
                self.version = GUI_VERSION

        if not settings.has_key('show_ccpn_names'):
            old_settings = settings
            settings = GUISettings()
            settings.reset()
            settings.update(old_settings)
            
        self.setSettings(settings)

    def setup(self):

        import os

        settings = self.getSettings()

        ## initialize list of recent files

        recent_files = list(settings['recent_files'])
        recent_files.reverse()

        self.recent_files = []
        [self.add_to_recent_files(f) for f in recent_files]

        ## Set geometry of main window

        z = self.root.winfo_toplevel()
        z.geometry(settings['geometry'])

    def message(self, m):

        m = str(m)

        if m == '\n':
            return
        
        if m[-1] <> '\n':
            m += '\n'

        AriaBaseClass.message(self, m)

    def write(self, m):
        self.output_text.text.insert(Tix.END, m)
        self.output_text.text.see(Tix.END)

    def set_title(self, t = None):

        import aria.ariabase as ariabase

        z = self.root.winfo_toplevel()

        if t is not None:
            t = ' - %s' % t
        else:
            t = ''
        
        z.wm_title('ARIA %s GUI' % ariabase.AriaBaseClass().get_version_string()+t)

    def insert_panels(self, tree, project):
        
        import aria.DataContainer as DC
        import panels
        from aria.OrderedDict import OrderedDict

        NP = self.new_panel

        project.ccpn_project = None

        project_node, p = NP(tree, 'Project', panels.ProjectPanel,
                             name = 'project')
        
        p.set_settings(project)

        ## data
        
        data_node, p = NP(project_node, 'Data')

        d = OrderedDict()

        self.data = d

        d[DC.DATA_SEQUENCE] = data_node
        ADD = self.add_data_container
        
        ## add sequence

        node, panel = ADD(DC.DATA_SEQUENCE,
                          project.getData(DC.DATA_SEQUENCE)[0],
                          'Molecular system', id = 'sequence')

        panel.project = project
        panel.root = self.root

        ## add spectra
        
        d[DC.DATA_SPECTRUM] = NP(data_node, 'Spectra')[0]

        for s in project.getData(DC.DATA_SPECTRUM):
            node, panel = ADD(DC.DATA_SPECTRUM, s)
            node.close()
            panel.project = project

        ## BARDIAUX 2.2 add symmetry data
        d[DC.DATA_SYMMETRY], p = \
                             NP(data_node, 'Symmetry', panels.SymmetryPanel)
        p.set_settings(project.getData(DC.DATA_SYMMETRY)[0])        

        # templates
        
        name = 'Initial Structure Ensemble'
        node, panel = NP(data_node, name, panels.DataTemplatePanel)
        d[DC.DATA_TEMPLATE_STRUCTURE] = node
        panel.set_settings(project)
        panel.project = project
        
        ## add dihedral restraints

        d[DC.DATA_DIHEDRALS] = NP(data_node, 'Dihedral angles')[0]

        for dh in project.getData(DC.DATA_DIHEDRALS):
            ADD(DC.DATA_DIHEDRALS, dh)

        ## add hbond restraints

        d[DC.DATA_HBONDS] = NP(data_node, 'Hydrogen bonds')[0]
        
        for hb in project.getData(DC.DATA_HBONDS):
            ADD(DC.DATA_HBONDS, hb)

        ## add rdc restraints

        d[DC.DATA_RDCS] = NP(data_node, 'RDCs')[0]

        for x in project.getData(DC.DATA_RDCS):
            ADD(DC.DATA_RDCS, x)

        ## add karplus restraints

        d[DC.DATA_KARPLUS] = NP(data_node, 'Scalar couplings')[0]

        for kp in project.getData(DC.DATA_KARPLUS):
            ADD(DC.DATA_KARPLUS, kp)

        ## add (un)ambiguous distance restraints

        d[DC.DATA_AMBIGUOUS] = NP(data_node, 'Ambiguous distances')[0]

        for x in project.getData(DC.DATA_AMBIGUOUS):
            ADD(DC.DATA_AMBIGUOUS, x)
            
        d[DC.DATA_UNAMBIGUOUS] = NP(data_node, 'Unambiguous distances')[0]

        for x in project.getData(DC.DATA_UNAMBIGUOUS):
            ADD(DC.DATA_UNAMBIGUOUS, x)

        ## add ssbond restraints

        d[DC.DATA_SSBONDS] = NP(data_node, 'Disulfide bridges (restraints)')[0]

        for x in project.getData(DC.DATA_SSBONDS):
            ADD(DC.DATA_SSBONDS, x)

        ## add ss-bridges [covalent]

        d[DC.DATA_SSBRIDGE], p = \
                             NP(data_node, 'Disulfide bridges (covalent)',
                                            panels.SSBridgePanel)
        p.set_settings(project)

        ## add his-patches

        d[DC.DATA_HISPATCH], p = \
                             NP(data_node, 'HIS patches', panels.HisPatchPanel)
        p.set_settings(project)

        ## BARDIAUX 2.2
        ## add cispro-patches

        d[DC.DATA_CISPROPATCH], p = \
                             NP(data_node, 'Cis-Proline patches', panels.CisProPatchPanel)
        p.set_settings(project)
        
        ## BARDIAUX 2.2
        ## add ZN patches

        d[DC.DATA_ZNPATCH] = NP(data_node, 'Zinc')[0]

        for x in project.getData(DC.DATA_ZNPATCH):
            ADD(DC.DATA_ZNPATCH, x)

        ## Other data    
        d[DC.DATA_OTHER] = NP(data_node, 'Other')[0]

        for x in project.getData(DC.DATA_OTHER):
            ADD(DC.DATA_OTHER, x)
            
        
        ## CCPN

        ccpn_node, p = NP(project_node, "CCPN data model", panels.CCPNPanel,
                          name='ccpn')
        p.project = project

        p.set_settings({'ccpn': project.ccpn_model,
                        'report': project.getReporter()['ccpn']})

        ## protocol
        
        protocol, p = NP(project_node, 'Protocol',
                         panels.ProtocolPanel,
                         name = 'protocol')
        
        p.set_settings(project)
        p.project = project
        
        self.iterations, p = NP(protocol, 'Iterations')
        node, panel = NP(protocol, 'Water refinement',
                         panels.WaterRefinementPanel)
        protocol_settings = project.getProtocol().getSettings()
        panel.set_settings(protocol_settings['water_refinement'])

        ## add iterations

        its = project.getProtocol().getSettings()['iteration_settings']

        for number, settings in its.items():
            self.add_iteration(settings)

        ## Structure generation, job_manager

        engine = project.getStructureEngine()
        struct_gen, p = NP(project_node, 'Structure Generation')
	
        Job_gen, p = NP(struct_gen, 'Job Manager')
	
        Ms, Ms_panel = NP(Job_gen, 'Mode',
                          panels.ModePanel,
                          name = 'job_manager')
        Ms_panel.set_settings(engine.getJobScheduler().getSettings())
	
        js, js_panel = NP(Job_gen, 'Host Manager',
                          panels.JobManagerPanel,
                          name = 'host')
        js_panel.set_settings(engine.getJobScheduler().getSettings(), self.getSettings())

        ## cns, dynamics
        cns, p = NP(struct_gen, 'CNS', panels.CNSPanel, name = 'cns')
        p.set_settings(engine.getSettings())

        ## annealing

        annealing, p = NP(cns, 'Annealing Parameters')

        ## add entry point for annealing sub-panels
        
        d = OrderedDict()
        
        d[DC.DATA_ANNEALING_AMBIG] = annealing
        d[DC.DATA_ANNEALING_UNAMBIG] = annealing
        d[DC.DATA_ANNEALING_DIHEDRAL] = annealing
        d[DC.DATA_ANNEALING_HBOND] = annealing
        d[DC.DATA_ANNEALING_FBHW] = annealing
        ## BARDIAUX 2.2
        d[DC.DATA_ANNEALING_SYM] = annealing
        ## BERNARD 2.3
        d[DC.DATA_ANNEALING_LOGHARMONIC] = annealing
        
        self.annealing = d

        node, panel = NP(cns, 'Dynamics', panels.DynamicsPanel)
        panel.set_settings(engine.getMDParameters())

        ## add annealing parameters

        engine = project.getStructureEngine()
        sa_params = engine.getAnnealingParameters()

        ADD = self.add_annealing_sub_panel

        ## ambiguous restraints

        params = sa_params.getParameters(DC.DATA_ANNEALING_AMBIG)
        ADD(DC.DATA_ANNEALING_AMBIG, params, 'Ambiguous restraints')

        ## unambiguous restraints

        params = sa_params.getParameters(DC.DATA_ANNEALING_UNAMBIG)
        ADD(DC.DATA_ANNEALING_UNAMBIG, params, 'Unambiguous restraints')

        ## Symmetry restraints

        params = sa_params.getParameters(DC.DATA_ANNEALING_SYM)
        ADD(DC.DATA_ANNEALING_SYM, params, 'Symmetry restraints')

        ## FBHW

        params = sa_params.getParameters(DC.DATA_ANNEALING_FBHW)
        ADD(DC.DATA_ANNEALING_FBHW, params, 'Flat Bottom Harmonic Wall')
        
        ##  Potentiel Restraint  ## BERNARD 2.3

        params = sa_params.getParameters(DC.DATA_ANNEALING_LOGHARMONIC)
        ADD(DC.DATA_ANNEALING_LOGHARMONIC, params, 'Log-Harmonic Potential')

        ## hbonds

        params = sa_params.getParameters(DC.DATA_ANNEALING_HBOND)
        ADD(DC.DATA_ANNEALING_HBOND, params, 'Hydrogen bonds')
        
        ## dihedrals

        params = sa_params.getParameters(DC.DATA_ANNEALING_DIHEDRAL)
        ADD(DC.DATA_ANNEALING_DIHEDRAL, params, 'Dihedral angle restraints')

        ## RDCs

        rdcs, p = NP(annealing, 'RDCs')
        self.annealing[DC.DATA_ANNEALING_RDC] = rdcs
        
        params = sa_params.getParameters(DC.DATA_ANNEALING_RDC)

        for p in params:
            name = 'class ' + str(p['class'])
            ADD(DC.DATA_ANNEALING_RDC, p, name)

        ## karplus

        karplus, p = NP(annealing, 'Scalar couplings')
        self.annealing[DC.DATA_ANNEALING_KARPLUS] = karplus
        
        params = sa_params.getParameters(DC.DATA_ANNEALING_KARPLUS)

        for p in params:
            name = 'class ' + str(p['class'])
            ADD(DC.DATA_ANNEALING_KARPLUS, p, name)

        ## analysis
            
        node, panel = NP(project_node, 'Analyses',
                         panels.AnalysisPanel)
        panel.set_settings(project.getAnalyser().getSettings())

        ## report

        repor, p = NP(project_node, 'Report',
                                  panels.ReportPanel)
        p.set_settings(project)

        ## results
        import AriaViewer
        self.cmaps, p = NP(project_node, 'Peak Maps', panels.CmapOptionsPanel,
                           name = 'cmap_opt')
        set = AriaViewer.DisplaySettings()
        set.reset()
        p.set_settings(set)
        its = project.getProtocol().getSettings()['iteration_settings']
        for number, settings in its.items():
            self.add_cmaps(settings)
                   

        karplus.close()
        rdcs.close()
        data_node.open()
        protocol.close()
        struct_gen.open()
        Job_gen.close()
        annealing.close()
        project_node.open()

    def load_recent(self, filename):

        if not self.close_project():
            return 0

        return self.load_project(filename)

    def load_project(self, filename, relax = 1, to_recent_files = 1,
                     new_project = 0):

        from aria.AriaXML import AriaXMLPickler
        import widgets
        from aria.Settings import Path
        import os

        ## Disable path-checking temporarily when loading
        ## existing project.

        Path.global_mandatory = 0

        pickler = AriaXMLPickler()

        try:
            if relax:
                project = pickler.load_relaxed(filename)
            else:
                project = pickler.load(filename)
            
        except Exception, msg:
            import aria.tools as tools
            
            print tools.last_traceback()
            self.message('Could not load project: ' + str(msg))

            Path.global_mandatory = 1

            return 0

        ## Check whether loaded XML file is
        ## a project file

        if not is_type(project, 'Project'):
            s = 'Specified XML file is not an ARIA project file.'
            self.message(s)

            return

        Path.global_mandatory = 1

        ## close current project

        self.project['data'] = project
        self.project['file'] = filename
        
        if to_recent_files and not new_project:
            self.add_to_recent_files(filename)

        if not new_project:
            project_name = project.getSettings()['name']
            title = '%s (%s)' % (project_name, filename)
        else:
            title = ''

        self.set_title(title)
        
        ## create new tree
        b_path = os.path.join(AriaBaseClass.install_path, PATH_GUI_MODULES)
        new_tree = widgets.Tree(self.tree_panel, bitmap_path = b_path)
        
        ## remove old panels and insert new ones
        
        self.control_panel.destroy_panels()
        self.insert_panels(new_tree, project)
        
        ## remove old tree

        if self.tree is not None:
            self.tree.destroy()

        new_tree.pack(expand = 1, fill = Tix.BOTH, padx = 10, pady = 10,
                      side = Tix.LEFT)
        self.tree = new_tree

        self.control_panel.project = self.project['data']
        self.control_panel.set_modified(1)
        self.control_panel.find('project').show()

        ## TODO: Hack
        ## Check whether local structure engine exists
        
        engine = project.getStructureEngine()
        entity = engine.getSettings().getEntity('local_executable')
        if entity.is_initialized():
            valid = os.path.exists(entity.get())
        else:
            valid = 0

        self.control_panel.find('cns').set_modified(not valid and not new_project)

        self.message('Project loaded.')

    def reload_project(self):
        if self.project['file'] is None:
            return

        if not self.close_project():
            return
        
        self.load_project(self.project['file'], to_recent_files = 0)

    def panels_missing(self, panels):
        EntriesMissingDialog(title = 'Error')
        panels[0].show()
        
    def open_project(self):

        import tkFileDialog

        if not self.close_project():
            return 0
                
        file_types = [("ARIA XML files", "*.xml"),
                      ("All files", "*")]
        
        dialog = tkFileDialog.Open(master = self.root,
                                   filetypes = file_types)
        
        filename = dialog.show()#initialdir=dir, initialfile=base)

        ## cancel button was pressed

        if filename == '':
            return

        self.load_project(filename)
        
    def add_to_recent_files(self, filename):

        import os

        ## make absolute path

        filename = os.path.expanduser(filename)
        filename = os.path.abspath(filename)

        current = self.recent_files
        updated = list(current)

        n_recent_files = self.getSettings()['n_recent_files']

        if filename in updated:
            updated.remove(filename)
        
        updated.insert(0, filename)
        if len(updated) > n_recent_files:
            updated = updated[:n_recent_files]

        fm = self.file_menu

        for s in current:
            index = fm.index(s)
            fm.delete(index)

        for name in updated:
            
            f = lambda s = self, n = name: s.load_recent(n)

            ## TODO: hard-coded: 'index-1' since separator comes
            ## right before 'Exit' item

            index = fm.index('Exit')
            if current:
                index -= 1
                
            fm.insert_command(index, label = name, command = f)

        if not current:
            fm.insert_separator(fm.index('Exit'))

        self.recent_files = updated
            
    def close_project(self):
        
        if self.project['data'] is not None:
            
            if self.control_panel.modified():
                c = NotSavedDialog()

                if c == YES:
                    return self.save_project()

                elif c == CANCEL:
                    return 0

        return 1

    def new_project(self):

        import os, time

        if not self.close_project():
            return 0

        filename = os.path.join(AriaBaseClass.data_path,
                                PROJECT_TEMPLATE)

        self.load_project(filename, relax = 1, new_project = 1)
        self.project['file'] = None

        ## validate panels

        invalid_panels = ['sequence']
        
        ## set 'preferences' as default values

        gui_settings = self.getSettings()

        engine = self.project['data'].getStructureEngine()
        
        try:
            cns_exec = gui_settings['cns_executable']
            engine.getSettings()['local_executable'] = cns_exec

        except:
            invalid_panels.append('cns')

        js_settings = engine.getJobScheduler().getSettings()

        try:
            js_settings['host_list'] = gui_settings['host_list']
        except:
            invalid_panels.append('job_manager')
            
        ## for some settings, the validation mechanism 
        ## had been disabled during the unpickling process.
        ## unless no proper values have been set to these
        ## settings, they remain invalid (i.e. ARIA cannot be run).
        ## to mark the respective panels, we simply flag them
        ## as 'modified'

        invalid_panels.append('project')

        [self.control_panel.find(p).set_modified() for p in invalid_panels]

        ## set some default values

        project_settings = self.project['data'].getSettings()

        project_settings['date'] = time.ctime()
        project_settings['author'] = gui_settings['author']

        ## dict which maps GUISetting entity names to AnalyserSettings
        ## entity names

        d = {'prosa_executable': 'prosa_executable',
             'procheck_executable': 'procheck_executable',
             'whatif_executable': 'whatif_executable'}

        analyser_settings = self.project['data'].getAnalyser().getSettings()
        
        for gui_name, analyser_name in d.items():
            analyser_settings[analyser_name] = gui_settings[gui_name]
        
        self.control_panel.find('project').refresh()
        
        ## invalidate control panel

        self.control_panel.set_modified()

    def dump_project(self, filename):
        from aria.AriaXML import AriaXMLPickler
        
        p = AriaXMLPickler()

        try:
            p.dump(self.project['data'], filename)
            self.control_panel.set_modified(0)
            return 1
        
        except Exception, msg:
            import aria.tools as tools
            print tools.last_traceback()
            print 'Could not save project:', msg
            return 0
        
    def save_project(self):

        if self.project['file'] is None:
            return self.save_project_as()

        panels = self.control_panel.savePanels()
##        panels = self.control_panel.get_modifieda_panels()
        
        if panels:
            self.panels_missing(panels)
            return None

        val = self.dump_project(self.project['file'])

        self.message('Project saved.')

        return val
        
    def save_project_as(self):
        
        import tkFileDialog

        if self.project['data'] is None:
            self.message('No project loaded')
            return 

        panels = self.control_panel.savePanels()
##        panels = self.control_panel.get_modified_panels()
        if panels:
            self.panels_missing(panels)
            return
        
        file_types = [("ARIA project files", "*.xml"),
                      ("All files", "*")]
        
        dialog = tkFileDialog.SaveAs(master = self.root, filetypes=file_types)
        
        filename = dialog.show()#initialdir=dir, initialfile=base)

        if filename == '':
            return

        if not self.dump_project(filename):
            return
        
        self.project['file'] = filename

        self.add_to_recent_files(filename)

        ## set project name as window title

        project_name = self.project['data'].getSettings()['name']
        title = '%s (%s)' % (project_name, filename)
        self.set_title(title)

    def exit(self):


        """Quit our mainloop. It is up to you to call root.destroy() after."""

        if self.project['data'] is not None:

            if self.control_panel.modified():
                c = NotSavedDialog()

                if c == YES:
                    if self.save_project() is None:
                        return
                elif c == 'cancel':
                    return

        ## save preferences
                
        s = self.getSettings()

        s.update(self.control_panel.getSettings())

        s['recent_files'] = self.recent_files

        ## Get geometry of various windows / panes
        
        s['geometry'] = self.root.winfo_geometry()

        tree_size = self.horizontal_pane.panecget('tree', 'size')
        output_size = self.vertical_pane.panecget('list', 'size')

        s['tree_panel_size'] = int(tree_size)
        s['output_panel_size'] = int(output_size)
        s['use_tooltips'] = self.useBalloons.get()
        s['no_main_window_tooltips'] = self.no_tooltips_main_window.get()
        
        try:
            Dump((self.version, s), self.settings_file)
        except:
            print 'Could not not save configration file.'
            
        self.__exit = 0

    def new_panel(self, parent_node, label, constructor = None,
                  menu_setup = None, name = None, open = 0):

        import panels

        new_node = parent_node.new_node(label)

        if constructor is not None:
            panel = constructor(node = new_node)
        else:
            panel = panels.Panel(node = new_node)

        self.control_panel.add_panel(panel, name = name)
        
        new_node.bindings['browse'] = panel.show

        ## add context menu

        if menu_setup is not None:
            popup = menu_setup(self.root)
            popup.bind_widget(new_node.getTree().hlist)

        if open:
            new_node.open()
        else:
            new_node.close()

        ## set panel instance as user-data

        new_node.setUserData(panel)
            
        return new_node, panel

    def add_annealing_sub_panel(self, data_type, params, node_basename,
                                id = None, open = 0):

        import aria.DataContainer as DC

        constructors = {DC.DATA_ANNEALING_RDC: panels.AnnealingRDCPanel,
                        DC.DATA_ANNEALING_FBHW: panels.AnnealingFBHWPanel,
                        DC.DATA_ANNEALING_KARPLUS: panels.AnnealingKarplusPanel,
                        DC.DATA_ANNEALING_DIHEDRAL: \
                        panels.AnnealingDihedralPanel,
                        DC.DATA_ANNEALING_HBOND: panels.AnnealingHBondPanel,
                        DC.DATA_ANNEALING_AMBIG: panels.AnnealingAmbigPanel,
                        DC.DATA_ANNEALING_UNAMBIG: panels.AnnealingUnambigPanel,
                        DC.DATA_ANNEALING_SYM: panels.AnnealingSymmetryPanel,   # BARDIAUX 2.2
                        DC.DATA_ANNEALING_LOGHARMONIC: panels.LogHarmonicPanel} # BERNARD 2.3

        settings = {DC.DATA_ANNEALING_RDC: DC.RDCParameters,
                    DC.DATA_ANNEALING_FBHW: DC.FBHWParameters,
                    DC.DATA_ANNEALING_KARPLUS: DC.KarplusParameters,
                    DC.DATA_ANNEALING_DIHEDRAL: DC.DihedralParameters,
                    DC.DATA_ANNEALING_HBOND: DC.HBondParameters,
                    DC.DATA_ANNEALING_AMBIG: DC.AmbiguousParameters,
                    DC.DATA_ANNEALING_UNAMBIG: DC.UnambiguousParameters,
                    DC.DATA_ANNEALING_SYM: DC.SymmetryParameters,            # BARDIAUX 2.2
                    DC.DATA_ANNEALING_LOGHARMONIC: DC.LogHarmonicParameters} # BERNARD 2.3
        
        master = self.annealing[data_type]

        name = node_basename

        if id is None:
            id = str(data_type)

        node, panel = self.new_panel(master, name, constructors[data_type],
                                     name = id, open = open)
    
        if params is None:
            raise 'this error should not occur'
        
        ## attach settings to respective panel
        
        panel.set_settings(params)

        return panel

    def add_data_container(self, data_type, data, node_basename = '#',
                           id = None, open = 0, show = 0):
        
        if self.project['data'] is None:
            self.message('No project loaded')
            return

        import aria.DataContainer as DC

        constructors = {DC.DATA_SPECTRUM: panels.DataSpectrumPanel,
                        DC.DATA_KARPLUS: panels.DataKarplusPanel,
                        DC.DATA_HBONDS: panels.DataHBondPanel,
                        DC.DATA_DIHEDRALS: panels.DataDihedralPanel,
                        DC.DATA_RDCS: panels.DataRDCPanel,
                        DC.DATA_SSBONDS: panels.DataSSBondPanel,
                        DC.DATA_SEQUENCE: panels.DataSequencePanel,
                        DC.DATA_TEMPLATE_STRUCTURE: panels.DataTemplatePanel,
                        DC.DATA_UNAMBIGUOUS: panels.DataUnambigPanel,
                        DC.DATA_AMBIGUOUS: panels.DataAmbigPanel,
                        DC.DATA_SYMMETRY : panels.SymmetryPanel,
                        DC.DATA_ZNPATCH : panels.DataZnPatchPanel,
                        DC.DATA_OTHER : panels.DataOtherPanel }

        settings = {DC.DATA_SPECTRUM: DC.SpectrumData,
                    DC.DATA_KARPLUS: DC.KarplusData,
                    DC.DATA_HBONDS: DC.HBondData,
                    DC.DATA_DIHEDRALS: DC.DihedralData,
                    DC.DATA_RDCS: DC.RDCData,
                    DC.DATA_SSBONDS: DC.SSBondData,
                    DC.DATA_SEQUENCE: DC.SequenceData,
                    DC.DATA_TEMPLATE_STRUCTURE: DC.TemplateData,
                    DC.DATA_UNAMBIGUOUS: DC.UnambiguousDistanceData,
                    DC.DATA_AMBIGUOUS: DC.AmbiguousDistanceData,
                    DC.DATA_SYMMETRY : DC.Symmetry,
                    DC.DATA_ZNPATCH : DC.ZnPatch,
                    DC.DATA_OTHER : DC.OtherData}

        master = self.data[data_type]

        if data_type == DC.DATA_SEQUENCE:
            name = node_basename
        else:
            name = node_basename + str(len(master.get_children())+1)

        if id is None:
            id = str(data_type)

        node, panel = self.new_panel(master, name, constructors[data_type],
                                     name = id, open = open)


        panel.project = self.project['data']

        if data is None:
            data = settings[data_type]()
            data.reset()
            self.project['data'].addData(data)
            panel.set_modified()
##            self.control_panel.set_modified()
        
        panel.set_settings(data)

        if show:
            node.show()
            panel.show()
            self.control_panel.frame.focus_set()

        return node, panel

    def add_iteration(self, iteration_settings = None, show = 0):

        if self.project['data'] is None:
            self.message('No project loaded')
            return

        name = '' + str(len(self.iterations.get_children()))
        
        if iteration_settings is None:

            from aria.Protocol import IterationSettings

            iteration_settings = IterationSettings()
            iteration_settings.reset()

            if self.project['data'] is not None:
                s = self.project['data'].getProtocol().getSettings()

                it_settings = s['iteration_settings']
                numbers = map(lambda s: s['number'], it_settings.values())
		
                new_num = max(numbers) + 1
                iteration_settings['number'] = new_num

                s.addIterationSettings(iteration_settings)
                
                name = str(new_num)
        
##                self.control_panel.set_modified()
                
        node, panel = self.new_panel(self.iterations, name,
                                     panels.IterationPanel,
                                     name = name)

        panel.set_settings(iteration_settings)

        if show:
            node.show()
            panel.show()
            self.control_panel.frame.focus_set()
            
    #### BARDIAUX
    def add_cmaps(self, iteration_settings = None, show = 0):

        if self.project['data'] is None:
            self.message('No project loaded')
            return
        

        #n_all = 'iteration_' + str(len(self.cmaps.get_children()))
        
        if iteration_settings is None:

            from aria.Protocol import IterationSettings

            iteration_settings = IterationSettings()
            iteration_settings.reset()

            if self.project['data'] is not None:
                s = self.project['data'].getProtocol().getSettings()

                it_settings = s['iteration_settings']
                numbers = map(lambda s: s['number'], it_settings.values())

                new_num = max(numbers) + 1
                iteration_settings['number'] = new_num

                s.addIterationSettings(iteration_settings)
                
                n_all = 'Iteration %s' % str(new_num)
        else:
            n_all = 'Iteration %s' % str(iteration_settings['number'])

        from aria.Protocol import REPORT_NOE_RESTRAINTS
        
        path =  self.project['data'].getInfrastructure().get_iteration_path(iteration_settings['number'])
        
        filename = os.path.join(path, REPORT_NOE_RESTRAINTS + '.pickle')
        
        if not os.path.exists(filename):
##             self.message('No Peak Maps available for iteration ' + str(iteration_settings['number']))
            return
        
        subnode, subpanel= self.new_panel(self.cmaps, n_all,
                                     panels.ContactMapPanel,
                                     name = n_all)

        import AriaViewer
        cmap_settings = AriaViewer.CmapSettings()


        cmap_settings['file'] = filename
        cmap_settings['it'] = iteration_settings['number']
        cmap_settings['selection'] = ""
            
        subpanel.set_settings(cmap_settings)  
        
        
        subnode.close()
        
    def make_main_menu(self):
        
        top = self.root
        w = Tix.Frame(top, bd=2, relief=RAISED)

        ## file menu

        file = Tix.Menubutton(w, text='Project', underline=0, takefocus=0)
        file.pack(side=LEFT)
        fm = Tix.Menu(file, tearoff=0)
        file['menu'] = fm
        
        fm.add_command(label = 'New', command = self.new_project,
                       underline = 0)
        fm.add_command(label = 'Open...', command = self.open_project,
                       underline=0)
        fm.add_command(label = 'Reload', command = self.reload_project,
                       underline = 0)
        fm.add_command(label = 'Save', command = self.save_project,
                       underline = 0)
        fm.add_command(label = 'Save as...', command = self.save_project_as,
                       underline = 1)

        fm.add_separator()
                       
        fm.add_command(label = 'Exit', underline = 1,
                     command = lambda s = self: s.exit())

        self.file_menu = fm

        ## edit menu

        edit = self.create_edit_menu(w)

        ## add menu

        add = self.create_add_menu(w)

        ## help menu

        help = Tix.Menubutton(w, text='Help', underline=0, takefocus=0)
        help.pack(side=RIGHT)
        hm = Tix.Menu(help, tearoff=0)
        help['menu'] = hm
        
        hm.add_checkbutton(label = 'No tooltips for tree window',
                           underline = 0,
                           command = lambda s = self: s.toggleMainTooltips(),
                           variable = self.no_tooltips_main_window)

        hm.add_checkbutton(label = 'Show tooltips', underline = 0,
                           command = lambda s = self: s.ToggleHelp(),
                           variable = self.useBalloons)

        hm.add_separator()
        hm.add_command(label = 'About...', command = AboutBox)
                      
        return w

    def show_preferences(self):
        prefs = self.control_panel.find('preferences')
        
        if prefs is None:

            ## create preferences panel

            import panels

            prefs = panels.PreferencesPanel()

            prefs.set_settings(self.getSettings(), 
                               self.control_panel.getSettings())

            self.control_panel.add_panel(prefs, name = 'preferences')

        prefs.show()

    def create_edit_menu(self, p):
        
        edit = Tix.Menubutton(p, text = 'Edit', underline = 0, takefocus = 0)
        edit.pack(side = LEFT)
        pm = Tix.Menu(edit, tearoff = 0)
        edit['menu'] = pm
        
        pm.add_command(label = 'Preferences...', underline = 0,
                       command = self.show_preferences)

        return edit

    def create_add_menu(self, p):

        from aria.OrderedDict import OrderedDict
        import aria.DataContainer as DC

        add = Tix.Menubutton(p, text = 'Add...', underline = 0)
        add.pack(side = LEFT)
        am = Tix.Menu(add, tearoff = 0)
        add['menu'] = am

        f = lambda TYPE, s = self.add_data_container: \
            lambda : s(TYPE, data = None, node_basename = '#',
                       open = 1, show = 1)

        names = OrderedDict()
        names['Spectrum'] = (f(DC.DATA_SPECTRUM), 0)
        names['Hydrogen bonds restraints'] = (f(DC.DATA_HBONDS), 0)
        names['RDCs'] = (f(DC.DATA_RDCS), 0)
        names['Dihedral angle restraints'] = (f(DC.DATA_DIHEDRALS), 0)
        names['Disulfide bridges'] = (f(DC.DATA_SSBONDS), 0)
        names['Ambiguous distance restraints'] = (f(DC.DATA_AMBIGUOUS), 0)
        names['Unambiguous distance restraints'] = (f(DC.DATA_UNAMBIGUOUS), 0)
        names['Scalar couplings'] = (f(DC.DATA_KARPLUS), 0)
        names['Zinc Coordination'] = (f(DC.DATA_ZNPATCH), 0)
        names['Other data'] = (f(DC.DATA_OTHER), 0)        
        names['separator 1'] = None
        names['Iteration'] = (lambda s = self.add_iteration: s(show = 1), 0)
	
        for name, f in names.items():
            if name.find('separator') >= 0:
                am.add_separator()
            else:
                f, pos = f
                am.add_command(label = name, command = f, underline = pos)

        return am

    def create_control_panel(self, p):
        
        import panels

        default_settings = self.getSettings()
        cp_settings = panels.ControlPanelSettings(default_settings)
        
        cp_settings['bitmap_path'] = os.path.join(AriaBaseClass.install_path,
                                                  PATH_GUI_MODULES)

        return panels.ControlPanel(cp_settings, p, self.balloon)

    def make_main_window(self):
        from widgets import MyPanedWindow

        settings = self.getSettings()

        parent = self.root
        pane = MyPanedWindow(parent, orientation = 'vertical')
        self.vertical_pane = pane

        size = settings['output_panel_size']
        p1 = pane.add('list', min = 10, size = size, allowresize = 1)
        p2 = pane.add('text', min = 10, allowresize = 1)

        p1['relief'] = 'flat'
        p2['relief'] = 'flat'

        import widgets
        self.output_text = widgets.MyScrolledText(p2, scrollbar = 'auto')
        self.output_text.text['wrap'] = 'none'
        self.output_text.pack(expand = 1, fill = Tix.BOTH,
                              padx = 4, pady = 6)

        pane2 = MyPanedWindow(p1, orientation = 'horizontal')
        self.horizontal_pane = pane2
        pane2.pack(expand=1, fill=Tix.BOTH, padx=1, pady=6)

        size = settings['tree_panel_size']
        self.tree_panel = pane2.add('tree', min = 250, size = size)
        self.control_panel = self.create_control_panel(pane2.add('main'))

        self.tree = None
        
        return pane

    def make_main_status(self):
        top = self.root
        w = Tix.Frame(top, relief=Tix.RAISED, bd=1)
        self.statusbar = Tix.Label(w, relief=Tix.SUNKEN, bd=1)
        self.statusbar.form(pady = 3, padx = 3, left=0, right='%70')
        
        return w

    def patch_tk_options(self):

        op_filename = os.path.expanduser(GUI.tkoptions_filename)

        lines = TK_OPTIONS_TEMPLATE % self.getSettings()
        
        f = open(op_filename, 'w')
        f.write(lines)
        f.close()

        self.root.option_readfile(op_filename)
        self.root.tk_strictMotif(0)

    def bindTooltips(self, widget, msg):
        from aria.tools import wrap_string

        msg = wrap_string(msg, AriaBaseClass.description_length)
        self.balloon.bind_widget(widget, msg = msg)

    def unbindTooltips(self, widget):
        try:
            self.balloon.unbind_widget(widget)
        except:
            pass

    def build(self):

        z = self.root.winfo_toplevel()
        self.patch_tk_options()

        self.balloon = Tix.Balloon(self.root)
        main_menu = self.make_main_menu()
        status_line = self.make_main_status()
        main_window = self.make_main_window()
        main_menu.pack(side=TOP, fill=X)

        status_line.pack(side=BOTTOM, fill=X)
        main_window.pack(side=TOP, expand=1, fill=BOTH, padx=4, pady=4)
##        self.balloon['statusbar'] = self.statusbar
        z.wm_protocol("WM_DELETE_WINDOW", lambda self=self: self.exit())
 
        ## enable tooltips

        use_balloons = self.getSettings()['use_tooltips']
        self.useBalloons.set(use_balloons)
        self.balloon['state'] = {'0': 'none', '1': 'both'}[use_balloons]

        ## enable toolstips for main windows

        b = self.getSettings()['no_main_window_tooltips']

        self.no_tooltips_main_window.set(b)
        self.toggleMainTooltips()

        self.set_title()
        
    def loop(self):
        import tkMessageBox, traceback
        while self.__exit < 0:
            try:
                self.root.tk.dooneevent(TCL_ALL_EVENTS)
            except SystemExit:
                self.__exit = 1
                break
            except KeyboardInterrupt:
                if tkMessageBox.askquestion ('Interrupt',
                                             'Really Quit?') == YES:
                    return
                else:
                    pass
                continue
            except:
                t, v, tb = sys.exc_info()
                text = ""
                for line in traceback.format_exception(t,v,tb):
                    text = text + line + '\n'
                try: tkMessageBox.showerror ('Error', text)
                except: pass
                tkinspect_quit (1)

    def destroy (self):
        self.root.destroy()
    
    def ToggleHelp(self):
        if self.useBalloons.get() == '1':
            self.balloon['state'] = 'both'
        else:
            self.balloon['state'] = 'none'

    def toggleMainTooltips(self):

        main_windows = {'log_window': self.output_text,
                        'tree_panel': self.tree_panel}

        [self.unbindTooltips(w) for w in main_windows.values()]
        
        if self.no_tooltips_main_window.get() == '0':
            for key, w in main_windows.items():
                msg = self.messages[key]
                self.bindTooltips(w, msg)

def go(project_file = None):
    
    gui_path = os.path.join(AriaBaseClass.install_path,
                            PATH_GUI_MODULES)

    # BARDIAUX extendNmr
    # if a tk already exists
    # open ARIA gui as Toplevel (with Tix)
    import Tkinter
    if Tkinter._default_root:
        import widgets
        root = widgets.ariaPopup()
    else:
        root = Tix.Tk()

    app = GUI(root)
    AriaBaseClass.display_deprecated = 0
    app.load_settings()

    app.build()
    app.setup()

    AriaBaseClass.log_stdout = None
    AriaBaseClass.log_gui = app.output_text
    
    if project_file is not None:
        app.load_project(project_file)

    ## redirect stdout to app-window

    #sys.stdout = app
    

    app.loop()
    app.destroy()
    
    AriaBaseClass.log_gui = None
    
if __name__ == '__main__':

    aria_root = AriaBaseClass.install_path

    modules = (PATH_MODULES, PATH_GUI_MODULES)

    for m in modules:
        
        p = os.path.join(aria_root, m)

        if not p in sys.path:
            sys.path.insert(0, p)

    if len(sys.argv) > 1:
        project_file = sys.argv[1]
    else:
        project_file = None
        
    go(project_file)



