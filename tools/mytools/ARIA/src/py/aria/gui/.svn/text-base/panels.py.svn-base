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


import Tix, Tkinter
from decorators import *
from aria.ariabase import *
from aria.Settings import Settings

DATA_TEXT = 'text'
DATA_PDB = 'pdb'

YES_NO_DICT = [{'label': 'Yes', 'value': YES}, \
               {'label': 'No', 'value': NO}]

YES_NO_GZIP_DICT = YES_NO_DICT + [{'label': 'GZip', 'value': 'gzip'}]

def file_and_format(p, format_choices=None, row=0, file_types = None,
                    label_filename = 'Filename: ', label_format = 'Format: ',
                    label_ccpn_id = 'CCPN id: ', with_ccpn_id=0,
                    with_format=1, ccpn_select_func=None):
    """
    file_types format: list of tuples: [(description, filename wildcard),...]

    example: [('Text files', '*.txt')]
    """

    rel_row = row
    result = []

    frame = Tix.Frame(p)
    
    label = Tix.Label(frame, text = label_filename)
    label.grid(row=rel_row, column=0, sticky = Tix.W)

    entry = DFileEntry(frame, file_types = file_types, width = 10, \
                       text = '')
    
    entry.entry.configure(width = 40)
    entry.grid(row=rel_row, column=1, columnspan=2, sticky = Tix.E)

    rel_row += 1
    result.append(entry)

    ## format

    if with_format:

        l = Tix.Label(frame, text = label_format)
        l.grid(row=rel_row, column=0, sticky = Tix.W)

        format = DOptionMenu(frame, label = '')

        for format_choice in format_choices:
            value = format_choice['value']
            label = format_choice['label']
            format.add_command(value, label = label)

        format.grid(row=rel_row, column=2, sticky = Tix.E)

        rel_row += 1
        result.append(format)

    ## widgets for CCPN access

    if with_ccpn_id:

        from widgets import DefaultButton as Button

        l = Tix.Label(frame, text = label_ccpn_id)
        l.grid(row=rel_row, sticky=Tix.W)

        ccpn_id = DEntry(frame, width = 25)
        ccpn_id.grid(row=rel_row, column=1, sticky=Tix.W+Tix.E)
        
        select = Button(frame, text = 'Select...', command=ccpn_select_func)
        select.grid(row=rel_row, column=2, sticky=Tix.E+Tix.W)

        rel_row += 1
        
        label_name = Tix.Label(frame, text = "CCPN Name: ")
        label_name.grid(row=rel_row, sticky=Tix.W)

        name = Tix.Label(frame, text = "-")
        name.grid(row=rel_row,  column=1, columnspan=2, sticky=Tix.W)

        rel_row += 1
        result.append(ccpn_id)

        result.append(name)

    frame.grid(row=row, sticky = Tix.W, columnspan=2)

    return result

class Panel(AriaBaseClass):

    def __init__(self, node = None):
        self.__widgets = {}
        self.__owner = None
        self.__initialized = 0
        self.__node = node
        self.__frame = None
        self.set_modified(0)
	
	## BARDIAUX 2.2
	## access message and error functions
	self._name = "GUI"
	AriaBaseClass.__init__(self)
        
    def setup(self):
        
        ## TODO: hard-coded __owner.frame

        self.__frame = Tix.Frame(self.__owner.frame)
        self.__widgets = self.create(self.__frame)
        
        return self.__frame

    def create(self, parent):
        pass

    def set_owner(self, p):
        self.__owner = p

    def get_owner(self):
        return self.__owner

    def get_node(self):
        return self.__node

    def _get_widgets(self):
        if self.__widgets is None:
            return []
        else:
            return reduce(lambda x, y: x + y, self.__widgets.values())

    def show(self):
        """
        check_panel: if non-zero, the panels entry fields are
        checked for validity. 
        """
        if not self.__owner.change_panel(self):
            return

        ## create widgets if necessary and
        ## attach entites
        if not self.__initialized:
            self.__frame = self.setup()
            self.update()

##            self.validate()
            
            self.__initialized = 1

        self.__owner.set_panel(self, self.__frame)
        
        ## disable control-panels delete button

        self.__owner.button_delete.config(state = Tix.DISABLED)

        if self.get_node() is not None:
            self.get_node().select()

        ## set focus to first widget in list

        w = self._get_widgets()

        if w:
            w[0]['widget'].focus_set()

    def validate(self):
        """
        checks the validity of the panel, i.e.
        whether the values stored in the widgets are
        returns 0 if panel had been saved, 1 otherwise
        """

        ## compile list of all widgets

        widgets = self._get_widgets()
        
        if not widgets:
            return self.set_modified(0)

        try:
            map(lambda w: w['widget'].save(), widgets)
            return self.set_modified(0)
            
        except:
            return self.set_modified()

    def save(self):

        from widgets import OKDialog

        ## compile list of all widgets

        widgets = self._get_widgets()
        
        if not widgets:
            return 1

        invalid_widget = None

        for widget in widgets:
            widget = widget['widget']

            try:
                widget.save()
                
            except Exception, msg:
                if invalid_widget is None:
                    exc_save, msg_save = Exception, msg
                    invalid_widget = widget

        if invalid_widget is not None:
            print invalid_widget
            invalid_widget.highlight()
            
            descr = invalid_widget.get_entity().getDescription()

            text = '%s:\n%s' % (descr, msg_save)
            
            dialog = OKDialog(title = 'Value error', text = text,
                              modal = 1)
            
            invalid_widget.highlight()
            invalid_widget.focus_set()

            return self.set_modified()

        ## Panel is valid
        self.set_modified(0)
        
        ## Return 1 if panel has been saved
        return 1

    def undo(self):

        ## compile list of all widgets

        widgets = self._get_widgets()
        if not widgets:
            return

        for widget in widgets:

            widget = widget['widget']

            try:
                widget.undo()
            except Exception, msg:
                dialog = OKDialog(title = 'Should not occur: Value error',
                                  text = msg, modal = 1)
                widget.highlight()
                break

    def focus_get(self):
        
        widgets = self._get_widgets()
        if not widgets:
            return

        widgets = [w['widget'] for w in widgets]
        
        current_tk_widget = self.__frame.focus_get()

        current_widget = None

        for widget in widgets:
            if widget.get_tk_widget() is current_tk_widget:
                current_widget = widget
                break

        return current_widget

    def reset_widget(self):
        """
        resets the widgets which has the current focus to
        its default value.
        """

        try:
            self.focus_get().reset()
        except:
            pass

    def reset(self):

        ## compile list of all widgets

        widgets = self._get_widgets()
        
        if not widgets:
            return
        
        for widget in widgets:

            try:
                widget['widget'].reset()
            except:
                pass

        ## some widgets may not have a default value.
        ## if some widgets are out of synch with their
        ## attached entities, the panel is marked as
        ## modified.

        self.validate()

    def refresh(self, overwrite=1):
        """
        loops through all widgets and attempts to
        load the data from the entity into the widget.
        data cannot be loaded, if the respective entity
        has never been initialized.
        """
        widgets = self._get_widgets()

        for widget in widgets:

            widget = widget['widget']

            if not overwrite and widget.has_changed():
                continue
            
            try:
                widget.load()
            except:
                pass

    def set_settings(self):
        pass

    def connect(self, widget, entity):

        check_type(entity, 'Entity')

        widget.set_entity(entity)
        widget.load()

    def get_settings(self):
        """
        used by the method 'update' to implement the
        default updating procedure. the method must return
        a dict of (name, settings) items which shall connected to the
        panels widgets.
        """

        return {}

    def update_widgets(self, s, use_name = None):
        """
        by default, all widgets are connected to the
        settings which are return from get_settings.
        """

        if self.__widgets is None:
            return

        if use_name is None:
            widgets = reduce(lambda x, y: x + y, self.__widgets.values())
        else:
            widgets = self.__widgets[use_name]

        for d in widgets:
            widget = d['widget']
            entity = d['entity']
            entity = s.getEntity(entity)
            self.connect(widget, entity)
            self.__owner.bind_tooltips(widget, entity.getDescription())

    def update(self):

        if self.__initialized:
            return

        settings = self.get_settings()
        for name, s in settings.items():
            self.update_widgets(s, name)

    def isInitialized(self):
        return self.__initialized

    def set_modified(self, m = 1):

        import os

        self.__modified = m
        
        if self.__owner is not None:
            s = self.__owner.getSettings()            
            if m:
                self.__owner.set_modified()
        else:
            s = None

        if self.get_node() is not None and s is not None:
            if m:
                name = os.path.join(s['bitmap_path'], 'warning')
            else:
                name = os.path.join(s['bitmap_path'], 'shaded_folder')

            self.get_node().setIcon(name)

        return not m

    def modified(self):
        return self.__modified

    def delete(self, project):
        pass
        
    def _delete(self, project):

        self.delete(project)
        prev = self.__node.get_previous()

        ## select prev. node

        if prev is not None:
            prev.select()
            prev.getUserData().show()

        self.__frame.forget()
        
        ## remove attached node

        self.__node.remove()

    def find_widget(self, entity_name):
        l = [w['widget'] for w in self._get_widgets() \
             if w['entity'] == entity_name]

        if len(l) == 1:
            return l[0]

        if not l:
            s = 'Could not find widget connected to entity %s.'
            self.error(KeyError, s % entity_name)

        else:
            s = 'More than one widget found for entity %s: %s'
            self.error(KeyError, s % (entity_name, str(l)))
        
class PanelEx(Panel):

    from aria.TypeChecking import STRING, INT

    widget_constructors = {'entry': DEntry,
                           'choice': DOptionMenu,
                           STRING: DOptionMenu,
                           INT: DIntOptionMenu,
                           'switch': DCheckbutton,
                           'file': DFileEntry,
                           'path': DPathEntry,
                           'ccpn_browser' : CCPNProjectBrowser,
                           'listbox': DListBox,
                           'text': DText,
                           'combo': DComboBox}

    def _create_widget(self, p, name, d, row, column, columnspan):

        if 'columnspan' in d:
            column_span = d['columnspan']
        else:
            column_span = columnspan

        widget_name = d['widget']

        ## Label and widget live on a separate frame

        if widget_name not in ('switch',):
            label = Tix.Label(p, text = name)
            label.grid(row = row, column = column - 1, sticky = Tix.W)
        else:
            label = None

        if widget_name is None:
            return None, label

        if not widget_name in self.widget_constructors:
            raise ValueError, 'Unknown widget "%s"' % widget_name

        constructor = self.widget_constructors[widget_name]

        if widget_name == 'entry':
            widget = constructor(p, width = 10)
            widget.grid(row = row, column = column, sticky = Tix.EW,
                        columnspan = column_span)

        elif widget_name == 'choice':

            if 'type' in d:
                constructor = self.widget_constructors[d['type']]

            widget = constructor(p, label = '')

            for choice in d['choices']:
                widget.add_command(choice['value'],
                                   label = choice['label'])
            widget.grid(row = row, column = column, sticky = Tix.E,
                        columnspan=column_span)

        elif widget_name == 'switch':
            choices = d['choices']
            widget = constructor(p, text = name,
                                 onvalue = choices['yes'],
                                 offvalue = choices['no'])
            widget.grid(row = row, column = column - 1, sticky = Tix.W,
                        columnspan = column_span)

        elif widget_name == 'combo':
            widget = constructor(p)
            widget.set_items(d['items'])
            widget.grid(row = row, column = column, sticky = Tix.W,
                        columnspan = column_span)

        elif widget_name == 'file':
            if 'file_types' in d:
                f_types = d['file_types']
            else:
                f_types = None

            widget = constructor(p, file_types = f_types, label = '')
            widget.entry.configure(width = 40)
            widget.grid(row = row, column = column, sticky = Tix.EW,
                        columnspan = column_span)

        elif widget_name == 'path':
            widget = constructor(p)
            widget.entry.configure(width = 40)
            widget.grid(row = row, column = column, sticky = Tix.EW,
                        columnspan = column_span)

        elif widget_name == 'listbox':
            widget = constructor(p)
            widget.grid(row = row, column = 1, sticky = Tix.EW)

        elif widget_name == 'text':
            widget = constructor(p, height = d['height'], width = 40,
                                 wrap='word')
            widget.grid(row = row, column = column, sticky = Tix.EW,
                        columnspan = column_span)
            
        elif widget_name == 'ccpn_browser':
            widget = constructor(p)
            widget.entry.configure(width = 40)
            widget.grid(row = row, column = column, sticky = Tix.EW,
                        columnspan = column_span)
            
        return widget, label

    def create_widgets(self, p, definitions, row = 0, columnspan = 1):
        
        controls = []
        i = row

        for name, d in definitions.items():
            if type(d) <> type(()):
                d = (d,)

            column = 1

            for widget_descr in d:

                if 'name' in widget_descr:
                    widget_name = widget_descr['name']
                else:
                    widget_name = name

                widget, label = self._create_widget(p, widget_name,
                                                    widget_descr,
                                                    i, column, columnspan)

                if widget is not None:
                
                    entity = widget_descr['entity']

                    descr = {'widget': widget, 'entity': entity,
                             'label': label}

                    controls.append(descr)

                column += 2
                
            i += 1

        return controls
    
    ############ bardiaux #############

    def create_text_item(self, who, pos, col, char, tag, bind=1):

        x = who.create_text(pos[0],pos[1],text = char,tags = tag,fill = col)

        if bind:
            who.tag_bind(tag,'<Button-1> ',func=self.circPeak)
            who.tag_bind(tag,'<Shift-Button-1> ',func=self.circAddPeak)
            who.tag_bind(tag,'<Button-3> ',func=self.circDel)        

        return x

class ProjectPanel(PanelEx):

    def get_project_defs(self):

        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['Name:'] = {'widget': 'entry', 'entity': 'name'}
        widgets['Version:'] = {'widget': 'entry', 'entity': 'version'}
        widgets['Author:'] = {'widget': 'entry', 'entity': 'author'}
        widgets['Date:'] =  {'widget': 'entry', 'entity': 'date',
                             'columnspan': 1}
        widgets['Description:'] = {'widget': 'text', 'height': 5,
                                   'entity': 'description'}
        widgets['Comment:'] = {'widget': 'text', 'height': 5,
                               'entity': 'comment'}
        widgets['References:'] = {'widget': 'text', 'height': 5,
                                  'entity': 'references'}

        widgets['Working directory:'] = {'widget': 'path',
                                         'entity': 'working_directory'}
        widgets['File root:'] = {'widget': 'entry', 'entity': 'file_root'}
        widgets['Temporary path:'] = {'widget': 'path', 'entity': 'temp_root'}
        widgets['Run nickname:'] = {'widget': 'entry', 'entity': 'run'}
        widgets['Cache files:'] = {'widget': 'choice', 'entity': 'cache',
                                   'choices': YES_NO_DICT}
        
        widgets['Cleanup:'] = {'widget': 'choice', 'entity': 'cleanup',
                               'choices': YES_NO_DICT}
        return widgets

    def set_date(self):
        import time

        self.get_settings()['project']['date'] = time.ctime()
        widget = self.find_widget('date')
        widget.load()
        self.set_modified()

    def create(self, p):

        from widgets import DefaultButton as Button

        CW = self.create_widgets

        widgets = {}

        generic = Tix.LabelFrame(p, label='Generic')
        defs = self.get_project_defs()
        widgets['project'] = CW(generic.frame, defs, columnspan = 2)

        Button(generic.frame, text = 'Today',
               command = self.set_date).grid(row = 3, column = 2,
                                             sticky = Tix.EW)

        generic.form(top = 0, left = 0, right = -1, fill = Tix.X)
        
        return widgets

    def set_settings(self, project):
        
        d = {}
        d['project'] = project.getSettings()
        
        self.__settings = d

    def get_settings(self):
        return self.__settings

class ReportPanel(PanelEx):

    def get_noe_restraint_list_defs(self):

        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['XML output:'] = {'widget': 'choice',
                                  'choices': YES_NO_GZIP_DICT,
                                  'entity': 'xml_output'}
        
        widgets['Text output:'] = {'widget': 'choice',
                                  'choices': YES_NO_GZIP_DICT,
                                  'entity': 'text_output'}
        
        widgets['Python pickle output:'] = \
                        {'widget': 'choice',
                         'choices': YES_NO_GZIP_DICT,
                         'entity': 'pickle_output'}

        return widgets

    def get_ccpn_defs(self):
        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['CCPN data-model output:'] = \
                      {'widget': 'choice',
                       'choices': YES_NO_GZIP_DICT,
                       'entity': 'enabled'}
        
        return widgets

    def get_molmol_defs(self):
        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['Restraint files:'] = \
                      {'widget': 'choice',
                       'choices': YES_NO_GZIP_DICT,
                       'entity': 'enabled'}
        
        return widgets
    
    ## BARDIAUX
    
    def get_updated_spectrum_defs(self):
        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['Export to ARIA XML:'] = \
                      {'widget': 'choice',
                       'choices': YES_NO_DICT,
                       'entity': 'write_assigned'}

        widgets['Force export:'] = \
                      {'widget': 'choice',
                       'choices': YES_NO_DICT,
                       'entity': 'write_assigned_force'}
        
        widgets['Iteration:'] = \
                              {'widget': 'choice',
                               'choices': ({'label': 'All', 'value': 'all'},
                                           {'label': 'Last', 'value': 'last'}),
                               'entity': 'iteration'}

        widgets['Unambiguous peaks only:'] = \
                              {'widget': 'choice',
                               'choices': YES_NO_DICT,
                               'entity': 'write_unambiguous_only'}
        return widgets

    def create(self, p):

        CW = self.create_widgets

        widgets = {}

        ## NOE restraints

        noe_restraint_list = Tix.LabelFrame(p, label='NOE restraint lists')
        defs = self.get_noe_restraint_list_defs()
        widgets['noe_restraint_list'] = CW(noe_restraint_list.frame, defs)

        ## BARDIAUX Updated spectrum

        updated_spectrum = Tix.LabelFrame(p, label='Assigned spectra')
        defs = self.get_updated_spectrum_defs()
        widgets['spectra'] = CW(updated_spectrum.frame, defs)

##         ## CCPN

##         f_ccpn = Tix.LabelFrame(p, label='Export to CCP data-model')
##         defs = self.get_ccpn_defs()
##         widgets['ccpn'] = CW(f_ccpn.frame, defs)

        ## MOLMOL

        f_molmol = Tix.LabelFrame(p, label='MOLMOL restraints')
        defs = self.get_molmol_defs()
        widgets['molmol'] = CW(f_molmol.frame, defs)

        noe_restraint_list.form(top = 0, left = 0, right = -1, fill=Tix.X)
        # f_ccpn.form(top = noe_restraint_list, left = 0, right = -1, fill=Tix.X)
##        f_ccpn.form(top = updated_spectrum, left = 0, right = -1, fill=Tix.X)
##        f_molmol.form(top = f_ccpn, left = 0, right = -1, fill = Tix.X)
        f_molmol.form(top = noe_restraint_list, left = 0, right = -1, fill = Tix.X)
        updated_spectrum.form(top = f_molmol, left = 0, right = -1, fill=Tix.X)

        return widgets

    def set_settings(self, project):

        from aria.Report import UpSpecSettings

        s = project.getReporter()

        d = {}
        d['noe_restraint_list'] = s['noe_restraint_list']
##         d['ccpn'] = s['ccpn']
        d['molmol'] = s['molmol']
        
        ## BARDIAUX

        if 'spectra' in s:
            val = s['spectra']
        else:
            val = UpSpecSettings()
            val.reset()

        d['spectra'] = val
        
        self.__settings = d

    def get_settings(self):
        return self.__settings

class DataPanel(PanelEx):

    def __init__(self, node = None, enable_delete_button = 1):
        PanelEx.__init__(self, node)
        self.__enable_button = enable_delete_button
    
    def read_ccpn_model(self):

        import os
        from widgets import MessageBox

        filename = self.project.ccpn_model['filename']

        err_msg = None
        
        if not os.path.exists(filename):

            if filename == '':
                err_msg = 'No CCPN project specified (see Node "CCPN data model").'
            else:
                err_msg = 'File "%s" not found or cannot be read.' % filename

        else:

            msg = 'Accessing CCPN data model "%s".\n\nPlease be patient ...' % \
                  self.project.ccpn_model['filename']

            dialog = MessageBox(msg)

            if self.project.ccpn_project is None:

                from memops.general.Io import loadProject

                try:
                    self.project.ccpn_project = loadProject(self.project.ccpn_model['filename'])


                except Exception, msg:
                    err_msg='Could not access CCPN project "%s".\n\nError message was:\n\n %s' % (filename, msg)

            dialog.destroy()

            if err_msg is None:
                return self.project.ccpn_project

        from widgets import OKDialog

        OKDialog(self.parent, 'Error', err_msg, modal=1)

    def edit_source(self, name, data_type):

        try:
            filename = self.get_settings()[name]['filename']
            
        except ValueError, msg:
            s = 'Could not access source-file. Please save ' + \
                'your data ("Commit" Button) and try again.'
            #print str(ValueError), msg
            #print s
	    self.warning(str(ValueError) + ' ' + str(msg))
	    self.warning(s)
            return
            
        self.get_owner().view(filename, data_type)

    def add_edit_button(self, p, row, name, text = 'Edit source file...',
                        data_type = DATA_TEXT):
        """
        'name' is used to select the correct settings from
        the instance' settings dict (get_settings())
        """

        from widgets import DefaultButton as Button

        b = Button(p, text = text,
                   command = lambda s = self, n = name,
                   dt = data_type: s.edit_source(n, dt))
        
        b.grid(row = row, columnspan = 2, sticky = Tix.W)
        
    def delete(self, project):

        for s in self.get_settings().values():
            project.delData(s)

    def show(self):
        PanelEx.show(self)

        ## enable control-panels delete button

        if self.__enable_button and self.get_owner() is not None:
            self.get_owner().button_delete.config(state = Tix.NORMAL)

        if self.get_owner().getSettings()['show_ccpn_names'] == YES:
            self.add_ccpn_names()

    def add_ccpn_names(self):
        pass

    def update_ccpn_name(self, widget, name):

        owner_settings = self.get_owner().getSettings()
        
        if name is None or owner_settings['show_ccpn_names'] <> YES:
            return

        font = widget.option_get('font', "Font")
        font += ' bold'
        
        widget.config(text=name, font=font)
            
class DataSpectrumPanel(DataPanel):

    def add_ccpn_names(self):

        s = self.get_settings()

        if not s['shifts']['ccpn_id'] or not s['peaks']['ccpn_id']:
            return
        
        from aria.importFromCcpn import getCcpnShiftList, getCcpnPeakList, getKeysFromString

        ccpn_project = self.project.ccpn_project
        
        if not ccpn_project:            
            ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return
        
        if s['shifts']['ccpn_id']:
            shift_list = getCcpnShiftList(ccpn_project, getKeysFromString(s['shifts']['ccpn_id']))
            self.update_ccpn_name(self.ccpn_names['shifts'], shift_list.name)
            
        if s['peaks']['ccpn_id']:
            peak_list = getCcpnPeakList(ccpn_project, getKeysFromString(s['peaks']['ccpn_id']))
            exp_name = peak_list.dataSource.experiment.name
            name = "%s:%s:%d" % (exp_name, peak_list.dataSource.name, peak_list.serial)
            self.update_ccpn_name(self.ccpn_names['peaks'], name)


    def select_shifts(self):

        from ccpnGui import gui_select_shiftList
        from aria.importFromCcpn import getShiftLists, getObjectKeyString
        from ccpnmr.format.general.Util import createSelection

        ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return

        shift_lists = getShiftLists(ccpn_project)

        selection_list, selection_dict = createSelection(shift_lists)

        shift_list = gui_select_shiftList(selection_list,
                                          selection_dict, gui=self.parent)

        if shift_list is None:
            return

        m = 'Chemical shift list "%s" is valid for the following experiments: %s' % (shift_list.name, ', '.join([a.name for a in shift_list.experiments]))
	self.message(m)

        settings = self.get_settings()['shifts']

        settings['ccpn_id'] = getObjectKeyString(shift_list)

        self.update_ccpn_name(self.ccpn_names['shifts'], shift_list.name)
        
        self.refresh(overwrite=0)


    def select_peaks(self):
        
        from ccpnGui import gui_select_peakList
        from aria.importFromCcpn import getKeysFromString, getObjectKeyString, getNoesyPeakLists, getCcpnChain, getCcpnChains
        from ccpnmr.format.general.Util import createSelection
        import aria.DataContainer as DC

        ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return

        mol_data = self.project.getData(DC.DATA_SEQUENCE)[0]

        molsystem_id = mol_data['ccpn_id']

        if molsystem_id == '':

            from widgets import OKDialog

            msg = 'Name of molecular system has not yet been specified (see entry field "CCPN name" in node "Sequence").'

            OKDialog(title='Error', text=msg)

            return

        chain = getCcpnChains(ccpn_project, getKeysFromString(molsystem_id))[0]

        peakLists = getNoesyPeakLists(ccpn_project, chain.molSystem)

        if not peakLists:

            from widgets import OKDialog

            msg = 'No CCPN peak or shift list present for\nmolecular system "%s".' % molsystem_id

            OKDialog(title='Error', text=msg)

            return


        peak_list_labels, peak_list_dict = createSelection(peakLists)
        peak_list = gui_select_peakList(peak_list_labels,
                                        peak_list_dict,
                                        gui=self.parent)

        if peak_list is None:
            return

        # BARDIAUX 2.2
        self.setCcpnExperimentData(peak_list)
        
        exp_name = peak_list.dataSource.experiment.name

        key = getObjectKeyString(peak_list)

        settings = self.get_settings()['peaks']

        settings['ccpn_id'] = key

        m = 'Using experiment=%s, spectrum=%s, peaklist=%d from key %s.' % (exp_name, peak_list.dataSource.name, peak_list.serial, key)
	self.message(m)

        ccpn_name = "%s:%s:%d" % (exp_name, peak_list.dataSource.name, peak_list.serial)
        self.update_ccpn_name(self.ccpn_names['peaks'], ccpn_name)

        self.refresh(overwrite=0)

    # BARDIAUX 2.2
    def setCcpnExperimentData(self, peakList):
        
        from aria.importFromCcpn import getCcpnExperimentData
        freq, mixing  = getCcpnExperimentData(peakList)

        settings = self.get_settings()['experiment_data']
        if freq is not None and settings['spectrometer_frequency'] == 0.:
            settings['spectrometer_frequency'] = freq
            
        if mixing is not None and settings['spectrum_mixing_time'] == 0.:
            settings['spectrum_mixing_time'] = mixing

        self.refresh(overwrite=0)
    

    def delete(self, project):
        project.delData(self.get_settings()['generic'])


    def set_settings(self, s):

        check_type(s, 'SpectrumData')

        self.__settings = {'shifts': s['shifts'],
                           'peaks': s['peaks'],
                           'generic': s,
                           'lower_bound_correction': \
                           s['peaks']['lower_bound_correction'],
                           'upper_bound_correction': \
                           s['peaks']['upper_bound_correction'],
                           'experiment_data' : s['experiment_data']} # BARDIAUX rMat

    def create(self, p):


        from aria.OrderedDict import OrderedDict

        self.parent = p

        CW = self.create_widgets

        data_formats = ({'label': 'XML', 'value': 'xml'},
                        {'label': 'CCPN', 'value': 'ccpn'})

        file_types = [('XML files', '*.xml'),
                      ('All files', '*')]

        ## generic:

        generic = Tix.LabelFrame(p, label = 'Generic')

        defs = OrderedDict()

        defs['Use Spectrum:'] = {'widget': 'choice', 'entity': 'enabled',
                                 'choices': YES_NO_DICT}
        
        defs['Use manual Assignments:'] = {'widget': 'choice',
                                           'entity': 'use_assignments',
                                           'choices': YES_NO_DICT}
        
        defs['Trust assigned peaks:'] = {'widget': 'choice',
                                         'entity': 'trust_assigned_peaks',
                                         'choices': YES_NO_DICT}
        
        defs['Filter diagonal peaks:'] = {'widget': 'choice',
                                         'entity': 'filter_diagonal_peaks',
                                         'choices': YES_NO_DICT}
        
        defs['Only fully assigned peaks:'] = {'widget': 'choice',
                                            'entity': 'filter_unassigned_peaks',
                                            'choices': YES_NO_DICT}
        
        ## BARDIAUX 2.2 : Structural rules for multimers
 
        defs['Enable structural rules :'] =  {'widget': 'choice', 'entity': 'structural_rules_enabled',
                                              'choices': YES_NO_DICT}
        
        d = {}
        d['generic'] = CW(generic.frame, defs)

        n_widgets = len(d['generic'])

        defs = OrderedDict()
        defs['Lower bound correction:'] = \
                    ({'widget': 'entry', 'entity': 'value'},
                     {'widget': 'choice', 'entity': 'enabled',
                      'choices': YES_NO_DICT, 'name': 'Enabled:'})

        d['lower_bound_correction'] = CW(generic.frame, defs,
                                         row = n_widgets)
        
        n_widgets += len(d['lower_bound_correction'])

        defs = OrderedDict()
        defs['Upper bound correction:'] = \
                    ({'widget': 'entry', 'entity': 'value'},
                     {'widget': 'choice', 'entity': 'enabled',
                      'choices': YES_NO_DICT, 'name': 'Enabled:'})

        d['upper_bound_correction'] = CW(generic.frame, defs,
                                         row = n_widgets)

        ## shifts

        shifts = Tix.LabelFrame(p, label = 'Chemical Shift List')
        file, format, ccpn_id, ccpn_name = file_and_format(shifts.frame, data_formats,
                                                file_types = file_types,
                                                with_ccpn_id=1, ccpn_select_func=self.select_shifts)

        self.ccpn_names = {}
        self.ccpn_names['shifts'] =  ccpn_name

        self.add_edit_button(shifts.frame, 1, 'shifts')
        
        defs = OrderedDict()
        
        defs['Default chemical-shift error [ppm]:'] = \
                      {'widget': 'entry', 
                       'entity': 'default_shift_error'}

        d['shifts'] = CW(shifts.frame, defs, row = 2)
        
        d['shifts'] += [{'widget': file, 'entity': 'filename'},
                        {'widget': format, 'entity': 'format'},
                        {'widget': ccpn_id, 'entity': 'ccpn_id'}]
        ## peaks

        peaks = Tix.LabelFrame(p, label = 'Spectrum')

        file, format, ccpn_id, ccpn_name = file_and_format(peaks.frame, data_formats,
                                                file_types = file_types,
                                                with_ccpn_id=1, ccpn_select_func=self.select_peaks)

        self.add_edit_button(peaks.frame, 1, 'peaks')

        d['peaks'] = [{'widget': file, 'entity': 'filename'},
                      {'widget': format, 'entity': 'format'},
                      {'widget': ccpn_id, 'entity': 'ccpn_id'}]
        
        defs = OrderedDict()
        defs['Proton1 freq. window [ppm]:'] = {'widget': 'entry',
                                              'entity': 'proton1_shift_err'}
        defs['Hetero1 freq. window [ppm]:'] = {'widget': 'entry',
                                              'entity': 'hetero1_shift_err'}
        defs['Proton2 freq. window [ppm]:'] = {'widget': 'entry',
                                              'entity': 'proton2_shift_err'}
        defs['Hetero2 freq. window [ppm]:'] = {'widget': 'entry',
                                              'entity': 'hetero2_shift_err'}
        defs['Peak type:'] = {'widget': 'choice',
                              'entity': 'volume_or_intensity',
                              'choices': ({'label': 'Volume',
                                           'value': 'volume'},
                                          {'label': 'Intensity',
                                           'value': 'intensity'})}

        self.ccpn_names['peaks'] =  ccpn_name

        d['peaks'] += CW(peaks.frame, defs, row = 2)

        # BARDIAUX rMat
        ## Experiment data
        experiment = Tix.LabelFrame(p, label = 'Experiment')
        
        defs = OrderedDict()
        defs['Molecule correlation time [ns]:'] = {'widget': 'entry',
                                                   'entity': 'molecule_correlation_time'}
        defs['Spectrometer frequency [MHz]:'] = {'widget': 'entry',
                                                 'entity': 'spectrometer_frequency'}
        defs['Mixing time [ms]:'] = {'widget': 'entry',
                                     'entity': 'spectrum_mixing_time'}

        defs['Ambiguity level (for multimers):'] = {'widget': 'choice',
                                     'entity': 'ambiguity_type',
                                    'choices' : ({'label' : 'Intra-molecular only',
                                                  'value' : 'intra'},
                                                 {'label' : 'Inter-molecular only',
                                                  'value' : 'inter'},
                                                 {'label' : 'Unknown',
                                                  'value' : 'all'})}
        
        d['experiment_data'] = CW(experiment.frame, defs, row = 2)
        
        generic.form(top = 0, left = 0, right = -1)
        shifts.form(top = generic, left = 0, right = -1)
        #peaks.form(top = shifts, bottom = -1)
        peaks.form(top = shifts, left = 0, right = -1)

        # BARDIAUX rMat
        experiment.form(top = peaks, left= 0,right = -1)

        return d

    def get_settings(self):

        return self.__settings

class CCPNPanel(DataPanel):

    def set_settings(self, settings):

        check_dict(settings)
        
        self.__settings = settings

        E = settings['ccpn'].getEntity('filename')
        E.set_callback(self.callback)
        
    def create(self, p):

        from aria.OrderedDict import OrderedDict

        f_project = Tix.LabelFrame(p, label='CCPN project folder')


        defs = OrderedDict()
        defs['Project'] = {'widget': 'ccpn_browser', 'entity': 'filename', 'columnspan':2}
        
        controls = self.create_widgets(f_project.frame, defs, row=0)
        

        defs = OrderedDict()

        choices = ({'label': 'Last iteration', 'value': 'last'},
                   {'label': 'All iterations', 'value': 'all'},
                   {'label': 'Disable', 'value': 'no'})

        defs['Export restraint lists:'] = {'widget': 'choice',
                                           'entity': 'export_noe_restraint_list',
                                           'type': STRING,
                                           'choices': choices}

        defs['Export assignments:'] = {'widget': 'choice',
                                      'entity': 'export_assignments',
                                      'type': STRING,
                                      'choices': YES_NO_DICT}

        defs['Export structures:'] = {'widget': 'choice',
                                      'entity': 'export_structures',
                                      'type': STRING,
                                      'choices': YES_NO_DICT}

        f_options =  Tix.LabelFrame(p, label='Options')

        controls_report = self.create_widgets(f_options.frame, defs, row=0)

        f_project.form(top=0, left=0, right=-1, fill=Tix.X)
        f_options.form(top=f_project, left=0, right=-1, fill=Tix.X)
        
        return {'ccpn': controls,
                'report': controls_report}

    def get_settings(self):
        return self.__settings

    def callback(self, E):
        self.project.ccpn_project = None

class DataKarplusPanel(DataPanel):

    def set_settings(self, settings):

        check_type(settings, 'KarplusData')
        self.__settings = {'data_karplus': settings}

    def add_ccpn_names(self):

        s = self.get_settings()['data_karplus']
        
        if not s['ccpn_id']:
            return
         
        from aria.importFromCcpn import getCcpnConstraintList, getKeysFromString

        ccpn_project = self.project.ccpn_project
        if not ccpn_project:            
            ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return
        
        if s['ccpn_id']:
            clist = getCcpnConstraintList(ccpn_project, getKeysFromString(s['ccpn_id']))
            name = "%d:%d:%s" % (clist.nmrConstraintStore.serial, cList.serial, cList.name)
            self.update_ccpn_name(self.ccpn_name, name)


    ## BARDIAUX
    def select_constraints(self):

        from ccpnGui import selectConstraintList

        ccpn_project = self.read_ccpn_model()

        kw = {'restraint_type' : 'JCouplingConstraintList',
              'msg' : 'Scalar Couplings'}
        
        key, msg, ccpn_name = selectConstraintList(ccpn_project, self.parent, **kw)
	self.message(msg)

        if key:
            settings = self.get_settings()['data_karplus']

            settings['ccpn_id'] = key

            self.update_ccpn_name(self.ccpn_name, ccpn_name)

            self.refresh(overwrite=0)

    def create(self, p):

        from aria.TypeChecking import INT
        from aria.OrderedDict import OrderedDict

        self.parent = p

        controls = []

        #file, format = file_and_format(p, ({'value': 'tbl', 'label': 'TBL'}, ))
        
        data_formats = ({'label': 'TBL', 'value': 'tbl'},
                        {'label': 'CCPN', 'value': 'ccpn'})

        file_types = [('TBL files', '*.tbl')]
        
        file, format, ccpn_id, ccpn_name = file_and_format(p, data_formats,
                                                file_types = file_types,
                                                with_ccpn_id=1, ccpn_select_func=self.select_constraints)

        self.ccpn_name = ccpn_name
        
        self.add_edit_button(p, 1, 'data_karplus')
        
        controls.append({'widget': file, 'entity': 'filename'})
        controls.append({'widget': format, 'entity': 'format'})
        controls.append({'widget': ccpn_id, 'entity': 'ccpn_id'})
        
        defs = OrderedDict()
        defs['Class:'] = {'widget': 'choice',
                          'entity': 'class',
                          'type': INT, 
                          'choices': ({'label': '1', 'value': '1'},
                                      {'label': '2', 'value': '2'},
                                      {'label': '3', 'value': '3'},
                                      {'label': '4', 'value': '4'},
                                      {'label': '5', 'value': '5'})}
        
        defs['Use restraints:'] = {'widget': 'choice', 'entity': 'enabled',
                                   'choices': YES_NO_DICT}
        
        controls += self.create_widgets(p, defs, row = 2)
        
        return {'data_karplus': controls}

    def get_settings(self):

        return self.__settings

class DataHBondPanel(DataPanel):

    def set_settings(self, settings):

        check_type(settings, 'HBondData')

        self.__settings = {'data_hbonds': settings}

    def add_ccpn_names(self):

        s = self.get_settings()['data_hbonds']
        
        if not s['ccpn_id']:
            return
         
        from aria.importFromCcpn import getCcpnConstraintList, getKeysFromString

        ccpn_project = self.project.ccpn_project
        if not ccpn_project:            
            ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return
        
        if s['ccpn_id']:
            cList = getCcpnConstraintList(ccpn_project, getKeysFromString(s['ccpn_id']))
            name = "%d:%d:%s" % (cList.nmrConstraintStore.serial, cList.serial, cList.name)
            self.update_ccpn_name(self.ccpn_name, name)


    ## BARDIAUX        
    def select_constraints(self):

        from ccpnGui import selectConstraintList

        ccpn_project = self.read_ccpn_model()

        kw = {'restraint_type' : 'HBondConstraintList',
              'msg' : 'Hydrogen Bonds'}
        
        key, msg, ccpn_name = selectConstraintList(ccpn_project, self.parent, **kw)
	self.message(msg)

        if key:
            settings = self.get_settings()['data_hbonds']

            settings['ccpn_id'] = key

            self.update_ccpn_name(self.ccpn_name, ccpn_name)

            self.refresh(overwrite=0)
            
    
    ## BARDIAUX        
    def create(self, p):
        from aria.OrderedDict import OrderedDict

        controls = []

        self.parent = p
        
        data_formats = ({'label': 'TBL', 'value': 'tbl'},
                        {'label': 'CCPN', 'value': 'ccpn'})

        file_types = [('TBL files', '*.tbl')]
        
        file, format, ccpn_id, ccpn_name = file_and_format(p, data_formats,
                                                file_types = file_types,
                                                with_ccpn_id=1, ccpn_select_func=self.select_constraints)

        self.ccpn_name  = ccpn_name
        
        self.add_edit_button(p, 1, 'data_hbonds')
        
        controls.append({'widget': file, 'entity': 'filename'})
        controls.append({'widget': format, 'entity': 'format'})
        controls.append({'widget': ccpn_id, 'entity': 'ccpn_id'})

        defs = OrderedDict()
        defs['Type:'] = {'widget': 'choice',
                         'entity': 'type',
                         'choices': ({'label': 'Standard',
                                      'value': 'standard'},
                                     {'label': 'CSI',
                                      'value': 'csi'})}
        
        defs['Use restraints:'] = {'widget': 'choice', 'entity': 'enabled',
                                   'choices': YES_NO_DICT}
        
        controls += self.create_widgets(p, defs, row = 2)
        
        return {'data_hbonds': controls}
    

    def get_settings(self):

        return self.__settings

class DataTemplatePanel(DataPanel):

    table_descr = \
'''
Template structures comprise the initial ensemble. If given, the initial structure ensemble is used by ARIA to calibrate all spectra for the first itertion and to create a seed NOE-assignment by analyzing the restraint-list with respect to the template structures. During setup, all template structures are copied from their source locations into the local data directory PROJECT_PATH/RUNxxx/data/templates.

To modify an entry field, <double-click> the cell. To keep the changes, press <return> or <click> another cell; to discard the changed, press <escape>.
'''
    def __init__(self, node = None):
        DataPanel.__init__(self, node, enable_delete_button = 0)

    def add_ccpn_model(self):

        from aria.importFromCcpn import getModels, getObjectKeyString, getKeysFromString, getCcpnChains
        from ccpnGui import gui_select_model
        from ccpnmr.format.general.Util import createSelection
        import aria.DataContainer as DC
        
        ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return

        mol_data = self.project.getData(DC.DATA_SEQUENCE)[0]
        
        molsystem_id = mol_data['ccpn_id']

        if molsystem_id == '':

            from widgets import OKDialog

            msg = 'Name of molecular system has not yet been specified (see entry field "CCPN name" in node "Sequence").'

            OKDialog(title='Error', text=msg)

            return
        
        chains = getCcpnChains(ccpn_project, getKeysFromString(molsystem_id))
        
        models = getModels(ccpn_project, chains)

        if not models:

            from widgets import OKDialog

            msg = 'No CCPN models present for\nmolecular system "%s".' % molsystem_id

            OKDialog(title='Error', text=msg)

            return        

        model_labels, model_dict = createSelection(models)
        
        model = gui_select_model(model_labels,
                                 model_dict,
                                 gui=self.parent)

        if model is None:
            return

        key = getObjectKeyString(model)

        row = self._insert_ccpn_model(key)

        

    def add_ccpn_ensemble(self):

        from aria.importFromCcpn import getStructureEnsembles, getObjectKeyString, getKeysFromString, getCcpnChains
        from ccpnGui import gui_select_ensemble
        from ccpnmr.format.general.Util import createSelection
        import aria.DataContainer as DC
        
        ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return

        mol_data = self.project.getData(DC.DATA_SEQUENCE)[0]
        
        molsystem_id = mol_data['ccpn_id']

        if molsystem_id == '':

            from widgets import OKDialog

            msg = 'Name of molecular system has not yet been specified (see entry field "CCPN name" in node "Sequence").'

            OKDialog(title='Error', text=msg)

            return
        
        chains = getCcpnChains(ccpn_project, getKeysFromString(molsystem_id))
        
        ensembles = getStructureEnsembles(ccpn_project, chains)

        if not ensembles:

            from widgets import OKDialog

            msg = 'No CCPN Structure ensemble present for\nmolecular system "%s".' % molsystem_id

            OKDialog(title='Error', text=msg)

            return        

        ensemble_labels, ensemble_dict = createSelection(ensembles)
        
        ensemble = gui_select_ensemble(ensemble_labels,
                                       ensemble_dict,
                                       gui=self.parent)

        if ensemble is None:
            return
        
        if not ensemble.models:
            return

        for model in ensemble.sortedModels():
            key = getObjectKeyString(model)
            self._insert_ccpn_model(key)

    
    def _insert_ccpn_model(self, key):

        n_rows = len(self.table.rows)
        
        # find a free row
        free = self.first_free_row
        if free >= n_rows:
            self.more_rows(n=1)
            
        self.table.rows[free][3].set(key)
        self.table.rows[free][2].set('ccpn')
        
        self.first_free_row += 1

        return free
        
    def save(self):
        import aria.DataContainer as DC

        templates = []

        d = {'enabled': 0,
             'filename': 1,
             'format': 2,
             'ccpn_id': 3}

        for i, row in enumerate(self.table.rows):

            valid = {'enabled': 0,
                     'filename': 0,
                     'format': 0,
                     'ccpn_id': 0}

            t = DC.TemplateData()
            ok = 0
            #t['ccpn_id'] = ''

            for entity_name, column in d.items():

                cell = row[column]

                try:
                    val = cell.get()
                except:
                    self.not_saved(cell.editor)
                    val = None
                    return self.set_modified(1)
                
                if val is not None:
                    cell.save_now()
                    t[entity_name] = row[column].get()
                    ok += 1
                    valid[entity_name] = 1
                else:
                    break
                
            if ok == len(d):
                templates.append(t)
                self.first_free_row+=1
            else:
                if ok == (len(d)-1) and (valid['ccpn_id'] + valid['filename']) == 1:
                    if valid['ccpn_id']:
                        cell = row[d['filename']]
                        cell.save_now()
                        t['filename'] = cell.get()
                    templates.append(t)
                    self.first_free_row+=1


        project = self.get_settings()['project']

        current = project.getData(DC.DATA_TEMPLATE_STRUCTURE)
        [project.delData(t) for t in current]
        
        [project.addData(t) for t in templates]

        self.set_modified(0)
        
        return 1

    def validate(self):
        """
        Checks whether all cells are valid.
        """
        for row in self.table.rows:
            for r in row:
                is_valid = r.isValid()
                if is_valid == 0:
                    return self.set_modified()
                
        return self.set_modified(0)

    def not_saved(self, widget):
        from tkMessageBox import showerror

        if widget is not None:
            try:
                widget.save()
            except Exception, msg:
                showerror('Value error',msg)
        
        self.set_modified()
    
    def more_rows(self, n = 1):

        l = len(self.table.rows)
        indices = range(l, l + n)
        
        self.table.append_row(n)
        self.init_cells(indices)

    def del_row(self):
        current_row = self.table.get_current_row_index()

        if current_row is None:
            return
        self.table.remove_row(current_row)
        self.first_free_row-=1

    def init_cells(self, indices):

        from aria.DataContainer import TemplateData
        from aria.Settings import YesNoChoice, Path

        table = self.table

        for i in indices:

            t = TemplateData()

            ## 0th column: enabled

            e = t.getEntity('enabled')
            table.set_widget(i, 0, DOptionMenu, e)
            editor = table.rows[i][0].editor
            editor.add_command(YES, label = 'Yes')
            editor.add_command(NO, label = 'No')

            ## By default, hosts are enabled
            table.rows[i][0].set(YES)

            ## 1st column: filename

            e = t.getEntity('filename')
            table.set_widget(i, 1, DFileEntry, e)

            ## 2nd column: format

            d = {'iupac': 'IUPAC',
                 'cns': 'CNS',
                 'dyana': 'DYANA',
                 'ccpn' : 'CCPN'}

            e = t.getEntity('format')
            table.set_widget(i, 2, DOptionMenu, e)
            editor = table.rows[i][2].editor

            for key, label in d.items():
                editor.add_command(key, label = label)

            table.rows[i][2].set('iupac')

            ## BARDIAUX
            ## 3rd column: ccpn_id
            e = t.getEntity('ccpn_id')
            table.set_widget(i, 3, DEntry, e)

            for j in range(4):
                table.rows[i][j].setCallback(self.not_saved)

            ## BARDIAUX
            ## 4th column not editable
##             old_label = self.table.rows[i][4].label
##             info = old_label.grid_info()
##             l = Tix.Label(self.table.window)
##             l.config(width = old_label['width'], relief = Tix.FLAT,
##                      bd = 1, bg = 'white')
##             old_label.forget()
##             self.table.rows[i][4].label = l


                
    def create(self, p, title = 'Initial structure ensemble'):

        import aria.tools as tools
        import aria.DataContainer as DC

        self.parent = p
        
        frame = Tix.LabelFrame(p, label = title)

        header = ('Enabled', 'Filename', 'Format', 'CCPN Id')
        
        c_widths = (10, 50, 10, 10)

        ## Initially, we have 10 rows
        
        self.table = DTable(frame.frame, (0, 4), header = header,
                            column_widths = c_widths, relief = Tix.FLAT)

        table_size = 10
        self.more_rows(table_size)

        ## fill table with host-list data

        project = self.get_settings()['project']
        templates = project.getData(DC.DATA_TEMPLATE_STRUCTURE)

        n_rows = len(templates)
        
        self.more_rows(max(n_rows - table_size, 0))

        kw = ('enabled', 'filename', 'format', 'ccpn_id')

        for i in range(n_rows):
            for j in range(len(kw)):
                self.table[i,j] = templates[i][kw[j]]

        self.first_free_row = n_rows

        self.table.grid(row = 0, column = 0, columnspan = 2)
##        self.table.grid_propagate(0)

        box = Tix.ButtonBox(frame.frame, orientation=Tix.HORIZONTAL,
                            relief = Tix.FLAT)

        descr = tools.wrap_string(self.table_descr,
                                  AriaBaseClass.description_length)

        self.get_owner().bind_tooltips(box, descr)
        
        box.add('add', text = 'Add row', width = 10,
                command = self.more_rows, borderwidth = 1)
        
        box.add('del', text = 'Delete row', width = 10,
                command = self.del_row, borderwidth = 1)

        box.add('add_ccpn1', text = 'Add CCPN Model',
                command = self.add_ccpn_model, borderwidth = 1)

        box.add('add_ccpn2', text = 'Add CCPN Structure Ensemble',
                command = self.add_ccpn_ensemble, borderwidth = 1)        

        box.grid(row = 1, column = 0, columnspan = 2,
                 sticky = Tix.EW, pady = 5)

        frame.form(left = 0, top = 0)

        return None

    def set_settings(self, project):
        
        self.__settings = {'project': project}
        
    def get_settings(self):
        return self.__settings
    
class DataDihedralPanel(DataPanel):

    def set_settings(self, settings):

        check_type(settings, 'DihedralData')

        self.__settings = {'data_dihedrals': settings}

    def add_ccpn_names(self):

        s = self.get_settings()['data_dihedrals']
        
        if not s['ccpn_id']:
            return
         
        from aria.importFromCcpn import getCcpnConstraintList, getKeysFromString

        ccpn_project = self.project.ccpn_project
        if not ccpn_project:            
            ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return
        
        if s['ccpn_id']:
            cList = getCcpnConstraintList(ccpn_project, getKeysFromString(s['ccpn_id']))
            name = "%d:%d:%s" % (cList.nmrConstraintStore.serial, cList.serial, cList.name)
            self.update_ccpn_name(self.ccpn_name, name)
            

    ## BARDIAUX        
    def select_constraints(self):

        from ccpnGui import selectConstraintList

        ccpn_project = self.read_ccpn_model()

        kw = {'restraint_type' : 'DihedralConstraintList',
              'msg' : 'Dihedral Angle'}
        
        key, msg, ccpn_name = selectConstraintList(ccpn_project, self.parent, **kw)
	self.message(msg)

        if key:
            settings = self.get_settings()['data_dihedrals']

            settings['ccpn_id'] = key

            self.update_ccpn_name(self.ccpn_name, ccpn_name)

            self.refresh(overwrite=0)
            

    def create(self, p):
        from aria.OrderedDict import OrderedDict

        self.parent = p

        controls = []

        #file, format = file_and_format(p, ({'value': 'tbl', 'label': 'TBL'}, ))
        
        data_formats = ({'label': 'TBL', 'value': 'tbl'},
                        {'label': 'CCPN', 'value': 'ccpn'})

        file_types = [('TBL files', '*.tbl')]
        
        file, format, ccpn_id, ccpn_name = file_and_format(p, data_formats,
                                                file_types = file_types,
                                                with_ccpn_id=1, ccpn_select_func=self.select_constraints)

        self.ccpn_name = ccpn_name
        
        self.add_edit_button(p, 1, 'data_dihedrals')
        
        controls.append({'widget': file, 'entity': 'filename'})
        controls.append({'widget': format, 'entity': 'format'})
        controls.append({'widget': ccpn_id, 'entity': 'ccpn_id'})
        
        defs = OrderedDict()
        defs['Type:'] = {'widget': 'choice',
                         'entity': 'type',
                         'choices': ({'label': 'Standard',
                                      'value': 'standard'},
                                     {'label': 'CSI',
                                      'value': 'csi'},
                                     {'label': 'TALOS',
                                      'value': 'talos'})}
        
        defs['Use restraints:'] = {'widget': 'choice', 'entity': 'enabled',
                                   'choices': YES_NO_DICT}
        
        controls += self.create_widgets(p, defs, row = 2)
        
        return {'data_dihedrals': controls}


    def get_settings(self):

        return self.__settings

class DataOtherPanel(DataPanel):

    def set_settings(self, settings):

        check_type(settings, 'OtherData')

        self.__settings = {'data_other_data': settings}
        
##     ## BARDIAUX        
##     def select_constraints(self):

##         from ccpnGui import selectConstraintList

##         ccpn_project = self.read_ccpn_model()

##         kw = {'restraint_type' : 'DihedralConstraintList',
##               'msg' : 'Dihedral Angle'}
        
##         key, msg = selectConstraintList(ccpn_project, self.parent, **kw)
## 	self.message(msg)

##         if key:
##             settings = self.get_settings()['data_dihedrals']

##             settings['ccpn_id'] = key

##             self.refresh(overwrite=0)
            

    def create(self, p):
        from aria.OrderedDict import OrderedDict

        self.parent = p

        controls = []

        #file, format = file_and_format(p, ({'value': 'tbl', 'label': 'TBL'}, ))
        
        data_formats = ({'label': 'TBL', 'value': 'tbl'},)
##                         {'label': 'CCPN', 'value': 'ccpn'})

        file_types = [('TBL files', '*.tbl')]
        
        file, format, ccpn_id, ccpn_name = file_and_format(p, data_formats,
                                                file_types = file_types,
                                                with_ccpn_id=1, ccpn_select_func=None)

        
        self.add_edit_button(p, 1, 'data_other_data')
        
        controls.append({'widget': file, 'entity': 'filename'})
        controls.append({'widget': format, 'entity': 'format'})
        controls.append({'widget': ccpn_id, 'entity': 'ccpn_id'})
        
        defs = OrderedDict()
        defs['Type:'] = {'widget': 'choice',
                         'entity': 'type',
                         'choices': ({'label': 'DNA/RNA planarity',
                                      'value': 'planarity'},)}

        
        defs['Use restraints:'] = {'widget': 'choice', 'entity': 'enabled',
                                   'choices': YES_NO_DICT}
        
        controls += self.create_widgets(p, defs, row = 2)
        
        return {'data_other_data': controls}


    def get_settings(self):

        return self.__settings

class DataRDCPanel(DataPanel):

    def set_settings(self, settings):

        check_type(settings, 'RDCData')

        self.__settings = {'data_rdc': settings}

    def add_ccpn_names(self):

        s = self.get_settings()['data_rdc']
        
        if not s['ccpn_id']:
            return
         
        from aria.importFromCcpn import getCcpnConstraintList, getKeysFromString

        ccpn_project = self.project.ccpn_project
        if not ccpn_project:            
            ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return
        
        if s['ccpn_id']:
            cList = getCcpnConstraintList(ccpn_project, getKeysFromString(s['ccpn_id']))
            name = "%d:%d:%s" % (cList.nmrConstraintStore.serial, cList.serial, cList.name)
            self.update_ccpn_name(self.ccpn_name, name)


    ## BARDIAUX        
    def select_constraints(self):

        from ccpnGui import selectConstraintList

        ccpn_project = self.read_ccpn_model()

        kw = {'restraint_type' : 'RdcConstraintList',
              'msg' : 'RDC'}
        
        key, msg, ccpn_name = selectConstraintList(ccpn_project, self.parent, **kw)
	self.message(msg)

        if key:
            settings = self.get_settings()['data_rdc']

            settings['ccpn_id'] = key

            self.update_ccpn_name(self.ccpn_name, ccpn_name)
            
            self.refresh(overwrite=0)
            

    def create(self, p):

        from aria.TypeChecking import INT
        from aria.OrderedDict import OrderedDict

        self.parent = p

        controls = []

        #file, format = file_and_format(p, ({'value': 'tbl', 'label': 'TBL'}, ))
        
        data_formats = ({'label': 'TBL', 'value': 'tbl'},
                        {'label': 'CCPN', 'value': 'ccpn'})

        file_types = [('TBL files', '*.tbl')]
        
        file, format, ccpn_id, ccpn_name = file_and_format(p, data_formats,
                                                file_types = file_types,
                                                with_ccpn_id=1, ccpn_select_func=self.select_constraints)

        self.ccpn_name = ccpn_name
        
        self.add_edit_button(p, 1, 'data_rdc')
        
        controls.append({'widget': file, 'entity': 'filename'})
        controls.append({'widget': format, 'entity': 'format'})
        controls.append({'widget': ccpn_id, 'entity': 'ccpn_id'})
        
        defs = OrderedDict()
        defs['Parameter class:'] = {'widget': 'choice',
                                    'entity': 'class',
                                    'type': INT, 
                                    'choices': ({'label': '1', 'value': '1'},
                                                {'label': '2', 'value': '2'},
                                                {'label': '3', 'value': '3'},
                                                {'label': '4', 'value': '4'},
                                                {'label': '5', 'value': '5'})}
        
        defs['Use restraints:'] = {'widget': 'choice', 'entity': 'enabled',
                                   'choices': YES_NO_DICT}
        
        controls += self.create_widgets(p, defs, row = 2)
        
        return {'data_rdc': controls}


    def get_settings(self):

        return self.__settings

    
## BARDIAUX 2.2
class SymmetryPanel(DataPanel):

    def __init__(self, node = None):
        DataPanel.__init__(self, node, enable_delete_button = 0)

    def set_settings(self, settings):

        check_type(settings, 'Symmetry')

        self.__settings = {'symmetry' : settings}
        
    def create(self, p):

        controls = []

        LF = Tix.LabelFrame

        CW = self.create_widgets
        
        generic = LF(p, label = 'Symmetry')

        from aria.OrderedDict import OrderedDict
        
        defs = OrderedDict()

        defs['Symmetry enabled:'] = {'widget': 'choice', 'entity': 'enabled',
                                     'choices': YES_NO_DICT}

        defs['Symmetry method:'] = {'widget': 'choice',
                                    'entity': 'method',
                                    'choices': ({'label': 'Standard', 'value': 'standard'},)}
        
        defs['Symmetry type:'] = {'widget': 'choice',
                                  'entity': 'symmetry_type',
                                  'choices': ({'label': 'None', 'value': 'None'},
                                              {'label': 'C2', 'value': 'C2'},
                                              {'label': 'C3', 'value': 'C3'},
                                              {'label': 'C5', 'value': 'C5'},
                                              {'label': 'D2', 'value': 'D2'})}
        
        defs['Number of monomers:'] = {'widget': 'entry',
                                         'entity': 'n_monomers'}
        
        d = []
        d += CW(generic.frame, defs)
        


        rest = Tix.LabelFrame(p, label = 'Symmetry related restraints')

        defs = OrderedDict()

        defs['NCS restraints enabled:'] = {'widget': 'choice', 'entity': 'ncs_enabled',
                                           'choices': YES_NO_DICT}

        defs['Packing restraints enabled:'] = {'widget': 'choice', 'entity': 'packing_enabled',
                                               'choices': YES_NO_DICT}
        

        d += CW(rest.frame, defs)

        generic.form(top = 0, left = 0, right = -1)
        rest.form(top = generic, left = 0, right = -1)

 

        return {'symmetry' : d}
    
    def get_settings(self):
        return self.__settings



from widgets import Dialog, DefaultButtonBox, BUTTON_WIDTH

class ChainCodeSelector(Dialog):

    def __init__(self, parent, mol_system_name, codes):

        self.codes = codes
        self.var = StringVar()
        self.mol_system_name = mol_system_name

        Dialog.__init__(self, parent, 'Select chain code', 1)

    def _iterall(self, iterables):
        if iterables:
            for head in iterables[0]:
                for remainder in self._iterall(iterables[1:]):
                    if head not in remainder:
                        yield [head] + remainder
        else:
            yield []
        
    def setup(self):

        frame = Tix.Frame(self,   bd=1, relief=Tix.RAISED) # padx=20, pady=10,

        om = Tix.OptionMenu(frame, label = 'Please select chain code for molecular system "%s":  ' % self.mol_system_name, variable=self.var)


        pairs = []    
        for i in range(1, len(self.codes)+1):
            for thing in self._iterall([self.codes] * i):
                thing.sort()
                if thing not in pairs:
                    pairs.append(thing)
                    c = ",".join(thing)
                    om.add_command(c, label='"%s"' % str(c))
                    

        om.pack(side=Tix.TOP)

        box = DefaultButtonBox(frame, orientation=Tix.HORIZONTAL)
        
        box.add('select', text = 'Select', width = BUTTON_WIDTH,
                command = lambda s = self: s.select(0))
        
        box.add('cancel', text = 'Cancel', width = BUTTON_WIDTH,
                command = lambda s = self: s.select(1))
        
        box.pack(side=Tix.BOTTOM,  anchor=Tix.CENTER)

        frame.pack()
        
    def select(self, cancel=0):

        if not cancel:
            self.result = self.var.get()
        else:
            self.result = None
                
        self.destroy()

class DataSequencePanel(DataPanel):

    def select_from_ccpn(self):

        from ccpnGui import gui_select_molSystem
        from aria.importFromCcpn import getObjectKeyString, getObjectsKeyString
        #from memops.general.Io import loadXmlProjectFile

        ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return
        
        mol_system = gui_select_molSystem(ccpn_project, gui=self.parent)

        if mol_system is None:
            return

        chain_codes = [chain.code for chain in mol_system.chains]

        if len(chain_codes) > 1:
            dialog = ChainCodeSelector(self.parent, mol_system.code, chain_codes)

            chain_code = dialog.result

            if chain_code is None:
                return

        else:
            chain_code = chain_codes[0]

        # BARDIAUX 
        chain_code = chain_code.split(",")
        
        if len(chain_code) > 1:
            chains = map(lambda c: mol_system.findFirstChain(code=c), chain_code)
            ccpn_id = getObjectsKeyString(chains)
            
        else:
            chain_code = chain_code[0]
            chain = mol_system.findFirstChain(code=chain_code)
            ccpn_id = getObjectKeyString(chain)
            
        #chain = mol_system.findFirstChain(code=chain_code)

        settings = self.get_settings()['data_sequence']

        #ccpn_id = getObjectKeyString(chain)
        
        settings['ccpn_id'] = ccpn_id

        self.find_widget('ccpn_id').load()
        

    def set_settings(self, settings):

        check_type(settings, 'SequenceData')

        self.__settings = {'data_sequence': settings}

    def __change_state(self, widget, *args):

        name = widget.get_entity().getName()
        base = name[:name.find('_')]
        
        val = widget.widget_get()
        if val == 'user_defined':
            state = Tix.NORMAL
        else:
            state = Tix.DISABLED

        entry = self.entries[base]
        
        for key in ('widget', 'label'):
            entry[key].configure(state = state)
            
    def create_file_selector(self, parent, label, items, entity_name,
                             entity_filename):
        
        from aria.OrderedDict import OrderedDict

        lf = Tix.LabelFrame(parent, label = label)
        
        defs = OrderedDict()

        choices = []

        ## values are labels, keys are values
        ## we want to have it the other way round
        ## reverse dict

        d = {}

        for key, value in items.items():
            d[value] = key

        labels = d.keys()
        labels.sort()
        
        for label in labels:
            value = d[label]
            choices.append({'label': label, 'value': value})

        defs['Definition Name:'] = {'widget': 'choice',
                                    'entity': entity_name,
                                    'choices': tuple(choices)}
        
        defs['Filename:'] = {'widget': 'file',
                             'entity': entity_filename}

        controls = self.create_widgets(lf.frame, defs)

        return controls, lf

    def create(self, p):

        from aria.DataContainer import SequenceData

        self.parent = p

        data_formats = ({'label': 'XML', 'value': 'xml'},
                        {'label': 'CCPN', 'value': 'ccpn'})
        
        controls = []

        LF = Tix.LabelFrame
        FUNC = self.create_file_selector
        
        f_seq = LF(p, label = 'Sequence file')
        
        file, format, ccpn_id, ccpn_name = file_and_format(f_seq.frame, data_formats,
                                                with_ccpn_id=1,
                                                ccpn_select_func=self.select_from_ccpn)

        self.add_edit_button(f_seq.frame, 1, 'data_sequence')

        controls.append({'widget': file, 'entity': 'filename'})
        controls.append({'widget': format, 'entity': 'format'})
        controls.append({'widget': ccpn_id, 'entity': 'ccpn_id'})

        c_link, f_link = FUNC(p, 'Linkage Definition',
                              SequenceData.linkage_names,
                              'linkage_name',
                              'linkage_filename')

        c_top, f_top = FUNC(p, 'Topology Definition',
                            SequenceData.topology_names,
                            'topology_name',
                            'topology_filename')

        c_param, f_param = FUNC(p, 'Parameter Definition',
                                SequenceData.parameter_names,
                                'parameter_name',
                                'parameter_filename')

        controls += c_link + c_top + c_param

        ## Set callbacks

        base_names = ('linkage', 'topology', 'parameter')
        
        ## dict that stores entry fields
        
        self.entries = {}

        for base in base_names:

            name = base + '_name'
            filename = base + '_filename'

            menu = [d['widget'] for d in controls \
                    if d['entity'] == name][0]

            entry = [d for d in controls \
                    if d['entity'] == filename][0]

            self.entries[base] = entry

            menu.set_trace_callback('rw', self.__change_state)

        f_seq.form(top = 0, left = 0, right = -1, fill = Tix.X)
        f_link.form(top = f_seq, left = 0, right = -1, fill = Tix.X)
        f_top.form(top = f_link, left = 0, right = -1, fill = Tix.X)
        f_param.form(top = f_top, left = 0, right = -1, fill = Tix.X)

        return {'data_sequence': controls}

    def get_settings(self):
        return self.__settings
    
class SSBridgePanel(PanelEx):

    table_descr = ''

    ## do not change order of tuple below!
    entity_names = ('residue1', 'segid1', 'residue2', 'segid2')
    
    def save(self):

        import aria.DataContainer as DC

        bridges = []

        for row in self.table.rows:

            bridge = DC.SSBridge()

            ok = 0

            for column in range(len(self.entity_names)):
                entity_name = self.entity_names[column]

                cell = row[column]

                try:
                    val = cell.get()
                except:
                    self.not_saved(cell.editor)
                    val = None
                    return self.set_modified(1)

                if val is not None:
                    
                    cell.save_now()
                    val = row[column].get()

                    ## cast values

                    entity = cell.editor.get_entity()
                    val = entity.cast(val)
                    bridge[entity_name] = val
                    ok += 1
                else:
                    break

            if ok == len(self.entity_names):
                bridges.append(bridge)

        ## remove old bridges
                
        current = self.project.getData(DC.DATA_SSBRIDGE)
        [self.project.delData(d) for d in current]

        ## add new bridges

        [self.project.addData(b) for b in bridges]

        self.set_modified(0)

        return 1
    
    def not_saved(self, widget):
        from tkMessageBox import showerror

        if widget is not None:
            try:
                widget.save()
            except Exception, msg:
                showerror('Value error',msg)
        
        self.set_modified()
    
    def validate(self):
        """
        Checks whether all cells are valid.
        """
        for row in self.table.rows:
            for r in row:
                is_valid = r.isValid()
                if is_valid == 0:
                    return self.set_modified()

        ## all cells are valid
                
        return self.set_modified(0)

    def more_rows(self, n = 1):

        l = len(self.table.rows)
        indices = range(l, l + n)
        
        self.table.append_row(n)
        self.init_cells(indices)

    def del_row(self):
        current_row = self.table.get_current_row_index()
        if current_row is None:
            return
        self.table.remove_row(current_row)

    def init_cells(self, indices):

        from aria.DataContainer import SSBridge


        for i in indices:

            bridge = SSBridge()
            bridge.reset()

            for j in range(len(self.entity_names)):
                name = self.entity_names[j]
                e = bridge.getEntity(name)
                self.table.set_widget(i, j, DEntry, e)

    def create(self, p):

        import aria.tools as tools
        import aria.DataContainer as DC

        frame = Tix.LabelFrame(p, label = 'SS-brigde list')

        header = ('Residue 1', 'Segid 1', 'Residue 2', 'Segid 2')
        
        self.table = DTable(frame.frame, (0, 4), header = header,
                            relief = Tix.FLAT)
        
        table_size = 10
        self.more_rows(table_size)
        
        ## fill table 

        bridges = self.project.getData(DC.DATA_SSBRIDGE)
        n_rows = len(bridges)
        self.more_rows(max(n_rows - table_size, 0))

        for i in range(n_rows):
            for j in range(len(self.entity_names)):
                name = self.entity_names[j]
                self.table[i,j] = bridges[i][name]

        self.table.grid(row = 0, column = 0, columnspan = 2)
##        self.table.grid_propagate(0)

        box = Tix.ButtonBox(frame.frame, orientation=Tix.HORIZONTAL,
                            relief = Tix.FLAT)

        descr = tools.wrap_string(self.table_descr,
                                  AriaBaseClass.description_length)
        
        self.get_owner().bind_tooltips(box, descr)
        
        box.add('add', text = 'Add row', width = 10,
                command = self.more_rows, borderwidth = 1)
        
        box.add('del', text = 'Delete row', width = 10,
                command = self.del_row, borderwidth = 1)

        box.grid(row = 1, column = 0, columnspan = 2,
                 sticky = Tix.EW, pady = 5)

        frame.form(left = 0, top = 0)

        return None

    def set_settings(self, project):
        self.project = project
        self.__settings = {}
        
    def get_settings(self):
        return self.__settings
    
class SSBridgePanel2(PanelEx):

    def save(self):

        from aria.DataContainer import DATA_SSBRIDGE

        valid_data = []

        for row in self.table.rows:
            try:
                [row[i].get() for i in range(4)]
                valid_data.append(row[0].editor.get_entity()._parent)
            except:
                pass

        ## remove old data from project

        p = self.project

        [p.delData(d) for d in p.getData(DATA_SSBRIDGE)]

        ## add new data

        [p.addData(d) for d in valid_data]

        return 1

    def _new_rows(self, n):
        from aria.DataContainer import SSBridge

        offset = len(self.table.rows)
        s = {}

        for i in range(n):
            s[offset + i] = SSBridge()
            s[offset + i].reset()
            
        return s

    def more_rows(self, n = 4):
        d = self._new_rows(n)
        self.table.append_row(n)
        self.set_widgets(d)

    def del_row(self):
        current_row = self.table.get_current_row_index()
        if current_row is None:
            return
        self.table.remove_row(current_row)

    def set_widgets(self, s):
        for i in s.keys():

            res1 = s[i].getEntity('residue1')
            segid1 = s[i].getEntity('segid1')
            res2 = s[i].getEntity('residue2')
            segid2 = s[i].getEntity('segid2')
            

            ## TODO: hacked

            res1._parent = s[i]
            res2._parent = s[i]

            self.table.set_widget(i, 0, DEntry, res1)
            self.table.set_widget(i, 1, DEntry, segid1)
            self.table.set_widget(i, 2, DEntry, res2)
            self.table.set_widget(i, 3, DEntry, segid2)
    
    def create(self, p):

        from aria.DataContainer import DATA_SSBRIDGE

        frame = Tix.LabelFrame(p, label = 'SS-brigde list')

        header = ('Residue 1', 'Segid 1', 'Residue 2', 'Segid 2')
        
        self.table = DTable(frame.frame, (0, 4), header = header,
                            relief = Tix.FLAT)
        
        data = self.project.getData(DATA_SSBRIDGE)
        n_rows = max(10, len(data))

        self.more_rows(n_rows)

        ## fill table with data stored in current project

        for i in range(len(data)):
            self.table[i,0] = data[i]['residue1']
            self.table[i,1] = data[i]['segid1']
            self.table[i,2] = data[i]['residue2']
            self.table[i,3] = data[i]['segid2']

        self.table.grid(row = 0, column = 0, columnspan = 2)
##        self.table.propagate(0)

        box = Tix.ButtonBox(frame.frame, orientation=Tix.HORIZONTAL,
                            relief = Tix.FLAT)
        box.add('add', text = 'Add rows', width = 10,
                command = self.more_rows)
        
        box.add('del', text = 'Delete row', width = 10,
                command = self.del_row)

        box.grid(row = 1, column = 0, columnspan = 2,
                 sticky = Tix.EW, pady = 5)

        frame.form(left = 0, top = 0)

        return None

    def set_settings(self, project):
        from aria.DataContainer import SSBridge

        self.project = project

        self.__settings = {}
        
    def get_settings(self):
        return self.__settings

class HisPatchPanel(PanelEx):

    table_descr = ''
    entity_names = ('residue', 'segid', 'proton')

    def save(self):

        import aria.DataContainer as DC

        patches = []

        for row in self.table.rows:

            patch = DC.HisPatch()
            ok = 0

            for column in range(len(self.entity_names)):
                entity_name = self.entity_names[column]

                cell = row[column]

                try:
                    val = cell.get()
                except:
                    self.not_saved(cell.editor)
                    val = None
                    return self.set_modified(1)

                if val is not None:

                    cell.save_now()
                    val = row[column].get()

                    ## cast values

                    entity = cell.editor.get_entity()
                    val = entity.cast(val)
                    patch[entity_name] = val
                    ok += 1
                else:
                    break

            if ok == len(self.entity_names):
                patches.append(patch)

        ## remove old data from project

        current = self.project.getData(DC.DATA_HISPATCH)
        [self.project.delData(d) for d in current]

        ## add new data

        [self.project.addData(d) for d in patches]

        self.set_modified(0)
        
        return 1

    def validate(self):
        """
        Checks whether all cells are valid.
        """
        for row in self.table.rows:
            for r in row:
                is_valid = r.isValid()
                if is_valid == 0:
                    return self.set_modified()

        ## all cells are valid
                
        return self.set_modified(0)

    def more_rows(self, n = 4):

        l = len(self.table.rows)
        indices = range(l, l + n)
        
        self.table.append_row(n)
        self.init_cells(indices)

    def del_row(self):
        current_row = self.table.get_current_row_index()
        if current_row is None:
            return
        self.table.remove_row(current_row)

    def init_cells(self, indices):

        from aria.DataContainer import HisPatch
        
        for i in indices:

            patch = HisPatch()
            patch.reset()
                        
            ## 0th column: residue

            e = patch.getEntity('residue')
            self.table.set_widget(i, 0, DEntry, e)

            ## 1st column: segid

            e = patch.getEntity('segid')
            self.table.set_widget(i, 1, DEntry, e)

            ## 2nd column: proton

            e = patch.getEntity('proton')
            self.table.set_widget(i, 2, DOptionMenu, e)
            editor = self.table.rows[i][2].editor
            editor.add_command('HISD', label = 'HISD')
            editor.add_command('HISE', label = 'HISE')            

    def create(self, p):

        import aria.tools as tools
        import aria.DataContainer as DC

        frame = Tix.LabelFrame(p, label = 'HIS-patches')

        header = ('Residue', 'Segid', 'Proton')
        
        self.table = DTable(frame.frame, (0, 3), header = header,
                            relief = Tix.FLAT)
        
        patches = self.project.getData(DC.DATA_HISPATCH)
        n_rows = max(10, len(patches))

        self.more_rows(n_rows)

        ## fill table with data stored in current project

        for i in range(len(patches)):
            for j in range(len(self.entity_names)):
                name = self.entity_names[j]
                self.table[i,j] = patches[i][name]

        self.table.grid(row = 0, column = 0, columnspan = 2)
##        self.table.propagate(0)

        box = Tix.ButtonBox(frame.frame, orientation=Tix.HORIZONTAL,
                            relief = Tix.FLAT)
        
        descr = tools.wrap_string(self.table_descr,
                                  AriaBaseClass.description_length)
        
        self.get_owner().bind_tooltips(box, descr)
        
        box.add('add', text = 'Add rows', width = 10,
                command = self.more_rows)
        
        box.add('del', text = 'Delete row', width = 10,
                command = self.del_row)

        box.grid(row = 1, column = 0, columnspan = 2,
                 sticky = Tix.EW, pady = 5)

        frame.form(left = 0, top = 0)

        return None

    def set_settings(self, project):
        self.project = project

        self.__settings = {}
        
    def get_settings(self):
        return self.__settings
    
## BARDIAUX 2.2    
class CisProPatchPanel(HisPatchPanel):

    table_descr = ''
    entity_names = ('residue', 'segid')

    def save(self):

        import aria.DataContainer as DC

        patches = []

        for row in self.table.rows:

            patch = DC.CisProPatch()
            ok = 0

            for column in range(len(self.entity_names)):
                entity_name = self.entity_names[column]

                cell = row[column]

                try:
                    val = cell.get()
                except:
                    self.not_saved(cell.editor)
                    val = None
                    return self.set_modified(1)

                if val is not None:

                    cell.save_now()
                    val = row[column].get()

                    ## cast values

                    entity = cell.editor.get_entity()
                    val = entity.cast(val)
                    patch[entity_name] = val
                    ok += 1
                else:
                    break

            if ok == len(self.entity_names):
                patches.append(patch)

        ## remove old data from project

        current = self.project.getData(DC.DATA_CISPROPATCH)
        [self.project.delData(d) for d in current]

        ## add new data

        [self.project.addData(d) for d in patches]

        self.set_modified(0)
        
        return 1

    def init_cells(self, indices):

        from aria.DataContainer import CisProPatch
        
        for i in indices:

            patch = CisProPatch()
            patch.reset()
                        
            ## 0th column: residue

            e = patch.getEntity('residue')
            self.table.set_widget(i, 0, DEntry, e)

            ## 1st column: segid

            e = patch.getEntity('segid')
            self.table.set_widget(i, 1, DEntry, e)          

    def create(self, p):

        import aria.tools as tools
        import aria.DataContainer as DC

        frame = Tix.LabelFrame(p, label = 'Cis-Proline patches')

        header = ('Residue', 'Segid')
        
        self.table = DTable(frame.frame, (0, 2), header = header,
                            relief = Tix.FLAT)
        
        patches = self.project.getData(DC.DATA_CISPROPATCH)
        n_rows = max(10, len(patches))

        self.more_rows(n_rows)

        ## fill table with data stored in current project

        for i in range(len(patches)):
            for j in range(len(self.entity_names)):
                name = self.entity_names[j]
                self.table[i,j] = patches[i][name]

        self.table.grid(row = 0, column = 0, columnspan = 2)
##        self.table.propagate(0)

        box = Tix.ButtonBox(frame.frame, orientation=Tix.HORIZONTAL,
                            relief = Tix.FLAT)
        
        descr = tools.wrap_string(self.table_descr,
                                  AriaBaseClass.description_length)
        
        self.get_owner().bind_tooltips(box, descr)
        
        box.add('add', text = 'Add rows', width = 10,
                command = self.more_rows)
        
        box.add('del', text = 'Delete row', width = 10,
                command = self.del_row)

        box.grid(row = 1, column = 0, columnspan = 2,
                 sticky = Tix.EW, pady = 5)

        frame.form(left = 0, top = 0)

        return None

    
class DataSSBondPanel(DataPanel, HisPatchPanel):
    
    table_descr = 'List the cysteines involved in ambigous disulfide bonds.\nThe CNS DISU patch will be applied on them.'
    entity_names = ('residue', 'segid')

    def save(self):

        import aria.DataContainer as DC

        patches = []

        for row in self.table.rows:

            patch = DC.CysPatch()
            ok = 0

            for column in range(len(self.entity_names)):
                entity_name = self.entity_names[column]

                cell = row[column]

                try:
                    val = cell.get()
                except:
                    self.not_saved(cell.editor)
                    val = None
                    return self.set_modified(1)

                if val is not None:

                    cell.save_now()
                    val = row[column].get()

                    ## cast values

                    entity = cell.editor.get_entity()
                    val = entity.cast(val)
                    patch[entity_name] = val
                    ok += 1
                else:
                    break

            if ok == len(self.entity_names):
                patches.append(patch)


        self.get_settings()['data_ssbonds']['cyspatch'] = patches

        self.set_modified(0)
        
        return 1

    def set_settings(self, settings):

        check_type(settings, 'SSBondData')

        self.__settings = {'data_ssbonds': settings}

    def init_cells(self, indices):

        from aria.DataContainer import CysPatch
        
        for i in indices:

            patch = CysPatch()
            patch.reset()
                        
            ## 0th column: residue

            e = patch.getEntity('residue')
            self.table.set_widget(i, 0, DEntry, e)

            ## 1st column: segid

            e = patch.getEntity('segid')
            self.table.set_widget(i, 1, DEntry, e)          

    def add_ccpn_names(self):

        s = self.get_settings()['data_ssbonds']
        
        if not s['ccpn_id']:
            return
         
        from aria.importFromCcpn import getCcpnConstraintList, getKeysFromString

        ccpn_project = self.project.ccpn_project
        if not ccpn_project:            
            ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return
        
        if s['ccpn_id']:
            cList = getCcpnConstraintList(ccpn_project, getKeysFromString(s['ccpn_id']))
            name = "%d:%d:%s" % (cList.nmrConstraintStore.serial, cList.serial, cList.name)
            self.update_ccpn_name(self.ccpn_name, name)
            
    ## BARDIAUX        
    def select_constraints(self):

        from ccpnGui import selectConstraintList

        ccpn_project = self.read_ccpn_model()

        kw = {'restraint_type' : 'DistanceConstraintList',
              'msg' : 'Disulfide Bonds'}
        
        key, msg, ccpn_name = selectConstraintList(ccpn_project, self.parent, 'DistanceConstraintList', 'Disulfide Bonds')
	self.message(msg)

        if key:
            settings = self.get_settings()['data_ssbonds']

            settings['ccpn_id'] = key

            self.update_ccpn_name(self.ccpn_name, ccpn_name)
            
            self.refresh(overwrite=0)
            

    def create(self, p):
        
        # BARDIAUX 2.2 modified to introduce CysPatch
        
        from aria.OrderedDict import OrderedDict

        self.parent = p

        controls = []

        generic = Tix.LabelFrame(p, label = 'Generic')
        
        data_formats = ({'label': 'TBL', 'value': 'tbl'},
                        {'label': 'CCPN', 'value': 'ccpn'})

        file_types = [('TBL files', '*.tbl')]
        
        
        file, format, ccpn_id, ccpn_name = file_and_format(generic.frame, data_formats,
                                                file_types = file_types,
                                                with_ccpn_id=1, ccpn_select_func=self.select_constraints)

        self.ccpn_name = ccpn_name
        
        self.add_edit_button(generic.frame, 1, 'data_ssbonds')
        
        controls.append({'widget': file, 'entity': 'filename'})
        controls.append({'widget': format, 'entity': 'format'})
        controls.append({'widget': ccpn_id, 'entity': 'ccpn_id'})
        
        defs = OrderedDict()
        
        defs['Use restraints:'] = {'widget': 'choice', 'entity': 'enabled',
                                   'choices': YES_NO_DICT}
        
        controls += self.create_widgets(generic.frame, defs, row = 2)
        
        generic.form(top = 0, left = 0, right = -1)# BARDIAUX
        
        ## BARDIAUX 2.2 Table
        import aria.tools as tools
        import aria.DataContainer as DC

        cys = Tix.LabelFrame(p, label = 'Cysteines involved')

        header = ('Residue', 'Segid')
        
        self.table = DTable(cys.frame, (0, 2), header = header,
                            relief = Tix.FLAT)
        
        patches = self.get_settings()['data_ssbonds']['cyspatch']
        n_rows = max(10, len(patches))

        self.more_rows(n_rows)

        ## fill table with data stored in current project

        for i in range(len(patches)):
            for j in range(len(self.entity_names)):
                name = self.entity_names[j]
                self.table[i,j] = patches[i][name]

        self.table.grid(row = 0, column = 0, columnspan = 2)

        box = Tix.ButtonBox(cys.frame, orientation=Tix.HORIZONTAL,
                            relief = Tix.FLAT)
        
        descr = tools.wrap_string(self.table_descr,
                                  AriaBaseClass.description_length)
        
        self.get_owner().bind_tooltips(cys.frame, descr)
        
        box.add('add', text = 'Add rows', width = 10,
                command = self.more_rows)
        
        box.add('del', text = 'Delete row', width = 10,
                command = self.del_row)

        box.grid(row = 1, column = 0, columnspan = 2,
                 sticky = Tix.EW, pady = 5)

        cys.form(left = 0, top = generic)
    
        return {'data_ssbonds': controls}


    def get_settings(self):

        return self.__settings

# BARDIAUX 2.2
# ZN patches
class DataZnPatchPanel(DataPanel):
    
##     def delete(self, project):
##         project.delData(self.get_settings()['generic'])


    def set_settings(self, settings):

        check_type(settings, 'ZnPatch')

        self.__settings = {'data': settings} # BARDIAUX rMat

    def create(self, p):

        from aria.OrderedDict import OrderedDict
        import aria.tools as tools

        self.parent = p

        CW = self.create_widgets

        generic = Tix.LabelFrame(p, label = 'Zinc tetrahedral coordination')

        defs = OrderedDict()

        controls = []

        defs['Type of coordination:'] = {'widget': 'choice',
                                         'entity': 'type',
                                         'choices': ({'label': 'Cys4', 'value': 'SSSS'},
                                                     {'label': 'Cys3Hisd', 'value': 'SSSD'},
                                                     {'label': 'Cys3Hise', 'value': 'SSSE'})}
        
        controls += CW(generic.frame, defs, row=2)
        generic.form(top = 0, left = 0, right = -1)

        resid_zn = Tix.LabelFrame(p, label = 'Zinc Ion')
        defs = OrderedDict()
        
        defs['Residue Number:']  = {'widget': 'entry', 'entity': 'residue_zn'}
        defs['Residue Segid:']  = {'widget': 'entry', 'entity': 'segid_zn'}
        controls += CW(resid_zn.frame, defs)

        resid_zn.form(top = generic, left = 0, right = -1)
        
        resid1 = Tix.LabelFrame(p, label = 'Cysteine 1')
        defs = OrderedDict()
        
        defs['Residue Number:']  = {'widget': 'entry', 'entity': 'residue1'}
        defs['Residue Segid:']  = {'widget': 'entry', 'entity': 'segid1'}
        controls += CW(resid1.frame, defs, row=3)

        resid1.form(top = resid_zn, left = 0, right = -1)

        msg = "The 1st Cysteine bound to the Zinc ion"
        descr = tools.wrap_string(msg, AriaBaseClass.description_length)
        self.get_owner().bind_tooltips(resid1, descr)    

        resid2 = Tix.LabelFrame(p, label = 'Cysteine 2')
        defs = OrderedDict()
        
        defs['Residue Number:']  = {'widget': 'entry', 'entity': 'residue2'}
        defs['Residue Segid:']  = {'widget': 'entry', 'entity': 'segid2'}
        controls += CW(resid2.frame, defs, row=3)

        resid2.form(top = resid1, left = 0, right = -1)
        
        msg = "The 2nd Cysteine bound to the Zinc ion."
        descr = tools.wrap_string(msg, AriaBaseClass.description_length)
        self.get_owner().bind_tooltips(resid2, descr)
        
        resid3 = Tix.LabelFrame(p, label = 'Cysteine 3')
        defs = OrderedDict()
        
        defs['Residue Number:']  = {'widget': 'entry', 'entity': 'residue3'}
        defs['Residue Segid:']  = {'widget': 'entry', 'entity': 'segid3'}
        controls += CW(resid3.frame, defs,row=3)

        resid3.form(top = resid2, left = 0, right = -1)
        
        msg = "The 3rd Cysteine bound to the Zinc ion."
        descr = tools.wrap_string(msg, AriaBaseClass.description_length)
        self.get_owner().bind_tooltips(resid3, descr)

        resid4 = Tix.LabelFrame(p, label = 'Cysteine/Histidine 4')
        defs = OrderedDict()
        
        defs['Residue Number:']  = {'widget': 'entry', 'entity': 'residue4'}
        defs['Residue Segid:']  = {'widget': 'entry', 'entity': 'segid4'}
        controls += CW(resid4.frame, defs, row=3)

        resid4.form(top = resid3, left = 0, right = -1)
        
        msg = "Either the 4th Cysteine or an Histidine (HISE/HISD) bound to the Zinc ion."
        descr = tools.wrap_string(msg, AriaBaseClass.description_length)
        self.get_owner().bind_tooltips(resid4, descr)



        return {'data' : controls}
        
    def get_settings(self):

        return self.__settings
    
class DataUnambigPanel(DataPanel):

    def set_settings(self, settings):
        check_type(settings, 'UnambiguousDistanceData')
        self.__settings = {'data': settings}

    def add_ccpn_names(self):

        s = self.get_settings()['data']

        if not s['ccpn_id']:
            return
         
        from aria.importFromCcpn import getCcpnConstraintList, getKeysFromString

        ccpn_project = self.project.ccpn_project
        if not ccpn_project:            
            ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return
        
        if s['ccpn_id']:
            cList = getCcpnConstraintList(ccpn_project, getKeysFromString(s['ccpn_id']))
            name = "%d:%d:%s" % (cList.nmrConstraintStore.serial, cList.serial, cList.name)
            self.update_ccpn_name(self.ccpn_name, name)
                        

    ## BARDIAUX
    def select_constraints(self):

        from ccpnGui import selectConstraintList

        ccpn_project = self.read_ccpn_model()
        
        kw = {'restraint_type' : 'DistanceConstraintList',
              'msg' : 'Unambiguous Distance'}
        
        key, msg, ccpn_name = selectConstraintList(ccpn_project, self.parent, **kw)
	self.message(msg)

        if key:
            settings = self.get_settings()['data']

            settings['ccpn_id'] = key

            self.update_ccpn_name(self.ccpn_name, ccpn_name)
            
            self.refresh(overwrite=0)
            

    def create(self, p):

        from aria.OrderedDict import OrderedDict

        self.parent = p

        controls = []
        
        generic = Tix.LabelFrame(p, label = 'Generic')
        
        data_formats = ({'label': 'TBL', 'value': 'tbl'},
                        {'label': 'CCPN', 'value': 'ccpn'})

        file_types = [('TBL files', '*.tbl')]
        
        file, format, ccpn_id, ccpn_name = file_and_format(generic.frame, data_formats,
                                                file_types = file_types,
                                                with_ccpn_id=1, ccpn_select_func=self.select_constraints)
        self.ccpn_name = ccpn_name
        
        self.add_edit_button(generic.frame, 1, 'data')
        
        controls.append({'widget': file, 'entity': 'filename'})
        controls.append({'widget': format, 'entity': 'format'})
        controls.append({'widget': ccpn_id, 'entity': 'ccpn_id'})

        defs = OrderedDict()
        
        defs['Use restraints:'] = {'widget': 'choice', 'entity': 'enabled',
                                   'choices': YES_NO_DICT}

        controls += self.create_widgets(generic.frame, defs, row=2)
        
        ## BARDIAUX 2.2    
        ## CCPN Specific:

        ccpn_specific = Tix.LabelFrame(p, label = 'CCPN specific')

        defs = OrderedDict()

        ## BARDIAUX 2.2
        defs['Add to Network-Anchoring:'] = {'widget': 'choice', 'entity': 'add_to_network',
                                             'choices': YES_NO_DICT}

        choices = ({'label': 'All iterations', 'value': 'all_iterations'},
                   {'label': 'Not 1st iteration', 'value': 'all_iterations_except_first'},
                   {'label': 'Disable', 'value': NO})
        
        defs['Calibration:'] = {'widget': 'choice', 'entity': 'calibrate',
                                'choices': choices}
        
        defs['Filter contributions:'] = {'widget': 'choice', 'entity': 'filter_contributions',
                                         'choices': YES_NO_DICT}

        defs['Run Network-Anchoring:'] = {'widget': 'choice', 'entity': 'run_network_anchoring',
                                         'choices': YES_NO_DICT}
        
        controls += self.create_widgets(ccpn_specific.frame, defs)#, row = 2)
        
        generic.form(top = 0, left = 0, right = -1)
        ccpn_specific.form(top = generic, left = 0, right = -1)

##         #
        
##         #file, format = file_and_format(p, ({'value': 'tbl', 'label': 'TBL'}, ))

##         data_formats = ({'label': 'TBL', 'value': 'tbl'},
##                         {'label': 'CCPN', 'value': 'ccpn'})

##         file_types = [('TBL files', '*.tbl')]
        
##         file, format, ccpn_id = file_and_format(p, data_formats,
##                                                 file_types = file_types,
##                                                 with_ccpn_id=1, ccpn_select_func=self.select_constraints)
      
##         self.add_edit_button(p, 1, 'data')
        
##         controls.append({'widget': file, 'entity': 'filename'})
##         controls.append({'widget': format, 'entity': 'format'})
##         controls.append({'widget': ccpn_id, 'entity': 'ccpn_id'})

##         defs = OrderedDict()
        
##         defs['Use restraints:'] = {'widget': 'choice', 'entity': 'enabled',
##                                    'choices': YES_NO_DICT}

##         #controls += self.create_widgets(p, defs, row = 2)

##         ## BARDIAUX 2.2
##         defs['Add to Network-Anchoring (CCPN format required):'] = {'widget': 'choice', 'entity': 'add_to_network',
##                                                                     'choices': YES_NO_DICT}

##         controls += self.create_widgets(p, defs, row = 2)        
       
        return {'data': controls}


    def get_settings(self):
        return self.__settings

class DataAmbigPanel(DataUnambigPanel):
    def set_settings(self, settings):
        check_type(settings, 'AmbiguousDistanceData')
        self.__settings = {'data': settings}

    ## BARDIAUX
    def select_constraints(self):

        from ccpnGui import selectConstraintList

        ccpn_project = self.read_ccpn_model()
        
        kw = {'restraint_type' : 'DistanceConstraintList',
              'msg' : 'Ambiguous Distance'}
        
        key, msg, ccpn_name = selectConstraintList(ccpn_project, self.parent, **kw)
	self.message(msg)

        if key:
            settings = self.get_settings()['data']

            settings['ccpn_id'] = key

            self.update_ccpn_name(self.ccpn_name, ccpn_name)
            
            self.refresh(overwrite=0)
            
            
    def get_settings(self):
        return self.__settings

class CNSPanel(PanelEx):
    
    def get_defs(self):
        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        DICT = YES_NO_DICT + [{'label': 'Always', 'value': ALWAYS}]

        d['Local CNS executable:'] = {'widget': 'file',
                                      'entity': 'local_executable'}
        
        d['Keep CNS output files:'] = {'widget': 'choice',
                                       'entity': 'keep_output',
                                       'choices': YES_NO_GZIP_DICT}
        
        d['Keep ambig.tbl / unambig.tbl'] = {'widget': 'choice',
                                             'entity': 'keep_restraint_files',
                                             'choices': YES_NO_GZIP_DICT}
        
        d['Create PSF file:'] = {'widget': 'choice',
                                 'entity': 'create_psf_file',
                                 'choices': DICT}

        d['Generate template PDB file:'] = {'widget': 'choice',
                                            'entity': 'generate_template',
                                            'choices': DICT}

        choices = ({'label': 'PROLSQ', 'value': 'PROLSQ'},
                   {'label': 'PARMALLH6', 'value': 'PARMALLH6'},
                   {'label': 'PARALLHDG', 'value': 'PARALLHDG'},
                   {'label': 'OPLSX', 'value': 'OPLSX'})
        
        d['Non-bonded parameters:'] = {'widget': 'choice',
                                       'entity': 'nonbonded_parameters',
                                       'choices': choices}
        return d

    def create(self, p):
        widgets = self.create_widgets(p, self.get_defs())
        return {'cns': widgets}

    def set_settings(self, s):
        self.__settings = {'cns': s}

    def get_settings(self):
        return self.__settings

## < Mareuil
class ModePanel(PanelEx):
    def Mode_defs(self):
        from aria.OrderedDict import OrderedDict 
        widgets =  OrderedDict()
        widgets['Mode'] = {'widget': 'choice',
                'entity': 'default_command',
                'choices': ({'label': 'LOCAL', 'value': 'LOCAL'},
                            {'label': 'CLUSTER', 'value': 'CLUSTER'},
                            {'label': 'GRID', 'value': 'GRID'})}
			    
        widgets['Job management'] = {'widget': 'choice',
                'entity': 'job_management',
                'choices': ({'label': 'Synchro', 'value': 'Synchro'},
                            {'label': 'Asynchro', 'value': 'Asynchro'})}
        return widgets    
	
    def create(self, p):
        widgets = self.create_widgets(p, self.Mode_defs())
        return {'Mode': widgets}

    def set_settings(self, s):
        self.__settings = {'Mode': s}

    def get_settings(self):
        return self.__settings
## Mareuil >

class JobManagerPanel(PanelEx):

    table_descr = \
'''
List of machines available for structure calculation. Every machine must meet all of the following requirements:

1. CNS binary must be accessible (cf. column "CNS executable")
2. The temporary directory, specified in the panel "Project", must be read/write accessible .
3. The machine must accessible without prompting for a password.

To start a remote job, ARIA invokes the following command: [command] [script]

The [command] is used to launch a csh [script] which is generated on the fly. The user can neither modify its name nor its contents. When running ARIA on a single machine, you may leave the command blank ("csh" is the default). For a parallel setup, you have to tell ARIA how to establish a connection to either of the machines - usually via "rsh" or "ssh".

Example: if you want to run calculations on the machine "my_machine", using "ssh", the command should be: "ssh my_machine csh". In that case, ARIA starts the following remote job:

         ssh my_machine csh [script]     (1)
         
When using multi-processor machines or queuing systems (cf. column #CPUs), ARIA runs <#CPUs> jobs in parallel by using command (1) several times in a row. Some systems must be supplied with the local name of the csh-script rather than the absolute name (i.e. refine.csh instead of /some/directory/refine.csh). In that case, set the field "use absolute path" to "no".

To append a new row to the table, press the button "Add row". To delete a row, first select any cell in that row and press "Delete row". To modify an entry field, <double-click> the cell. To keep the changes, press <return> or <click> another cell; to discard changed, press <escape>.
'''
    
    def save(self, verbose = 1):

        from aria.JobManager import HostSettings

        settings = self.get_settings()['job_manager']['host_list']

        hosts = []

        d = {'enabled': 0,
             'submit_command': 1,
             'state_command': 2,
             'output_command': 3,
             'n_cpu': 4,
             'executable': 5,
             'use_absolute_path': 6}

        for row in self.table.rows:

            hs = HostSettings()

            ok = 0

            for entity_name, column in d.items():
            
                cell = row[column]
                if entity_name == "state_command" or entity_name == "output_command":
                    try:
                        val = cell.get()
                        if val == None:
                            val = ''
                    except:
                        val = ''
                if entity_name != "state_command" and entity_name != "output_command":
                    try:
                        val = cell.get()
                    except:
                        self.not_saved(cell.editor)
                        val = None
                        return self.set_modified(1)

                if val is not None:
                    if val.strip() == '' and entity_name != 'state_command' and entity_name != 'output_command':
                        continue
                    
                    cell.save_now()
                    #val = row[column].get()
                    if entity_name == 'n_cpu':
                        val = int(val)
                    hs[entity_name] = val
                    ok += 1
                else:
                    break

            if ok == len(d):
                hosts.append(hs)

        settings = self.get_settings()['job_manager']
        settings['host_list'] = hosts

        if not len(hosts) and verbose:
            s = 'List of available hosts is void or incomplete. If you start a calculation, the structure calculation engine will be run on the local host.'
            self.message(s)

         
        return 1
    
    def validate(self):
        """
        Checks whether all cells are valid.
        """
        for row in self.table.rows:
            for r in row:
                is_valid = r.isValid()
                if is_valid == 0:
                    return self.set_modified()

        ## all cells are valid
                
        return self.set_modified(0)

    def not_saved(self, widget):
        from tkMessageBox import showerror

        if widget is not None:
            try:
                widget.save()
            except Exception, msg:
                showerror('Value error',msg)
        
        self.set_modified()
    
    def more_rows(self, n = 1):

        l = len(self.table.rows)
        indices = range(l, l + n)
        
        self.table.append_row(n)
        self.init_cells(indices)

    def del_row(self):
        current_row = self.table.get_current_row_index()
        if current_row is None:
            return
        self.table.remove_row(current_row)

    def init_cells(self, indices):

        from aria.JobManager import HostSettings

        gui_settings = self.get_settings()['gui']
        default_cns_executable = gui_settings['cns_executable']
        
        for i in indices:

            hs = HostSettings()

            ## 0th column: enabled

            e = hs.getEntity('enabled')
            self.table.set_widget(i, 0, DOptionMenu, e)
            editor = self.table.rows[i][0].editor
            editor.add_command(YES, label = 'Yes')
            editor.add_command(NO, label = 'No')
            self.table.rows[i][0].set(e.get_default_value())

            e = hs.getEntity('submit_command')
            self.table.set_widget(i, 1, DEntry, e)
            
            e = hs.getEntity('state_command')
            self.table.set_widget(i, 2, DEntry, e)
            
            e = hs.getEntity('output_command')
            self.table.set_widget(i, 3, DEntry, e)
  
            ## 2nd column: #CPUs

            e = hs.getEntity('n_cpu')
            self.table.set_widget(i, 4, DEntry, e)
            self.table[i,4] = 1

            ## 3rd column consists of DFileEntry widgets

            e = hs.getEntity('executable')
            self.table.set_widget(i, 5, DFileEntry, e)
            self.table[i,5] = default_cns_executable

            ## 4th column: use_absolute_path

            e = hs.getEntity('use_absolute_path')
            self.table.set_widget(i, 6, DOptionMenu, e)
            editor = self.table.rows[i][6].editor
            editor.add_command(YES, label = 'Yes')
            editor.add_command(NO, label = 'No')
            self.table.rows[i][6].set(e.get_default_value())

    def create(self, p):

        import aria.tools as tools
        frame = Tix.LabelFrame(p, label = 'Host-list')
        header = ('Enabled', 'Submit Command', 'State Command', 'Output Command', '#CPUs', 'CNS executable', 'Use absolute path')
        
        c_widths = (8, 20, 20, 20, 7, 20, 8)
        self.table = DTable(frame.frame, (0, len(c_widths)), header = header,
                            column_widths = c_widths, relief = Tix.FLAT)

        table_size = 15
        self.more_rows(table_size)
        
        ## fill table with host-list data

        hosts = self.get_settings()['job_manager']['host_list']
        n_rows = len(hosts)
        self.more_rows(max(n_rows - table_size, 0))

        kw = ('enabled', 'submit_command', 'state_command', 'output_command', 'n_cpu', 'executable', 'use_absolute_path')

        for i in range(n_rows):
            for j in range(len(kw)):
                self.table[i,j] = hosts[i][kw[j]]

        self.table.grid(row = 0, column = 1, columnspan = 2)
##        self.table.grid_propagate(0)

        box = Tix.ButtonBox(frame.frame, orientation=Tix.VERTICAL,
                            relief = Tix.FLAT)

        descr = tools.wrap_string(self.table_descr,
                                  AriaBaseClass.description_length)
        
        self.get_owner().bind_tooltips(box, descr)
        
        box.add('add', text = 'Add row', width = 10,
                command = self.more_rows, borderwidth = 1)
        
        box.add('del', text = 'Delete row', width = 10,
                command = self.del_row, borderwidth = 1)

        box.grid(row = 0, column = 0, columnspan = 1,
                 sticky = Tix.N, pady = 5)
        
        frame.form(left = 0, top = 0)

        return None

    def set_settings(self, job_manager_settings, gui_settings):
        self.__settings = {'job_manager': job_manager_settings,
                           'gui': gui_settings}
        
    def get_settings(self):
        return self.__settings
    ## Mareuil >
    
class AnnealingAmbigPanel(PanelEx):

    def get_defs(self):
        
        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['Use ambiguous restraint starting at iteration:'] = \
               {'widget': 'entry', 'entity': 'first_iteration'}
        d['Force constant high temp.:'] = \
                 {'widget': 'entry', 'entity': 'k_hot'}
        d['Initial force constant cool1:'] = \
                   {'widget': 'entry', 'entity': 'k_cool1_initial'}
        d['Final force constant cool1:'] = \
                 {'widget': 'entry', 'entity': 'k_cool1_final'}
        d['Force constant cool2:'] = \
                 {'widget': 'entry', 'entity': 'k_cool2'}
        
        return d

    def create(self, p):
        widgets = self.create_widgets(p, self.get_defs())
        return {'annealing_ambig': widgets}

    def set_settings(self, s):
        self.__settings = {'annealing_ambig': s}

    def get_settings(self):
        return self.__settings

class AnnealingUnambigPanel(PanelEx):

    def get_defs(self):
        
        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['Use unambiguous restraints starting at iteration:'] = \
               {'widget': 'entry', 'entity': 'first_iteration'}
        d['Force constant high temp.:'] = \
                 {'widget': 'entry', 'entity': 'k_hot'}
        d['Initial force constant cool1:'] = \
                   {'widget': 'entry', 'entity': 'k_cool1_initial'}
        d['Final force constant cool1:'] = \
                 {'widget': 'entry', 'entity': 'k_cool1_final'}
        d['Force constant cool2:'] = \
                 {'widget': 'entry', 'entity': 'k_cool2'}
        
        return d

    def create(self, p):
        widgets = self.create_widgets(p, self.get_defs())
        return {'annealing_unambig': widgets}

    def set_settings(self, s):
        self.__settings = {'annealing_unambig': s}

    def get_settings(self):
        return self.__settings

class AnnealingHBondPanel(PanelEx):

    def get_defs(self):
        
        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['Use H-Bonds starting at iteration:'] = {'widget': 'entry',
                                                   'entity': 'first_iteration'}
        d['Force constant high temp.:'] = \
                 {'widget': 'entry', 'entity': 'k_hot'}
        d['Initial force constant cool1:'] = \
                   {'widget': 'entry', 'entity': 'k_cool1_initial'}
        d['Final force constant cool1:'] = \
                 {'widget': 'entry', 'entity': 'k_cool1_final'}
        d['Force constant cool2:'] = \
                 {'widget': 'entry', 'entity': 'k_cool2'}
        
        return d

    def create(self, p):
        widgets = self.create_widgets(p, self.get_defs())
        return {'annealing_dihedrals': widgets}

    def set_settings(self, s):
        self.__settings = {'annealing_dihedrals': s}

    def get_settings(self):
        return self.__settings

## BARDIAUX 2.2
class AnnealingSymmetryPanel(PanelEx):

    def get_defs(self):
        
        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['NCS restraints force constant:'] = {'widget': 'entry',
                                               'entity': 'k_ncs'}
        d['Packing Force constant hot:'] = \
                   {'widget': 'entry', 'entity': 'k_packing_hot'}
        
        d['Packing Force constant cool1:'] = \
                   {'widget': 'entry', 'entity': 'k_packing_cool1'}
        
        d['Packing Force constant cool2:'] = \
                   {'widget': 'entry', 'entity': 'k_packing_cool2'}
        
        d['Use Packing until iteration:'] = \
                   {'widget': 'entry', 'entity': 'last_iteration_packing'}
        return d

    def create(self, p):
        widgets = self.create_widgets(p, self.get_defs())
        return {'annealing_symmetry': widgets}

    def set_settings(self, s):
        self.__settings = {'annealing_symmetry': s}

    def get_settings(self):
        return self.__settings
    
class AnnealingDihedralPanel(PanelEx):

    def get_defs(self):
        
        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['Force constant high temp.:'] = \
                 {'widget': 'entry', 'entity': 'k_hot'}
        d['Force constant cool1:'] = \
                 {'widget': 'entry', 'entity': 'k_cool1'}
        d['Force constant cool2:'] = \
                 {'widget': 'entry', 'entity': 'k_cool2'}
        
        return d

    def create(self, p):
        widgets = self.create_widgets(p, self.get_defs())
        return {'annealing_dihedrals': widgets}

    def set_settings(self, s):
        self.__settings = {'annealing_dihedrals': s}

    def get_settings(self):
        return self.__settings

class AnnealingKarplusPanel(PanelEx):

    def get_defs(self):
        
        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['a:'] = {'widget': 'entry', 'entity': 'a'}
        d['b:'] = {'widget': 'entry', 'entity': 'b'}
        d['c:'] = {'widget': 'entry', 'entity': 'c'}
        d['d:'] = {'widget': 'entry', 'entity': 'd'}
        d['Force constant high temp.:'] = \
                 {'widget': 'entry', 'entity': 'k_hot'}
        d['Force constant cool1:'] = \
                 {'widget': 'entry', 'entity': 'k_cool1'}
        d['Force constant cool2:'] = \
                 {'widget': 'entry', 'entity': 'k_cool2'}
        
        return d

    def create(self, p):
        widgets = self.create_widgets(p, self.get_defs())
        return {'annealing_karplus': widgets}

    def set_settings(self, s):
        self.__settings = {'annealing_karplus': s}

    def get_settings(self):
        return self.__settings

class AnnealingFBHWPanel(PanelEx):

    def get_defs(self):
        
        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['Left switch-point high temp.:'] = {'widget': 'entry',
                                              'entity': 'm_rswitch_hot'}
        d['Left switch-point cool1:'] = {'widget': 'entry',
                                         'entity': 'm_rswitch_cool1'}
        d['Left switch-point cool2:'] = {'widget': 'entry',
                                         'entity': 'm_rswitch_cool2'}

        d['Left asymptote high temp.:'] = {'widget': 'entry',
                                           'entity': 'm_asymptote_hot'}
        d['Left asymptote cool1:'] = {'widget': 'entry',
                                      'entity': 'm_asymptote_cool1'}
        d['Left asymptote cool2:'] = {'widget': 'entry',
                                      'entity': 'm_asymptote_cool2'}

        d['Right switch-point high temp.:'] = {'widget': 'entry',
                                              'entity': 'rswitch_hot'}
        d['Right switch-point cool1:'] = {'widget': 'entry',
                                         'entity': 'rswitch_cool1'}
        d['Right switch-point cool2:'] = {'widget': 'entry',
                                         'entity': 'rswitch_cool2'}

        d['Right asymptote high temp.:'] = {'widget': 'entry',
                                            'entity': 'asymptote_hot'}
        d['Right asymptote cool1:'] = {'widget': 'entry',
                                       'entity': 'asymptote_cool1'}
        d['Right asymptote cool2:'] = {'widget': 'entry',
                                       'entity': 'asymptote_cool2'}
         
        return d

    def create(self, p):
        widgets = self.create_widgets(p, self.get_defs())
        return {'annealing_fbhw': widgets}

    def set_settings(self, s):
        self.__settings = {'annealing_fbhw': s}

    def get_settings(self):
        return self.__settings

class AnnealingRDCPanel(PanelEx):

    def get_sa_defs(self):

        from aria.OrderedDict import OrderedDict

        d = OrderedDict()
        
        d['High temperature:'] = \
                {'widget': 'entry', 'entity': 'k_hot'}
        d['Cool1 temperature:'] = \
                 {'widget': 'entry', 'entity': 'k_cool1'}
        d['Cool2 temperature:'] = \
                 {'widget': 'entry', 'entity': 'k_cool2'}
         
        return d

    def get_potential_shape_defs(self):

        from aria.OrderedDict import OrderedDict

        d = OrderedDict()
        
        d['Border initial high temp.:'] = \
                   {'widget': 'entry', 'entity': 'border_hot_initial'}
        d['Border final high temp.:'] = \
                 {'widget': 'entry', 'entity': 'border_hot_final'}
        d['Center initial high temp.:'] = \
                 {'widget': 'entry', 'entity': 'center_hot_initial'}
        d['Center final high temp.:'] = \
                 {'widget': 'entry', 'entity': 'center_hot_final'}
        d['Border initial cool1:'] = \
                   {'widget': 'entry', 'entity': 'border_cool1_initial'}
        d['Border final cool1:'] = \
                 {'widget': 'entry', 'entity': 'border_cool1_final'}
        d['Center intiial cool1:'] = \
                 {'widget': 'entry', 'entity': 'center_cool1_initial'}
        d['Center final cool1:'] = \
                 {'widget': 'entry', 'entity': 'center_cool1_final'}
        d['Border initial cool2:'] = \
                   {'widget': 'entry', 'entity': 'border_cool2_initial'}
        d['Border final cool2:'] = \
                 {'widget': 'entry', 'entity': 'border_cool2_final'}
        d['Center initial cool2:'] = \
                 {'widget': 'entry', 'entity': 'center_cool2_initial'}
        d['Center final cool2:'] = \
                 {'widget': 'entry', 'entity': 'center_cool2_final'}
         
        return d

    def get_generic_defs(self):

        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['Type:'] = {'widget': 'choice',
                      'choices': ({'label': 'SANI', 'value': 'SANI'},
                                  {'label': 'VANGLE', 'value': 'VANGLE'}),
                      'entity': 'method'}

        d['Use RDCs from iteration:'] = \
               {'widget': 'entry', 'entity': 'first_iteration'}
         
        return d

    def get_tensor_defs(self):

        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['Rhombicity:'] = {'widget': 'entry', 'entity': 'r'}
        d['Magnitude:'] = {'widget': 'entry', 'entity': 'd'}
         
        return d

    def create(self, p):

        CW = self.create_widgets

        widgets = {}

        ## general

        label = 'Generic'
        f_generic = Tix.LabelFrame(p, label = label)
        defs = self.get_generic_defs()
        widgets['generic'] = CW(f_generic.frame, defs)

        ## alignment tensor

        label = 'Alignment tensor'
        f_tensor = Tix.LabelFrame(p, label = label)
        defs = self.get_tensor_defs()
        widgets['tensor'] = CW(f_tensor.frame, defs)

        ## Force constants for SA protocol

        label = 'Force constants: SA protocol'
        f_sa = Tix.LabelFrame(p, label = label)
        defs = self.get_sa_defs()
        widgets['sa_protocol'] = CW(f_sa.frame, defs)

        ## potential shape

        label_potential = 'Potential shape (VEAN)'
        f_potential = Tix.LabelFrame(p, label = label_potential)
        defs = self.get_potential_shape_defs()
        widgets['potential_shape'] = CW(f_potential.frame, defs)

        f_generic.form(top = 0, left = 0, right = -1, fill=Tix.X)
        f_tensor.form(top = f_generic, left = 0, right = -1, fill = Tix.X)
        f_sa.form(top = f_tensor, left = 0, right = -1, fill = Tix.X)
        f_potential.form(top = f_sa, left = 0, right = -1, fill = Tix.X)
        
        return widgets

    def set_settings(self, s):
        self.__settings = {'generic': s,
                           'tensor': s,
                           'sa_protocol': s,
                           'potential_shape': s}

    def get_settings(self):
        return self.__settings

class DynamicsPanel(Panel):

    def create(self, p):
        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['Type:'] = {'widget': 'choice',
                            'choices':
                            ({'label': 'Torsion angle', 'value': 'torsion'},
                             {'label': 'Cartesian',
                              'value': 'cartesian'}),
                            'entity': 'md_type'}

        widgets['Random seed:'] = \
                        {'widget': 'entry', 'entity': 'random_seed'}
        widgets['TAD high temperature:'] = \
                     {'widget': 'entry', 'entity': 'tad_temp_high'}
        widgets['TAD time-step factor:'] = \
                     {'widget': 'entry', 'entity': 'tad_timestep_factor'}
        widgets['Cartesian High temperature:'] = \
                           {'widget': 'entry',
                            'entity': 'cartesian_temp_high'}
        widgets['Cartesian 1st iteration:'] = \
                           {'widget': 'entry',
                            'entity': 'cartesian_first_iteration'}
        widgets['Time-step:'] = {'widget': 'entry', 'entity': 'timestep'}
        widgets['Cool1 final temperature:'] = \
                       {'widget': 'entry', 'entity': 'temp_cool1_final'}
        widgets['Cool2 final temperature:'] = \
                       {'widget': 'entry', 'entity': 'temp_cool2_final'}
        widgets['High-temp steps:'] = {'widget': 'entry',
                                       'entity': 'steps_high'}
        # not used any longer
        # BARDIAUX 2.3: used with kept_structures
        widgets['Refine steps:'] = {'widget': 'entry',
                                    'entity': 'steps_refine'}
        
        widgets['Cool1 steps:'] = {'widget': 'entry', 'entity': 'steps_cool1'}
        widgets['Cool2 steps:'] = {'widget': 'entry', 'entity': 'steps_cool2'}

        i = 0
        controls = []

        for name, d in widgets.items():

            label = Tix.Label(p, text = name)
            label.grid(row = i, sticky = Tix.W)

            widget = d['widget']
            entity = d['entity']

            if widget == 'entry':
                entry = DEntry(p, width = 40)
                controls.append({'widget': entry, 'entity': entity})

            elif widget == 'choice':
                entry = DOptionMenu(p, label = '')
                for choice in d['choices']:
                    entry.add_command(choice['value'], label = choice['label'])

                controls.append({'widget': entry, 'entity': entity})
            else:
                raise ValueError, 'Unknown widget %s' % widget

            entry.grid(row = i, column = 1, sticky = Tix.W)

            i += 1

        return {'dynamics': controls}

    def get_settings(self):
        return self.__settings

    def set_settings(self, s):
        self.__settings = {'dynamics': s}

## BERNARD 2.3
"""
Nilges M., Bernard A., Bardiaux B., Malliavin T.E., Habeck M., 
Rieping W.(2008) Accurate NMR structures through minimization of
an extended hybrid energy. Structure 16,9:1305-1312
"""

class LogHarmonicPanel(PanelEx):
    
    def get_defs(self):

        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['Use log-harmonic potential (SA cool2):'] = \
                {'widget': 'choice',
                 'entity': 'enabled',
                 'choices': YES_NO_DICT}
        
        d['Use automatic restraints weight:'] =\
                {'widget': 'choice',
                 'entity': 'use_auto_weight',
                 'choices': YES_NO_DICT}
        
        d['Weight for unambig restraints:'] = \
                {'widget': 'entry', 'entity': 'weight_unambig'}

        d['Weight for ambig restraints:'] = \
                {'widget': 'entry', 'entity': 'weight_ambig'}

        d['Weight for hbond restraints:'] = \
                {'widget': 'entry', 'entity': 'weight_hbond'}

        return d

    def create(self, p):

        
        widgets = self.create_widgets(p, self.get_defs(), row = 0)

        frame = Tix.Frame(p)

        text = '''\nIf you use the Log-Harmonic potential, please quote the following 
reference:

Nilges M., Bernard A., Bardiaux B., Malliavin T.E., Habeck M., 
Rieping W.(2008) Accurate NMR structures through minimization of
an extended hybrid energy. Structure 16,9:1305-1312
'''        
        label = Tix.Label(p, text = text, justify=Tix.LEFT)
        label.grid(sticky=Tix.W, row=len(widgets), columnspan=10)
        
        return {'logharmonic_potential': widgets}

    def set_settings(self, s):
        self.__settings = {'logharmonic_potential': s}

    def get_settings(self):
        return self.__settings


class ProtocolPanel(DataPanel):

    def __init__(self, node = None):
        DataPanel.__init__(self, node, enable_delete_button = 0)

    def add_ccpn_names(self):

        s = self.get_settings()['initial_structure']
        
        if not s['ccpn_id']:
            return
         
        from aria.importFromCcpn import getCcpnModel, getKeysFromString

        ccpn_project = self.project.ccpn_project
        
        if not ccpn_project:            
            ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return
        
        if s['ccpn_id']:
            model = getCcpnModel(ccpn_project, getKeysFromString(s['ccpn_id']))
            structureEnsemble = model.structureEnsemble
            molSystem = structureEnsemble.molSystem
            name = "%s:%d:model_%d" % (molSystem.code, structureEnsemble.ensembleId, model.serial)
            self.update_ccpn_name(self.ccpn_name, name)

    def select_ccpn_model(self):

        from aria.importFromCcpn import getModels, getObjectKeyString, getKeysFromString, getCcpnChains
        from ccpnGui import gui_select_model
        from ccpnmr.format.general.Util import createSelection
        import aria.DataContainer as DC
        
        ccpn_project = self.read_ccpn_model()

        if ccpn_project is None:
            return

        mol_data = self.project.getData(DC.DATA_SEQUENCE)[0]
        
        molsystem_id = mol_data['ccpn_id']

        if molsystem_id == '':

            from widgets import OKDialog

            msg = 'Name of molecular system has not yet been specified (see entry field "CCPN name" in node "Sequence").'

            OKDialog(title='Error', text=msg)

            return
        
        chains = getCcpnChains(ccpn_project, getKeysFromString(molsystem_id))#[0]
        
        models = getModels(ccpn_project, chains)

        if not models:

            from widgets import OKDialog

            msg = 'No CCPN models present for\nmolecular system "%s".' % molsystem_id

            OKDialog(title='Error', text=msg)

            return        

        model_labels, model_dict = createSelection(models)
        
        model = gui_select_model(model_labels,
                                 model_dict,
                                 gui=self.parent)

        if model is None:
            return

        key = getObjectKeyString(model)

        settings = self.get_settings()['initial_structure']

        settings['ccpn_id'] = key

        structureEnsemble = model.structureEnsemble
        molSystem = structureEnsemble.molSystem
        name = "%s:%d:model_%d" % (molSystem.code, structureEnsemble.ensembleId, model.serial)
        self.update_ccpn_name(self.ccpn_name, name)        

        self.refresh(overwrite=0)
        
    def create(self, p):
        from aria.OrderedDict import OrderedDict

        self.parent = p
        
        LF = Tix.LabelFrame

        widgets = {}
        d = OrderedDict()

        d['Use floating chirality assignment:'] = \
               {'widget': 'choice',
                'entity': 'floating_assignment',
                'choices': YES_NO_DICT}

        f_float = LF(p, label = 'Floating Chirality Assignment')

        widgets['protocol'] = self.create_widgets(f_float.frame, d)

        label = 'Initial Structure For Minimization Protocol'
        f_initial = LF(p, label = label)
        
        controls = []

        label_format = 'PDB format: '

##         file, format = file_and_format(f_initial.frame,
##                                        ({'label': 'IUPAC',
##                                          'value': 'iupac'},
##                                         {'label': 'CNS',
##                                          'value': 'cns'},
##                                         {'label': 'Dyana',
##                                          'value': 'dyana'}),
##                                        label_format = label_format)

        data_formats = ({'label': 'IUPAC',
                         'value': 'iupac'},
                        {'label': 'CNS',
                         'value': 'cns'},
                        {'label': 'Dyana',
                         'value': 'dyana'},
                        {'label' : 'CCPN',
                         'value' : 'ccpn'})
        
        file, format, ccpn_id, ccpn_name = file_and_format(f_initial.frame, data_formats,
                                                           with_ccpn_id=1, label_format = label_format,
                                                           ccpn_select_func=self.select_ccpn_model)

        self.ccpn_name = ccpn_name
        
        controls.append({'widget': file, 'entity': 'filename'})
        controls.append({'widget': format, 'entity': 'format'})
        controls.append({'widget': ccpn_id, 'entity': 'ccpn_id'})

        n_widgets = len(controls)

        d = OrderedDict()
        
        d['Use as initial structure:'] = \
               {'widget': 'choice', 'entity': 'enabled',
                'choices': YES_NO_DICT}
        
        self.add_edit_button(f_initial.frame, n_widgets, 'initial_structure',
                             text = 'View...', data_type = DATA_PDB)

        controls += self.create_widgets(f_initial.frame, d,
                                        row = n_widgets + 1)

        widgets['initial_structure'] = controls

        f_float.form(top = 0, left = 0, right = -1, fill = Tix.X)
        f_initial.form(top = f_float, left = 0, right = -1, fill = Tix.X)
        
        return widgets

    def get_settings(self):
        return self.__settings

    def set_settings(self, project):
        import aria.DataContainer as DC
        
        s = project.getProtocol().getSettings()
        initial_structure = project.getData(DC.DATA_INITIAL_STRUCTURE)[0]
        
        self.__settings = {'protocol': s,
                           'initial_structure': initial_structure}

class WaterRefinementPanel(PanelEx):

    def set_settings(self, settings):

        check_type(settings, 'WaterRefinementParameters')

        self.__settings = {'water_refinement_params': settings}

    def create(self, p):

        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        d['Refine last iteration:'] = {'widget': 'choice',
                                       'entity': 'enabled',
                                       'choices': YES_NO_DICT}
        
        d['Solvent:'] = {'widget': 'choice',
                         'entity': 'solvent',
                         'choices': ({'label': 'Water', 'value': 'water'},
                                     {'label': 'DMSO', 'value': 'dmso'})}
        d['Number of structures:'] = {'widget': 'entry',
                                      'entity': 'n_structures'}
        d['Write solvent molecules into PDB-files:'] = \
                 {'widget': 'choice', 'entity': 'write_solvent_molecules',
                  'choices': YES_NO_DICT}

        widgets = self.create_widgets(p, d)

        return {'water_refinement_params': widgets}

    def get_settings(self):

        return self.__settings

class IterationPanel(PanelEx):

    def delete(self, project):
        s = self.__settings['iteration']
        project.getProtocol().getSettings().delIterationSettings(s)

    def getSettings(self):
        try:
            return self.__settings['iteration']
        except:
            return None
    
    def iteration_defs(self):

        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['Number of structures:'] = {'widget': 'entry',
                                            'entity': 'number_of_structures'}
        
        widgets['Sort criterion:'] = {'widget': 'choice',
                                     'choices': ({'label': 'Total Energy',
                                                 'value': 'total_energy'},
                                                 {'label': 'Restraint Energy',
                                                 'value': 'restraint_energy'},
                                                 {'label': 'NOE violations',
                                                 'value': 'noe_violations'},
                                                 {'label': 'Restraint violations',
                                                 'value': 'restraint_violations'},),
                                     'entity': 'sort_criterion'}
                   
        widgets['Use n best structures:'] = \
                     {'widget': 'entry',
                      'entity': 'number_of_best_structures'}

        widgets['Keep n previous structures:'] = \
                     {'widget': 'entry',
                      'entity': 'number_of_kept_structures'}

        return widgets

    def assignment_defs(self):
        return {}

    def merger_defs(self):

        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['Method:'] = {'widget': 'choice',
                              'choices': ({'label': 'Standard',
                                           'value': 'standard'},
                                          {'label': 'No merging',
                                           'value': 'no_merging'},
                                          {'label': 'Combination',
                                           'value': 'combination'}),
                              'entity': 'method'}

        return widgets

    def calibrator_defs(self):
        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['Distance Cutoff [A]:'] = {'widget': 'entry',
                                       'entity': 'distance_cutoff'}
        
        widgets['Estimator:'] = {'widget': 'choice',
                                 'choices': ({'label': 'Ratio of averages',
                                              'value': 'ratio_of_averages'},
                                              {'label': 'Fixed Bound',
                                              'value': 'fixed_bound'}),
                                 'entity': 'estimator'}

        widgets['Spin Diffusion Correction:'] = {'widget': 'choice',
                                                 'entity': 'relaxation_matrix',
                                                 'choices': tuple(YES_NO_DICT)}

        widgets['Error Estimator:'] = {'widget': 'choice',
                                       'entity': 'error_estimator',
                                       'choices': ({'label': 'Distance',
                                                    'value': 'distance'},
                                                   {'label': 'Volume',
                                                    'value': 'intensity'},)}        

        return widgets

    def va_defs(self):

        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()
	
	#Mareuil
	widgets['Sigma Mode:'] = {'widget': 'choice',
                                      'entity': 'sigma_mode',
                                      'choices': ({'label': 'auto',
                                                    'value': 'auto'},
                                                   {'label': 'fix',
                                                    'value': 'fix'})} 
	#Mareuil

        widgets['Violation tolerance:'] = {'widget': 'entry',
                                          'entity': 'violation_tolerance'}
        
        widgets['Violation threshold:'] = {'widget': 'entry',
                                           'entity': 'violation_threshold'}

        return widgets

    def pa_defs(self):
        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['Ambiguity cutoff:'] = {'widget': 'entry',
                                        'entity': 'weight_cutoff'}

        widgets['Max. number of contributions:'] = \
                      {'widget': 'entry',
                       'entity': 'max_contributions'}
        return widgets

    ## BARDIAUX : Network Anchoring
    def net_defs(self):
        
        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['Enabled:'] =  {'widget': 'choice', 'entity': 'enabled',
                                   'choices': tuple(YES_NO_DICT)}

        widgets['High NA score per residue threshold:'] = \
                      {'widget': 'entry',
                       'entity': 'high_residue_threshold'}


        widgets['Minimal NA score per residue threshold:'] = \
                      {'widget': 'entry',
                       'entity': 'min_residue_threshold'}
        
        widgets['Minimal NA score per atom threshold:'] = \
                      {'widget': 'entry',
                       'entity': 'min_atom_threshold'}
        return widgets
    
    def create(self, p):

        LF = Tix.LabelFrame
        CW = self.create_widgets

        widgets = {}

        it_frame = LF(p, label = 'Generic')
        defs = self.iteration_defs()
        widgets['iteration'] = CW(it_frame.frame, defs)

        assigner_frame = LF(p, label = 'Assignment')
        defs = self.assignment_defs()
        widgets['assignment'] = CW(assigner_frame.frame, defs)

        merger_frame = LF(p, label = 'Merging')
        defs = self.merger_defs()
        widgets['merging'] = CW(merger_frame.frame, defs)

        calib_frame = LF(p, label = 'Calibration')
        defs = self.calibrator_defs()
        widgets['calibration'] = CW(calib_frame.frame, defs)

        va_frame = LF(p, label = 'Violation Analysis')
        defs = self.va_defs()
        widgets['violation_analysis'] = CW(va_frame.frame, defs)

        pa_frame = LF(p, label = 'Partial Assignment')
        defs = self.pa_defs()
        widgets['partial_assignment'] = CW(pa_frame.frame, defs)

        ## BARDIAUX : Network Anchoring
        net_frame = LF(p, label = 'Network Anchoring (NA)')
        defs = self.net_defs()
        widgets['network_anchoring'] = CW(net_frame.frame, defs)

        it_frame.form(top = 0, left = 0, right = -1)
        assigner_frame.form(top = it_frame, left = 0, right = -1)
        merger_frame.form(top = assigner_frame, left = 0, right = -1)
        calib_frame.form(top = merger_frame, left = 0, right = -1)
        va_frame.form(top = calib_frame, left = 0, right = -1)
        pa_frame.form(top = va_frame, left = 0, right = -1)
        ## BARDIAUX 2.2
        net_frame.form(top = pa_frame, left = 0, right = -1)
        
        return widgets

    def set_settings(self, iteration_settings):

        d = {'assignment': 'peak_assigner_settings',
             'merging': 'merger_settings',
             'calibration': 'calibrator_settings',
             'violation_analysis': 'violation_analyser_settings',
             'partial_assignment': 'contribution_assigner_settings',
             'network_anchoring' : 'network_anchoring_settings'}

        self.__settings = {}

        for name, setting in d.items():
            self.__settings[name] = iteration_settings[setting]

        self.__settings['iteration'] = iteration_settings

    def get_settings(self):
        return self.__settings

    def show(self):

        PanelEx.show(self)

        ## enable control-panels delete button

        if self.get_owner() is not None:
            self.get_owner().button_delete.config(state = Tix.NORMAL)

class AnalysisPanel(PanelEx):

    def create(self, p):

        from aria.OrderedDict import OrderedDict

        LF = Tix.LabelFrame
        CW = self.create_widgets

        widgets = []

        ## procheck

        frame1 = LF(p, label = 'Procheck')
        defs = OrderedDict()
        defs['Enabled:'] = {'widget': 'choice',
                            'entity': 'procheck_enabled',
                            'choices': YES_NO_DICT}

        defs['Executable:'] = {'widget': 'file',
                               'entity': 'procheck_executable'}
                                        
        widgets += CW(frame1.frame, defs)

        ## whatif

        frame3 = LF(p, label = 'Whatif')
        defs = OrderedDict()
        defs['Enabled:'] = {'widget': 'choice',
                            'entity': 'whatif_enabled',
                            'choices': YES_NO_DICT}
        
        defs['Executable:'] = {'widget': 'file',
                               'entity': 'whatif_executable'}
                                        
        widgets += CW(frame3.frame, defs)

        ## prosa

        frame4 = LF(p, label = 'Prosa')
        defs = OrderedDict()
        defs['Enabled:'] = {'widget': 'choice',
                            'entity': 'prosa_enabled',
                            'choices': YES_NO_DICT}

        defs['Executable:'] = {'widget': 'file',
                               'entity': 'prosa_executable'}
                                        
        widgets += CW(frame4.frame, defs)

        ## prosa

        frame5 = LF(p, label = 'Molprobity Clashlist')
        defs = OrderedDict()
        defs['Enabled:'] = {'widget': 'choice',
                            'entity': 'clashlist_enabled',
                            'choices': YES_NO_DICT}

        defs['Executable:'] = {'widget': 'file',
                               'entity': 'clashlist_executable'}
                                        
        widgets += CW(frame5.frame, defs)
        

        ## CNS analysis scriptse

        f_cns = LF(p, label = 'CNS analysis scripts')
        defs = OrderedDict()

        defs['Enabled:'] = {'widget': 'choice',
                            'entity': 'cns_analyses',
                            'choices': YES_NO_DICT}

        widgets += CW(f_cns.frame, defs)
        
        frame1.form(top = 0, left = 0, right = -1)
        frame3.form(top = frame1, left = 0, right = -1)
        frame4.form(top = frame3, left = 0, right = -1)
        frame5.form(top = frame4, left = 0, right = -1)        
        f_cns.form(top = frame5, left = 0, right = -1)

        return {'analysis': widgets}

    def set_settings(self, settings):
        self.__settings = {'analysis': settings}

    def get_settings(self):
        return self.__settings

class PreferencesPanel(PanelEx):

    def save(self):

        val = PanelEx.save(self)
        val = val and self.jm_panel.save(verbose = 0)

        return val

    def gui_defs(self):
        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['Default author:'] = {'widget': 'entry',
                                      'entity': 'author'}

        widgets['Number of recent files:'] = {'widget': 'entry',
                                              'entity': 'n_recent_files'}

        widgets['Font size:'] = {'widget': 'entry',
                                 'entity': 'font_size'}

        widgets['CNS default executable:'] = {'widget': 'file',
                                              'entity': 'cns_executable'}

        widgets['PROSA II executable:'] = {'widget': 'file',
                                           'entity': 'prosa_executable'}
        widgets['Procheck executable:'] = {'widget': 'file',
                                           'entity': 'procheck_executable'}
        widgets['Whatif executable:'] = {'widget': 'file',
                                         'entity': 'whatif_executable'}
        
        return widgets

    def control_panel_defs(self):
        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()
        
        widgets['Text editor:'] = {'widget': 'file',
                                   'entity': 'text_editor'}
        
        widgets['PDB viewer:'] = {'widget': 'file',
                                  'entity': 'pdb_viewer'}

        widgets['Show CCPN Names:'] = {'widget': 'choice',
                                       'entity': 'show_ccpn_names',
                                       'choices' : YES_NO_DICT}        

        return widgets

    def create(self, p):

        d = {}

        CW = self.create_widgets
        
        prefs = Tix.LabelFrame(p, label = 'Preferences')

        gui_defs = self.gui_defs()
        
        d['gui'] = CW(prefs.frame, gui_defs)
        d['cp'] = CW(prefs.frame, self.control_panel_defs(),
                     row = len(gui_defs))

        hosts = Tix.Frame(p)

        s = self.get_settings()

        self.jm_panel = JobManagerPanel()
        self.jm_panel.set_owner(self.get_owner())
        self.jm_panel.set_settings(s['gui'], s['gui'])
        self.jm_panel.create(hosts, title = 'Default host list')
        
        prefs.form(top = 0, left = 0, right = -1)
        hosts.form(top = prefs, left = 0, right = -1)

        return d

    def set_settings(self, gui_settings, control_panel_settings):

        self.__settings = {'gui': gui_settings,
                           'cp': control_panel_settings}

    def get_settings(self):
        return self.__settings

class ControlPanelSettings(Settings):

    def __init__(self, default_settings = None):

        from aria.Settings import Path, String, YesNoChoice

        ccpn_names_descr = 'Show detailed names for CCPN Data (this can make the display of Data panels slower)'

        d = {'text_editor': String(),
             'pdb_viewer': String(),
             'bitmap_path': Path(),
             'show_ccpn_names' : YesNoChoice(description = ccpn_names_descr)}

        Settings.__init__(self, d, default_settings)

    def create_default_values(self):

        d = {'text_editor': 'emacs',
             'pdb_viewer': 'rasmol',
             'show_ccpn_names' : YES}

        return d

class ControlPanel(AriaBaseClass):

    messages = {}

    messages['button_save'] = \
"""                            
If pressed, values in all edit fields are first checked for validity and then, if valid, stored. Note: Switching to another panel (by clicking another node in the tree window), has the same effect.
"""

    messages['button_undo'] = \
"""
All entry fields are set to values that are currently stored.
"""

    messages['button_reset'] = \
"""
Resets the focused entry field to its default value. If you also want to store that value, you either have to press the 'Commit' button or switch to another panel. Note: For some entry fields it does not make sense to define a default value.
"""

    messages['button_reset_panel'] = \
"""
All entry fields on the panel are reset to their default values. If you also want to store these values, you either have to press the 'Commit' button or switch to another panel. Note: For some entry fields it does not make sense to define a default value.
"""

    messages['button_delete'] = \
"""
Sometimes one might want to delete the current panel (e.g. for removing <Iterations>, or <Spectra> etc). 
"""

    from aria.tools import wrap_string

    for name, msg in messages.items():
        descr = wrap_string(msg, AriaBaseClass.description_length)
        messages[name] = descr

    del wrap_string

    def __init__(self, settings, parent, tooltips):

        check_type(settings, 'ControlPanelSettings')

        AriaBaseClass.__init__(self)
        self.setSettings(settings)
        
        self.tooltips = tooltips
        tooltips.config(bg='black')

        self.frames = {}
        self.current = None
        self.panels = {}
        self.frame = self.create(parent)
        self.set_modified(0)

        ## to assign a unique name to those panels
        ## with None as name

        self.id = 0

    def bind_tooltips(self, widget, msg):
        if msg not in ('', ' ', None):
            self.tooltips.bind_widget(widget, msg = msg)

    def create(self, parent):

        from widgets import DefaultButton as Button

        scrolled = Tix.ScrolledWindow(parent, scrollbar = 'auto')
        controls = Tix.Frame(scrolled.window)
        button_frame = Tix.Frame(parent)
        
        save = Button(button_frame, text = 'Commit', command =
                      lambda s = self: s.save())
        self.bind_tooltips(save, self.messages['button_save'])
        
        undo = Button(button_frame, text = 'Undo', command =
                      lambda s = self: s.undo())
        self.bind_tooltips(undo, self.messages['button_undo'])

        reset = Button(button_frame, text = 'Reset', command =
                       lambda s = self: s.reset_widget())
        self.bind_tooltips(reset, self.messages['button_reset'])

        reset_panel = Button(button_frame, text = 'Reset all',
                             command = lambda s = self: s.reset_panel())
        self.bind_tooltips(reset_panel, self.messages['button_reset_panel'])
        
        delete = Button(button_frame, text = 'Delete panel', command =
                        lambda s = self: s.delete())
        self.bind_tooltips(delete, self.messages['button_delete'])

        buttons = (save, undo, reset, reset_panel, delete)

        [b.pack(side = Tix.LEFT) for b in buttons]

        controls.pack(fill = Tix.BOTH, padx = 10)
        button_frame.pack(side = Tix.BOTTOM, pady = 10)

        scrolled.pack(expand = 1, fill = Tix.BOTH)

        self.button_save = save
        self.button_delete = delete

        return controls

    def add_panel(self, panel, name = None):
        panel.set_owner(self)

        if name is None:
            name = str(self.id)
            self.id += 1

        self.panels[name] = panel

    def find(self, name):
        if name in self.panels:
            return self.panels[name]
        else:
            return None

    def save(self):
        if self.current is not None:
            self.current.save()
            
        self.set_modified()

    def reset_widget(self):
        if self.current is not None:
            self.current.reset_widget()

    def reset_panel(self):
        if self.current is not None:
            self.current.reset()

    def undo(self):
        if self.current is not None:
            self.current.undo()

    def go(self):

        from aria.Project import ProjectThread

        if not hasattr(self, 'project'):
            return

        project = self.project

        project.getStructureEngine().getSettings()['prepare'] = 0
        project.setup()
        project.load_and_preprocess_data()

        thread = ProjectThread(project)
        thread.start()

    def delete(self):

        import widgets

        if self.current is None:
            return

        node = self.current.get_node()
        parent = node.getParent()
        
        if parent is not None:
            node_name = '%s.%s' % (parent.getName(), node.getName())
        else:
            node_name = node.getName()

        text = 'Do you want to delete\n' + \
               'node "%s"? [No Undo!]' % node_name

        dialog = widgets.YesNoDialog(self.frame, 'Warning',
                                     text = text, modal = 1)
        if dialog.get() <> 'yes':
            return

        current_panel = self.current
        self.current = None
        current_panel._delete(self.project)
        self.set_modified()

        ## remove delete panel from panel-list

        found = 0

        for k, v in self.panels.items():
            if v == current_panel:
                found = 1
                break

        if found:
            del self.panels[k]

    def set_panel(self, panel, frame):

        if panel <> self.current:
            
            self.frames[panel] = frame
            
            if self.current is not None:
                self.frames[self.current].forget()
            
            frame.pack(side = Tix.LEFT)
            self.current = panel

    def destroy_panels(self):

        map(lambda p: p.destroy(), self.frames.values())
        
        self.panels = {}
        self.frames = {}
        self.current = None

    def set_modified(self, m = 1):
        self.__modified = m

    def modified(self):
        if self.current is not None:
            ok = self.current.validate()
            if ok:
                self.current.save()
            self.__modified = self.__modified or not ok
            
        return self.__modified

    def savePanels(self):
        """
        Saves all panels
        """

        panels = self.panels.values()

        m = filter(lambda x: x.modified(), panels)

        for panel in panels:
            if panel.isInitialized() and panel.validate():
                panel.save()

        return m

    def get_modified_panels(self):

        ## check current panel since it might have been
        ## changed
        
        if self.current is not None:
            ok = self.current.validate()
            if ok:
                self.current.save()
                
        m = filter(lambda x: x.modified(), self.panels.values())

        return m

    def view(self, filename, data_type):
        import os

        s = self.getSettings()

        if data_type == DATA_TEXT:
            cmd = s['text_editor']
        elif data_type == DATA_PDB:
            cmd = s['pdb_viewer']

        os.system('%s %s &' % (cmd,filename))

    def change_panel(self, new_panel):

        allow_change = 1
        
        if self.current is None or new_panel is self.current:
            return allow_change

        self.current.validate()

##         if not self.current.save():

##             from widgets import YesNoDialog as Dialog

##             text = \
## """
## Some entry fields are missing or invalid.
## Switch to other panel anyway?
## """
##             dialog = Dialog(title = 'Error', text = text, modal = 1)
            
##             choices = {'yes': 1, 'no': 0}
##             allow_change = choices[dialog.get()]

##         if not allow_change:
##             node = self.current.get_node()
##             if node is not None:
##                 node.select()

        return allow_change

############################
######  BARDIAUX ############
    
## import threading
## class LoadThread(threading.Thread):

##     def __init__(self, filename):
##         threading.Thread.__init__(self)
##         self.file = filename
##         self.out = None

##     def run(self):
##         from tools import Load
##         self.out = Load(self.file)

class PlotPanel(AriaBaseClass):


    from aria.tools import wrap_string


    del wrap_string

    def __init__(self, settings, parent, tooltips):

        check_type(settings, 'ControlPanelSettings')

        AriaBaseClass.__init__(self)
        self.setSettings(settings)
        
        self.tooltips = tooltips
        tooltips.config(bg='black')

        self.frames = {}
        self.current = None
        self.panels = {}
        self.frame = self.create(parent)
        self.set_modified(0)

        ## to assign a unique name to those panels
        ## with None as name

        self.id = 1

    def bind_tooltips(self, widget, msg):
        if msg not in ('', ' ', None):
            self.tooltips.bind_widget(widget, msg = msg)

    def create(self, parent):

        from widgets import DefaultButton as Button

        scrolled = Tix.ScrolledWindow(parent, scrollbar = 'auto')
        controls = Tix.Frame(scrolled.window)
        button_frame = Tix.Frame(parent)

        close = Button(button_frame, text = 'Dismiss', command =
                       lambda p = parent: parent.destroy())
        self.bind_tooltips(close, "Close window")

        save  = Button(button_frame, text = 'Save as Text', command =
                      lambda s = self: s.save_as_text())
        
##         undo = Button(button_frame, text = 'Undo', command =
##                       lambda s = self: s.undo())
##         self.bind_tooltips(undo, self.messages['button_undo'])

##         reset = Button(button_frame, text = 'Reset', command =
##                        lambda s = self: s.reset_widget())
##         self.bind_tooltips(reset, self.messages['button_reset'])

##         reset_panel = Button(button_frame, text = 'Reset all',
##                              command = lambda s = self: s.reset_panel())
##         self.bind_tooltips(reset_panel, self.messages['button_reset_panel'])
        
        delete = Button(button_frame, text = 'Delete panel', command =
                        lambda s = self: s.delete())
##         self.bind_tooltips(delete, self.messages['button_delete'])

##         buttons = (save, undo, reset, reset_panel, delete)
        buttons = (close,save,)
        [b.pack(side = Tix.LEFT) for b in buttons]

        controls.pack(fill = Tix.BOTH, padx = 10)
        button_frame.pack(side = Tix.BOTTOM, pady = 10)

        scrolled.pack(expand = 1, fill = Tix.BOTH)

##         self.button_save = save
        self.button_delete = delete

        return controls

    def add_panel(self, panel, name = None):
        panel.set_owner(self)

        if name is None:
            name = str(self.id)
            self.id += 1

        self.panels[name] = panel

    def find(self, name):
        if name in self.panels:
            return self.panels[name]
        else:
            return None

    def save(self):
        if self.current is not None:
            self.current.save()
            
        self.set_modified()

    def reset_widget(self):
        if self.current is not None:
            self.current.reset_widget()

    def reset_panel(self):
        if self.current is not None:
            self.current.reset()

    def undo(self):
        if self.current is not None:
            self.current.undo()

    def go(self):

        from aria.Project import ProjectThread

        if not hasattr(self, 'project'):
            return

        project = self.project

        project.getStructureEngine().getSettings()['prepare'] = 0
        project.setup()
        project.load_and_preprocess_data()

        thread = ProjectThread(project)
        thread.start()

    def delete(self):

        import widgets

        if self.current is None:
            return

        node = self.current.get_node()
        parent = node.getParent()
        
        if parent is not None:
            node_name = '%s.%s' % (parent.getName(), node.getName())
        else:
            node_name = node.getName()

        text = 'Do you want to delete\n' + \
               'node "%s"? [No Undo!]' % node_name

        dialog = widgets.YesNoDialog(self.frame, 'Warning',
                                     text = text, modal = 1)
        if dialog.get() <> 'yes':
            return

        current_panel = self.current
        self.current = None
        current_panel._delete(self.project)
        self.set_modified()

        ## remove delete panel from panel-list

        found = 0

        for k, v in self.panels.items():
            if v == current_panel:
                found = 1
                break

        if found:
            del self.panels[k]

    def set_panel(self, panel, frame):

        if panel <> self.current:
            
            self.frames[panel] = frame
            
            if self.current is not None:
                self.frames[self.current].forget()
            
            frame.pack(side = Tix.LEFT)
            self.current = panel

    def destroy_panels(self):

        map(lambda p: p.destroy(), self.frames.values())
        
        self.panels = {}
        self.frames = {}
        self.current = None

    def set_modified(self, m = 1):
        self.__modified = m

    def modified(self):
        if self.current is not None:
            ok = self.current.validate()
            if ok:
                self.current.save()
            self.__modified = self.__modified or not ok
            
        return self.__modified

    def savePanels(self):
        """
        Saves all panels
        """

        panels = self.panels.values()

        m = filter(lambda x: x.modified(), panels)

        for panel in panels:
            if panel.isInitialized() and panel.validate():
                panel.save()

        return m

    def get_modified_panels(self):

        ## check current panel since it might have been
        ## changed
        
        if self.current is not None:
            ok = self.current.validate()
            if ok:
                self.current.save()
                
        m = filter(lambda x: x.modified(), self.panels.values())

        return m

    def view(self, filename, data_type):
        import os

        s = self.getSettings()

        if data_type == DATA_TEXT:
            cmd = s['text_editor']
        elif data_type == DATA_PDB:
            cmd = s['pdb_viewer']

        os.system('%s %s &' % (cmd,filename))

    def change_panel(self, new_panel):

        allow_change = 1
        
        if self.current is None or new_panel is self.current:
            return allow_change

        self.current.validate()

##         if not self.current.save():

##             from widgets import YesNoDialog as Dialog

##             text = \
## """
## Some entry fields are missing or invalid.
## Switch to other panel anyway?
## """
##             dialog = Dialog(title = 'Error', text = text, modal = 1)
            
##             choices = {'yes': 1, 'no': 0}
##             allow_change = choices[dialog.get()]

##         if not allow_change:
##             node = self.current.get_node()
##             if node is not None:
##                 node.select()

        return allow_change
    
    ## bardiaux
    def get_parent(self):
        return self.parent

    def save_as_text(self):

        c = self.find('contribs')
        c.save_as_text()

class ContactMapPanel(PanelEx):

    L_offset = 30.
    R_offset = 30.
    MAX_SIZE = 600

    

    def delete(self, project):
        s = self.__settings['cmap']
        #project.getProtocol().getSettings().delIterationSettings(s)

    def getSettings(self):
        try:
            return self.__settings['cmap']
        except:
            return None
    

    def test_defs(self):

        from aria.OrderedDict import OrderedDict

        widgets = OrderedDict()

        widgets['x:'] = {'widget': 'entry',
                           'entity': 'x'}

        widgets['y:'] = {'widget': 'entry',
                           'entity': 'y'}
        

        return widgets
    
    def scale_it(self, x, s):

        return x * (((self.MAX_SIZE - self.L_offset - self.R_offset)) / float(s)) + self.L_offset

        
    def create(self, p):
        
        self.parent = p

        LF = Tix.LabelFrame

        self.selected = []
        self.plot = []

       
        can = Tix.Canvas(p, height = self.MAX_SIZE, width = self.MAX_SIZE,bg='white', relief = Tix.RIDGE )
        self.canvas = can

        widgets = {}
    

        self.canvas.form(top = 0, left = 0, right = -1)

        
        from widgets import DefaultButton as Button

        button_frame = Tix.Frame(p)
        
        print_button = Button(button_frame, text = 'Save as Postscript', command = self._print )    
        print_button.grid()

        button_frame.form(top = self.canvas, left = 0, right = -1)
        
        return None

    def _print(self):

        # Print the Canvas as postcript
        uplevel =  self.get_owner() 
        cmap_settings = self.get_settings()
        selection = uplevel.find('cmap_opt').get_settings()['display']['display']

        base = "peakmap_it%s_%s.ps" % (cmap_settings['it'], selection)
            
        import tkFileDialog
        
        file_types = [("Postscript", "*.ps"),
                      ("Encapsulated Postscript", "*.eps"),
                      ("All files", "*")]
        
        dialog = tkFileDialog.SaveAs(master = self.parent, filetypes=file_types)

        filename = dialog.show(initialfile=base)

        if not filename:
            return
        
        data = self.canvas.postscript()
        
        if not data:
            return

        # patch the Title
        data = data.split('\n')
        data[3] = "%%%%Title: %s" % base
        f = open(filename, 'w')
        f.write('\n'.join(data))
        f.close()
        
        self.message('%s saved.' % filename)

    
    def plot_axes(self):

        self.axes_widgets = []

        if self.get_settings()['posit']:


            ###################################
            #  axes
            axis_max = int(self.size) / 5 + 1
            for step in range(0,axis_max):

                xx = step * 5

                val = self.scale_it(xx,self.size)
                #->
                self.canvas.create_text(val, self.MAX_SIZE - (self.L_offset - 15) ,text = str(xx+self.lower), tag = "axe")
                self.canvas.create_line(val, self.MAX_SIZE - (self.L_offset -5) ,val, self.MAX_SIZE - (self.L_offset - 10), tag = "axe")
                # |
                # V
                self.canvas.create_text(self.L_offset - 15, self.MAX_SIZE - val,text = str(xx+self.lower), tag = "axe")
                self.canvas.create_line(self.L_offset - 5, self.MAX_SIZE - val,self.L_offset - 10, self.MAX_SIZE - val, tag = "axe")


            self.canvas.create_rectangle(self.L_offset - 5, self.R_offset - 5, self.MAX_SIZE - (self.L_offset - 5), self.MAX_SIZE - (self.R_offset -5), tag = "axe")

            self.axes_widgets = self.canvas.find_withtag("axe")

            self.plot += self.axes_widgets
            
        
    def update_plot(self):
        
        nb = 0
        peaks = self.get_settings()['posit']

        all = self.canvas.find_all()
        [self.canvas.delete(i) for i in all if i in self.plot]

        self.plot = []
        coord = {}
        self.wid_pk = {}


        if not peaks:
            
            x = self.create_text_item(self.canvas,
                                          [self.MAX_SIZE/2.,
                                           self.MAX_SIZE/2.], "black",
                                          "No Peaks",str(nb), bind=0)
            coord.setdefault(x,[self.MAX_SIZE/2.,self.MAX_SIZE/2.])
            self.wid_pk.setdefault(x,nb)
            self.plot.append(x)
            nb+=1

        else:          

            self.lower = min(min([pk[0] for pk in peaks]),min([pk[1] for pk in peaks])) -1
            self.size = max(max([pk[0] for pk in peaks]),max([pk[1] for pk in peaks]))

            self.size = self.size - self.lower
            
            peaks = [ [self.scale_it(pk[0]-self.lower, self.size),self.MAX_SIZE - self.scale_it(pk[1]-self.lower, self.size)] for pk in peaks ]


            for i in peaks:
            
                x = self.create_text_item(self.canvas, i, "red", "+","point")

                coord.setdefault(x,i)
                self.wid_pk.setdefault(x,nb)
                self.plot.append(x)
                nb+=1

        settings = self.get_settings()
        settings['xy'] = coord
        self.set_settings(settings)

        return None

    def set_settings(self, cmap_settings):

        self.__settings = {}

        self.__settings = cmap_settings

    def get_settings(self):
        return self.__settings

##     def show(self):

##         PanelEx.show(self)

##         ## enable control-panels delete button

## ##         if self.get_owner() is not None:
## ##             self.get_owner().button_delete.config(state = Tix.NORMAL)

    def show(self):

        file = self.get_settings()['file']
        self.selection = self.get_settings()['selection']
        
        from aria.tools import Load
        from AriaViewer import get_cmap_list

        try:
            x = self.pickle
            
        except:
            
            self.message('Loading restraint data. This might take a while ...')

            self.pickle = Load(file)
            
##         if len(self.get_settings()['posit']) < 1:
            
            
## ##             popup = Tix.Tk()
## ##             popup.wm_title('Loading iteration data')
            
## ##             m = Meter(popup,relief='ridge', bd=3)
## ##             m.pack(fill='x')
## ##             m.set(0.0, 'Loading...')

## ##             loading = LoadThread(file)
## ##             loading.start()
## ##             value = 0.
## ##             while loading.isAlive():
## ##                 value += 0.00005
## ##                 m.set(value)                
## ##             m.set(1., 'Pickled file loaded')
## ##             popup.destroy()
## ##             self.pickle = loading.out

##             print 'Loading restraint data. This might take a while ...'

##             self.pickle = Load(file)
            
            
        selection = self.get_owner().find('cmap_opt').get_settings()['display']['display']

        if self.selection != selection:
        
            xx, pk, ct = get_cmap_list(self.pickle,selection)
             
            cmap_settings = self.get_settings()
            cmap_settings['posit'] = xx
            cmap_settings['plist'] = pk
            cmap_settings['clist'] = ct

            self.set_settings(cmap_settings)
            cmap_settings['selection'] = selection
           

        PanelEx.show(self)

        self.update_plot()            

        self.plot_axes()


        if self.get_settings()['posit']:
            self.showPopUp()

    def showPopUp(self):

        uplevel =  self.get_owner() 

        cmap_settings = self.get_settings()
        selection = uplevel.find('cmap_opt').get_settings()['display']['display']
    
        
        try:
            uplevel.pk_win.wm_title('Peaks and Contributions lists for iteration ' + str(cmap_settings['it']) + ' - ' + selection)
            
        except:

            uplevel.pk_win = Tix.Toplevel(self.parent)
            uplevel.pk_win.wm_title('Peaks and Contributions lists for iteration ' + str(cmap_settings['it']) + ' - ' + selection)
            z = uplevel.pk_win.winfo_toplevel()
            z.wm_protocol("WM_DELETE_WINDOW", lambda p= uplevel.pk_win : p.destroy())

            cp_settings = ControlPanelSettings()
            balloon = Tix.Balloon(uplevel.pk_win)

            uplevel.win = PlotPanel(cp_settings, uplevel.pk_win, balloon)

            contrib_panel = ContribPanel()

            uplevel.win.add_panel(contrib_panel, name = 'contribs')

            contribs = uplevel.win.find('contribs')
            from AriaViewer import ContribListSettings, PeakListSettings
            contribs.set_settings(ContribListSettings(),PeakListSettings())
            contribs.show()

            
        

    ########## for CANVAS #########


    def getCoor(self,event):
        match  = event.widget
        no = match.find_closest(event.x, event.y)[0]
        pos = self.get_settings()['xy'][no]

        return pos[0],pos[1],no

    
    def circAddPeak(self,event):


        x , y, npk = self.getCoor(event)
        if npk not in self.selected:
            match  = event.widget
            match.create_oval(x-4,y-4,x+4,y+4, tag =str(x)+"-"+str(y))
            self.selected.append(npk)
 
        self.showSelected()

    def circDel(self, event):

        x, y, npk = self.getCoor(event)
        match = event.widget
        no = match.find_withtag(str(x)+"-"+str(y))
        try:
            match.delete(no[0])
        except:
            pass

        self.selected.remove(npk)
        self.showSelected()


    def circPeak(self,event):

        #SC = self.scale
        match  = event.widget

        x, y, npk= self.getCoor(event)
        
        hit = match.find_all()
        [match.delete(i) for i in hit if \
         i not in self.plot and \
         i not in self.axes_widgets]

        
        ## TODO : hack w/ getter/setter
        self.selected = []
        
        ov = match.create_oval(x-4,y-4,x+4,y+4, tag = str(x)+"-"+str(y))
 
        ## TODO : hack w/ getter/setter
        self.selected.append(npk)
        self.showSelected()

    def showSelected(self):
        
        from AriaViewer import setPeak, setContrib
        from AriaViewer import ContribListSettings, PeakListSettings

        self.showPopUp()
        uplevel = self.get_owner()
        contribs = uplevel.win.find('contribs')
        
        p_list = []
        c_list = []
        for n_wid in self.selected:
            npk = self.wid_pk[n_wid]
            p_list = p_list + self.get_settings()['plist'][npk]
            c_list = c_list + self.get_settings()['clist'][npk]
            
        #p_list = [setPeak(pk) for pk in p_list]
        pp = PeakListSettings()
        pp['pk_list'] = p_list

        #c_list = [setContrib(c) for c in c_list]        
        cc = ContribListSettings()
        cc['c_list'] = c_list
        
        contribs.set_settings(cc,pp)
        contribs.update_tab()
#        self.pk_win.focus_force()
#        self.pk_win.transient()
#        self.pk_win.grab_set()
        uplevel.pk_win.tkraise()


################ ContribPanle
        
class ContribPanel(PanelEx):

    def save_as_text(self):
        
        import tkFileDialog
        
        file_types = [("Text", "*.dat"),
                      ("All files", "*")]
        
        dialog = tkFileDialog.SaveAs(master = self.parent, filetypes=file_types)
        
        filename = dialog.show()

        if not filename:
            return

        out = open(filename, 'w')
        
        kw = ('id', 'ref_peak', 'spec', 'dist', 'low','up', 'weight', 'deg_viol','dist_avg','state','viol','type')

        out.write("# Restraints\n")
        out.write("#%s\n\n" % " ".join(["%7s" % x for x in kw]))

        l = []
        for p in self.get_settings()['peaks']['pk_list']:
                s = ["%7s" % p[x] for x in kw]
                l.append(" " + " ".join(s))

        out.write( "\n".join(l) )

        kw = ('ariapeak','id', 'dist', 'weight', 'res1', 'at1','seg1', 'res2', 'at2','seg2')

        out.write("\n\n# Contributions\n")
        out.write("#%s\n\n" % " ".join(["%7s" % x for x in kw]))
        
        l = []
        for p in self.get_settings()['contrib']['c_list']:
                s = ["%7s" % p[x] for x in kw]
                l.append(" " + " ".join(s))
                
        out.write( "\n".join(l) )                
        out.close()

        self.message('%s saved.' % filename)
        
    def validate(self):
        """
        Checks whether all cells are valid.
        """
        for row in self.table.rows:
            for r in row:
                is_valid = r.isValid()
                if is_valid == 0:
                    return self.set_modified()

        ## all cells are valid
                
        return self.set_modified(0)

    def not_saved(self, widget):
        from tkMessageBox import showerror

        if widget is not None:
            try:
                widget.save()
            except Exception, msg:
                showerror('Value error',msg)
        
        self.set_modified()
    
    def more_rows_a(self, n = 1):

        l = len(self.table_a.rows)
        indices = range(l, l + n)
        
        self.table_a.append_row(n)
        self.init_cells_a(indices)
        
    def more_rows_c(self, n = 1):

        l = len(self.table_c.rows)
        indices = range(l, l + n)
        
        self.table_c.append_row(n)
        self.init_cells_c(indices)

    def del_row(self):
        current_row = self.table.get_current_row_index()
        if current_row is None:
            return
        self.table.remove_row(current_row)

    def init_cells_a(self, indices):

        from AriaViewer import PeakSettings

##         gui_settings = self.get_settings()['gui']
##         default_cns_executable = gui_settings['cns_executable']
        
        for i in indices:

            ps = PeakSettings()

            ## 0th column: enabled

            e = ps.getEntity('id')
            self.table_a.set_widget(i, 0, DEntry, e)
            
            e = ps.getEntity('ref_peak')
            self.table_a.set_widget(i, 1, DEntry, e)
            
            e = ps.getEntity('spec')
            self.table_a.set_widget(i, 2, DEntry, e)
            
            e = ps.getEntity('dist')
            self.table_a.set_widget(i, 3, DEntry, e)

            e = ps.getEntity('low')
            self.table_a.set_widget(i, 4, DEntry, e)

            e = ps.getEntity('up')
            self.table_a.set_widget(i, 5, DEntry, e) 
            
            e = ps.getEntity('weight')
            self.table_a.set_widget(i, 6, DEntry, e)

            e = ps.getEntity('deg_viol')
            self.table_a.set_widget(i, 7, DEntry, e)

            e = ps.getEntity('dist_avg')
            self.table_a.set_widget(i, 8, DEntry, e)
            
            e = ps.getEntity('state')
            self.table_a.set_widget(i, 9, DEntry, e)

            e = ps.getEntity('viol')
            self.table_a.set_widget(i, 10, DEntry, e)            

            e = ps.getEntity('type')
            self.table_a.set_widget(i, 11, DEntry, e)

            
    def init_cells_c(self, indices):

        from AriaViewer import ContribSettings

##         gui_settings = self.get_settings()['gui']
##         default_cns_executable = gui_settings['cns_executable']
        
        for i in indices:

            cs = ContribSettings()

            ## 0th column: enabled
            e = cs.getEntity('ariapeak')
            self.table_c.set_widget(i, 0, DEntry, e)
            
            e = cs.getEntity('id')
            self.table_c.set_widget(i, 1, DEntry, e)
            
            e = cs.getEntity('dist')
            self.table_c.set_widget(i, 2, DEntry, e)
            
            e = cs.getEntity('weight')
            self.table_c.set_widget(i, 3, DEntry, e)
            
            e = cs.getEntity('res1')
            self.table_c.set_widget(i, 4, DEntry, e)

            e = cs.getEntity('at1')
            self.table_c.set_widget(i, 5, DEntry, e)

            e = cs.getEntity('seg1')
            self.table_c.set_widget(i, 6, DEntry, e) 
            
            e = cs.getEntity('res2')
            self.table_c.set_widget(i, 7, DEntry, e)

            e = cs.getEntity('at2')
            self.table_c.set_widget(i, 8, DEntry, e)

            e = cs.getEntity('seg2')
            self.table_c.set_widget(i, 9, DEntry, e)
            

    def create(self, p, title = 'Aria peaks'):

        import aria.tools as tools

        self.parent = p
        #######" Peaks
        frame1 = Tix.LabelFrame(p, label = 'Restraints')

        header = ('Id','Ref peak', 'Spectrum','Dist','Lower','Upper','Weight','Violation','Avg dist','State','Violated','Type')

        
        c_widths = (6,6,12,6,6,6,6,6,6,6,6,10)
        
        self.table_a = DTable(frame1.frame, (0, len(c_widths)), header = header,
                            column_widths = c_widths, relief = Tix.FLAT)

        table_size = 0
        self.more_rows_a(table_size)
        
        ## fill table with host-list data
        
        kw = ('id', 'ref_peak', 'spec', 'dist', 'low','up', 'weight', 'deg_viol','dist_avg','state','viol','type')
        
        pk = self.get_settings()['peaks']['pk_list']
        n_rows = len(pk)
        
        self.more_rows_a( max(n_rows - table_size, 0))
       
        for i in range(n_rows):
            for j in range(len(kw)):
                self.table_a[i,j] = pk[i][kw[j]]

        self.table_a.grid(row = 0, column = 1, columnspan = 2)

        
        ################ Contrib
        frame2 = Tix.LabelFrame(p, label = 'Contributions')

        header = ('ariapeak','id', 'dist','weight','res 1','at 1', 'seg 1', 'res 2', 'at 2', 'seg 2')
        
        
        c_widths = (6,6,6,6,6,12,6,6,12,6)
        
        self.table_c = DTable(frame2.frame, (0, len(c_widths)), header = header,
                            column_widths = c_widths, relief = Tix.FLAT)

        table_size = 0
        self.more_rows_c(table_size)
        
        kw = ('ariapeak','id', 'dist', 'weight', 'res1', 'at1','seg1', 'res2', 'at2','seg2')
        

        contribs = self.get_settings()['contrib']['c_list'] 
        n_rows = len(contribs)
        
        self.more_rows_c(max(n_rows - table_size, 0))
       
        for i in range(n_rows):
            for j in range(len(kw)):
                self.table_c[i,j] = contribs[i][kw[j]]

        self.table_c.grid(row = 0, column = 1, columnspan = 2)

        
        frame1.form(left = 0, top = 0)
        frame2.form(left = 0, top = frame1)

        return None

    def set_settings(self, cmap_set, peak_set):
        
        self.__settings = {'contrib': cmap_set,
                           'peaks' : peak_set}
        
        
    def get_settings(self):
        return self.__settings

    def update_tab(self):

        for i in range(0,len(self.table_a.rows)):
            self.table_a.remove_row(0)

        for i in range(0,len(self.table_c.rows)):
            self.table_c.remove_row(0)
        
        kw = ('id', 'ref_peak', 'spec', 'dist', 'low','up', 'weight', 'deg_viol','dist_avg','state','viol','type')        
        pk = self.get_settings()['peaks']['pk_list']
        n_rows = len(pk)
        self.more_rows_a(n_rows)
       
        for i in range(n_rows):
            for j in range(len(kw)):
                self.table_a[i,j] = pk[i][kw[j]]
                # color row if violated
                if kw[j] == 'viol' and pk[i][kw[j]] == 'yes':
                    for cell in self.table_a.rows[i]:
                        cell.label.configure(bg ='#FFFF99')
        
        kw = ('ariapeak','id', 'dist', 'weight', 'res1', 'at1','seg1', 'res2', 'at2','seg2')        
        contribs = self.get_settings()['contrib']['c_list']
        n_rows = len(contribs)
        
        self.more_rows_c(n_rows)
       
        for i in range(n_rows):
            for j in range(len(kw)):
                self.table_c[i,j] = contribs[i][kw[j]]
                
##                 # color row if inter
##                 if kw[j] == 'type' and contribs[i][kw[j]] == 'INTER':
##                     for cell in self.table_c.rows[i]:
##                         cell.label.configure(bg ='#CCFFCC')



###### Cmap options panel
                        
class CmapOptionsPanel(DataPanel):

    def __init__(self, node = None):
        DataPanel.__init__(self, node, enable_delete_button = 0)
    
    
    def get_defs(self):
        from aria.OrderedDict import OrderedDict

        d = OrderedDict()

        DICT =[{'label': 'All', 'value': 'all'},{'label': 'Ambiguous', 'value': 'ambig'},{'label': 'Unambiguous', 'value': 'unambig'}]

        d['Show distance restraints:'] = {'widget': 'choice',
                                          'entity': 'display',
                                          'choices': DICT}

        return d

    def create(self, p):

        frame = Tix.Frame(p)

        text = '''\nThe tool "peak map analysis" provides a graphical representation
of the restraint lists used to calculate the structures.

In order to use this tool, the option "Python pickle output"
in node "Report" needs to be enabled.

The maps are organised on an per-iteration basis. Dots in
a graph (which looks similar to a contact map) indicate pairs
of residues that interact via a distance restraint.

By clicking the dots, a window pops up giving additional
information on the respective reference cross peaks as well as on
the contributions that make up the restraint.
'''        
        label = Tix.Label(p, text = text, justify=Tix.LEFT)
        label.grid(sticky=Tix.W, row=0, columnspan=10)

        widgets = self.create_widgets(p, self.get_defs(), row=1)
        return {'display': widgets}
    
    def set_settings(self, s):
        
        self.__settings = {'display': s}  

    def get_settings(self):
        return self.__settings

