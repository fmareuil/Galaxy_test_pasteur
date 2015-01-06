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


import Tix, widgets
from aria.ariabase import *
from aria.TypeChecking import check_type

class WidgetSaveError(Exception):
    def __str__(self):
        return Exception.__str__(self)

class WidgetDecorator:

    def __init__(self):
        self.set_variable(None)
        self.is_modified = 0

    def set_entity(self, e):
        check_type(e, 'Entity', NONE)
        self.__entity = e
        self.__prev_value = None

    def get_entity(self):
        return self.__entity

    def widget_get(self):
        """
        get widgets state
        """

        return self.get()

    def widget_set(self, that):
        """
        set widgets state
        """

        self.set(that)
        self.is_modified = 1

    def set_variable(self, v):
        self.__var = v

    def get_variable(self):
        return self.__var

    def is_valid(self):
        """
        Returns 1 if the entity would be able to
        store widgets current value
        """

        entity = self.get_entity()

        if entity is None:
            return 0

##         if not entity.is_initialized:
##             return 1

        value = self.widget_get()
        
        try:
            value = entity.cast(value)
        except Exception, msg:
            return 0

        return entity.is_valid(value)

    def validate(self):
        """
        Validates the editors value against
        the entity. if the value is not valid
        the entity raises an exception.
        """

        entity = self.get_entity()

        if entity is None:
            return 0

        if not entity.is_initialized():
            return 1

        value = self.widget_get()
        entity.validate(value)
        
        return 1

    def save(self):
        """
        attempts to store widgets value into entity
        """

        from aria.Settings import EntityCastError

        entity = self.get_entity()

        if entity is None:
            raise WidgetSaveError, 'Cannot save. Entity not set.'

        value = self.widget_get()
        value = entity.cast(value)
        entity.set(value)

        ## save current value
        self.__prev_value = value

    def undo(self):
        if self.__prev_value is not None:
            self.widget_set(self.__prev_value)

    def reset(self):
        entity = self.get_entity()
        self.__prev_value = entity.get()
        entity.reset()
        self.load()

    def load(self):
        """
        sets widgets value according to entity's value.
        """

        if self.get_entity() is None:
            raise ValueError, 'Cannot load. Entity not set.'

        try:
            val = self.get_entity().get()
        except:
            return

        self.widget_set(val)
        if self.__prev_value is None:
            self.__prev_value = val

    def highlight(self):
        pass

    def set_trace_callback(self, mode, callback):
        v = self.get_variable()
        
        if v is not None:
            self.__callback = callback
            v.trace(mode, self.callback)
        else:
            name = str(self.__class__.__name__)
            raise Exception, '%s: no variable set.' % name

    def callback(self, *args, **kw):
        args = [self] + list(args)
        self.__callback(*args, **kw)

    def has_changed(self):
        try:
            self.__prev_value
            has_attr = 1
        except:
            has_attr = 0

        if not has_attr:
            return 1

        try:
            curr_val = self.get_entity().cast(self.widget_get())
        except:
            return 1
        
        return self.__prev_value <> curr_val

    def get_tk_widget(self):
        return self

    def __str__(self):
        return self.__class__.__name__
    
## decorated widgets

class DEntry(Tix.Entry, WidgetDecorator):

    def __init__(self, *args, **kw):

        WidgetDecorator.__init__(self)
        
        if 'textvariable' not in kw:
            var = Tix.StringVar()
            kw['textvariable'] = var
        else:
            var = None

        self.set_variable(var)
        
        Tix.Entry.__init__(self, *args, **kw)

    def widget_set(self, v):
        self.delete(0, Tix.END)
        self.insert(0, str(v))

    def highlight(self):
        self.select_range(0, Tix.END)

class DText(Tix.Text, WidgetDecorator):

    def __init__(self, *args, **kw):

        WidgetDecorator.__init__(self)

        Tix.Text.__init__(self, *args, **kw)

        self.callback = None
        self.bind('<Key>', self.key)

    def key(self, *args, **kw):
        if self.callback is not None:
            self.callback(*args, **kw)

    def widget_get(self):
        return self.get(1.0, Tix.END)[:-1]

    def widget_set(self, v):
        self.delete(1.0, Tix.END)
        self.insert(1.0, str(v))

    def highlight(self):

        l = len(self.widget_get())
        
        self.tag_remove(Tix.SEL, 1.0, 1.0)
        self.tag_add(Tix.SEL, 1.0, l)
        self.tag_remove(Tix.SEL, l, Tix.END)

    def set_trace_callback(self, mode, callback):
        self.callback = callback

class DLabelEntry(Tix.LabelEntry, WidgetDecorator):
    def __init__(self, *args, **kw):

        Tix.LabelEntry.__init__(self, *args, **kw)
        WidgetDecorator.__init__(self)

        if 'textvariable' not in kw:
            var = Tix.StringVar()
            self.entry.configure(textvariable = var)
        else:
            var = None

        self.set_variable(var)

    def widget_get(self):
        return self.entry.get()

    def widget_set(self, v):
        self.entry.delete(0, Tix.END)
        self.entry.insert(0, str(v))

    def get_tk_widget(self):
        return self.entry

## class DFileEntry(Tix.FileEntry, WidgetDecorator):

##     def __init__(self, *args, **kw):
##         Tix.FileEntry.__init__(self, *args, **kw)
##         WidgetDecorator.__init__(self)

##         if 'textvariable' not in kw:
##             var = Tix.StringVar()
##             self.entry.configure(textvariable = var)
##         else:
##             var = None

##         self.set_variable(var)

##     def widget_get(self):
##         return self.entry.get()

##     def widget_set(self, v):
##         self.entry.delete(0, Tix.END)
##         self.entry.insert(0, str(v))

##     def highlight(self):
##         self.entry.select_range(0, Tix.END)

from Tix import StringVar

## TODO: for testing, remove later

class MyStringVar(StringVar):
    """Value holder for strings variables."""

    def get(self):
        """Return value of variable as string."""
        value = self._tk.globalgetvar(self._name)
        print 'StringVar:', self._name, 'get', value, type(value)
        return StringVar.get(self)

MyStringVar = StringVar

class DOptionMenu(Tix.OptionMenu, WidgetDecorator):

    def __init__(self, *args, **kw):

        WidgetDecorator.__init__(self)
        
        var = MyStringVar() 

        kw['variable'] = var

        self._translator = {'True': YES, 'False': NO}

        self.set_variable(var)
        Tix.OptionMenu.__init__(self, *args, **kw)

    def widget_get(self):
        name = self.get_variable().get()

        if name == 'yes': name = name.replace('', '')

        try:
            name = self._translator[name]
        except:
            pass

        return name

    def widget_set(self, name):
        if name == 'yes': name = '' + name

        self.get_variable().set(name)

class DIntOptionMenu(DOptionMenu):

    def widget_get(self):
        return int(self.get_variable().get())

class DCheckbutton(Tix.Checkbutton, WidgetDecorator):

    def __init__(self, *args, **kw):

        WidgetDecorator.__init__(self)
        
        var = Tix.StringVar()
        kw['variable'] = var
        
        self.set_variable(var)

        Tix.Checkbutton.__init__(self, *args, **kw)

    def widget_get(self):
        return self.get_variable().get()

    def widget_set(self, n):
        self.get_variable().set(n)

class DListBox(Tix.ScrolledListBox, WidgetDecorator):

    def __init__(self, *args, **kw):
        """
        argument for callbacks:

        <Button-1>: self, event, current_selection
        """
        
        Tix.ScrolledListBox.__init__(self, *args, **kw)
        WidgetDecorator.__init__(self)
        
        self.listbox.bind('<Button-1>', self._button1)

        self.__callbacks = {'<Button-1>': None}

    def insert(self, pos, item):
        self.listbox.insert(pos, item)

    def delete(self, start, end):
        self.listbox.delete(start, end)

    def widget_get(self):
        return list(self.listbox.get(0, Tix.END))

    def widget_set(self, l):
        lb = self.listbox
        lb.delete(0, Tix.END)
        map(lambda ll, lb = lb: lb.insert(Tix.END, ll), l)

    def set_callback(self, event_name, callback):
        self.__callbacks[event_name] = callback

    def _button1(self, event):
        ## TODO: does not work at the moment
        callback = self.__callbacks['<Button-1>']
        if callback is None:
            return
        
        sel = self.listbox.curselection()
        callback(self, sel)

class DComboBox(Tix.ComboBox, WidgetDecorator):

    def __init__(self, *args, **kw):

        var = Tix.StringVar()
        kw['variable'] = var

        Tix.ComboBox.__init__(self, *args, **kw)
        WidgetDecorator.__init__(self)

        self.listbox = self.slistbox.listbox
        self.set_variable(var)
        
    def delete(self, start, end):
        self.listbox.delete(start, end)

    def index(self, index):
        return self.listbox.index(index)

    def set_items(self, items):
        self.delete(0, Tix.END)
        [self.insert(Tix.END, i) for i in items]

    def widget_get(self):
        return self.get_variable().get()

    def widget_set(self, l):
        items = self.listbox.get(0, Tix.END)
        if l in items:
            self.get_variable().set(l)

class DBrowseEntry(widgets.BrowseEntry, WidgetDecorator):
    
    def __init__(self, *args, **kw):
        widgets.BrowseEntry.__init__(self, *args, **kw)
        WidgetDecorator.__init__(self)

    def widget_get(self):
        return self.entry.get()

    def widget_set(self, v):
        self.entry.delete(0, Tix.END)
        self.entry.insert(0, str(v))

    def highlight(self):
        self.entry.select_range(0, Tix.END)

    def get_tk_widget(self):
        return self.entry

class DPathEntry(DBrowseEntry):

    def browse(self):
        from tkFileDialog import Directory
        import os

        current_dir = self.widget_get()
        if current_dir is None or current_dir == '':
            current_dir = os.getcwd()

        return Directory(self, title = 'Select a directory',
                         initialdir = current_dir).show()

class CCPNProjectBrowser(DBrowseEntry):

    def browse(self):
        from tkFileDialog import Directory
        import os

        current_dir = self.widget_get()
        if current_dir is None or current_dir == '':
            current_dir = os.getcwd()

        try:
            from memops.editor.OpenProjectPopup import OpenProjectPopup
            return OpenProjectPopup(self, transient=True, modal=True).\
                   file_select.getDirectory()
        except:
            return Directory(self, title = 'Select a CCPN project (directory)',
                             initialdir = current_dir).show()
        
class DFileEntry(DBrowseEntry):

    def __init__(self, parent, file_types = None,
                 title = 'Select a file', text = None, *args, **kw):
        DBrowseEntry.__init__(self, parent, text = text)

        if 'textvariable' not in kw:
            var = Tix.StringVar()
            self.entry.configure(textvariable = var)
        else:
            var = None

        self.set_variable(var)

        self.__file_types = file_types
        self.__title = title

    def browse(self):
        import tkFileDialog, os

        current_file = self.widget_get()

        if current_file is not None:
            path, name = os.path.split(current_file)
        else:
            path = os.getcwd()
            name = None

        kw={}
        kw['title'] = self.__title

        if self.__file_types is not None:
            kw['filetypes'] = self.__file_types

        dialog = tkFileDialog.Open(**kw)
        val = dialog.show(initialdir = path, initialfile = name)

        if val == '':
            val = None
        
        return val

class DCell(widgets.Cell):

    def __init__(self, parent, row, column, width = 10):
        widgets.Cell.__init__(self, parent, row, column, width)
        self.setCallback(None)

    def setCallback(self, f):
        """
        Permits to install a callback which is invoked
        when the entity could not save widget's value.
        """
        
        self.__callback = f

    def getCallback(self):
        return self.__callback

    def bind_events(self):

        if is_type(self.editor, 'DEntry') or \
           is_type(self.editor, 'DOptionMenu'):
            widgets.Cell.bind_events(self)
            
        elif is_type(self.editor, 'DFileEntry') or \
             is_type(self.editor, 'DPathEntry'):
            self._bind(self.editor.entry)

    def create_editor(self, constructor, entity):
        editor = widgets.Cell.create_editor(self, constructor)
        editor.set_entity(entity)
        
        try:
            val = str(entity.get())
            self.label.configure(text = val)
            self._label_value = val
            self.editor_set(val)
        except:
            pass

        return editor

    def editor_get(self):
        if self.editor is None:
            val = None
        else:
            val = self.editor.widget_get()

        return val

    def editor_set(self, v):
        self.editor.widget_set(v)

    def save(self):

        import aria.Settings as Settings
        
        try:
            self.editor.save()
            saved = widgets.Cell.save(self)
            
        except Settings.EntityCastError, msg:
            entity = self.editor.get_entity()
            
            if entity is not None:
                err_msg = entity.getErrorMessage()
            else:
                err_msg = ''
                
            msg = ('Could not cast value: %s\nTo leave the ' + \
                   'cell without saving, press <esc>') % err_msg
            
            saved = 0

        except Settings.EntityValueError, msg:
            
            entity = self.editor.get_entity()

            if entity is not None:
                err_msg = entity.getErrorMessage()
            else:
                err_msg = ''
                
            msg = ('Could not save value: %s\nTo leave the ' + \
                   'cell without saving, press <esc>') % err_msg
            
            saved = 0

        if not saved:
            self.__invoke_callback(self.editor)

        return saved

    def __invoke_callback(self, editor):
        f = self.getCallback()
        if f is not None:
            f(editor)

    def get(self):
        if self.editor is None:
            val = None
        elif self._in_edit_mode:
            val = self.editor_get()
            if val is not None:
                self.editor.validate()
        else:
            val = self._label_value

        return val

    def set(self, v):
        if self.editor is None:
            return

        self.editor.get_entity().set(v)
        widgets.Cell.set(self, v)

    def isValid(self):
        if self.editor is None:
            return 0

        if self._in_edit_mode:
            return self.editor.is_valid()
        else:
            return 1

class DTable(widgets.Table):

    def create_cell(self, row, column):
        from aria.Settings import String
        
        width = self.column_widths[column]
        
        cell = DCell(self.window, row, column, width)

        ## default cell: entry-field, string entity

        entity = String()
        cell.create_editor(DEntry, entity)

        return cell

    def set_widget(self, i, j, constructor, entity):
        check_type(entity, 'Entity')

        cell = self.rows[i][j]

        widget = cell.create_editor(constructor, entity)

if __name__ == '__main__':


    from aria.Settings import *

    root = Tix.Tk()
    table = DTable(root, (10,3), column_widths = (35, 15, 30),
                   header = ('A', 'Command', 'CNS executable'))

    entity = Path(exists=1)
    i = Integer()
#    i.set(4)
    table.set_widget(0,0,DFileEntry, entity)
    table.set_widget(1,0,DFileEntry, entity)
    table.set_widget(0,1,DEntry,i)
    table.pack()
