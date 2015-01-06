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


import Tkinter, Tix
from aria.ariabase import YES, NO
BUTTON_WIDTH = 8

class MyScrolledText(Tix.ScrolledText):

    def write(self, m):
        self.text.insert(Tix.END, m)
	self.text.see(Tix.END)

    def flush(self):
        pass

            
class MyPanedWindow(Tix.PanedWindow):
    """
    Some extensions that are missing in the Python
    Tix wrappers
    """

    def panecget(self, pane_name, key):
        args = (self._w, 'panecget', pane_name, '-'+key)
        return self.tk.call(args)

    def setsize(self, pane_name, size, direction = ''):
        args = (self._w, 'setsize', pane_name, str(size),
                direction)
        self.tk.call(args)

class Menu(Tix.Menu):

    ## TODO: does not work at the moment.

    def add_command(self, *args, **kw):
        Tix.Menu.add_command(self, *args, **kw)
        if 'accelerator' in kw:
            self.tk.bind('<%s>' % kw['accelerator'], kw['command'])
            print 'bound'

class Tree(Tix.Tree):

    def __init__(self, parent, bitmap_path = None):

        Tix.Tree.__init__(self, parent, options = 'separator "/"')
        
        self.__nodes = {}
        
        self['ignoreinvoke'] = '1'
        self['browsecmd'] = self.__cmd_browse
        self['opencmd'] = self.__cmd_open
        self['command'] = self.__cmd_command
        self['closecmd'] = self.__cmd_close

        self.hlist.configure(selectmode = 'single', bg='white')
        
        # BARDAIUX 2.2
        # Tcl8.4 fix
        if Tix.TkVersion >= 8.4:
            self.hlist.configure(selectbackground = '#000080',\
                                 selectforeground = 'white')

        self.bitmap_path = bitmap_path

    def _add(self, node, node_id, *args, **kw):
        node.setId(node_id)
        node.setTree(self)
        node._setParent(self)

        if self.bitmap_path is None:
            name = 'folder'
        else:
            import os
            name = os.path.join(self.bitmap_path, 'shaded_folder')

        self.hlist.add(node_id, itemtype=Tix.IMAGETEXT, text = node.getName())

        self.__nodes[node_id] = node

        node.setIcon(name)
        
        if 'mode' in kw:
            self.setmode(node_id, kw[mode])

    def _del(self, node):
        id = node.getId()
        self.hlist.delete_entry(id)

        ## remove node from parent's list of children

        p = node.getParent()
        
        if p is not None:
            c = p.get_children()
            del c[id]
            if not c:
                p.getTree().setmode(p.getId(), 'none')
            

        ## remove node from all-nodes list

        del self.__nodes[id]

    def find(self, id):
        if id in self.__nodes:
            return self.__nodes[id]
        else:
            return None

    def new_node(self, *args, **kw):
        new_node = Node(*args, **kw)
        self.add_node(new_node)

        return new_node

    def add_node(self, node):
        node_id = str(node.getId())
        self._add(node, node_id)

    def __cmd_open(self, node_id):
##         image = self.tk.call('tix', 'getimage', 'act_fold')
##         node = self.__nodes[node_id]

##         tree = self.hlist
        
##         tree.tk.call(tree._w, 'item', 'configure', node_id, 0,
##                      '-image', image)
      
        self.__nodes[node_id]._cmd_open()
        
    def __cmd_command(self, node_id):
        self.__nodes[node_id]._cmd_command()

    def __cmd_browse(self, node_id):
        self.__nodes[node_id]._cmd_browse()

    def __cmd_close(self, node_id):

##         image = self.tk.call('tix', 'getimage', 'folder')

##         node = self.__nodes[node_id]
        
##         tree = self.hlist
##         tree.tk.call(tree._w, 'item', 'configure', node_id, 0,
##                      '-image', image)

        self.__nodes[node_id]._cmd_close()

    def open(self, node_id):
        
##         image = self.tk.call('tix', 'getimage', 'act_fold')
##         node = self.__nodes[node_id]

##         tree = self.hlist
        
##         tree.tk.call(tree._w, 'item', 'configure', node_id, 0,
##                      '-image', image)

        Tix.Tree.open(self, node_id)

    def close(self, node_id):
        
##         image = self.tk.call('tix', 'getimage', 'folder')

##         node = self.__nodes[node_id]
        
##         tree = self.hlist
##         tree.tk.call(tree._w, 'item', 'configure', node_id, 0,
##                      '-image', image)

        Tix.Tree.close(self, node_id)
        
class Node:

    id = 0

    def __init__(self, name, open = None, command = None, browse = None,
                 close = None):
        
        self.__name = name
        self.__tree = None
        self.__children = {}
        self.__id = Node.id
        
        Node.id += 1

        if open is None:
            open = self.cmd_open

        if close is None:
            close = self.cmd_close

        if command is None:
            command = self.cmd_command

        if browse is None:
            browse = self.cmd_browse

        self.bindings = {'open': open,
                         'close': close,
                         'command': command,
                         'browse': browse}

        self._setParent(None)
        self.is_open = 0

        self.__user_data = None

    def setUserData(self, d):
        self.__user_data = d

    def getUserData(self):
        return self.__user_data

    def _setParent(self, p):
        self.__parent = p

    def getParent(self):
        return self.__parent

    def getName(self):
        return self.__name

    def setTree(self, t):
        self.__tree = t

    def getTree(self):
        return self.__tree

    def setId(self, id):
        self.__id = id

    def getId(self):
        return self.__id

    def setIcon(self, name):
        tree = self.getTree()
        if tree is None:
            return

        tree = tree.hlist
        
        image = tree.tk.call('tix', 'getimage', name)
        tree.tk.call(tree._w, 'item', 'configure', 
                     self.getId(), 0, '-image', image)

    def add_node(self, node, mode = 'close'):
        
        if self.__tree is None:
            s = 'Node "%s" is not connected yet.\n' + \
                'Nodes have to be connected before ' + \
                'children can be added.' % self.getName()

            raise StandardError, s
        
        node_id = str(self.getId()) + '/%d' % node.getId()
        self.__tree._add(node, node_id)
        self.__children[node_id] = node

        self.getTree().setmode(self.getId(), mode)
        node._setParent(self)

    def new_node(self, *args, **kw):
        new_node = Node(*args, **kw)
        
        self.add_node(new_node)
        return new_node

    def get_children(self):
        return self.__children

    def show(self):
        self.getParent().show_entry()
        self.select()
        self.getTree().hlist.see(self.getId())

    def _cmd_open(self):
        tree = self.getTree()
        entries = self.__children.keys()

        if entries:
            map(tree.hlist.show_entry, entries)

            self.is_open = 1

        return self.bindings['open']()

    def _cmd_close(self):
        tree = self.getTree()
        entries = self.__children.keys()

        if entries:
            map(tree.hlist.hide_entry, entries)

            self.is_open = 0

        return self.bindings['close']()

    def _cmd_command(self):
        
        tree = self.getTree()
        entries = self.__children.keys()

        if entries:
            if not self.is_open:
                map(tree.hlist.show_entry, entries)
                tree.setmode(self.getId(), 'close')
            else:
                map(tree.hlist.hide_entry, entries)
                tree.setmode(self.getId(), 'open')

            self.is_open = not self.is_open
            
        return self.bindings['command']()

    def _cmd_browse(self):
        return self.bindings['browse']()

    def cmd_open(self):
        pass

    def cmd_close(self):
        pass

    def cmd_command(self):
        pass

    def cmd_browse(self):
        pass

    def open(self):
        self.is_open = 1
        self.getTree().open(self.getId())

    def close(self):
        self.is_open = 0
        self.getTree().close(self.getId())

    def show_entry(self):
        self.getTree().hlist.show_entry(self.getId())

    def see(self):
        self.getTree().hlist.see(self.getId())

    def select(self):
        tree = self.getTree()
        tree.hlist.selection_clear()
        tree.hlist.selection_set(self.getId(), self.getId())
        self.see()

    def get_previous(self):
        tree = self.getTree()
        prev = tree.hlist.info_prev(self.getId())
        return tree.find(prev)

    def get_next(self):
        tree = self.getTree()
        next = tree.hlist.info.next(self.getId())
        return tree.find(next)

    def remove(self):
        self.getTree()._del(self)

class Dialog(Tkinter.Toplevel):

    def __init__(self, parent, title, modal = 0, resizable = 0):

        import Tix, Tkinter

        if parent is None:
            parent = Tkinter._default_root

        Tkinter.Toplevel.__init__(self, parent)

        if not resizable:
            self.resizable(0,0)

        place_x = parent.winfo_rootx() + 100
        place_y = parent.winfo_rooty() + 100

##        self.geometry("1x1+0+0")
        self.geometry("+%d+%d" % (place_x, place_y))

        self.title(title)

##         self.grid_rowconfigure(0, weight=1)
##         self.grid_columnconfigure(0, weight=1)
        
        self.setup()

        if modal:
            self.set_modal()
        else:
            self.modal = 0
            self.show()

    def set_modal(self):
        self.modal = 1
        self.focus_set()
        try:
            self.grab_set()
        except:
            pass
        self.mainloop()

    def destroy(self):
        
        if self.modal:
            self.quit()

        Tkinter.Toplevel.destroy(self)

    def setup(self):
        raise NotImplementedError

    def show(self):
        return
        self.list()

import tkSimpleDialog

class OKDialog(Dialog):

    def __init__(self, parent = None, title = '', text = '', modal = 0, image = None):

        self.text = text
        self.image = image
        Dialog.__init__(self, parent, title, modal)

    def setup(self):

        import Tix, Tkinter

        if self.image:
            img = Tix.Label(self, padx=20, pady=10, bd=1, relief=Tix.RAISED,
                          anchor=Tix.CENTER, image = self.image, bg= 'white')
            img.image =self.image
            img.pack(fill=Tix.X)
            
        label = Tix.Label(self, padx=20, pady=10, bd=1, relief=Tix.RAISED,
                          anchor=Tix.CENTER, text=self.text)

        if self.image:
            label.image =self.image
    
        box = DefaultButtonBox(self, orientation=Tix.HORIZONTAL)
        
        box.add('ok', text = 'OK', width = BUTTON_WIDTH,
                command = lambda s = self: s.__ok())
        
        box.pack(side=Tix.BOTTOM, fill=Tix.X)#, anchor=Tix.CENTER)
        label.pack(side = Tix.TOP, fill = Tix.X)

    def __ok(self):
        self.ok()
        self.destroy()

    def ok(self):
        pass

class MessageBox(Dialog):

    def __init__(self, text='', parent=None, title = 'Message', modal = 0):

        self.__text = '\n%s\n' % text

        Dialog.__init__(self, parent, title, modal)

        self.update_idletasks()

        import time
        time.sleep(.1)
        
    def setup(self):

        import Tkinter

        frame = Tkinter.Frame(self)
        frame.grid(row=0, column=0, sticky=Tkinter.NSEW)
        frame.update_idletasks()
        
        label = Tkinter.Label(frame, padx=20, pady=10, bd=1, relief=Tix.RAISED,
                              anchor=Tix.CENTER, text=self.__text)

        label.pack(side = Tix.TOP, fill = Tix.X)
        frame.update_idletasks()
##         label.grid(row=0, column=0, sticky=Tkinter.NSEW)

class YesNoDialog(Dialog):

    def __init__(self, parent = None, title = '', text = '', modal = 0):
        
        self.text = text
        self.__val = None
        Dialog.__init__(self, parent, title, modal)

    def setup(self):

        import Tix, Tkinter
        
        label = Tix.Label(self, padx=20, pady=10, bd=1, relief=Tix.RAISED,
                          anchor=Tix.CENTER, text=self.text)
    
        box = DefaultButtonBox(self, orientation=Tix.HORIZONTAL)
        box.add(YES, text = 'Yes', width = BUTTON_WIDTH,
                command = lambda s = self: s.destroy(YES))
        
        box.add(NO, text = 'No', width = BUTTON_WIDTH,
                command = lambda s = self: s.destroy(NO))
        
        box.pack(side=Tix.BOTTOM, fill=Tix.X)
        label.pack(side = Tix.TOP, fill = Tix.X)

    def get(self):
        return self.__val

    def destroy(self, v):
        self.__val = v
        Dialog.destroy(self)
        
class YesNoCancelDialog(Dialog):

    def __init__(self, parent = None, title = '', text = '', modal = 0):

        self.text = text
        self.__val = None
        Dialog.__init__(self, parent, title, modal)

    def setup(self):

        import Tix

        label = Tix.Label(self, padx=20, pady=10, bd=1, relief=Tix.RAISED,
                          anchor=Tix.CENTER, text=self.text)
    
        box = DefaultButtonBox(self, orientation=Tix.HORIZONTAL)
        box.add(YES, text = 'Yes', width = BUTTON_WIDTH,
                command = lambda s = self: s.destroy(YES))
        
        box.add(NO, text = 'No', width = BUTTON_WIDTH,
                command = lambda s = self: s.destroy(NO))
        
        box.add('cancel', text = 'Cancel', width = BUTTON_WIDTH,
                command = lambda s = self: s.destroy('cancel'))
        
        box.pack(side=Tix.BOTTOM, fill=Tix.X)
        label.pack(side = Tix.TOP, fill = Tix.X)

    def get(self):
        return self.__val

    def destroy(self, v = None):
        self.__val = v
        Dialog.destroy(self)

## class ButtonDialog(Dialog):

##     def __init__(self, title = '', text = '', modal = 0, labels = ('Ok',)):

##         self.text = text
##         self.labels = labels
##         self.__val = None
##         Dialog.__init__(self, title, modal)

##     def setup(self):

##         import Tix, Tkinter
        
##         label = Tix.Label(self, padx=20, pady=10, bd=1, relief=Tix.RAISED,
##                           anchor=Tix.CENTER, text=self.text)
    
##         box = DefaultButtonBox(self, orientation=Tix.HORIZONTAL)

##         self.values = {}
##         self.funcs = []

##         for i in range(len(self.labels)):

##             value = self.labels[i].lower()
##             self.values[value] = self.labels[i]

##             f = lambda s = self: s.destroy(self.labels[i])
##             self.funcs.append(f)
            
##             box.add(value, text = self.values[value], command = self.funcs[-1])

##         box.add(value+'2', text = self.values[value], 
##                 command = lambda s = self: s.destroy('shit'))
        
##         box.pack(side=Tix.BOTTOM, fill=Tix.X)
##         label.pack(side = Tix.TOP, fill = Tix.X)

##     def get(self):
##         return self.__val

##     def destroy(self, v = None):
##         self.__val = v
##         Dialog.destroy(self)

class DirChooser(Dialog):
    
    def __init__(self, parent, title, current_dir = None, modal = 1):

        import os
        
        if current_dir is None:
            self.current_dir = os.getcwd()
        else:
            self.current_dir = current_dir

        Dialog.__init__(self, parent, title, modal, resizable = 1)

    def setup(self):

        dir = Tix.DirTree(self)
        dir.hlist['width'] = 60
        dir.hlist['height'] = 20

        try:
            dir.chdir(self.current_dir)
        except:
            pass

        box = DefaultButtonBox(self, orientation=Tix.HORIZONTAL)
        
        box.add('ok', text = 'Select', width = BUTTON_WIDTH,
                command = lambda s = self: s.copy_name(dir))
        
        box.add('cancel', text = 'Cancel', width = BUTTON_WIDTH,
                command = lambda s = self: s.destroy())
        
        dir.pack(expand=1, fill=Tkinter.BOTH, padx=4, pady=4,
                 side = Tkinter.TOP)
        box.pack(side = Tix.BOTTOM, fill=Tix.X)

        self.val = None
        
    def copy_name (self, dir):
        self.val = dir.cget('value')
        self.destroy()

    def get(self):
        return self.val

class BrowseEntry(Tix.Frame):

    def __init__(self, parent, text = None):

        Tix.Frame.__init__(self, parent)

        self.entry = Tix.Entry(self)
        self.button = DefaultButton(self, text = 'browse...',
                                    command = lambda s = self: s.__browse(),
                                    width=10)


        if text is not None:
            self.label = Tix.Label(self, text = text)
            self.label.pack(side = Tix.LEFT, anchor = Tix.E, fill = Tix.X)
            
        self.entry.pack(side = Tix.LEFT, anchor = Tix.E, fill = Tix.X,
                        expand = 1)
        self.button.pack(side = Tix.RIGHT, anchor = Tix.E)

    def __browse(self):
        val = self.browse()
        self.focus_set()
        if val is None:
            return
            
        self.entry.delete(0, Tix.END)
        self.entry.insert(0, str(val))

    def browse(self):
        """
        returned value will be copied into entry field
        """
        pass

    def focus_set(self):
        self.entry.focus_set()

    def configure(self, *args, **kw):
        self.entry.configure(*args, **kw)

        ## some options also apply to the button

        if 'state' in kw:
            self.button.configure(*args, **kw)

class Cell:

    _last_cell = None
    _current_cell = None

    def __init__(self, parent, row, column, width = 10):

        self.parent = parent
        
        self.label = Tix.Label(parent, width = width, relief = Tix.FLAT,
                               bd = 1, bg = 'white',
                               highlightbackground='grey',
                               highlightthickness=1)
        
        self.label.bind('<Button-1>',
                        func = lambda x, s = self: s.__select(x))
        self.label.bind('<Double-1>',
                        func = lambda x, s = self: s.__edit_mode(x))
        # BARDIAUX 2.2
        # <Double-1> seems not bound with Mac os X quartz-wm
        self.label.bind('<Shift-Button-1>',
                        func = lambda x, s = self: s.__edit_mode(x), add="+")
        
        self.label.bind('<Enter>',
                        func = lambda x, s = self: s.__enter(x))

        self.label.bind('<Leave>',
                        func = lambda x, s = self: s.__leave(x))

        self.label.grid(row = row, column = column, padx = 0, pady = 0,
                        sticky = Tix.EW)
        self.info = self.label.grid_info()

        self.editor = None
        self._in_edit_mode = 0
        self._label_value = None

    def bind_events(self):
        self._bind(self.editor)

    def _bind(self, tk_widget):
         tk_widget.bind('<Return>', func = lambda x,
                        s = self: s.__save(x))
         tk_widget.bind('<Escape>', func = lambda x,
                        s = self: s.__restore(x))
         
    def create_editor(self, constructor):
        width = self.label.cget('width')
        self.editor = constructor(self.parent, relief = Tix.FLAT, bd=1,
                                  width = width)
        self.bind_events()

        return self.editor

    def editor_get(self):
        return self.editor.get()

    def editor_set(self, v):
        self.editor.delete(0, Tix.END)
        self.editor.insert(0, str(v))

    def __enter(self, event):
        self.label.configure(relief = Tix.SUNKEN)
        self.label.configure(highlightbackground='black')

    def __leave(self, event):
        self.label.configure(relief = Tix.FLAT)
        self.label.configure(highlightbackground='grey80')

    def __save(self, event):
        """
        Called when the user presses <return> within
        edit-mode
        """
        
        if self._in_edit_mode:
            saved = self.save()
        else:
            saved = 1

        if saved:
            self._restore()

    def __restore(self, event):
        self._restore()
        self._in_edit_mode = 0

    def _restore(self):
        """
        Restore cell's style to its default.
        If cell was in edit-mode, its value
        is saved.
        """

        ## Set labels background to white
        
        self.label.configure(bg = 'white')

        ## - pack label at its original position
        
        self.label.grid(**self.info)

        if self.editor is not None:
            self.editor.grid_forget()

    def _set_active(self):
        """
        Called when instance has become the active cell
        """

        cell = self._last_cell

        saved = 1
        
        if cell in [self, None]:
            return saved

        if cell._in_edit_mode:
            saved = cell.save()
            
        if saved:
            cell._restore()

        return saved
        
    def save(self):
        """
        Reads cell's editable widget and
        store its value in cell's label.
        Returns 1 if editor's value has been saved.
        """

        ## - unpack widget
        ## - save cell's value

        if self.editor is not None:
            self._label_value = self.editor_get()
            self.label.configure(text = self._label_value)

        ## We left edit-mode

        self._in_edit_mode = 0

        return 1

    def __edit_mode(self, event):
        """
        Called when a cell is double-clicked
        """

        last_cell = self.__class__._last_cell
        if last_cell._in_edit_mode:
            return

        ## Current cell has be become the active cell

        self._set_active()
            
        if self.editor is not None:

            ## Pack editable widget
            
            self.editor.grid(row = self.info['row'],
                             column = self.info['column'], sticky = Tix.EW)

##             if self._label_value is not None:
##                 self.editor_set(self._label_value)

            ## Set focus
            
            self.editor.focus_set()

            ## Temporarily unpack label
            
            self.label.grid_forget()

        ## Indicate, that cell is being edited

        self._in_edit_mode = 1

        ## Now, current cell is last cell

        last_cell = self
        
    def __select(self, event):
       
        """
        Called when a cell is clicked
        """

        ## We are current cell.

        self.__class__._current_cell = self

        ## Change style of last cell to
        ## its default.

        ok = self._set_active()

        if not ok:
            return

        ## Set background of currently
        ## selected cell to grey
        
        self.label.configure(bg = 'grey80')

        ## We are last cell
        
        self.__class__._last_cell = self

    def save_now(self):
        """
        Leaves edit-mode an attempts to store editors
        value.
        """

        saved = self.save()
        if saved:
            self._restore()

    def get(self):
        ## If we try to get the cell's value
        ## when the cell is being edited,
        ## store its value first.

        if self._in_edit_mode:
            return self.editor_get()
        else:
            if self._label_value is not None:
                val = self.editor_get()
            else:
                val = None
            return val

    def set(self, v):
        self._label_value = str(v)
        self.label.configure(text = self._label_value)

        if self.editor is not None:
            self.editor_set(v)
        
    def destroy(self):
        self.label.grid_forget()
        self.label.destroy()

        if self.editor is not None:
            self.editor.grid_forget()
            self.editor.destroy()

        last_cell = self.__class__._last_cell

        if last_cell is self:
            self.__class__._last_cell = None

class Table(Tix.ScrolledWindow):

    def __init__(self, parent, dims, column_widths = None, header = None,
                 *args, **kw):

        Tix.ScrolledWindow.__init__(self, parent, scrollbar='auto',
                                    *args, **kw)

        self.n_columns = dims[1]
        self.rows = []
        self.__n_rows = 0
        
        if header is not None and len(header) > self.n_columns:
            s = 'Length of header must be <= no. of columns.'
            raise ValueError, s

        if column_widths is not None and len(column_widths) <> self.n_columns:
            s = 'Number of column widths must match no. of columns.'
            raise ValueError, s

        if column_widths is None:
            column_widths = [10] * self.n_columns

        if header is not None:
            self.has_header = 1
            self.column_widths = self.create_header(column_widths,
                                                    header, self.n_columns)
            
        else:
            self.column_widths = column_widths
            self.has_header = 0

        self.append_row(dims[0])
        self.grid_propagate()

    def create_cell(self, row, column):
        width = self.column_widths[column]
        
        cell = Cell(self.window, row, column, width)
        cell.create_editor(Tix.Entry)

        return cell

    def __add_row(self):

        offset = self.__n_rows
        if self.has_header:
            offset += 1
##        if self.has_header:
##            offset = len(self.rows) + 1
##        else:
            
##            offset = len(self.rows)

        cells = []

        for i in range(self.n_columns):
            cell = self.create_cell(offset, i)
            cells.append(cell)

        self.rows.append(cells)

        self.__n_rows += 1

    def append_row(self, n = 1):
        new_cells = map(lambda i, f = self.__add_row: f(), range(n))

    def remove_row(self, index):
        for cell in self.rows[index]:
            cell.destroy()

        del self.rows[index]

    def create_header(self, column_widths, header, n_columns):
        header = list(header)
        header += [''] * (n_columns - len(header))

        widths = []
        
        for i in range(len(header)):
            width = max(column_widths[i], len(header[i]))
            label = Tix.Label(self.window, width = width,
                              relief = Tix.GROOVE, bd = 1,
                              bg = 'grey80', text = header[i])
            label.grid(row = 0, column = i, sticky = Tix.EW)
            widths.append(width)

        return widths

    def destroy(self, *args, **kw):
        self.forget()
        
        for row in self.rows:
            [cell.destroy() for cell in row]
            
        Tix.ScrolledWindow.destroy(self, *args, **kw)

    def __getitem__(self, a):
        if type(a) == type(()):
            if len(a) > 2:
                raise IndexError, 'Table is 2-dimensional.'

            cell = self.rows[a[0]][a[1]]

            return cell.get()

        else:
            return [row.get() for row in self.rows[a]]

    def __setitem__(self, a, val):
        if type(a) == type(()):
            if len(a) > 2:
                raise IndexError, 'Table is 2-dimensional.'

            cell = self.rows[a[0]][a[1]]
            cell.set(val)

        else:
            return [cell.set(val) for cell in self.rows[a]]

    def get_current_cell(self):
        if not self.rows:
            return None
        else:
            return self.rows[0][0].__class__._current_cell

    def get_current_row_index(self):
        cs = self.get_current_cell()

        for i in range(len(self.rows)):
            if cs in self.rows[i]:
                return i

        return None

class DefaultButton(Tkinter.Button):

    def __init__(self, *args, **kw):

        Tkinter.Button.__init__(self, *args, **kw)

        if not 'borderwidth' in kw:
            self.configure(borderwidth = 1)

class DefaultButtonBox(Tix.ButtonBox):

    def add(self, *args, **kw):

        if not 'borderwidth' in kw:
            kw['borderwidth'] = 1

        Tix.ButtonBox.add(self, *args, **kw)

# BARDIAUX extendNmr
# inspired by Tix.Tk()
class ariaPopup(Tkinter.Toplevel, Tix.tixCommand):
    """Toplevel widget of Tix which represents mostly the main window
    of an application. It has an associated Tcl interpreter."""
    def __init__(self, master=None, cnf={}, **kw):
        Tkinter.Toplevel.__init__(self, master=master, cnf=cnf, **kw)
        import os
        tixlib = os.environ.get('TIX_LIBRARY')
        self.tk.eval('global auto_path; lappend auto_path [file dir [info nameof]]')
        if tixlib is not None:
            self.tk.eval('global auto_path; lappend auto_path {%s}' % tixlib)
            self.tk.eval('global tcl_pkgPath; lappend tcl_pkgPath {%s}' % tixlib)
        # Load Tix - this should work dynamically or statically
        # If it's static, tcl/tix8.1/pkgIndex.tcl should have
        #               'load {} Tix'
        # If it's dynamic under Unix, tcl/tix8.1/pkgIndex.tcl should have
        #               'load libtix8.1.8.3.so Tix'
        self.tk.eval('package require Tix')

    def destroy(self):
        # For safety, remove an delete_window binding before destroy
        self.protocol("WM_DELETE_WINDOW", "")
        Tkinter.Toplevel.destroy(self)
        

if __name__ == '__main__':

    root = Tix.Tk()
    table = Table(root, (3,4), header=('column1',))
    table.pack()
    table.grid_propagate(0)
