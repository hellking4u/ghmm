# EditObj
# Copyright (C) 2001-2002 Jean-Baptiste LAMY -- jiba@tuxfamily
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from Tkinter import *
import os, sys, types

sys.path.append('/home/rungsari/EditObj-0.4/lib/python2.2/site-packages')
from editobj import *
import editobj, editobj.eventobj as eventobj, editobj.editor as editor, editobj.treewidget as treewidget


class EditWindow(Toplevel):
  def __init__(self, o = None, view_hidden = 0, right_menu = 1, dialog = 0, command = None, preamble = None, cancel = None, grab = 0):
    # Get, or create if needed, the root Tk.
    try:
      from Tkinter import _default_root
      tkroot = _default_root
    except ImportError: pass
    if not tkroot:
      tkroot = Tk(className = 'EditObj')
      tkroot.withdraw()
      
    self.history = []
    self.dialog  = dialog
    
    if not cancel:
      from cancel import Context
      cancel = Context()
      
    if grab and cancel:
      #import cancel as cancel_module
      #
      #self.cancel_context = cancel
      #self.cancel = cancel_module.CancelmentStack()
      #
      #def _call():
      #  self.cancel()
      #  
      #self.cancel.__call__ = _call
      #
      #cancel.add_cancel(self.cancel)
      
      # If grab (= modal), push the cancel stack so as all the modification are cancelled at the same time.
      cancel.push()
    self.cancel = cancel
    
    if not dialog:
      menubar = Menu(tkroot)
      self.picklemenu = picklemenu = Menu(menubar)
      picklemenu.add_command(label = 'Load...'   , command = self.pickle_load)
      picklemenu.add_command(label = 'Save'      , command = self.pickle_save)
      picklemenu.add_command(label = 'Save as...', command = self.pickle_save_as)
      menubar.add_cascade(label='Pickle', menu = picklemenu)
      
      self.copymenu = copymenu = Menu(menubar)
      copymenu.add_command(label = 'Copy...'     , command = self.copy_copy)
      copymenu.add_command(label = 'Deep copy...', command = self.copy_deepcopy)
      menubar.add_cascade(label='Copy', menu = copymenu)
      
      self.editmenu = editmenu = Menu(menubar)
      if self.cancel:
        editmenu.add_command(label = 'Undo'                 , command = self.edit_undo)
        editmenu.add_command(label = 'Redo'                 , command = self.edit_redo)
        editmenu.add_separator()
      editmenu.add_command(label = 'New hierarchy view...', command = self.edit_newview)
      editmenu.add_command(label = 'New property view...' , command = self.edit_newpropview)
      editmenu.add_command(label = 'Back'                 , command = self.edit_back)
      editmenu.add_command(label = 'Console...'           , command = self.display_console)
      menubar.add_cascade(label='Edit', menu = editmenu)
    else:
      menubar = None
      
    self.filename = None
    
    Toplevel.__init__(self, tkroot, menu = menubar)
    self.bind("<Destroy>", self.__del)
    
    self.withdraw()
    
    if preamble:
      t = Text(self, width = 0, height = 10, wrap = "word", highlightthickness = 0, font = "Helvetica -12", selectbackground = "#CCCCFF")
      t.pack(fill = BOTH, side = TOP)
      t.insert("end", preamble)
      
    if dialog:
      def _command(event = None):
        # Focus the OK button so as the currently active editor will be unfocused (as some editors validate themselves when they are unfocused).
        ok.focus_set()
        self.update() # Ensure the Tk unfocus event will be sent
        self.destroy()
        if command:
          def do_it():
            command()
            if cancel: cancel.add_post_cancel(cancel_it) # Do it AFTER the other cancellable operation !!!
            
          def cancel_it():
            command()
            if cancel: cancel.add_post_redo(do_it)
            
          do_it()
          
        if grab and cancel: cancel.pop()
        
      ok = Tkinter.Button(self, text = "OK", command = _command)
      ok.pack(expand = 0, fill = X, side = BOTTOM)
      self.wm_protocol("WM_DELETE_WINDOW", _command)
      
    self.scrollpane = treewidget.ScrollPane(self, 0, 1)
    self.scrollpane.pack(expand = 1, fill = BOTH, side = BOTTOM)
    self.propframe = EditPropertyFrame(self, view_hidden, right_menu, self.cancel, bd = 0)
    self.propframe.edit = self.edit
    self.scrollpane.setContent(self.propframe)
    
    self.edited          = None
    self.hierarchy       = None
    self.hierarchyedited = None
    self.buttonsframe    = None
    self.console         = None
    self.edit(o)
    
    if dialog: ok.tkraise()
    
    self.deiconify()
    if grab: self.grab_set()
    
  def createButtonsFrame(self):
    self.buttonsframe = Frame(self, borderwidth = 2, relief = "raised")
    self.buttonsframe.pack(expand = 0, fill = BOTH, side = TOP)
    self.buttonsframe.columnconfigure(5, weight = 1)
    self.addbutton = Button(self.buttonsframe, padx = 5, text = editobj.TRANSLATOR("Add..."), command = self.button_add, relief = "flat")
    self.addbutton.grid(row = 0, column = 0, sticky = E + W)
    self.removebutton = Button(self.buttonsframe, padx = 5, text = editobj.TRANSLATOR("Remove"), command = self.button_remove, relief = "flat")
    self.removebutton.grid(row = 0, column = 1, sticky = E + W)
    self.upbutton = Button(self.buttonsframe, padx = 5, text = editobj.TRANSLATOR("Up"), command = self.button_up, relief = "flat")
    self.upbutton.grid(row = 0, column = 2, sticky = E + W)
    self.downbutton = Button(self.buttonsframe, padx = 5, text = editobj.TRANSLATOR("Down"), command = self.button_down, relief = "flat")
    self.downbutton.grid(row = 0, column = 3, sticky = E + W)
    if self.propframe.right_menu:
      self.clipbutton = Button(self.buttonsframe, padx = 5, text = editobj.TRANSLATOR("Copy"), command = self.button_clip, relief = "flat")
      self.clipbutton.grid(row = 0, column = 4, sticky = E + W)
      
  #def __del__(self):
  #  print "__del__ de", self
    
  def __del(self, event):
    if (type(event.widget) is StringType and event.widget[1:] == self.winfo_name()) or event.widget == self:
      if not self.dialog:
        self.picklemenu.destroy()
        self.copymenu  .destroy()
        self.editmenu  .destroy()
        
      self.unsetevent()
      #if self.hierarchy: self.hierarchy.node.unseteventtree()
      editor.saveconfig() # Save the config file.
      
  def setevent  (self): pass
  def unsetevent(self): pass
  
  def edit(self, o, newview = 0):
    if newview: edit(o, self.propframe.view_hidden)
    else:
      locked = self.winfo_ismapped()
      if locked: self.pack_propagate(0)
      
      o = editable(o)
      
      self.history.append(o)
      self.edited = o
      if hasattr(o, "filename"): self.filename = o.filename
      EditPropertyFrame.edit(self.propframe, o)
      self.scrollpane.updateContentSize()
      
      if self.console:
        self.console.dict["obj"] = o
        
      # If the new edited object is in the same hierarchy than the former, keep the same hierarchy. Else, create one.
      if (self.hierarchyedited is None) or not ((self.hierarchyedited is o) or in_hierarchy(o, self.hierarchyedited)):
        if hasattr(o, "__edit__"): # Only if it is a new hierarchy
          if o.__edit__(self):
            self.destroy()
            return
          
        if is_hierarchy(o):
          self.unsetevent()
          self.hierarchyedited = self.edited
          self.setevent()
          if self.hierarchy is None:
            self.createButtonsFrame()
            self.hierarchy = treewidget.Tree(self)
            self.hierarchy.pack(expand = 1, fill = BOTH)
          ItemNode(self.hierarchy, o, editor = self.edit)
        else:
          if not self.hierarchyedited is None:
            self.hierarchy.destroy()
            self.unsetevent()
            self.hierarchy = None
            self.hierarchyedited = None
            self.buttonsframe.destroy()
            self.buttonsframe = None
      else:
        if hasattr(o, "__subedit__"): # Else, a new edition as a sub component
          o.__subedit__(self)
          
      if locked: self.pack_propagate(1)
      
  def button_add(self):
    if not self.hierarchy.selection: return
    
    o = self.hierarchy.selection[0]
    if not is_hierarchy(o.item): return
    
    items = items_of(o.item)
    if not o.event_ok: # eventobj doesn't work => update manually.
      def callback(into, added):
        o.update()
        o.updatechildren()
    else: callback = None
    d = editor.AddDialog(o.item, callback, self.cancel)
    
  def button_remove(self):
    if not self.hierarchy.selection: return
    
    o = self.hierarchy.selection[0]
    if   hasattr(o.itemparent, "remove"):
      def _do_it    (): o.itemparent.remove(o.item)
      def _cancel_it():
        if hasattr(o.itemparent, "append"): o.itemparent.append(o.item)
        else:                               o.itemparent.add   (o.item)
    elif is_dict(o.itemparent):
      for key, value in o.itemparent.items():
        if value is o.item:
          def _do_it    (): del o.itemparent[key]
          def _cancel_it(): o.itemparent[key] = o.item
          break
    else:
      items = items_of(o.itemparent)
      index = items.index(o.item)
      def _do_it    (): del items[index]
      def _cancel_it(): items.insert(index, o.item)
      
    def do_it():
      _do_it()
      
      if not o.parent.event_ok: # eventobj doesn't work => update manually.
        o.parent.update()
        o.parent.updatechildren()
        
      if self.cancel: self.cancel.add_cancel(cancel_it)
      
    def cancel_it():
      _cancel_it()
      
      if not o.parent.event_ok: # eventobj doesn't work => update manually.
        o.parent.update()
        o.parent.updatechildren()
        
      if self.cancel: self.cancel.add_redo(do_it)
      
    do_it()
    
  def button_up(self):
    o = self.hierarchy.selection[0]
    items = items_of(o.itemparent)
    i = items.index(o.item)
    
    def do_it():
      items[i - 1], items[i] = items[i], items[i - 1]
      
      if not o.parent.event_ok: # eventobj doesn't work => update manually.
        o.parent.update()
        o.parent.updatechildren()
        
      if self.cancel: self.cancel.add_cancel(cancel_it)
      
    def cancel_it():
      items[i - 1], items[i] = items[i], items[i - 1]
      
      if not o.parent.event_ok: # eventobj doesn't work => update manually.
        o.parent.update()
        o.parent.updatechildren()
        
      if self.cancel: self.cancel.add_redo(do_it)
      
    do_it()
    
  def button_down(self):
    o = self.hierarchy.selection[0]
    items = items_of(o.itemparent)
    i = items.index(o.item)
    
    if not o.parent.event_ok: # eventobj doesn't work => update manually.
      o.parent.update()
      o.parent.updatechildren()
      
    def do_it():
      if i + 1 == len(items): items[0    ], items[i] = items[i], items[0]
      else:                   items[i + 1], items[i] = items[i], items[i + 1]
      
      if not o.parent.event_ok: # eventobj doesn't work => update manually.
        o.parent.update()
        o.parent.updatechildren()
        
      if self.cancel: self.cancel.add_cancel(cancel_it)
      
    def cancel_it():
      if i + 1 == len(items): items[0    ], items[i] = items[i], items[0]
      else:                   items[i + 1], items[i] = items[i], items[i + 1]
      
      if not o.parent.event_ok: # eventobj doesn't work => update manually.
        o.parent.update()
        o.parent.updatechildren()
        
      if self.cancel: self.cancel.add_redo(do_it)
      
    do_it()
    
  def button_clip(self):
    editobj.clipboard = self.hierarchy.selection[0].item
    
  def pickle_save(self):
    obj = self.hierarchyedited or self.edited
    if hasattr(obj, "save"):
      try:
        obj.save()
        return
      except TypeError: sys.excepthook(*sys.exc_info())
    if self.filename:
      import cPickle as pickle
      pickle.dump(obj, open(self.filename, "w"))
    else: self.pickle_save_as()
  def pickle_save_as(self):
    obj = self.hierarchyedited or self.edited
    import cPickle as pickle, tkFileDialog
    self.filename = tkFileDialog.asksaveasfilename()
    if hasattr(obj, "save"):
      try:
        obj.save(self.filename)
        return
      except TypeError: pass
    pickle.dump(obj, open(self.filename, "w"))
  def pickle_load(self):
    import cPickle as pickle, tkFileDialog
    self.filename = tkFileDialog.askopenfilename()
    self.edit(pickle.load(open(self.filename, "r")))
    
  def copy_copy(self):
    import copy
    edit(copy.copy(self.edited), self.propframe.view_hidden)
  def copy_deepcopy(self):
    import copy
    edit(copy.deepcopy(self.edited), self.propframe.view_hidden)
    
  def edit_undo(self): self.cancel.cancel()
  def edit_redo(self): self.cancel.redo()
  def edit_newview(self):
    if self.hierarchyedited is None: self.edit_newpropview()
    else: edit(self.hierarchyedited, self.propframe.view_hidden)
  def edit_newpropview(self):
    edit(self.edited, self.propframe.view_hidden)
  def edit_back(self):
    if len(self.history) >= 2:
      self.history.pop()
      self.edit(self.history.pop())
  def display_console(self):
    if not self.console:
      import editobj.console
      dict = { "root" : self.hierarchyedited or self.edited, "obj" : self.edited }
      #dict.update(editobj.EVAL_ENV)
      self.console = editobj.console.Console(self, dict = dict, globals = editobj.EVAL_ENV)
      
      self.console.text.insert("end", """\nYou can access the currently edited obj as "obj", and the root of the hierarchy as "root".\n""")
      self.console.text.insert("end", sys.ps1)
      self.console.text.configure(width = 10, height = 10)
      
      self.console.pack(fill = BOTH, expand = 1, side = BOTTOM)
    else:
      self.console.destroy()
      self.console = None
      
      
class EditPropertyFrame(Frame):
  def __init__(self, master, view_hidden = 0, right_menu = 1, cancel = None, **opts):
    self.view_hidden   = view_hidden
    self.edited        = None
    self.event_ok      = 0
    self.update_cancel = None
    self.right_menu    = right_menu
    self.cancel        = cancel
    
    Frame.__init__(self, master, opts)
    self.columnconfigure(0, weight = 0)
    self.columnconfigure(1, weight = 1)
    self.columnconfigure(2, weight = 0)
    self.bind("<Destroy>", self.__del)
    
  def __del(self, event = None):
    self.unsetevent()
    
  def setevent(self):
    try:
      eventobj.addevent(self.edited, self.edited_changed)
      self.event_ok = 1
    except eventobj.NonEventableError: 
      self.event_ok = 0
  def unsetevent(self):
    if self.event_ok:
      eventobj.removeevent(self.edited, self.edited_changed)
      self.event_ok = 0
      
  def edit(self, o):
    self.unsetevent()
    
    self.fields = {}
    
    # Unfocus the current widget -- avoid change in this widget to be lost !
    focused = self.focus_lastfor()
    if hasattr(focused, "focus_out"): focused.focus_out() # Hack !!! The "else" case should work with ANY editor widget, but it doesn't...?? In particular, it works with TextEditor, but not with StringEditor !
    else:
      focused.tk_focusNext().focus_set()
      focused.update() # Ensure the Tk unfocus event will be sent
      
    for widget in self.children.values(): widget.destroy()
    
    self.edited = o
    if not o is None:
      self.edited_attrs = attrs_of(o)
      attrs = map(lambda attr: (editobj.TRANSLATOR(attr), attr), self.edited_attrs)
      attrs.sort()
      
      line = 1
      for (translation, attr) in attrs:
        if self.view_hidden or attr[0] != "_":
          fieldclass = editor.find(o, attr)
          if not fieldclass: continue
          
          label = Label(self, text = translation)
          label.grid(row = line, column = 0, sticky = W + N + S)
          
          field = fieldclass(self, o, attr)
          self.fields[attr] = field
          
          if self.right_menu and field.require_right_menu:
            field.grid(row = line, column = 1, sticky = E + W + N + S)
            
            change = editor.TypeEditor(self, o, attr)
            change.grid(row = line, column = 2, sticky = E + W + N + S)
          else:
            field.grid(row = line, column = 1, columnspan = 1 + field.expand_right, sticky = E + W)
            
          line = line + 1
          
      self.setevent()
      
  def update_field(self, name):
    field = self.fields.get(name)
    if field is None:
      if (name[0] == "_") and not self.view_hidden: return
      self.edit(self.edited) # An unknown field ?? Rebuild all.
      
  def edited_changed(self, o, attr, value, oldvalue):
    if not(attr is eventobj.ADDITION or attr is eventobj.REMOVAL):
      if (value is None) and (not attr in self.edited_attrs): self.edit(self.edited) # A new attribute
      
      if self.update_cancel: self.after_cancel(self.update_cancel)
      self.update_cancel = self.after(100, self.update)
      
  def update(self):
    for field in self.fields.values(): field.update()
    self.update_cancel = None

class ItemNode(treewidget.Node):
  def __init__(self, parent, item, itemparent = None, editor = edit, key = None, show_key = 0):
    self.item          = item
    self.itemparent    = itemparent
    self.edit          = editor
    self.key           = key
    self.show_key      = show_key
    self.event_ok      = 0
    self.update_cancel = None
    treewidget.Node.__init__(self, parent)
    self.setevent()
    
  def setevent(self):
    try:
      eventobj.addevent(self.item, self.content_changed)
      self.event_ok = 1
    except eventobj.NonEventableError:
      self.event_ok = 0
  def unsetevent(self):
    if self.event_ok:
      eventobj.removeevent(self.item, self.content_changed)
      self.event_ok = 0
  def destroy(self):
    treewidget.Node.destroy(self)
    self.unsetevent()
    
  def unseteventtree(self):
    self.unsetevent()
    for child in self.children: child.unseteventtree()
    
  def __unicode__(self):
    if (not self.show_key) or (self.key is None): return unicode(self.item)
    return unicode(self.key) + ": " + unicode(self.item)
  
  def createchildren(self, oldchildren = ()):
    children = []
    
    the_list = None
    
    if   is_list(self.item): the_list = self.item
    elif is_dict(self.item): the_dict_items = self.item
    else:
      if hasattr(self.item, "children"): items = self.item.children
      else:                              items = self.item.items
      if callable(items): items = items()
      the_list = items
      
    if not the_list is None:
      i = 0
      for subitem in the_list:
        for oldchild in oldchildren: # Re-use the child if possible.
          if oldchild.item is subitem:
            oldchild.key = i # Index may have changed...
            oldchild.update() # ... so we need to update.
            children.append(oldchild)
            oldchildren.remove(oldchild) # Do not re-use it twice !
            break
        else: children.append(ItemNode(self, subitem, self.item, self.edit, i, 0))
        i = i + 1
      for child in oldchildren:
        child._undrawtree()
        child.unseteventtree()
        
    else:
      for child in oldchildren:
        for key, item in the_dict_items:
          if item is child.item and key == child.key:
            children.append(child) # Re-use this one.
            items.remove((key, item))
            break
        else: # Delete this one.
          child._undrawtree()
          child.unseteventtree()
      for key, item in the_dict_items: # The new ones stay left.
        children.append(ItemNode(self, item, self.item, self.edit, key, 1))
        
    return children
  
  def iseditable(self): return 0
  def isexpandable(self): return is_hierarchy(self.item) and (len(items_of(self.item)) > 0)
  def settext(self, text): pass
  
  def select(self, event = None):
    treewidget.Node.select(self)
    typ = type(self.item)
    if (typ is types.IntType) or (typ is types.FloatType) or (typ is types.StringType) or (typ is types.ComplexType) or (typ is types.LongType):
      # self.item is not an editable object ! => Wrap it !
      self.edit(StringOrNumberWrapper(self.key, self.item, self.itemparent, self))
    else: self.edit(self.item)
    
  def expand(self, event = None):
    treewidget.Node.expand(self, event)
    if hasattr(self.item, "__editchildren__"): self.item.__editchildren__(1)
    
  def collapse(self, event = None):
    treewidget.Node.collapse(self, event)
    if hasattr(self.item, "__editchildren__"): self.item.__editchildren__(0)
    
  def content_changed(self, o, name, value, oldvalue):
    if self.update_cancel: self.tree.after_cancel(self.update_cancel)
    if name is eventobj.ADDITION or name is eventobj.REMOVAL:
      self.update_cancel = self.tree.after(100, self._delay_update_all)
    else:
      self.update_cancel = self.tree.after(100, self._delay_update)
      
  def _delay_update(self):
    self.update_cancel = None
    self.update()
    
  def _delay_update_all(self):
    self.update_cancel = None
    self.update()
    self.updatechildren()
    
class StringOrNumberWrapper:
  def __init__(self, key, value, parent, itemnode):
    self._key = key
    self.value = value
    self._parent = parent
    self._itemnode = itemnode
    
  def __setattr__(self, name, value):
    if name == "value" and hasattr(self, "_parent"):
      self._itemnode.item = value
      self._itemnode.update()
      self._parent[self._key] = value
    self.__dict__[name] = value
    
  def __cmp__(self, other): return cmp(self.value, other)
  def __eq__ (self, other): return self.value == other


