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

import sys, types, copy, bisect, Tkinter

sys.path.append('/home/rungsari/EditObj-0.4/lib/python2.2/site-packages') 
import editobj, editobj.eventobj as eventobj

from Tkinter import *
from UserList import UserList
from UserDict import UserDict
import string, os, os.path
import cPickle as pickle

class Editor:
  require_right_menu = 1
  expand_right       = 0
  
  def __init__(self, master, obj, attr):
    self.obj  = obj
    self.attr = attr
    
  def get_value(self): return getattr(self.obj, self.attr)
  
  def set_value(self, value):
    old_value = getattr(self.obj, self.attr, None)
    
    def do_it():
      if hasattr(self.obj, "set_" + self.attr): # implicit setter
        getattr(self.obj, "set_" + self.attr)(value)
      else: setattr(self.obj, self.attr, value)
      
      if self.master.cancel: self.master.cancel.add_cancel(cancel_it)
      
    def cancel_it():
      if hasattr(self.obj, "set_" + self.attr): # implicit setter
        getattr(self.obj, "set_" + self.attr)(old_value)
      else: setattr(self.obj, self.attr, old_value)
      
      if self.master.cancel: self.master.cancel.add_redo(do_it)
      
    do_it()
    
  def update(self): return 0
  

class BoolEditor(Editor, Tkinter.Checkbutton):
  require_right_menu = 0
  expand_right       = 0
  
  def __init__(self, master, obj, attr):
    Editor.__init__(self, master, obj, attr)
    
    self.value     = Tkinter.BooleanVar()
    self.old_val   = self.get_value()
    self.value.set(self.old_val)
    
    Tkinter.Checkbutton.__init__(self, master, variable = self.value, command = self.validate)
    
  def validate(self, event = None):
    self.old_val = value = self.value.get()
    self.set_value(value)
    
  def get_value(self):
    if getattr(self.obj, self.attr, 0): return 1
    return 0
  
  def set_value(self, value):
    if value: Editor.set_value(self, 1)
    else:     Editor.set_value(self, 0)
    
  def update(self):
    val = self.get_value()
    
    if not (val == self.old_val):
      self.old_val   = val
      self.value.set(self.old_val)
      
    return 1
  

class EntryEditor(Editor, Tkinter.Entry):
  def __init__(self, master, obj, attr):
    Editor.__init__(self, master, obj, attr)
    
    self.unicode   = 0 # Will be overriden by subclasses' get_value()
    self.value     = Tkinter.StringVar()
    self.old_val   = getattr(obj, attr, "")
    self.old_value = self.get_value()
    
    if self.unicode: 
      self.value.set(self.old_value.encode("latin")) # Entry doesn't support Unicode ???
    else:
      self.value.set(self.old_value)
      
    Tkinter.Entry.__init__(self, master, width = 25, textvariable = self.value, selectbackground = "#CCCCFF")
    self.bind("<FocusOut>"       , self.focus_out)
    self.bind("<KeyPress-Return>", self.validate)
    
  def focus_out(self, event = None):
    value = self.value.get()
    
    if self.unicode: value = unicodify(value)
    else:            value = encode   (value)
    
    if value != self.old_value: self.validate()
      
  def validate(self, event = None):
    self.old_value = value = self.value.get()
    
    if self.unicode: value = unicodify(value)
    else:            value = encode   (value)
    
    choices = config.available_attrs.setdefault(self.attr, [])
    if (not value in choices) and (not value in defaultchoices):
      bisect.insort(choices, value)
      
    self.set_value(value)
    
  def get_value(self):
    return repr(getattr(self.obj, self.attr, ""))
  
  def set_value(self, value):
    if value and (value[0] != "<"):
      try: value = editobj.eval(value)
      except: sys.excepthook(*sys.exc_info())
      
      Editor.set_value(self, value)
      
  def update(self):
    try:
      val = getattr(self.obj, self.attr)
      if not val is self.old_val:
        #print "update !!!"
        
        self.old_val   = val
        self.old_value = self.get_value()
        if self.unicode: 
          self.value.set(self.old_value.encode("latin")) # Entry doesn't support Unicode ???
        else:
          self.value.set(self.old_value)
          
    except AttributeError: pass # Not readable => do not update !
    return 1
  
class StringEditor(EntryEditor):
  def get_value(self):
    value = getattr(self.obj, self.attr, "")
    self.unicode = isinstance(value, unicode)
    return value
  
  def set_value(self, value):
    Editor.set_value(self, value)
    
class IntEditor(EntryEditor):
  require_right_menu = 0
  
  def get_value(self): return str(getattr(self.obj, self.attr))
  
  def set_value(self, value):
    Editor.set_value(self, int(editobj.eval(value)))
    
class FloatEditor(EntryEditor):
  require_right_menu = 0
  
  def get_value(self): return str(getattr(self.obj, self.attr, ""))
  
  def set_value(self, value):
    Editor.set_value(self, float(editobj.eval(value)))
    
    
class TextEditor(Editor, Tkinter.Frame):
  """A muti-line text editor."""
  def __init__(self, master, obj, attr):
    Editor.__init__(self, master, obj, attr)
    Tkinter.Frame.__init__(self, master)
    
    self.columnconfigure(0, weight = 1)
    
    value = getattr(obj, attr, "")
    self.unicode = value.__class__ is unicode
    
    self.text = Tkinter.Text(self, width = 25, height = 5, wrap = "word", font = "Helvetica -12", selectbackground = "#CCCCFF")
    self.text.insert("end", value)
    self.text.bind("<FocusOut>", self.validate)
    self.text.grid(sticky = "nsew")
    
    bar = Tkinter.Scrollbar(self, orient = VERTICAL)
    bar.grid(row = 0, column = 1, sticky = "nsew")
    bar['command'] = self.text.yview
    self.text['yscrollcommand'] = bar.set
    
  def validate(self, event = None):
    value = self.text.get("0.0", "end")
    if value[-1] == "\n": value = value[0:-1]
    if self.unicode: value = unicodify(value)
    else:            value = encode   (value)
    
    choices = config.available_attrs.setdefault(self.attr, [])
    
    if (not value in choices) and (not value in defaultchoices):
      bisect.insort(choices, value)
      
    self.set_value(value)
    
  def get_value(self):
    return getattr(self.obj, self.attr, "")
  
  def update(self):
    try:
      val = getattr(self.obj, self.attr)
      if not val is self.old_val:
        self.old_val   = val
        self.old_value = self.get_value()
        self.value.set(self.old_value)
    except AttributeError: pass # Not readable => do not update !
    return 1
  

def RangeEditor(min, max):
  """A slider-based editor."""
  class _RangeEditor(Editor, Tkinter.Scale):
    def __init__(self, master, obj, attr):
      Editor.__init__(self, master, obj, attr)
      Tkinter.Scale.__init__(self, master, orient = Tkinter.HORIZONTAL, from_ = min, to = max)
      self.bind("<ButtonRelease>", self.validate)
      self.update()
      
    def validate(self, event): self.set_value(self.get())
      
    def get_value(self): return int(getattr(self.obj, self.attr, 0))
    
    def update(self):
      try: self.set(getattr(self.obj, self.attr))
      except AttributeError: pass # Not readable => do not update !
      return 1
    
  return _RangeEditor

class _ListEditor(Editor, Tkinter.OptionMenu):
  require_right_menu = 0
  expand_right       = 1
  
  def __init__(self, master, obj, attr, choices):
    Editor.__init__(self, master, obj, attr)
    
    self.value = Tkinter.StringVar()
    self.value.set(self.get_value())
    
    Tkinter.OptionMenu.__init__(self, master, self.value, command = self.validate, *choices)
    
  def validate(self, event = None):
    self.set_value(self.value.get())
    
  def update(self):
    self.value.set(self.get_value())
    return 1
  
class _LongListEditor(Editor, Tkinter.Frame):
  require_right_menu = 0
  expand_right       = 1
  
  def __init__(self, master, obj, attr, choices):
    Editor.__init__(self, master, obj, attr)
    Tkinter.Frame.__init__(self, master)
    
    self.choices = choices
    
    self.columnconfigure(0, weight = 1)
    self.listbox = Tkinter.Listbox(self, exportselection = 0, selectbackground = "#CCCCFF")
    i = 0
    for choice in choices:
      self.listbox.insert(i, choice)
      i = i + 1
    self.listbox.grid(sticky = "nsew")
    
    bar = Tkinter.Scrollbar(self, orient = VERTICAL)
    bar.grid(row = 0, column = 1, sticky = "nsew")
    bar['command'] = self.listbox.yview
    self.listbox['yscrollcommand'] = bar.set
    self.listbox.bind("<ButtonRelease>", self.validate)
    
    self.update()
    
  def validate(self, event = None):
    self.set_value(self.choices[int(self.listbox.curselection()[0])])
    
  def update(self):
    i     = 0
    value = self.get_value()
    while 1:
      if self.choices[i] == value: break
      i = i + 1
    self.listbox.activate(i)
    self.listbox.selection_set(i)
    self.listbox.see(i)
    return 1
  
def ListEditor(values):
  """Create and return a new option-menu like Editor class. The values available are given by VALUES."""
  str2values = {}
  for value in values: str2values[str(value)] = value
  
  if len(values) < 30: super = _ListEditor
  else:                super = _LongListEditor
  
  class ListEditor(super):
    def __init__(self, master, obj, attr):
      values_str = str2values.keys()
      values_str.sort()
      super.__init__(self, master, obj, attr, values_str)
      
    def get_value(self): return str(getattr(self.obj, self.attr))
    
    def set_value(self, value):
      try: value = str2values[value]
      except KeyError:
        for val in values:
          text = str(val)
          if text == value:
            str2values[text] = val
            value = val
            break
            
      Editor.set_value(self, value)
      
  return ListEditor

def EnumEditor(values):
  """Create and return a new option-menu like Editor class. The values available are given by VALUES, which is intented to be a dictionary mapping values to their string representations (e.g. { 0 : "disabled", 1 : "enabled" })."""
  items = values.items()
  items.sort()
  values_str = zip(*items)[1]
  
  if len(values) < 30: super = _ListEditor
  else:                super = _LongListEditor
  
  class ListEditor(super):
    def __init__(self, master, obj, attr):
      super.__init__(self, master, obj, attr, values_str)
      
    def get_value(self): return values.get(getattr(self.obj, self.attr), "")
    
    def set_value(self, value):
      for key, val in values.items():
        if val == value:
          Editor.set_value(self, key)
          break
        
  return ListEditor

def LambdaListEditor(lambd, value_tansformer = None):
  """Create and return a new option-menu like Editor class. The values available are the return values of LAMBD (which must be callable with one arg, the object edited), eventually transformed by VALUE_TRANSFORMER (which should be a callable with one arg)."""
  
  class ListEditor(_ListEditor):
    def __init__(self, master, obj, attr):
      _ListEditor.__init__(self, master, obj, attr, (getattr(obj, attr, ""),))
      
      self.values = None
      self.bind("<ButtonPress>", self.on_list)
      
    def on_list(self, event = None):
      if self.values: return
      
      self.values = lambd(self.obj)
      self.str2value = {}
      
      menu     = self["menu"]
      variable = self.value
      callback = self.validate
      
      menu.delete(0)
      
      for value in self.values:
        text = str(value)
        self.str2value[text] = value
        menu.add_command(label = text, command = Tkinter._setit(variable, text, callback))
        
    def get_value(self): return str(getattr(self.obj, self.attr))
    
    def set_value(self, value):
      value = self.str2value[value]
      if value_tansformer: value = value_tansformer(value)
      Editor.set_value(self, value)
      
  return ListEditor


class CallableEditor(Editor, Tkinter.Button):
  """An editor that displays a button; when it is clicked, the corresponding attribute is called."""
  require_right_menu = 0
  expand_right       = 0
  
  def __init__(self, master, obj, attr):
    Editor.__init__(self, master, obj, attr)
    Tkinter.Button.__init__(self, master, text = attr, command = self.button_click)
    
  def button_click(self, event = None):
    getattr(self.obj, self.attr)()


class WithButtonEditor(Editor, Tkinter.Frame, object):
  require_right_menu = 0
  expand_right       = 1
  
  def __init__(self, master, obj, attr, internal_editor_class, button_text):
    Editor.__init__(self, master, obj, attr)
    
    Tkinter.Frame.__init__(self, master)
    self.columnconfigure(0, weight = 1)
    
    self.internal_editor = internal_editor_class(self, obj, attr)
    #self.internal_editor.grid(row = 0, col = 0, **{"in" : self, "sticky" : "EW"})
    self.internal_editor.grid(row = 0, col = 0, sticky = "NSEW")
    
    self.button = Tkinter.Button(self, text = button_text, command = self.button_click)
    self.button.grid(row = 0, col = 1)
    
  def button_click(self): pass
  
  def get_value(self): return self.internal_editor.get_value()
  def set_value(self, value): self.internal_editor.set_value(value)
  def update(self): self.internal_editor.update()
  
  def get_cancel(self): return self.master.cancel
  cancel = property(get_cancel)
  
class FilenameEditor(WithButtonEditor):
  def __init__(self, master, obj, attr):
    WithButtonEditor.__init__(self, master, obj, attr, StringEditor, "...")
    
  def button_click(self):
    import tkFileDialog
    s = tkFileDialog.askopenfilename()
    if s: self.set_value(s)
    
class DirnameEditor(WithButtonEditor):
  def __init__(self, master, obj, attr):
    WithButtonEditor.__init__(self, master, obj, attr, StringEditor, "...")
    
  def button_click(self):
    import tkFileDialog
    s = tkFileDialog.askdirectory()
    if s: self.set_value(s)
    

class SubEditor(Editor, Tkinter.Frame, object):
  """An editor that edits its values in an inner property frame."""
  def __init__(self, master, obj, attr):
    Editor.__init__(self, master, obj, attr)
    
    Tkinter.Frame.__init__(self, master)
    self.columnconfigure(0, weight = 1)
    
    import main
    self.propframe = main.EditPropertyFrame(self)
    self.propframe.edit(self.get_value())
    self.propframe.pack()
    
  def get_cancel(self): return self.master.cancel
  cancel = property(get_cancel)


class TypeEditor(Tkinter.OptionMenu):
  def __init__(self, master, obj, attr):
    self.obj           = obj
    self.attr          = attr
    self.value         = Tkinter.StringVar()
    
    self.choices = config.available_attrs.setdefault(attr, [])
    apply(Tkinter.OptionMenu.__init__, [self, master, self.value, ""] + self.choices + ["(edit...)", "(set in clipboard)"] + defaultchoices, {"command" : self.validate})
    
  def validate(self, event = None):
    import editobj
    
    value = self.value.get()
    self.value.set("")
    
    if   value == "(edit...)":            editobj.edit(getattr(self.obj, self.attr))
    elif value == "(set in clipboard)":   editobj.clipboard = getattr(self.obj, self.attr)
    elif value == "(paste clipboard)":    setattr(self.obj, self.attr, editobj.clipboard); self.master.fields[self.attr].update()
    elif value == "(deepcopy clipboard)": setattr(self.obj, self.attr, copy.deepcopy(editobj.clipboard)); self.master.fields[self.attr].update()
    else:                                 self.master.fields[self.attr].set_value(value)



attrs_editors  = {
  # Hide items and children, since they are edited by the tree view.
  (None, "items")    : None,
  (None, "children") : None,
  }

DEFAULT_EDITOR = EntryEditor

def register_attr(attr, editor, clazz = None):
  attrs_editors[clazz, attr] = editor
  
def find(obj, attr):
  if isinstance(obj, editobj.MultiEdit): clazz = obj._objs[0].__class__
  else:                                  clazz = obj.__class__
  
  try:                   mro = clazz.__mro__
  except AttributeError: mro = _mro(clazz)
  
  NO_REPLY = 0 # 0 is not a valid editor !!! (None is one)
  
  for clazz in mro:
    editor = attrs_editors.get((clazz, attr), NO_REPLY)
    if not editor is NO_REPLY: return editor
    
  return attrs_editors.get((None, attr), DEFAULT_EDITOR)

def _mro(clazz):
  mro = [clazz]
  for c in clazz.__bases__: mro.extend(_mro(c))
  return mro







class _Config:
  def __init__(self):
    self.available_items = {}
    self.available_attrs = {}
    
CONFIG_FILE = os.path.expanduser("~" + os.sep + ".editobj.config")
try: config = pickle.load(open(CONFIG_FILE, "r"))
except:
  print "EditObj : no config file -- I create a new ~/.editobj.config"
  config = _Config()

def saveconfig():
  pickle.dump(config, open(CONFIG_FILE, "w"))

def encode(s):
  if type(s) is unicode:
    return s.encode("latin")
  return s

def unicodify(s):
  if type(s) is str:
    try:                 return unicode(s, "utf8")
    except UnicodeError: return unicode(s, "latin")
  return s

def register_children(clazz, code_expressions):
  code_expressions = map(unicodify, code_expressions)
  clazz = str(clazz)
  try:    config.available_items[clazz].extend(code_expressions)
  except: config.available_items[clazz] = list(code_expressions)

def register_auto_child(clazz, code_expression):
  code_expression = unicodify(code_expression)
  clazz = str(clazz)
  config.available_items[clazz] = code_expression

def register_values(attr, code_expressions):
  code_expressions = map(unicodify, code_expressions)
  try:    config.available_attrs[attr].extend(code_expressions)
  except: config.available_attrs[attr] = list(code_expressions)


defaultchoices = [u"(paste clipboard)", u"(deepcopy clipboard)"]

class AbstractDialog(Tkinter.Toplevel):
  def __init__(self, choices, showkey = 0):
    Tkinter.Toplevel.__init__(self)
    self.columnconfigure(0, weight = 1)
    self.columnconfigure(1, weight = 0)
    
    row = 0
    
    if showkey:
      self.entryk = Tkinter.Entry(self)
      self.entryk.insert(0, "key")
      self.entryk.grid(row = 0, columnspan = 2, sticky = "EW")
      row = row + 1
      
    self.var = Tkinter.StringVar()
    self.entry = Tkinter.Entry(self, textvariable = self.var)
    self.entry.grid(row = row, columnspan = 2, sticky = "EW")
    self.entry.focus_set()
    self.entry.bind("<KeyPress-Return>", self.validate)
    row = row + 1
    
    self.list = Tkinter.Listbox(self)
    if len(choices) > 0:
      apply(self.list.insert, [0] + choices)
    self.list.bind("<ButtonRelease-1>", self.list_selected)
    self.list.bind("<Double-1>", self.validate)
    self.list.bind("<KeyPress-Delete>", self.deletechoice)
    self.list.grid(row = row, column = 0, sticky = "NSEW")
    self.rowconfigure(row, weight = 1)
    
    self.bar = Tkinter.Scrollbar(self, orient = Tkinter.VERTICAL)
    self.bar.grid(row = row, column = 1, sticky = "NSEW")
    row = row + 1
    
    self.list['yscrollcommand'] = self.bar.set
    self.bar['command'] = self.list.yview
    
    
    self.ok = Tkinter.Button(self, text = "OK", command = self.validate)
    self.ok.grid(row = row, columnspan = 2, sticky = "EW")
    
  def list_selected(self, event = None):
    self.var.set(self.list.selection_get())
    
  def deletechoice(self, event = None):
    self.choices.remove(self.list.selection_get())
    self.list.delete(ACTIVE)
    
  def validate(self, event = None):
    self.choosen(self.var.get())
    self.destroy()
    

class AddDialog(AbstractDialog):
  def __init__(self, addinto, callback = None, cancel = None): # Node is the Node to update
    self.addinto  = addinto
    self.callback = callback
    self.cancel   = cancel
    
    if hasattr(addinto, "_EventObj_stuff"): clazz = addinto._EventObj_stuff.clazz
    else:                                   clazz = addinto.__class__
    if (clazz is UserList) or (clazz is UserDict): choicesfor = clazz
    else:
      choicesfor = list_or_dict_class_of(clazz)
      if choicesfor is None:
        choicesfor = clazz
        
    self.choices = config.available_items.setdefault(str(choicesfor), [])
    if not isinstance(self.choices, list):
      self.choosen(self.choices)
      return
    
    # Add the current item value or class to the list of possible item.
    items = editobj.items_of(addinto)
    for item in items:
      try:
        if isinstance(addinto, eventobj._EventObj): # EventObj change the class of the instance.
          clazz = item._EventObj_stuff.clazz
          self.addchoice(clazz.__module__ + "." + clazz.__name__ + "()")
        else: self.addchoice(item.__class__.__module__ + "." + item.__class__.__name__ + "()")
      except: pass
      
    AbstractDialog.__init__(self, self.choices + defaultchoices, editobj.is_dict(addinto))
    88
  def addchoice(self, choice):
    if (isinstance(self.choices, list)) and (not choice in self.choices) and (not choice in defaultchoices):
      self.choices.append(choice)
      self.choices.sort()
      
  def choosen(self, text):
    self.addchoice(text)
    
    if   text == "(paste clipboard)"   : added = editobj.clipboard
    elif text == "(deepcopy clipboard)": added = copy.deepcopy(editobj.clipboard)
    else:                                added = editobj.eval(text, locals = {"parent" : self.addinto})
    
    if isinstance(self.addinto, UserDict) or type(self.addinto) is DictType:
      key = self.entryk.get()
      def _do_it    (): self.addinto[key] = added
      def _cancel_it(): del self.addinto[key]
    elif hasattr(self.addinto, "append"):
      def _do_it    (): self.addinto.append(added)
      def _cancel_it(): self.addinto.remove(added)
    elif hasattr(self.addinto, "add"):
      def _do_it    (): self.addinto.add   (added)
      def _cancel_it(): self.addinto.remove(added)
    else:
      def _do_it    (): editobj.items_of(self.addinto).append(added)
      def _cancel_it(): editobj.items_of(self.addinto).remove(added)
      
    def do_it():
      _do_it()
      if self.callback: self.callback(self.addinto, added)
      
      if self.cancel: self.cancel.add_cancel(cancel_it)
      
    def cancel_it():
      _cancel_it()
      if self.callback: self.callback(self.addinto, added)
      
      if self.cancel: self.cancel.add_redo(do_it)
      
    do_it()
    
def list_or_dict_class_of(clazz):
  "Figure out what is the class that extends UserList or UserDict in the class inheritance tree..."
  
  # Guess that the the type of the item is associated with the class that inherit UserList/UserDict
  if (list in clazz.__bases__) or (dict in clazz.__bases__) or (UserList in clazz.__bases__) or (UserDict in clazz.__bases__):
    return clazz
  # Recursive
  for base in clazz.__bases__:
    answer = list_or_dict_class_of(base)
    if answer: return answer
    

def choices_for_attr(attrname):
  return config.available_attrs.setdefault(attrname, [])

