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

"""EditObj -- Display a Tkinter-based dialog box to edit any Python instance, list or dict.

This is what Java calls a "Bean Editor", but for Python :-)

The module only provide the function edit(obj), that returns the Tk top-level window.
EditObj will edit all attributes and, for list or dict, all items.
The upper part of the window display the list/dict items in a hierarchical tree view, and the lower part lists all the edited instance's attributes.

An object can provide an __edit__(window) method, that will be called on edition, with the EditObj's window as argument. If it returns true, it is assumed that the edition is entirely done by the object itself.
If an object is edited as a child / item / component of another, __subedit__ will be called, with the EditObj's window as argument.
When the children of an object are made visible or invisible in the tree view, the object __editchildren__ method is called (if it has one) with the visibility (1 or 0).
An object can also provide a __wrapedit__() method, that returns a wrapper for the object. The wrapper will be edited instead of the object.

editobj.TRANSLATOR can be set to a translator function (such as the ones from the gettext module), if you want to translate the names of the properties.

Quick example :
>>> import UserList, editobj
>>> class C(UserList.UserList):
...   def __init__(self):
...     self.x = 1
...     self.append("A string.")
...     self.append(UserList.UserList([1, 2, 3]))
>>> toplevel = editobj.edit(C())
>>> toplevel.mainloop()"""


TRANSLATOR = lambda key: key

import sys, Tkinter
from UserList import UserList
from UserDict import UserDict

def edit(o, view_hidden = 0, **kargs):
  from editobj.main import EditWindow
  
  return EditWindow(o, view_hidden, **kargs)


def debug_edit(o, view_hidden = 0):
  oldclass = o.__class__
  
  import main
  
  w = main.EditWindow(o, view_hidden)
  
  import gc
  #gc.set_debug(gc.DEBUG_COLLECTABLE | gc.DEBUG_SAVEALL | gc.DEBUG_INSTANCES)
  #gc.set_debug(gc.DEBUG_LEAK)
  
  gc.collect()
  weak = 0
  
  def call_gc():
    weak
    print "gc...", gc.collect(), o.__class__ is oldclass
    #print gc.garbage
    
    
    if gc.garbage:
      gc.collect()
      gc.collect()
      gc.collect()
      
      print "EditWindow"       , filter(lambda o: isinstance(o, main.EditWindow), gc.garbage)
      print "EditPropertyFrame", filter(lambda o: isinstance(o, main.EditPropertyFrame), gc.garbage)
      print "ItemNode"         , filter(lambda o: isinstance(o, main.ItemNode), gc.garbage)
      print "EditPropertyFrame utilisé par :"
      for user in gc.get_referrers(weak()):
        print user
        print
      
      print o.__class__, o.__class__ is oldclass
      try: print o._EventObj_stuff.events
      except: pass
      
      import sys
      sys.exit()
      
    Tkinter._default_root.after(5, call_gc)
  call_gc()
  
  import weakref
  def callback(o):
    print "La fénêtre a été fermée !!!"
  weak = weakref.ref(w.propframe, callback)
  print weak
  
  w = None
  
  return w

clipboard = None

EVAL_ENV = {}

_eval = eval
def eval(text, globals = EVAL_ENV, locals = EVAL_ENV):
  """An eval that automatically import modules."""
  i = text.find ("(")
  j = text.rfind(".", 0, i)
  
  if j != -1:
    modulename = text[:j]
    try: int(modulename)
    except ValueError:
      classname = text[j + 1:]
      
      #try:
        #module = sys.modules.get(modulename)
        #if not module:
        #  module = __import__(modulename, globals, locals, ())
        #  print module, id(module)
        #globals[modulename.split(".")[0]] = module # Add the loaded module.
      #except: pass
      
  try:
    return _eval(text, globals, locals)
  except:
    sys.excepthook(*sys.exc_info())
    print "Error ?", "-- consider input as a string."
    return text
  
def is_hierarchy(obj): return is_list(obj) or is_dict(obj) or hasattr(obj, "children") or hasattr(obj, "items")
def is_list     (obj): return isinstance(obj, list) or isinstance(obj, UserList)
def is_dict     (obj): return isinstance(obj, dict) or isinstance(obj, UserDict)
def suscriptable(obj): return hasattr(obj, "__getitem__")
def eventable   (obj): return hasattr(obj, "__dict__")
def items_of    (obj):
  if   is_list(obj): return obj
  elif is_dict(obj): return obj.values()
  else:
    if   hasattr(obj, "children"): items = obj.children
    elif hasattr(obj, "items"   ): items = obj.items
    else: return ()
    if callable(items): items = items()
    return items

def is_getset(o): return isinstance(o, property) or (type(o).__name__ == "getset_descriptor")
def attrs_of(obj):
  try:    attrs = obj.__dict__.keys()
  except: attrs = []
  
  if hasattr(obj, "__members__"): attrs.extend(obj.__members__)
  
  klass = obj.__class__
  attrs.extend(filter(lambda attr: is_getset(getattr(klass, attr)), dir(klass)))
  return attrs
  

def editable(o):
  """Wrap o if needed (try to call o.__wrapedit__()). Return o itself if no wrapper is provided."""
  try: return o.__wrapedit__()
  except:
    #if isinstance(o, Tkinter.Wm):
    #  import tk_wrapper
    #  return tk_wrapper.Tk_wrapper(o)
    return o

def in_hierarchy(item, hierarchy):
  """A recursive __contains__, for working with hierarchical tree.
The caller must assume that, if needed, item and hierarchy are wrapped (by _getEditedObject)."""
  for i in items_of(hierarchy):
    # Wrap i, if needed (already done for item and hierarchy).
    i = editable(i)
    if item == i: return 1
    if in_hierarchy(item, i): return 1
    
def remove(list, item):
  if hasattr(list, "remove"): list.remove(item)
  elif is_dict(list):
    for key, value in list.items():
      if value is item:
        del list[key]
        break
  else:
    items_of(list).remove(item)
    
def append(list, item):
  if hasattr(list, "append"): list.append(item)
  elif is_dict(list):
    for key, value in list.items():
      if value is item:
        del list[key]
        break
  else:
    items_of(list).remove(item)
    
def insert(list, index, item):
  if hasattr(list, "remove"): list.remove(item)
  elif is_dict(list):
    for key, value in list.items():
      if value is item:
        del list[key]
        break
  else:
    items_of(list).remove(item)


class MultiEdit(object):
  def __init__(self, objs):
    self.__dict__["_objs"] = objs
    
  def __addevent__(self, event):
    for obj in self._objs:
      eventobj.addevent(obj, event)
      
  def __removeevent__(self, event):
    for obj in self._objs:
      eventobj.removeevent(obj, event)
      
  def __hasevent__(self, event = None):
    return eventobj.hasevent(self._objs[0], event)
  
  def _get_members(self):
    attrs = []
    for obj in self._objs:
      for attr in attrs_of(obj):
        if not attr in attrs: attrs.append(attr)
    return attrs
  __members__ = property(_get_members)
  
  def __getattr__(self, attr):
    return getattr(self._objs[0], attr)
    
#     value = getattr(self._objs[0], attr)
#     for obj in self._objs:
#       val = getattr(obj, attr)
#       if val != value:
#         if   val.__class__ is str  : return ""
#         elif val.__class__ is int  : return 0
#         elif val.__class__ is float: return 0.0
#         else:                        return None
#     return value
  
  def __setattr__(self, attr, value):
    for obj in self._objs:
      setattr(obj, attr, value)
      
      
