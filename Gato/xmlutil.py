#!/usr/bin/env python
################################################################################
#
#       This file is part of Gato (Graph Algorithm Toolbox) 
#       version _VERSION_ from _BUILDDATE_. You can find more information at 
#       http://www.zpr.uni-koeln.de/~gato
#
#	file:   HMMEd.py
#	author: Alexander Schliep (schliep@zpr.uni-koeln.de)
#
#       Copyright (C) 1998-2002, Alexander Schliep, Winfried Hochstaettler and 
#       ZAIK/ZPR, Universitaet zu Koeln
#                                   
#       Contact: schliep@zpr.uni-koeln.de, wh@zpr.uni-koeln.de             
#
#       Information: http://gato.sf.net
#
#       This library is free software; you can redistribute it and/or
#       modify it under the terms of the GNU Library General Public
#       License as published by the Free Software Foundation; either
#       version 2 of the License, or (at your option) any later version.
#
#       This library is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#       Library General Public License for more details.
#
#       You should have received a copy of the GNU Library General Public
#       License along with this library; if not, write to the Free
#       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
#
#       This file is version $Revision$ 
#                       from $Date$
#             last change by $Author$.
#
################################################################################

from DataStructures import Point2D,EdgeWeight
from Graph import Graph, SubGraph

from xml.dom.minidom import *
import EditObjectAttributesDialog
from EditObjectAttributesDialog import EditObjectAttributesDialog, ValidatingString, ValidatingInt, ValidatingFloat, PopupableInt, Probability, DefaultedInt, DefaultedString

import whrandom
import string
import types
import copy

import sys 
import HMMXML

import logging
log = logging.getLogger("xmlutil.py")

def typed_assign(var, val):
    result = type(var)(val)
    result.__dict__ = var.__dict__
    #result.__dict__ = copy.copy(var.__dict__)
    return result

def listFromCSV(s, type):
    return map(type,string.split(s,','))

def csvFromList(list, perRow = None):
    if perRow == None:
        return string.join(map(str,list), ', ')
    else:
        result = ""
        for start in xrange(0, len(list), perRow):
            result += string.join(map(str,list[start:start+perRow]), ', ') + ',\n'
        return result[0:len(result)-2]

def writeContents(XMLDoc, XMLNode, data):
    contents = XMLDoc.createTextNode("%s" % data)
    XMLNode.appendChild(contents)

def writeData(XMLDoc, XMLNode, dataKey, dataValue):
    data = XMLDoc.createElement("data")
    data.setAttribute('key', "%s" % dataKey)
    contents = XMLDoc.createTextNode("%s" % dataValue)
    data.appendChild(contents)
    XMLNode.appendChild(data)

def writeXMLData(XMLDoc, XMLNode, dataKey, XMLData):
    data = XMLDoc.createElement("data")
    data.setAttribute('key', "%s" % dataKey)
    data.appendChild(XMLData)
    XMLNode.appendChild(data)

def writeXMLTextNode(XMLDoc, XMLNode, keyName, XMLData):
    data = XMLDoc.createElement(keyName)
    contents = XMLDoc.createTextNode("%s" % XMLData)
    data.appendChild(contents)
    XMLNode.appendChild(data)

class NamedDistributions:
    def __init__(self, itsHMM):
        self.initialize()
        self.itsHMM = itsHMM
        self.code2name = {-1:'None'}
        self.name2code = {'None':-1}
        self.maxCode = 0

    def initialize(self):
        self.dist = {}
        self.order = {}
       
    def addDistribution(self, name, order, p): 
        self.dist[name] = p
        self.order[name] = order
        self.code2name[self.maxCode] = name
        self.name2code[name] = self.maxCode
        self.maxCode += 1
        print self.dist
        
    def deleteDistribution(self, name):
        del self.dist[name]
        del self.order[name]
        del self.code2name[self.name2code[name]]
        del self.name2code[name]
   
    def fromDOM(self, XMLNode):
        self.initialize()
        datas = XMLNode.getElementsByTagName("hmm:background")
        for data in datas:
            dataKey = data.attributes['key'].nodeValue
            dataOrder = int(data.attributes['order'].nodeValue)
            dataValue = ""
            for child in data.childNodes:
                dataValue += child.nodeValue
            p = listFromCSV(dataValue, types.FloatType)
            self.addDistribution(dataKey, dataOrder, p)

    def toDOM(self, XMLDoc, XMLNode):
        for name in self.dist.keys():
            print "background: name = ", name, self.dist[name], self.order[name]
            background_elem = XMLDoc.createElement("hmm:background")
            background_elem.setAttribute('key', "%s" % name)
            background_elem.setAttribute('order', "%s" % self.order[name])
            if self.order[name] == 0:
                contents = XMLDoc.createTextNode(csvFromList(self.dist[name]))
            else:
                contents = XMLDoc.createTextNode(csvFromList(self.dist[name],
                                                             self.itsHMM.hmmAlphabet.size()))               
            background_elem.appendChild(contents)
            XMLNode.appendChild(background_elem)

    def names(self):
        return self.dist.keys()

   

class XMLElementWriter:
    
    def __init__(self):
        import string
        self._string = string
        self.element_keys = ['id', 'hmm:low', 'hmm:high', 'hmm:type', 'for', 'gd:type', 'code', 'key', 'order', 'x', 'y', 'source', 'target']
        
    def _write_data(self, writer, data):
        "Writes datachars to writer."
        replace = self._string.replace
        data = replace(data, "&", "&amp;")
        data = replace(data, "<", "&lt;")
        data = replace(data, "\"", "&quot;")
        data = replace(data, ">", "&gt;")
        writer.write(data)

    def writexml(self, XMLNode, writer, indent="", addindent="", newl=""): 
        """This version of write xml makes sure text nodes are 
        surrounded by tags for easier reading rather than being 
        on lines by themselves."""  

        # indent = current indentation
        # addindent = indentation to add to higher levels
        # newl = newline string
        
        if ( XMLNode.nodeType == XMLNode.TEXT_NODE ): 
            self._write_data(writer, "%s%s%s"%(indent, XMLNode.data, newl))
        else:
            writer.write(indent+"<" + XMLNode.nodeName)
            
            # build attribute list
            a_names = []
            try:
                for key in self.element_keys:
                    if ( XMLNode.getAttribute(key) != ""):
                        a_names.append( key )
                        a_names.sort()
            except AttributeError:
                a_names = []
        
            for a_name in a_names:
                writer.write(" %s=\"" % a_name)
                self._write_data(writer, XMLNode.getAttribute(a_name))
                writer.write("\"")
            if XMLNode.childNodes:
                writer.write(">%s"%(newl))
                for node in XMLNode.childNodes:
                    if node.nodeType!=node.TEXT_NODE: 
                        self.writexml(node,writer,indent+addindent,addindent,newl)
                    else:
                        writer.seek(writer.tell()-1)
                        self.writexml(node,writer,"",addindent,"")

                if XMLNode.childNodes[-1].nodeType!=node.TEXT_NODE:
                    writer.write("%s</%s>%s" % (indent,XMLNode.nodeName,newl))
                else:
                    writer.write("</%s>%s" % (XMLNode.nodeName,newl))
            else:
                writer.write("/>%s"%(newl))
                     

#
# Notice for this function:
# minidom.Element.toprettyxml is not so pretty and its output format is not compatible
# with XMLIO parser. The particular problem with minidom.toprettyxml is that
# it put the text data of a text node on a new line, instead of immediately after the element tag.
# Because XMLIO cannot parse this format, thus we need our own pretty print program 
#
def toprettyxml( XMLDoc ):
    # we can't use cStringIO since it doesn't support Unicode strings
    from StringIO import StringIO
    writer = StringIO()
    prettydoc = XMLElementWriter()
    writer.write('<?xml version="1.0" ?>\n')
    for node in XMLDoc.childNodes:
        prettydoc.writexml(node, writer, "","  ", "\n")
    return writer.getvalue();
    

class DOM_Map:
    def __init__(self):
        self.initialize()

    def initialize(self):
        self.name = {}
        self.desc = {}
        self.hasDesc = None
        self.name2code = {}
        
    def addCode(self, code, name, desc = None):
        self.name[code] = name
        if desc != None:
            self.desc[code] = desc
            self.hasDesc = 1
        self.name2code[name] = code

    def low(self):
        if len(self.name.keys()) > 0:
            return min(self.name.keys())
        else:
            return 0
                
    def high(self):
        if len(self.name.keys()) > 0:
            return max(self.name.keys())
        else:
            return 0
    
    def fromDOM(self, XMLNode):
        pass

    def symbolsFromDom(self, XMLNode):
        symbols = XMLNode.getElementsByTagName("symbol")
        
        for symbol in symbols:
            symbolCode = ValidatingInt(int(symbol.getAttribute("code")))
            symbolName = ValidatingString(symbol.firstChild.nodeValue)
            symbolDesc = symbol.getAttribute("desc")
            if symbolDesc != None:
                self.addCode(symbolCode, symbolName, ValidatingString(symbolDesc))
            else:
                self.addCode(symbolCode, symbolName)
                
    def toDOM(self, XMLDoc, XMLNode):
        XMLNode.setAttribute('hmm:low', "%s" % self.low())
        XMLNode.setAttribute('hmm:high', "%s" % self.high())
        map = XMLDoc.createElement("map")  
        for key in self.name.keys():
            symbol = XMLDoc.createElement("symbol")
            symbol.setAttribute('code', "%s" % key)
            if self.hasDesc and self.desc[key] != "":
                symbol.setAttribute('desc', "%s" % self.desc[key])
            writeContents(XMLDoc, symbol, "%s" % self.name[key])
            map.appendChild(symbol)
        XMLNode.appendChild(map)
   
    def buildList(self):
	return self.name.keys()
     
# -------------------------------------------
#  Exceptions

class HMMEdError(Exception):
    def __init__(self, message):
	print "\n\n Unknown error types. Please report \n\n"
	
class NotValidHMMType(HMMEdError):
    def __init__(self,message):
       print "\n\n Probabilities missing xception: " + str(message) + "\n"

class AlphabetErrorType(HMMEdError):
    def __init__(self,message):   
        print "\n\n Alphabet exception: " + str(message) + "\n"
	
	
class DiscreteHMMAlphabet(DOM_Map):
    def __init__(self):
        DOM_Map.__init__(self)
        self.hmm_type = 'discrete'

    def fromDOM(self, XMLNode):
        """Take dom subtree representing a <hmm:alphabet</hmm:alphabet> element"""
        self.initialize()
        # Not reading: hmm:low hmm:high
        if XMLNode.getAttribute("hmm:type") == self.hmm_type:
            self.symbolsFromDom(XMLNode)
        else:
            print "DiscreteHMMAlphabet wrong type %s" % XMLNode.getAttribute("hmm:type") 

    def toDOM(self, XMLDoc, XMLNode):
        hmmalphabet = XMLDoc.createElement("hmm:alphabet")
        hmmalphabet.setAttribute('hmm:type', 'discrete')
	hmmalphabet.setAttribute('hmm:low', "%s" % self.low())
        hmmalphabet.setAttribute('hmm:high', "%s" % self.high())
        map = XMLDoc.createElement("map")  
        for key in self.name.keys():
            symbol = XMLDoc.createElement("symbol")
            symbol.setAttribute('code', "%s" % key)
            if self.hasDesc and self.desc[key] != "":
                symbol.setAttribute('desc', "%s" % self.desc[key])
            writeContents(XMLDoc, symbol, "%s" % self.name[key])
            map.appendChild(symbol)
        hmmalphabet.appendChild(map)
	#  DOM_Map.toDOM(self, XMLDoc, hmmalphabet)
        XMLNode.appendChild(hmmalphabet)

    def toGHMM(self, XMLDoc, XMLNode):
        hmmalphabet = XMLDoc.createElement("alphabet")
        for key in self.name.keys():
            alphabet = XMLDoc.createElement("symbol")
            alphabet.setAttribute('id', "%s" % key)
            hmmalphabet.appendChild(alphabet)
        XMLNode.appendChild(hmmalphabet)
      
   
    def size(self):
        return len(self.name.keys())
    
    def buildList(self):
	return self.name.keys()

    def buildAlphabets(self, nrOfSymbols):
	alphas = range(nrOfSymbols)
	alphas = map( lambda x: 'a'+ str(x), alphas)
        for code in range(nrOfSymbols):
	    self.addCode( code, alphas[code], desc = None)
 
	
class HMMClass(DOM_Map):
    def __init__(self):
        DOM_Map.__init__(self)
        self.code2name = {-1:'None'}
        self.name2code = {'None':-1}
        self.maxCode = 0
    
    def fromDOM(self, XMLNode):
        """Take dom subtree representing a <hmm:class></hmm:class> element"""
        self.initialize()
        self.symbolsFromDom(XMLNode)
        # copy self.name to self.code2name
        for key in self.name.keys():            
            self.code2name[key] = self.name[key]
        
    def toDOM(self, XMLDoc, XMLNode):        
        hmmclass = XMLDoc.createElement("hmm:class")   
        DOM_Map.toDOM(self, XMLDoc, hmmclass)
        XMLNode.appendChild(hmmclass)

    def size(self):
        return len(self.name.keys())
            
class HMMState:

    def __init__(self, nodeIndex, itsHMM):

        self.initial  = Probability("0.0")
        self.label    =  ValidatingString("None")
        self.itsHMM   = itsHMM
	print type(self.label)
	
        self.index = nodeIndex # The node index in the underlying graph
        self.id    = DefaultedInt() # identification by the canvas, not always the same
	
	self.state_class = PopupableInt(-1)
        self.state_class.setPopup(itsHMM.hmmClass.code2name, itsHMM.hmmClass.name2code, 10)

	self.order = DefaultedInt()
        self.order.setDefault(1, 0)

        self.emissions = []

        self.tiedto = DefaultedString()
        self.tiedto.setDefault(1, '')
        self.desc = self.id

        self.reading_frame = PopupableInt(-1)
        code2name = {-1:'None', 0:'0', 1:'1', 2:'2'}
        name2code = {'None':-1, '0':0, '1':1, '2':2}
        self.reading_frame.setPopup(code2name, name2code, 4)

        self.duration = DefaultedInt()
        self.duration.setDefault(1, 0)

        self.background = PopupableInt(-1)
        self.background.setPopup(self.itsHMM.backgroundDistributions.code2name, self.itsHMM.backgroundDistributions.name2code, 10)

        self.editableAttr = ['label', 'state_class', 'initial', 'order', 'background']
        self.xmlAttr = self.editableAttr + ['ngeom', 'emissions']
        
    editableAttr = ['label', 'initial', 'order', 'background', 'tiedto', 'reading_frame']
    xmlAttr = editableAttr + ['ngeom', 'emissions']
    # ['id', 'state_class', 'label', 'order', 'initial', 'tiedto', 'reading_frame', 'duration', 'background']

    #
    # Set_<Attribute> Methods: for integration with class editobj.Editor.set_value()
    # the name of the method must be in the form of set_<attr name>(self,value).
    # Otherwise, EditObj cannot propogate new values. (See the base method editobj.EntryEditor.set_value())
    # When necessary, we also update the Graph here.
    #
    def set_label(self, value):
        # Get the id2index from the Graph
        oldv = self.itsHMM.id2index[self.id]
        self.label = typed_assign(self.label, value)
        
        # We only show the label out of the editable items
        self.itsHMM.G.labeling[oldv] = "%s\n%s" % (self.id, self.label) # XXX Hack Aaaargh!
        self.itsHMM.itsEditor.UpdateVertexLabel(oldv, 0)


    def set_initial(self, value):
        self.initial = typed_assign(self.initial, value)

    def fromDOM(self, XMLNode):

        self.id = typed_assign(self.id, int(XMLNode.attributes['id'].nodeValue)) # state's id
        
        self.index = self.itsHMM.G.AddVertex()
        
        datas = XMLNode.getElementsByTagName("data")
        for data in datas:
            dataKey = data.attributes['key'].nodeValue
            dataValue = data.firstChild.nodeValue

            #print dataValue
                       
            if dataKey == 'class':
                #if len(self.itsHMM.hmmClass.name2code.keys()) == 1:
                #    key = self.itsHMM.hmmClass.name2code.keys()
                #    self.state_class = typed_assign(self.state_class, self.itsHMM.hmmClass.name2code[key[0]])
                #else:
                self.state_class = typed_assign(self.state_class, int(dataValue)) # code for the state class

            elif  dataKey == 'label':
                self.label = type(self.label)(dataValue.encode('ascii', 'replace'))

            elif  dataKey == 'order':
                if dataValue == None: # use default value
                    self.order = typed_assign(self.order, self.order.defaultValue)
                    self.order.useDefault = 1
                else:
                    self.order = typed_assign(self.order, int(dataValue))
                    self.order.useDefault = 0

            elif  dataKey == 'initial':
                self.initial = typed_assign(self.initial, float(dataValue))
                
            elif  dataKey == 'tiedto':
                
                if dataValue == None: # use default value
                    self.tiedto = typed_assign(self.tiedto, self.tiedto.defaultValue)
                    self.tiedto.useDefault = 1
                else:
                    self.tiedto = typed_assign(self.tiedto, dataValue.encode('ascii', 'replace'))
                    self.tiedto.useDefault = 0

            elif dataKey == 'reading-frame':
                self.reading_frame = typed_assign(self.reading_frame, int(dataValue))

            elif dataKey == 'background':
                self.background = typed_assign(self.background, self.itsHMM.backgroundDistributions.name2code[dataValue])

            elif dataKey == 'duration':
                self.duration = typed_assign(self.duration, int(dataValue))
                self.duration.useDefault = 0
                                    
            elif dataKey == 'ngeom':
                # We only use pos
                pos = XMLNode.getElementsByTagName('pos')[0] # Just one pos ...                
                self.pos = Point2D(float(pos.attributes['x'].nodeValue),
                                   float(pos.attributes['y'].nodeValue))
                
            elif dataKey == 'emissions':
                # collect all strings from childnodes
                dataValue = ""
                for child in data.childNodes:
                    dataValue += child.nodeValue
                self.emissions = listFromCSV(dataValue, types.FloatType)
                #print self.emissions
                    
            else:
                print "HMMState.fromDOM: unknown key %s of value %s" % (dataKey, dataValue)
        

    def toDOM(self, XMLDoc, XMLNode, initial_sum):
        node = XMLDoc.createElement("node")
        node.setAttribute('id', "%s" % self.id)

        # Mandatory elems
        writeData(XMLDoc, node, 'label', self.label)
        writeData(XMLDoc, node, 'class', self.state_class)
        writeData(XMLDoc, node, 'initial', self.initial / initial_sum)
       
	pos = self.itsHMM.G.embedding[self.index] 
        pos_elem = XMLDoc.createElement("pos")
        pos_elem.setAttribute('x', "%s" % pos.x)
        pos_elem.setAttribute('y', "%s" % pos.y)
        writeXMLData(XMLDoc, node, 'ngeom', pos_elem)

        if not self.order.useDefault:
            writeData(XMLDoc, node, 'order', self.order)

        if self.reading_frame != -1:
            writeData(XMLDoc, node, 'reading-frame', self.reading_frame)

        if self.background != -1:
            writeData(XMLDoc, node, 'background', self.itsHMM.backgroundDistributions.code2name[self.background])

        if not self.duration.useDefault:
            writeData(XMLDoc, node, 'duration', self.duration)

        if not self.tiedto == '': # XXX tied emission
            writeData(XMLDoc, node, 'tiedto', self.tiedto)
            self.emissions = self.itsHMM.state[self.itsHMM.id2index[int(self.tiedto)]].emissions # XXX
            
        if self.order.useDefault:
            order = 0
        else:
            order = self.order

        # XXX Produce uniform emission probs, if we dont have the correct number of
        # parameters            
        size = self.itsHMM.hmmAlphabet.size()**(order+1)
        if len(self.emissions) != size:
            tmp = [1.0/self.itsHMM.hmmAlphabet.size()] * self.itsHMM.hmmAlphabet.size()
            if order == 0:
                self.emissions = tmp
            else:
                self.emissions = tmp * self.itsHMM.hmmAlphabet.size()**order

        if order > 0:
            writeData(XMLDoc, node, 'emissions', csvFromList(self.emissions,
                                                             self.itsHMM.hmmAlphabet.size()))
        else:
            writeData(XMLDoc, node, 'emissions', csvFromList(self.emissions))            
        XMLNode.appendChild(node)

    def toGHMM(self, XMLDoc, XMLNode, initial_sum):
        node = XMLDoc.createElement("state")
        node.setAttribute('id', "%s" % self.id)

        writeXMLTextNode(XMLDoc, node, 'initial', self.initial / initial_sum)
        # ignore order
        writeXMLTextNode(XMLDoc, node, 'emission', string.join(map(str,self.emissions),'\n'))
        XMLNode.appendChild(node)
        
        
class HMM:    
    def __init__(self, XMLFileName = None, G = None):

        # self.itsEditor = itsEditor
        if ( G is None ):
            self.G = Graph()
        else:
            self.G = G

        self.G.directed = 1
        self.G.euclidian = 0
        self.G.simple = 0
        self.Pi = {}
        self.id2index = {}
            
        self.hmmAlphabet = DiscreteHMMAlphabet()
        self.hmmClass    = HMMClass()
        
        self.editableAttr = {}
        self.editableAttr['HMM'] = ['desc']
        self.desc = ValidatingString()       

        self.state = {}

        self.backgroundDistributions = NamedDistributions(self)

        self.DocumentName = "graphml"
        if XMLFileName != None:
            self.OpenXML(XMLFileName)


    def Clear(self):
        self.G.Clear()
        self.Pi = {}
        self.id2index = {}
            
        self.hmmAlphabet = DiscreteHMMAlphabet()
        self.hmmClass    = HMMClass()
        self.backgroundDistributions = NamedDistributions(self)
        
        self.editableAttr = {}
        self.editableAttr['HMM'] = ['desc']
        self.desc = ValidatingString()       
        self.state = {}
        self.DocumentName = "graphml"        


    def AddState(self, index, label='None'):
        state = HMMState(-1, self)
        state.id = max(self.id2index.keys()) + 1
        state.index = index
        self.id2index[state.id] = state.index
        self.state[state.index] = state # XXX Use canvas id
        state.label = typed_assign(state.label, state.id)
        self.G.labeling[state.index] = "%s" % (state.label)
        return state.index
        
    def DeleteState(self, index):
	""" The method only deletes a map between index and its state object.
	    The caller must delete the corresponding vertex in the owner Graph self.G. """
	del self.id2index[self.state[index].id]
	del self.state[index]

    def fromDOM(self, XMLNode):
        
        # self.hmmClass.fromDOM(XMLNode.getElementsByTagName("hmm:class")[0])
        class_elements = XMLNode.getElementsByTagName("hmm:class")
        if class_elements != []:
            for tag in XMLNode.getElementsByTagName("hmm:class"):
                self.hmmClass.fromDOM(tag)

        # One "hmm:alphabet" XML element
        self.hmmAlphabet.fromDOM(XMLNode.getElementsByTagName("hmm:alphabet")[0]) 
        self.backgroundDistributions.fromDOM(XMLNode)

        nodes = XMLNode.getElementsByTagName("node")
        for n in nodes:
            state = HMMState(-1, self)
            state.fromDOM(n)
            self.state[state.index] = state # key must be string
            self.id2index[state.id] = state.index
            self.G.embedding[state.index] = state.pos
            self.G.labeling[state.index] = "%s\n%s" % (state.id, state.label) # XXX Hack Aaaargh!

        edges = XMLNode.getElementsByTagName("edge")
        nr_classes = int(self.hmmClass.high()-self.hmmClass.low())+1
        for i in range(nr_classes):
            self.G.edgeWeights[i] = EdgeWeight(self.G)

        for edge in edges:
            i = self.id2index[int(edge.attributes['source'].nodeValue)]
            j = self.id2index[int(edge.attributes['target'].nodeValue)]

            datas = edge.getElementsByTagName("data")
            for data in datas:
                dataKey = data.attributes['key'].nodeValue
                # dataValue = data.firstChild.nodeValue

            if dataKey == 'prob':
                #p = float(dataValue)
                # collect all strings from childnodes
                dataValue = ""
                for child in data.childNodes:
                    dataValue += child.nodeValue
                p = listFromCSV(dataValue, types.FloatType)
                self.G.AddEdge(i, j)
                if len(p) == 1: # only one class
                    for cl in range(nr_classes):
                        p.append(0.0)
                        
                for cl in range(nr_classes):
                    self.G.edgeWeights[cl][(i,j)] = p[cl]

    def modelCheck(self):
	
        # Compute sums of initial probabilities for renormalization 
        initial_sum = 0.0
        for s in self.state:
            initial_sum = initial_sum + self.state[s].initial

        print "Initial sum", initial_sum
        print "Alphabet size is",self.hmmAlphabet.size()
        
	if initial_sum == 0.0:
	    raise NotValidHMMType("Initial state is not specified.")
	    
	if self.hmmAlphabet.size() == 0.0:
	    raise AlphabetErrorType("Alphabet object is empty. You must create alphabet before saving.")
	
    def toDOM(self, XMLDoc, XMLNode):
        graphml = XMLDoc.createElement("graphml")
        XMLNode.appendChild(graphml)

        # Create key elements
        hmmtype = XMLDoc.createElement("key")
        hmmtype.setAttribute('id', 'emissions')
        hmmtype.setAttribute('gd:type', 'HigherDiscreteProbDist') # what's your type?
        hmmtype.setAttribute('for', 'node')
        graphml.appendChild(hmmtype)
        
        self.hmmClass.toDOM(XMLDoc, graphml)
        self.hmmAlphabet.toDOM(XMLDoc, graphml) 
        self.backgroundDistributions.toDOM(XMLDoc, graphml) 

        graph = XMLDoc.createElement("graph")

        # Compute sums of initial probabilities for renormalization 
        initial_sum = 0.0
        for s in self.state.keys():
            initial_sum = initial_sum + self.state[s].initial
        
        for s in self.state.keys():
            self.state[s].toDOM(XMLDoc, graph, initial_sum)
        
        # Compute sums of outgoing probabilities for renormalization of transition probabilities
        # NOTE: need dictionaries here
        out_sum = {}
        nr_classes = int(self.hmmClass.high())-int(self.hmmClass.low())+1
        for v in self.G.vertices:
            out_sum[v] = [0.0]*nr_classes

        for cl in range(1): # XXX Assuming one transition class
            for e in self.G.Edges():
                if self.G.edgeWeights[cl].has_key(e):
                    out_sum[e[0]][cl] = out_sum[e[0]][cl] + self.G.edgeWeights[cl][e]
                
        for e in self.G.Edges():
            transitions = []
            edge_elem = XMLDoc.createElement("edge")
            edge_elem.setAttribute('source', "%s" % self.state[e[0]].id)
            edge_elem.setAttribute('target', "%s" % self.state[e[1]].id)
            # writeData(XMLDoc, edge_elem, 'prob', self.G.edgeWeights[cl][e] / out_sum[e[0]])
            # XXX Assuming one transition class for cl in range(nr_classes):
            for cl in range(1):
                if self.G.edgeWeights[cl].has_key(e) and out_sum[e[0]][cl]:
                    transitions.append(self.G.edgeWeights[cl][e]/ out_sum[e[0]][cl])
                else:
                    transitions.append(0.0)
                
            writeData(XMLDoc, edge_elem, 'prob', csvFromList( transitions ))

            graph.appendChild(edge_elem)  
            
        graphml.appendChild(graph)

    def AlphabetType(self):
	""" return the type of emission domain 
	    XXX should call the method in HMMAlphabet
	"""
	return int
    
    def ClassType(self):
	pass
    
    def DistributionType(self):
	pass
    
    def buildMatrices(self):    
	""" return A, B, pi """
	pi = []
	B  = []
	A  = []
	nstates = len(self.state.keys())
	orders = {}
	k = 0 # C style index
	for s in self.state.values(): # ordering from XML
	    orders[s.index] = k
	    k = k + 1
	    
	for s in self.state.values(): # a list of indices
	    pi.append(s.initial)
	    
	    size = self.hmmAlphabet.size() # XXX first order only
	    if size != len(s.emissions):
		raise ValueError
	    else:
		B.append(s.emissions) # emission
	    
	    # transition probability
	    v = s.index
	    outprobs = [0.0] * nstates
	    for outid in self.G.OutNeighbors(v)[:]:
		myorder = orders[outid]
		outprobs[myorder] = self.G.edgeWeights[0][(v,outid)]
	    A.append(outprobs)
	return [A, B, pi]

    def OpenXML(self, fileName):
        dom = xml.dom.minidom.parse(fileName)
        if dom.documentElement.tagName == "ghmm":
            sys.stderr.write("Do not support ghmm format")
            raise FormatError
            dom.unlink()
            #self.DocumentName = "ghmm"
            #ghmmdom  = dom
            #ghmml = GHMMXML()
            #dom   = ghmml.GraphMLDOM(ghmmdom)
            #ghmmdom.unlink()
        else:
            assert dom.documentElement.tagName == "graphml"   
	    self.fromDOM(dom)
	    # dom.unlink()

    def WriteXML(self, fileName):
        try:
            self.modelCheck()   # raise exceptions here
            doc = xml.dom.minidom.Document()
            self.toDOM(doc, doc)
            file = open(fileName, 'w')
            # xml.dom.ext.PrettyPrint(doc, file)        
            file.write(toprettyxml(doc)) # problem with white spaces
            file.close()
            doc.unlink()
        except HMMEdError:
            print "HMMEdError: No file was written due to errors in the model."
            
    def WriteGHMM(self, fileName):
	self.modelCheck()   # raise exceptions here
        doc = xml.dom.minidom.Document()
        ghmm = doc.createElement("ghmm")
        doc.appendChild(ghmm)
        self.toGHMM(doc, ghmm)
        file = open(fileName, 'w')
        # xml.dom.ext.PrettyPrint(doc, file)        
        file.write(toprettyxml(doc)) # problem with white spaces
        file.close()
        doc.unlink()
        
    def SaveAs(self, fileName):
        if ( self.DocumentName == "graphml" ):
            self.WriteXML(fileName)
        else:
            self.WriteGHMM(fileName)
            
    def SaveAsGHMM(self, fileName):
        self.WriteGHMM(fileName)

            
        
################################################################################
if __name__ == '__main__':

    hmmobj = HMM()
    hmmobj.OpenXML(sys.argv[1])
    hmmobj.WriteXML("utz.xml")
    

