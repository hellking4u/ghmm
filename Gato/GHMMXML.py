#!/usr/bin/env python
#
#

#-----------------------------------
from sys import *
from xml.dom.minidom import *
import copy
import string

from HMMXML import *

from Graph import *
from Embedder import CircularCoords, RandomCoords


import whrandom
from math import log10

        
#------------------------------------------------------


class GHMMXML :
        
    def __init__(self):
        self.__adict = {}
        self.graphML  = GraphML()
        
    # XML Parsing
    # Ignore the "sequence" element

    def generateKeys(self):
        self.graphML.desc['initial']   = None
        self.graphML.type['initial']   = 'float'
        self.graphML.domain['initial'] = 'node'

        #self.graphML.desc['class']   = None
        #self.graphML.type['class']   = 'string'
        #self.graphML.domain['class'] = 'node'

        self.graphML.desc['label']   = None
        self.graphML.type['label']   = 'string'
        self.graphML.domain['label'] = 'node'

        #self.graphML.desc['order']   = None
        #self.graphML.type['order']   = 'int'
        #self.graphML.domain['order'] = 'node'

        self.graphML.desc['prob']   = None
        self.graphML.type['prob']   = 'float'
        self.graphML.domain['prob'] = 'edge'

        self.graphML.desc['emissions']   = None
        self.graphML.type['emissions']   = 'HigherDiscreteProbDist'
        self.graphML.domain['emissions'] = 'node'

        #
        # Default Visual representation
        #
        self.graphML.desc['npaint']  = 'npaint'
        self.graphML.type['npaint']  = 'paint'
        self.graphML.domain['npaint']= 'node'

        self.graphML.desc['epaint']  = 'npaint'
        self.graphML.type['epaint']  = 'paint'
        self.graphML.domain['epaint']= 'edge'

        self.graphML.desc['ngeom']  = 'npaint'
        self.graphML.type['ngeom']  = 'point'
        self.graphML.domain['ngeom']= 'node'

        self.graphML.desc['egeom']  = 'npaint'
        self.graphML.type['egeom']  = 'line'
        self.graphML.domain['egeom']= 'edge'

        # Add attributes to prototypes
        for keyName in self.graphML.desc.keys():
            keyValue  = None
            keyDomain = self.graphML.domain[keyName]
            keyType   = self.graphML.type[keyName]
            
            if not keyType in self.graphML.dataFactory.Types():
                keyValue = Empty()
                if keyName == 'epaint':
                    keyValue.__dict__['red']  = '0'
                    keyValue.__dict__['blue'] = '0'
                    keyValue.__dict__['green'] = '0'
                    keyValue.__dict__['style'] = 'solid'

                if keyName == 'npaint':
                    keyValue.__dict__['red']   = '255'
                    keyValue.__dict__['blue']  = '220'
                    keyValue.__dict__['green'] = '220'

                if keyName == 'ngeom':
                    keyValue.__dict__['shape']  = 'circle'
                    keyValue.__dict__['width']  = '50'
                    keyValue.__dict__['height'] = '50'

                if keyName == 'egeom':
                   keyValue.__dict__['shape'] = 'line'
                   keyValue.__dict__['width'] = '1'     
                    
            if keyDomain == 'all' or keyDomain == 'node':
                self.graphML.nodeProto.__dict__[keyName] = keyValue
           
            if keyDomain == 'all' or keyDomain == 'edge':
                self.graphML.edgeProto.__dict__[keyName] = keyValue
           
            if keyDomain == 'graph':
                self.graphML.graphProto.__dict__[keyName] = keyValue
        

    #-------------------------------------------------------
    
    def handleGhmmL(self, ghmml):
        self.generateKeys()

        hmms = ghmml.getElementsByTagName("hmm")
        self.handleHMM(hmms[0])  # There should only be one HMM 

    #-------------------------------------------------------
 
    def handleHMM(self, hmm):

        # required elements are alphabet, state
        # Should check for it

        self.graphML.graph = self.graphML.NewGraph()
            
        alphabetNode = hmm.getElementsByTagName("alphabet")
        assert len(alphabetNode) == 1  # There should only be one HMMAlphabet per HMM

        self.handleAlphabetNode(hmm)
        
        stateNodes   = hmm.getElementsByTagName("state")
        for state in stateNodes:
            self.handleStateNode(state)
        
        edgeNodes    = hmm.getElementsByTagName("transition")
        for edge in edgeNodes:
            self.handleEdge(edge)

        # Assign coordinates to vertices and edges
        G = self.generateRandomGraph()

    #-------------------------------------------------------
    
    def generateRandomGraph(self):
        G = Graph()
        G.directed=1  # directed graph
        
        for state in self.graphML.graph.nodes:
            vid = G.AddVertex()
            print state.id + '=' + '%s' % vid
            self.__adict[state.id] = G.vertices[vid-1]

        for e in self.graphML.graph.edges:
            print 'Edge:' + e.source + ' ' + e.target
            #a = G.vertices[self.__adict[e.source]]
            #b = G.vertices[self.__adict[e.target]]
            #if a != b:
            #    G.AddEdge(a,b)

        #   layout
        if RandomCoords(G):
            for state in self.graphML.graph.nodes:
                vid = self.__adict[state.id]
                state.__dict__['ngeom'].__dict__['x'] = G.xCoord[vid]
                state.__dict__['ngeom'].__dict__['y'] = G.yCoord[vid] 
            
        return G


    #-------------------------------------------------------

    
    def handleStateNode(self, stateNode):

        k         = stateNode.attributes['id'].nodeValue.strip()
        nNode     = self.graphML.AddNode(k)
        
        initials  = stateNode.getElementsByTagName("initial")
        for init in initials:
            # print getText(init.childNodes)
            dataInitial = getText(init.childNodes).strip()
            print dataInitial
            nNode.__dict__['initial'] = self.graphML.dataFactory(self.graphML.type['initial'], dataInitial)

            
        emissions = stateNode.getElementsByTagName("emission")
        for emiss in emissions:
            strEmission = getText( emiss.childNodes )
            strEmission = string.strip(strEmission)
            tmpSeq      = strEmission.splitlines(strEmission.count('\n'))
            print tmpSeq

            dataEmission = repr( map( (lambda x: float(x.strip())), tmpSeq)) 
            dataEmission = dataEmission[1:(len(dataEmission)-2)] # remove leading and trailing '[' ']'
            
            nNode.__dict__['emissions'] = self.graphML.dataFactory(self.graphML.type['emissions'], dataEmission)

        dataLabel = stateNode.attributes['id'].nodeValue.strip()
        nNode.__dict__['label']  = self.graphML.dataFactory(self.graphML.type['label'], dataLabel)
        # nNode.__dict__['class']  = self.graphML.dataFactory(self.graphML.type['class'], '0')

    #-------------------------------------------------------
    
    def handleEdge(self, edgeNode):
        edgeSource = edgeNode.attributes['source'].nodeValue
        edgeTarget = edgeNode.attributes['target'].nodeValue
        e = self.graphML.AddEdge(edgeSource, edgeTarget)

        for prob in edgeNode.getElementsByTagName("prob"):
            probValue = getText(prob.childNodes).strip()
            e.__dict__['prob'] = self.graphML.dataFactory(self.graphML.type['prob'], probValue)
        
    #-------------------------------------------------------
    
    def handleAlphabetNode(self, HMMNode):
        symbols = HMMNode.getElementsByTagName("symbol")
        
        self.graphML.hmmAlphabet      = DiscreteAlphabet()
        self.graphML.hmmAlphabet.low  = 0 
        self.graphML.hmmAlphabet.high = len(symbols)

        for i in range(len(symbols)):
            symbolRep = str(symbols[i].getAttribute("id"))
            self.graphML.hmmAlphabet.map[i] = symbolRep
            print symbols[i].attributes['id'].nodeValue

    #-------------------------------------------------------
    
    def __WriteData(self, doc, e, object, keys):
        #
        # skill visual element (npaint or epaint)
        # we use the global instead
        #
        print doc, e, object, keys
        for k in keys:
            if k == 'npaint' or k == 'epaint' :
                print ''
            else:
                elem = doc.createElement('data')
                elem.setAttribute('key', k)
                if k == 'ngeom':
                    poselem = doc.createElement('pos')
                    poselem.setAttribute('x',"%s" % object.__dict__[k].__dict__['x'])
                    poselem.setAttribute('y',"%s" % object.__dict__[k].__dict__['y'])
                    elem.appendChild(poselem)
                elif k == 'emissions':
                    strVal = ''
                    for p in object.__dict__[k]:
                        strVal = strVal + "%s" % p + ', '
                    strVal.strip()
                    strVal = strVal[0:len(strVal)-2]
                    contents = doc.createTextNode( strVal )
                    elem.appendChild(contents)                    
                else:
                    contents = doc.createTextNode("%s" % object.__dict__[k])
                    elem.appendChild(contents)

                e.appendChild(elem) # append this element to the DOM tree


    #-------------------------------------------------------

    def __FillVisualElement(self, doc, elem, key, object):

        if key == 'npaint' or key == 'epaint':
            paintElem = doc.createElement('paint')
            for k in ['red', 'green', 'blue']:
                paintElem.setAttribute( k, "%s" % object.__dict__[k] )
            elem.appendChild(paintElem)
        elif key == 'ngeom':
            pointElem = doc.createElement('point')
            for k in ['shape', 'width', 'height']:
                pointElem.setAttribute( k, "%s" % object.__dict__[k] )
            elem.appendChild(pointElem)
        elif key == 'egeom':
            lineElem = doc.createElement('line')
            for k in ['shape', 'width']:
                lineElem.setAttribute( k, "%s" % object.__dict__[k] )
            elem.appendChild(lineElem)
            
    #-------------------------------------------------------
    
    def GraphMLDOM(self, XMLDOM):
         self.handleGhmmL(XMLDOM)
         doc = Document()
         graphml = doc.createElement("graphml")
         doc.appendChild(graphml)
         
         node_keys = []
         edge_keys = []

         gml = self.graphML

         # hmm:alphabet
         elem = doc.createElement("hmm:alphabet")
         elem.setAttribute('hmm:type',"discrete")
         elem.setAttribute('hmm:low', str(self.graphML.hmmAlphabet.low))
         elem.setAttribute('hmm:high', str(self.graphML.hmmAlphabet.high))
         mapElem = doc.createElement('map')
         for i in range(len(self.graphML.hmmAlphabet.map)):
             symelem = doc.createElement('symbol')
             symelem.setAttribute('code', str(i))
             contents = doc.createTextNode(self.graphML.hmmAlphabet.map[i])
             symelem.appendChild(contents)
             mapElem.appendChild(symelem)
         elem.appendChild(mapElem)
         graphml.appendChild(elem)

         # hmm:class
         elem = doc.createElement("hmm:class")
         elem.setAttribute('hmm:low', '0')
         elem.setAttribute('hmm:high', '1')

         mapElem = doc.createElement("map")
         symbol  = doc.createElement("symbol")
         symbol.setAttribute('code', '0')
         symbol.setAttribute('desc', 'Simple')

         contents = doc.createTextNode("N")

         symbol.appendChild(contents)
         mapElem.appendChild(symbol)
         elem.appendChild(mapElem)
         graphml.appendChild(elem)
                                
         for k in gml.desc.keys():
             elem = doc.createElement("key")
             elem.setAttribute('id', k)
             
             if gml.domain[k] is not None:
                 elem.setAttribute('for', gml.domain[k] )          
                 elem.setAttribute('gd:type', gml.type[k])
                 graphml.appendChild(elem)
                 
             if gml.domain[k] == 'node' or gml.domain[k] == 'all':
                 node_keys.append(k)
                 self.__FillVisualElement(doc, elem, k, self.graphML.nodeProto.__dict__[k])

             if gml.domain[k] == 'edge' or gml.domain[k] == 'all':
                 edge_keys.append(k)
                 self.__FillVisualElement(doc, elem, k, self.graphML.edgeProto.__dict__[k])

                 
         graphelem = doc.createElement("graph")
         graphml.appendChild(graphelem)

                                          
         for n in gml.graph.nodes:
             elem = doc.createElement("node")        
             elem.setAttribute('id', n.id)
             self.__WriteData(doc,elem,n,node_keys)
             graphelem.appendChild(elem)

         for e in gml.graph.edges:
             elem = doc.createElement("edge")        
             elem.setAttribute('source', e.source)
             elem.setAttribute('target', e.target)
             self.__WriteData(doc,elem,e,edge_keys)
             graphelem.appendChild(elem)
    
         return doc
     
#
# Copied from HMMXML.py
#


def __WriteXML(gml):
    doc = Document()
    graphml = doc.createElement("graphml")
    doc.appendChild(graphml)

    node_keys = []
    edge_keys = []

    for k in gml.desc.keys():
        elem = doc.createElement("key")
        elem.setAttribute('id', k)
        if gml.domain[k] is not None:
            elem.setAttribute('for', gml.domain[k] )          
        elem.setAttribute('gd:type', gml.type[k])
        graphml.appendChild(elem)

        if gml.domain[k] == 'node' or gml.domain[k] == 'all':
            node_keys.append(k)
        if gml.domain[k] == 'edge' or gml.domain[k] == 'all':
            edge_keys.append(k)

    graphelem = doc.createElement("graph")
    graphml.appendChild(graphelem)

    for n in gml.graph.nodes:
        elem = doc.createElement("node")        
        elem.setAttribute('id', n.id)
        __WriteData(doc,elem,n,node_keys)
        graphelem.appendChild(elem)

    for e in gml.graph.edges:
        elem = doc.createElement("edge")        
        elem.setAttribute('source', e.source)
        elem.setAttribute('target', e.target)
        __WriteData(doc,elem,e,edge_keys)
        graphelem.appendChild(elem)
    
    print doc.toprettyxml()

#
#
# Main 
#
#
if __name__ == '__main__':
    dom = parse(argv[1])
    assert dom.documentElement.tagName == "ghmm"   

    ghmml = GHMMXML()
    gdom  = ghmml.GraphMLDOM(dom)
    print gdom.toprettyxml()
     
    dom.unlink()
