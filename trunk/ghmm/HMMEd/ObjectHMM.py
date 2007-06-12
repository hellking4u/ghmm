#!/usr/bin/env python
################################################################################
#
#       This file is part of Gato (Graph Algorithm Toolbox) 
#
#	file:   ObjectHMM.py
#	author: Janne Grunau
#
#       Copyright (C) 1998-2002, Alexander Schliep
#                                   
#       Contact: schliep@molgen.mpg.de
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
#       This file is version $Revision: 1.29 $ 
#                       from $Date: 2005/02/22 11:12:56 $
#             last change by $Author: schliep $.
#
################################################################################

from Gato.ObjectGraph import *
from Gato.EditObjectAttributesDialog import *
from Gato.MapEditor import NamedCollectionEditor
from Gato import ProbEditorBasics, ProbEditorDialogs, ProbEditorContinuous

import Tkinter

import ghmmwrapper, ghmmhelper, HMMEditor

import copy

class DiscreteHMMAlphabet:
    def __init__(self, names = [], description = "alphabet_1"):
        self.id = description
        self.name = {}
        self.name2code = {}
        for i,name in enumerate(names):
            self.name[i]         = name
            self.name2code[name] = i

    def size(self):
        return len(self.name.keys())

    def edit(self, master, name):
        pass

    def editDialog(self, master, hmm):
        self.hmm = hmm
        editor = NamedCollectionEditor(master, self)
        self.hmm = None

    def names(self):
        return self.name.values()

    def add(self, name):
        key = len(self.name)
        self.name[key] = name
        self.name2code[name] = key
        # update all states emissions
        for state in self.hmm.vertices.values():
            state.emission.grow()
        # also update background distributions if needed
        if self.hmm.modelType & ghmmwrapper.kBackgroundDistributions:
            self.hmm.backgroundDistributions.grow()

    def delete(self, name):
        key = self.name2code[name]
        for i in range(key, len(self.name)-1):
            self.name[i] = self.name[i+1]
            self.name2code[self.name[i]] = i

        del self.name2code[name]
        del self.name[len(self.name)-1]

        # update all states emissions
        for state in self.hmm.vertices.values():
            state.emission.shrink(key)
        # also update background distributions if needed
        if self.hmm.modelType & ghmmwrapper.kBackgroundDistributions:
            self.hmm.backgroundDistributions.shrink()
        # and tie groups
        if self.hmm.modelType & ghmmwrapper.kTiedEmissions:
            self.hmm.tie_groups.shrink()

    def GetKeys(self):
        return self.name.keys()

    def GetKey(self, symbol):
        return self.name2code[symbol]
    
    def GetSymbol(self, key):
        return self.name[key]

    def ReadCAlphabet(self, calphabet):
        self.id = calphabet.description
        if self.id is None:
            self.id = str(calphabet.size)
        for i in xrange(0, calphabet.size):
            name = calphabet.getSymbol(i)
            self.name[i] = name
            self.name2code[name] = i

    def WriteCAlphabet(self):
        calphabet = ghmmwrapper.ghmm_alphabet(len(self.name), self.id )
        #calphabet.description = self.id 
        for i in xrange(len(self.name)):
            calphabet.setSymbol(i, self.name[i])

        return calphabet


class DiscreteHMMBackground:
    def __init__(self, eclass):
        self.EmissionClass = eclass
        self.nextKey = 1
        self.val2pop = {0:"no background"}
        self.name = {}
        self.name2code = {}
        self.values = {}
        self.hmm = None

    def size(self):
        return len(self.name.keys())

    def edit(self, master, name):
        self.values[self.name2code[name]].edit(master, "backgound distribution \"%s\""%name)

    def editDialog(self, master, hmm):
        self.hmm = hmm
        editor = NamedCollectionEditor(master, self)
        self.hmm = None

    def names(self):
        return self.name.values()

    def add(self, name, alphabet=None):
        key = self.nextKey
        self.nextKey += 1
        
        if self.hmm is not None:
            e = self.EmissionClass(self.hmm.alphabet)
        elif alphabet is not None:
            e = self.EmissionClass(alphabet)
        else:
            e = self.EmissionClass()

        self.name[key] = name
        self.name2code[name] = key
        self.values[key] = e
        self.val2pop[key] = name

    def delete(self, name):
        key = self.name2code[name]
        del self.name[key]
        del self.name2code[name]
        del self.values[key]
        del self.val2pop[key]

    def grow(self):
        for emission in self.values.values():
            emission.grow()

    def shrink(self):
        for emission in self.values.values():
            emission.shrink()

    def getOrders(self):
        keys = self.name.keys()
        keys.sort()
        return [self.values[k].order for k in keys]

    def getWeights(self):
        keys = self.name.keys()
        keys.sort()
        return [self.values[k].weights for k in keys]

    def getNames(self):
        keys = self.name.keys()
        keys.sort()
        return [self.name[k] for k in keys]

    def ReadCBackground(self, alphabet, bp):
        number  = bp.n
        M       = bp.m
        for i in xrange(number):
            key = self.nextKey
            self.nextKey += 1

            e = self.EmissionClass(alphabet)
            e.order   = bp.getOrder(i)
            e.weights = ghmmhelper.double_array2list(bp.getWeights(i), M**(e.order+1))

            name = bp.getName(i)

            self.name[key] = name
            self.name2code[name] = key
            self.values[key] = e
            self.val2pop[key] = name

class TieGroups(DiscreteHMMBackground):

    def __init__(self, eclass):
        self.EmissionClass = eclass
        self.nextKey = 1
        self.val2pop = {0:"untied"}
        self.name = {}
        self.name2code = {}
        self.values = {}
        self.hmm = None

    def edit(self, master, name):
        self.values[self.name2code[name]].edit(master, "Emissions of tie group \"%s\""%name)

    def editEmissions(self, master, id):
        self.values[id].edit(master, "Emissions of tie group \"%s\""%self.name[id])

    def ReadCBackground(self, alphabet, bp):
        pass


class UniformDensity(ProbEditorContinuous.box_function):
    def getParameters(self):
        return (self.stop, self.start, 0, ghmmwrapper.uniform)

class NormalDensity(ProbEditorContinuous.gauss_function):
    def getParameters(self):
        return (self.mu, self.sigma, 0, ghmmwrapper.normal)

class NormalDensityTruncRight(ProbEditorContinuous.gauss_tail_function_right):
    def getParameters(self):
        return (self.mu, self.sigma, self.tail, ghmmwrapper.normal_right)
    
class NormalDensityTruncLeft(ProbEditorContinuous.gauss_tail_function_left):
    def getParameters(self):
        return (self.mu, self.sigma, self.tail, ghmmwrapper.normal_left)

class Emission:
    def __init__(self):
        pass

    def edit(self, master):
        pass

    def writeParameters(self, cstate):
        pass

    def ReadCState(self, cstate, M):
        pass

class DiscreteEmission(Emission):
    def __init__(self, alphabet):
        Emission.__init__(self)
        self.alphabet = alphabet
        if self.alphabet.size() > 0:
            self.weights = [1.0 / self.alphabet.size()] * self.alphabet.size()
        else:
            self.weights = []
        self.order = 0

    def grow(self):
        s = float(self.alphabet.size()-1)
        self.weights = [(x*s) for x in self.weights] + [1.0]
        self.weights = [x/sum(self.weights) for x in self.weights]
        
    def shrink(self, index):
        del self.weights[index]
        s = sum(self.weights)
        if s > 0.0:
            self.weights = [x / s for x in self.weights]
        elif self.alphabet.size() > 0:
            self.weights = [1.0 / self.alphabet.size()] * self.alphabet.size()
        else:
            self.weights = []

    def edit(self, master, description):
        transition_probabilities = ProbEditorBasics.ProbDict({})
        for key in self.alphabet.GetKeys():
            weight = self.weights[key]
            label  = self.alphabet.GetSymbol(key)
            transition_probabilities.update({label:weight})
            
        if transition_probabilities.sum == 0:
            key_list = transition_probabilities.keys()
            for key in key_list:
                transition_probabilities[key]=1.0/len(key_list)
        e = ProbEditorBasics.emission_data(transition_probabilities)
        d = ProbEditorDialogs.emission_dialog(master, e, description)

        if d.success():
            # write back normalized probabilities
            for label in transition_probabilities.keys():
                key = self.alphabet.GetKey(label)
                self.weights[key] = transition_probabilities[label]/transition_probabilities.sum

    def writeParameters(self, cstate):
        cstate.b = ghmmhelper.list2double_array(self.weights)

    def ReadCState(self, cstate, M):
        self.weights = ghmmhelper.double_array2list(cstate.b, M)


class DiscreteHigherOrderEmission(DiscreteEmission):
    def __init__(self, alphabet, order=0):
        self.alphabet = alphabet
        self.order    = ValidatingInt(order)
        M = self.alphabet.size()
        if M >  0:
            self.weights = [1.0 / M] * M**(order+1)
        else:
            self.weights = []

    def grow(self):
        print "not implemented"
        
    def shrink(self, index):
        print "not implemented"

    def edit(self, master, description):
        if self.order > 0:
            print "not implemented"
        else:
            DiscreteEmission.edit(self, master, description)

    def ChangeOrder(self, neworder):
        M = self.alphabet.size()
        if neworder < order:
            self.weights = self.weights[0:M**(neworder+1)]
        elif neworder > order:
            self.weights += [1.0 / M] * (M**(neworder+1) - M**neworder)

    def ReadCState(self, cstate, M):
        self.order   = ValidatingInt(cstate.order)
        self.weights = ghmmhelper.double_array2list(cstate.b, M**(self.order+1))


class ContinuousEmission(Emission):
    def __init__(self):
        Emission.__init__(self)
        self.weights  = [1.0]
        self.plotList = []

    def edit(self, master, stateId):
        if len(self.plotList) == 0:
            self.plotList.append(NormalDensity(mu=0,sigma=1,color="LightGreen"))

        tmp = [copy.copy(x) for x in self.plotList]
        top = Tkinter.Toplevel(master)
        d = HMMEditor.ContinuousEmissionEditor(top, tmp)
        d.pack(expand=1,fill=Tkinter.BOTH)
        top.withdraw()
        top.title( ("State %d emissions" % stateId))        
        top.update_idletasks()
        top.deiconify()
        top.wait_window(top)
        
        if d.success():
            self.plotList = d.plot_list
            self.weights  = [ d.dict[str(i)] for i in xrange(1,len(self.plotList)+1)]

    def writeParameters(self, cstate):
        cstate.M = len(self.plotList)
        
        muL=[]; uL=[]; aL=[]; density_tL=[]; cL=[]
        for d in xrange(len(self.plotList)): 
            (mu, u, a, density_t) = self.plotList[d].getParameters()
            muL.append(mu); uL.append(u); aL.append(a); density_tL.append(density_t);

        # write parameters in the state
        density_array  = ghmmwrapper.density_array_alloc(len(density_tL))
        for (i,density) in enumerate(density_tL):
            ghmmwrapper.density_array_setitem(density_array, i, density)
        cstate.density = density_array

        cstate.mue = ghmmhelper.list2double_array(muL)
        cstate.u   = ghmmhelper.list2double_array(uL)
        cstate.a   = ghmmhelper.list2double_array(aL)
        cstate.c   = ghmmhelper.list2double_array(self.weights)

        # editor doesn't supports fixed mixture components
        cstate.mixture_fix = ghmmhelper.list2int_array([0] * len(self.plotList))

    def ReadCState(self, cstate, M):
        M = cstate.M
        for i in xrange(M):
            if cstate.getDensity(i) == ghmmwrapper.normal:
                self.plotList.append(NormalDensity())
            elif cstate.getDensity(i) == ghmmwrapper.normal_left:
                self.plotList.append(NormalDensityTruncLeft())
            elif cstate.getDensity(i) == ghmmwrapper.normal_right:
                self.plotList.append(NormalDensityTruncRight())
            elif cstate.getDensity(i) == ghmmwrapper.uniform:
                self.plotList.append(UniformDensity())


class DiscretePairEmission(Emission):
    pass


class State(VertexObject):

    def __init__(self, emission=Emission(), hmm=None):
        VertexObject.__init__(self)
        self.num = -1   # continuous id (0..Order()-1) used for writing transitions
        self.editableAttr = {'labeling':"Name", 'initial':"Initial Probability", 'fixed':"fixed emissions"}
        self.initial = Probability()
        self.emission = emission
        self.itsHMM = hmm
        self.fixed = ValidatingBool(False)

    def update(self):
        pass

    def normalize(self):
        weights  = [e.GetWeight() for e in self.outEdges]
        s = float(sum(weights))
        for i,edge in enumerate(self.outEdges):
            edge.SetWeight(weights[i]/s)

    def editProperties(self, parent, attributes = None):
        self.update()
        self.desc = ('Properties of State %d (%s)' % (self.id, self.labeling))
        if self.itsHMM.modelType & ghmmwrapper.kHigherOrderEmissions:
            self.order = self.emission.order
            self.editableAttr['order'] = "Order"
        if attributes == None:
            editBox = EditObjectAttributesDialog(parent, self, self.editableAttr)
        else:
            editableAttr = {}
            for attr in attributes:
                editableAttr[attr] = self.editableAttr[attr]
            editBox = EditObjectAttributesDialog(parent, self, editableAttr)
        if self.itsHMM.modelType & ghmmwrapper.kHigherOrderEmissions:
            self.emission.ChangeOrder(self.order)
            del self.order
            del self.editableAttr['order']

    def editEmissions(self, master):
        self.emission.edit(master, self.id)

    def WriteCState(self, cstate):
        cstate.pi         = self.initial
        cstate.in_states  = len(self.inEdges)
        cstate.out_states = len(self.outEdges)

        cstate.desc      = str(self.labeling)
        cstate.xPosition = int(self.embedding.x)
        cstate.yPosition = int(self.embedding.y)

        cstate.fix = self.fixed

        self.WriteTransitions(cstate)
        self.emission.writeParameters(cstate)

    def WriteTransitions(self, cstate):
        inID = [edge.tail.num for edge in self.inEdges]
        inA  = [edge.GetEdgeWeight(0) for edge in self.inEdges]
        cstate.in_id = ghmmhelper.list2int_array(inID)
        cstate.in_a  = ghmmhelper.list2double_array(inA)

        outID = [edge.head.num for edge in self.outEdges]
        outA = [edge.GetEdgeWeight(0) for edge in self.outEdges]
        cstate.out_id = ghmmhelper.list2int_array(outID)
        cstate.out_a  = ghmmhelper.list2double_array(outA)

    def ReadCState(self, cmodel, cstate, i):
        self.initial  = ValidatingFloat(cstate.pi)
        self.fixed    = ValidatingBool(cstate.fix)
        self.labeling = ValidatingString(cstate.desc)

        if self.itsHMM is not None:
            self.itsHMM.SetEmbedding(self.id, cstate.xPosition, cstate.yPosition)

        self.emission.ReadCState(cstate, cmodel.M)


class ContinuousState(State):

    def WriteTransitions(self, cstate):
        inID = [edge.tail.id for edge in self.inEdges]
        inA  = [[edge.GetEdgeWeight(0) for edge in self.inEdges]]
        cstate.in_id = ghmmhelper.list2int_array(inID)
        (mat, lens)  = ghmmhelper.list2double_matrix(inA)
        cstate.in_a  = mat
        
        outID = [edge.head.id for edge in self.outEdges]
        outA  = [[edge.GetEdgeWeight(0) for edge in self.outEdges]]
        cstate.out_id = ghmmhelper.list2int_array(outID)
        (mat, lens)   = ghmmhelper.list2double_matrix(outA)
        cstate.out_a  = mat


class SilentState(State):
    def __init__(self, emission=Emission(), hmm=None, silent=ValidatingBool(False)):
        State.__init__(self, emission, hmm)
        self.init(silent)

    def init(self, silent):
        self.editableAttr['silent'] = "Silent"
        self.silent = silent

    def editEmissions(self, master):
        if not self.silent:
            State.editEmissions(self, master)

    def ReadCState(self, cmodel, cstate, i):
        State.ReadCState(self, cmodel, cstate, i)
        self.silent = ValidatingBool(cmodel.getSilent(i))

    
class BackgroundState(State):
    def __init__(self, emission=Emission(), hmm=None):
        State.__init__(self, emission, hmm)
        self.init()

    def init(self):
        self.background = PopupableInt()
        self.background.setPopup(self.itsHMM.backgroundDistributions.val2pop)

    def update(self):
        if len(self.itsHMM.backgroundDistributions.names()) > 0:
            self.editableAttr['background'] = "Background distribution"
            self.background.setPopup(self.itsHMM.backgroundDistributions.val2pop)

    def ReadCState(self, cmodel, cstate, i):
        State.ReadCState(self, cmodel, cstate, i)
        self.update()
        self.background = PopupableInt(cmodel.getBackgroundID(i)+1)


class LabeledState(State):
    def __init__(self, emission=Emission(), hmm=None):
        State.__init__(self, emission, hmm)
        self.init()

    def init(self):
        self.label = PopupableInt()
        self.label.setPopup(self.itsHMM.label_alphabet.name)

    def update(self):
        if len(self.itsHMM.label_alphabet.names()) > 0:
            self.editableAttr['label'] = "State label"
            self.label.setPopup(self.itsHMM.label_alphabet.name)

    def ReadCState(self, cmodel, cstate, i):
        State.ReadCState(self, cmodel, cstate, i)
        self.update
        self.label = PopupableInt(cmodel.getStateLabel(i))


class TiedState(State):
    def __init__(self, emission=Emission(), hmm=None):
        State.__init__(self, emission, hmm)
        self.init()

    def init(self):
        self.tiedto = PopupableInt(0)
        self.tiedto.setPopup(self.itsHMM.tie_groups.val2pop)

    def update(self):
        if len(self.itsHMM.tie_groups.names()) > 0:
            self.editableAttr['tiedto'] = "Emissions tied to state"
            self.tiedto.setPopup(self.itsHMM.tie_groups.val2pop)

    def editEmissions(self, master):
        if self.tiedto > 0:
            self.itsHMM.tie_groups.editEmissions(master, self.tiedto)
        else:
            State.editEmissions(self, master)

    def ReadCState(self, cmodel, cstate, i):
        State.ReadCState(self, cmodel, cstate, i)
        tied = cmodel.getTiedTo(i)
        if tied == i:
            self.itsHMM.tie_groups.add("tiedto_state%d"%tied, self.itsHMM.alphabet)
        self.update()
        self.tiedto = PopupableInt(self.itsHMM.tie_groups.name2code["tiedto_state%d"%tied])


class SilentBackgroundState(SilentState, BackgroundState):
        
    def __init__(self, emission=Emission(), hmm=None, silent=ValidatingBool(False)):
        State.__init__(self, emission, hmm)
        SilentState.init(self, silent)
        BackgroundState.init(self)

    def update(self):
        BackgroundState.update(self)

    def ReadCState(self, cmodel, cstate, i):
        SilentState.ReadCState(self, cmodel, cstate, i)
        BackgroundState.ReadCState(self, cmodel, cstate, i)


class SilentLabeledState(SilentState, LabeledState):

    def __init__(self, emission=Emission(), hmm=None, silent=ValidatingBool(False)):
        State.__init__(self, emission, hmm)
        SilentState.init(self, silent)
        LabeledState.init(self)

    def update(self):
        LabeledState.update(self)

    def ReadCState(self, cmodel, cstate, i):
        SilentState.ReadCState(self, cmodel, cstate, i)
        LabeledState.ReadCState(self, cmodel, cstate, i)


class SilentTiedState(SilentState, TiedState):
        
    def __init__(self, emission=Emission(), hmm=None, silent=ValidatingBool(False)):
        State.__init__(self, emission, hmm)
        SilentState.init(self, silent)
        TiedState.init(self)

    def update(self):
        TiedState.update(self)

    def ReadCState(self, cmodel, cstate, i):
        SilentState.ReadCState(self, cmodel, cstate, i)
        TiedState.ReadCState(self, cmodel, cstate, i)


class LabeledBackgroundState(LabeledState, BackgroundState):
    def __init__(self, emission=Emission(), hmm=None):
        State.__init__(self, emission, hmm)
        LabeledState.init(self)
        BackgroundState.init(self)

    def update(self):
        BackgroundState.update(self)
        LabeledState.update(self)

    def ReadCState(self, cmodel, cstate, i):
        LabeledState.ReadCState(self, cmodel, cstate, i)
        BackgroundState.ReadCState(self, cmodel, cstate, i)


class LabeledTiedState(LabeledState, TiedState):
    def __init__(self, emission=Emission(), hmm=None):
        State.__init__(self, emission, hmm)
        LabeledState.init(self)
        TiedState.init(self)

    def update(self):
        LabeledState.update(self)
        TiedState.update(self)

    def ReadCState(self, cmodel, cstate, i):
        LabeledState.ReadCState(self, cmodel, cstate, i)
        TiedState.ReadCState(self, cmodel, cstate, i)


class TiedBackgroundState(TiedState, BackgroundState):
    def __init__(self, emission=Emission(), hmm=None):
        State.__init__(self, emission, hmm)
        TiedState.init(self)
        BackgroundState.init(self)

    def update(self):
        BackgroundState.update(self)
        TiedState.update(self)

    def ReadCState(self, cmodel, cstate):
        TiedState.ReadCState(self, cmodel, cstate, i)
        BackgroundState.ReadCState(self, cmodel, cstate, i)


class SilentLabeledBackgroundState(SilentState, LabeledState, BackgroundState):
    def __init__(self, emission=Emission(), hmm=None, silent=ValidatingBool(False)):
        State.__init__(self, emission, hmm)
        SilentState.init(self, silent)
        LabeledState.init(self)
        BackgroundState.init(self)

    def update(self):
        BackgroundState.update(self)
        LabeledState.update(self)

    def ReadCState(self, cmodel, cstate, i):
        SilentState.ReadCState(self, cmodel, cstate, i)
        LabeledState.ReadCState(self, cmodel, cstate, i)
        BackgroundState.ReadCState(self, cmodel, cstate, i)


class SilentLabeledTiedState(SilentState, LabeledState, TiedState):
    def __init__(self, emission=Emission(), hmm=None, silent=ValidatingBool(False)):
        State.__init__(self, emission, hmm)
        SilentState.init(self, silent)
        LabeledState.init(self)
        TiedState.init(self)
        
    def update(self):
        LabeledState.update(self)
        TiedState.update(self)

    def ReadCState(self, cmodel, cstate, i):
        SilentState.ReadCState(self, cmodel, cstate, i)
        LabeledState.ReadCState(self, cmodel, cstate, i)
        TiedState.ReadCState(self, cmodel, cstate, i)


class SilentTiedBackgroundState(SilentState, TiedState, BackgroundState):
    def __init__(self, emission=Emission(), hmm=None, silent=ValidatingBool(False)):
        State.__init__(self, emission, hmm)
        SilentState.init(self, silent)
        TiedState.init(self)
        BackgroundState.init(self)

    def update(self):
        BackgroundState.update(self)
        TiedState.update(self)

    def ReadCState(self, cmodel, cstate):
        SilentState.ReadCState(self, cmodel, cstate, i)
        TiedState.ReadCState(self, cmodel, cstate, i)
        BackgroundState.ReadCState(self, cmodel, cstate, i)


class LabeledTiedBackgroundState(LabeledState, TiedState, BackgroundState):
    def __init__(self, emission=Emission(), hmm=None):
        State.__init__(self, emission, hmm)
        LabeledState.init(self)
        TiedState.init(self)
        BackgroundState.init(self)

    def update(self):
        BackgroundState.update(self)
        LabeledState.update(self)
        TiedState.update(self)

    def ReadCState(self, cmodel, cstate, i):
        LabeledState.ReadCState(self, cmodel, cstate, i)
        TiedState.ReadCState(self, cmodel, cstate, i)
        BackgroundState.ReadCState(self, cmodel, cstate, i)


class SilentLabeledTiedBackgroundState(SilentState, LabeledState, TiedState, BackgroundState):
    def __init__(self, emission=Emission(), hmm=None, silent=ValidatingBool(False)):
        State.__init__(self, emission, hmm)
        SilentState.init(self, silent)
        LabeledState.init(self)
        TiedState.init(self)
        BackgroundState.init(self)

    def update(self):
        BackgroundState.update(self)
        LabeledState.update(self)
        TiedState.update(self)

    def ReadCState(self, cmodel, cstate, i):
        SilentState.ReadCState(self, cmodel, cstate, i)
        LabeledState.ReadCState(self, cmodel, cstate, i)
        TiedState.ReadCState(self, cmodel, cstate, i)
        BackgroundState.ReadCState(self, cmodel, cstate, i)


class Transition(EdgeObject):

    def __init__(self, tail, head):
        EdgeObject.__init__(self, tail, head)
        #self.weight = Probability()
        self.editableAttr = {'weight':"Weight"}

    def GetWeight(self):
        return self.GetEdgeWeight(0)

    def SetWeight(self, value):
        self.SetEdgeWeight(0, value)

    def edit(self, parent, attributes = None):
        if attributes == None:
            editBox = EditObjectAttributesDialog(parent, self, self.editableAttr)
        else:
            editableAttr = {}
            for attr in attributes:
                editableAttr[attr] = self.editableAttr[attr]
            editBox = EditObjectAttributesDialog(parent, self, editableAttr)

    def ReadCTransition(self, state, cos, i):
        assert (cos == 1)
        self.weight = state.getOutProb(i)


class SwitchedTransition(Transition):
    
    def __init__(self, tail, head, noTransitionClasses=2):
        Transition.__init__(self, tail, head)
        self.weight = [Probability()] * noTransitionClasses
        self.noClasses = noTransitionClasses

    def GetEdgeWeight(self,i):
        if i < self.noClasses:
            return self.weight[i]

    def SetEdgeWeight(self,i,value):
        if i < self.noClasses:
            self.weight[i] = value

    def ReadCTransition(self, state, cos, i):
        self.weight = ghmmhelper.double_array2list(state.out_a[i], cos)



class ObjectHMM(ObjectGraph):
    """
    """

    def __init__(self, stateClass, transitionClass, emissionClass=Emission(), alphabet=None, type=0):
        ObjectGraph.__init__(self, stateClass, transitionClass)
        self.simple    = 0
        self.euclidian = 0
        self.directed  = 1
        self.modelType = 0

        self.alphabet = alphabet
        self.emissionClass = emissionClass
        self.type = type

        self.vertices_ids = {0:'untied'}

        # editable attributes per EditPropertiesDialog
        # common properties:
        self.desc         = "New HMM properties"
        self.editableAttr = {'name':"Name"}
        
        self.alphatype          = PopupableInt()
        self.alphatype_val2pop  = {0:"binary", 1:"dice", 2:"DNA",
                                   3:"amino acids", 4:"custom"}
        self.alphatype.setPopup(self.alphatype_val2pop)
        self.labels             = ValidatingBool()
        self.background         = ValidatingBool()
        self.maxOrder           = DefaultedInt(0)
        self.maxOrder.setDefault(1, 0)
        self.name               = ValidatingString("")
        self.silent             = ValidatingBool()
        self.switching          = DefaultedInt(1)
        self.switching.setDefault(1, 1)
        self.tied               = ValidatingBool()
        
        # discrete emissions only properties
        if type == 0:
            self.editableAttr['tied']       = "Tied emissions"
            self.editableAttr['silent']     = "Silent states"
            self.editableAttr['maxOrder']   = "Higher order emissions"
            self.editableAttr['background'] = "Background distributions"
            self.editableAttr['labels']     = "State labels"
            self.editableAttr['alphatype']  = "Alphabet"
        # continuous emissions only properties
        elif type == 1:
            self.editableAttr['switching'] = "No. of transition Classes"
        else:
            print "invalid type"

    def AddVertex(self):
        """ Add an isolated vertex. Returns the id of the new vertex """
        if self.alphabet is not None:
            e = self.emissionClass(self.alphabet)
        else:
            e = self.emissionClass()
        v = self.vertexClass(e, self)
        v.id = self.GetNextVertexID()
        self.vertices[v.id] = v
        self.vertices_ids[v.id+1] = str(v.id)
        return v.id

    def DeleteVertex(self, v):
        del self.vertices_ids[v+1]
        ObjectGraph.DeleteVertex(self, v)

    def AddEdge(self,tail,head):
        out_edges = len(self.vertices[tail].outEdges)
        if out_edges > 0:
            weight = 1.0/out_edges
        else:
            weight = 1.0
        ObjectGraph.AddEdge(self,tail,head)
        edge = self.edges[tail,head]
        edge.SetWeight(weight)
        self.vertices[tail].normalize()

    def DeleteEdge(self,tail,head):
        ObjectGraph.DeleteEdge(self,tail,head)
        self.vertices[tail].normalize()

    def SetLabeling(self,v, value):
        self.vertices[v].labeling = ValidatingString(value)

    def edit(self, parent, attributes = None):
        if attributes == None:
            editBox = EditObjectAttributesDialog(parent, self, self.editableAttr)
        else:
            editableAttr = {}
            for attr in attributes:
                editableAttr[attr] = self.editableAttr[attr]
            editBox = EditObjectAttributesDialog(parent, self, editableAttr)

        mt = self.computeModelType()
        if mt > 0:
            self.initHMM(mt)
        else:
            print "invalid model type:", mt


    def computeModelType(self):
        modelType = 0
        if self.type == 0:
            modelType += ghmmwrapper.kDiscreteHMM
            if self.maxOrder > 0:
                modelType += ghmmwrapper.kHigherOrderEmissions
        elif self.type == 1:
            modelType += ghmmwrapper.kContinuousHMM
        elif self.type == 2:
            modelType += ghmmwrapper.kDiscreteHMM
            modelType += ghmmwrapper.kPairHMM
        else:
            print "invalid type:", self.type

        if self.switching > 1:
            modelType += ghmmwrapper.kTransitionClasses
            
        if self.tied:
            modelType += ghmmwrapper.kTiedEmissions

        if self.silent:
            modelType += ghmmwrapper.kSilentStates

        if self.background:
            modelType += ghmmwrapper.kBackgroundDistributions

        if self.labels:
            modelType += ghmmwrapper.kLabeledStates

        return modelType


    def initHMM(self, modelType):
        # set the right emission type
        if modelType & ghmmwrapper.kDiscreteHMM:
            if modelType & ghmmwrapper.kPairHMM:
                emissionClass = DiscretePairEmission
                # alphabet missing
            else:
                if modelType & ghmmwrapper.kHigherOrderEmissions:
                    emissionClass = DiscreteHigherOrderEmission
                else:
                    emissionClass = DiscreteEmission
                alphabet = self.initAlphabet()
        elif modelType & ghmmwrapper.kContinuousHMM:
            emissionClass = ContinuousEmission
            alphabet = None
        else:
            print "not a valid model type"

        # set the right transition type
        if modelType & ghmmwrapper.kTransitionClasses:
            edgeClass = SwitchedTransition
        else:
            edgeClass = Transition
            
        # masking unnecessary model type flags out
        mt = modelType
        if modelType & ghmmwrapper.kDiscreteHMM:
            mt -= ghmmwrapper.kDiscreteHMM
        if modelType & ghmmwrapper.kContinuousHMM:
            mt -= ghmmwrapper.kContinuousHMM
        if modelType & ghmmwrapper.kPairHMM:
            mt -= ghmmwrapper.kPairHMM
        if modelType & (ghmmwrapper.kHigherOrderEmissions):
            mt -= ghmmwrapper.kHigherOrderEmissions

        # setting the right vertex type
        if mt == (ghmmwrapper.kSilentStates):
            vertexClass = SilentState
        elif mt == (ghmmwrapper.kTiedEmissions):
            vertexClass = TiedState
        elif mt == (ghmmwrapper.kBackgroundDistributions):
            vertexClass = BackgroundState
        elif mt == (ghmmwrapper.kLabeledStates):
            vertexClass = LabeledState
        # 2
        elif mt == (ghmmwrapper.kSilentStates + ghmmwrapper.kTiedEmissions):
            vertexClass = SilentTiedState
        elif mt == (ghmmwrapper.kSilentStates + ghmmwrapper.kLabeledStates):
            vertexClass = SilentLabeledState
        elif mt == (ghmmwrapper.kSilentStates + ghmmwrapper.kBackgroundDistributions):
            vertexClass = SilentBackgroundState
        elif mt == (ghmmwrapper.kTiedEmissions + ghmmwrapper.kBackgroundDistributions):
            vertexClass = TiedBackgroundState
        elif mt == (ghmmwrapper.kTiedEmissions + ghmmwrapper.kLabeledStates):
            vertexClass = LabeledTiedState
        elif mt == (ghmmwrapper.kLabeledStates + ghmmwrapper.kBackgroundDistributions):
            vertexClass = LabeledBackgroundState
        # 3
        elif mt == (ghmmwrapper.kSilentStates + ghmmwrapper.kTiedEmissions + ghmmwrapper.kBackgroundDistributions):
            vertexClass = SilentTiedBackgroundState
        elif mt == (ghmmwrapper.kSilentStates + ghmmwrapper.kTiedEmissions + ghmmwrapper.kLabeledStates):
            vertexClass = SilentLabeledTiedState
        elif mt == (ghmmwrapper.kSilentStates + ghmmwrapper.kLabeledStates + ghmmwrapper.kBackgroundDistributions):
            vertexClass = SilentLabeledBackgroundState
        elif mt == (ghmmwrapper.kTiedEmissions + ghmmwrapper.kLabeledStates + ghmmwrapper.kBackgroundDistributions):
            vertexClass = LabeledTiedBackgroundState
        # 4
        elif mt == (ghmmwrapper.kSilentStates + ghmmwrapper.kTiedEmissions + ghmmwrapper.kLabeledStates + ghmmwrapper.kBackgroundDistributions):
            vertexClass = SilentLabeledTiedBackgroundState
        else:
            vertexClass = (modelType & ghmmwrapper.kContinuousHMM) and ContinuousState or State

        # initialize state labels
        if mt & ghmmwrapper.kLabeledStates:
            self.label_alphabet = DiscreteHMMAlphabet(description = "state labels")

        # initialize background distributions
        if mt & ghmmwrapper.kBackgroundDistributions:
            self.backgroundDistributions = DiscreteHMMBackground(emissionClass)
            
        # initialize background distributions
        if mt & ghmmwrapper.kTiedEmissions:
            self.tie_groups = TieGroups(emissionClass)

        self.__init__(vertexClass, edgeClass, emissionClass, alphabet)
        self.modelType = modelType


    def initAlphabet(self):
        if self.alphatype == 0:
            return DiscreteHMMAlphabet(["0", "1"])
        elif self.alphatype == 1:
            return DiscreteHMMAlphabet(["1", "2", "3", "4", "5", "6"])
        elif self.alphatype == 2:
            return DiscreteHMMAlphabet(["A", "C", "G", "T"])
        elif self.alphatype == 3:
            return DiscreteHMMAlphabet(["ala", "arg", "asn", "asp", "asx",
                                        "cys", "glu", "gln", "glx", "gly",
                                        "his", "ile", "leu", "lys", "met",
                                        "phe", "pro", "ser", "thr", "try",
                                        "tyr", "val"])
        elif self.alphatype == 4:
            return DiscreteHMMAlphabet()
        else:
            print "invalid alphabet type"
            return None

    def openXML(self, filename="test.xml"):
        # simple check of file 
        filedata = ghmmwrapper.ghmm_xmlfile_parse(filename)
        if filedata == None:
            raise UnknownFileTypeException(filename)
        if filedata.noModels > 1:
            raise UnsupportedFileException(filename + "more than one HMM per file currently not supported")

        # initialize model and set auxiliary data accordingly to the model type
        self.initHMM(filedata.modelType)

        #cmodel = filedata.getModel(0)
        if self.modelType & ghmmwrapper.kContinuousHMM:
            cmodel = filedata.get_cmodel(0)
            cos = cmodel.cos
        elif self.modelType & ghmmwrapper.kDiscreteHMM:
            if self.modelType & ghmmwrapper.kPairHMM:
                cmodel = filedata.get_dpmodel(0)
                cos = 1
            elif self.modelType & ghmmwrapper.kTransitionClasses:
                cmodel = filedata.get_dsmodel(0)
                cos = cmodel.cos
            else:
                cmodel = filedata.get_dmodel(0)
                cos = 1

        # Add all states
        vdict = {}
        for i in xrange(cmodel.N):
            vdict[i] = self.AddVertex()

        # Add all transitions
        for i in xrange(cmodel.N):
            state = cmodel.getState(i)
            for j in xrange(state.out_states):
                outID = state.getOutState(j)
                tail = vdict[i]
                head = vdict[outID]
                self.AddEdge(tail, head)

        # Add alphabet if appropiate
        if self.modelType & ghmmwrapper.kDiscreteHMM:
            self.alphabet = DiscreteHMMAlphabet()
            self.alphabet.ReadCAlphabet(cmodel.alphabet)
        
        # Add label alphabet
        if self.modelType & ghmmwrapper.kLabeledStates:
            self.label_alphabet = DiscreteHMMAlphabet()
            self.label_alphabet.ReadCAlphabet(cmodel.label_alphabet)

        # Add background distributions if appropiate
        if self.modelType & ghmmwrapper.kBackgroundDistributions:
            self.backgroundDistributions = DiscreteHMMBackground(self.emissionClass)
            self.backgroundDistributions.ReadCBackground(self.alphabet, cmodel.bp)

        # Add switching functions if appropiate
        if self.modelType & ghmmwrapper.kTransitionClasses:
            print "TODO: transition classes???"

        # Set all states' values and set transition weights
        for i in xrange(cmodel.N):
            state = cmodel.getState(i)
            self.vertices[vdict[i]].ReadCState(cmodel, state, i)
            for j in xrange(state.out_states):
                outID = state.getOutState(j)
                tail = vdict[i]
                head = vdict[outID]
                self.edges[tail, head].ReadCTransition(state, cos, j)

    def writeXML(self, filename="test.xml"):

        if self.modelType & ghmmwrapper.kContinuousHMM:
            self.cmodel = ghmmwrapper.ghmm_cmodel()
            self.cmodel.s = ghmmwrapper.cstate_array_alloc(self.Order())
        elif self.modelType & ghmmwrapper.kDiscreteHMM:
            if self.modelType & ghmmwrapper.kPairHMM:
                self.cmodel = None
            elif self.modelType & ghmmwrapper.kTransitionClasses:
                self.cmodel = None
            else:
                self.cmodel   = ghmmwrapper.ghmm_dmodel()
                self.cmodel.s = ghmmwrapper.dstate_array_alloc(self.Order())
                self.cmodel.M = self.alphabet.size()
                self.cmodel.alphabet = self.alphabet.WriteCAlphabet()
                
        if self.modelType & ghmmwrapper.kTransitionClasses:
            self.cmodel.cos = maxcos()
        else:
            self.cmodel.cos = 1

        # sort state IDs
        sortedIDs = self.vertices.keys()
        sortedIDs.sort()

        # fill state property arrays according to the model type with default values
        self.cmodel.N = self.Order()

        # fill silent array
        if self.modelType & ghmmwrapper.kSilentStates:
            self.cmodel.silent = ghmmhelper.list2int_array([self.vertices[id].silent for id in sortedIDs])

        # fill tied to array
        if self.modelType & ghmmwrapper.kTiedEmissions:
            tied_list = [ghmmwrapper.kUntied] * self.Order()
            tieddict = {}
            # map python id to continious C array indeces
            for i, id in enumerate(sortedIDs):
                if self.vertices[id].tiedto > 0:
                    tiedto = self.vertices[id].tiedto-1
                    if tieddict.has_key(tiedto):
                        tieddict[tiedto].append(i)
                    else:
                        tieddict[tiedto] = [i]
            # tiedto has to be sorted, the first entry points to it self
            for k in tieddict.keys():
                temp = tieddict[k]
                temp.sort()
                first = temp[0]
                for index in temp:
                    tied_list[index] = first
            self.cmodel.tied_to = ghmmhelper.list2int_array(tied_list)

        # fill background id arrary
        if self.modelType & ghmmwrapper.kBackgroundDistributions:
            N = self.backgroundDistributions.size()
            M = self.alphabet.size()
            orders = ghmmhelper.list2int_array(self.backgroundDistributions.getOrders())
            (weights,lengths) = ghmmhelper.list2double_matrix(self.backgroundDistributions.getWeights())
            self.cmodel.bp = ghmmwrapper.ghmm_dbackground(N, M, orders, weights)
            for i,name in enumerate(self.backgroundDistributions.getNames()):
                self.cmodel.bp.setName(i, name)
            self.cmodel.background_id = ghmmhelper.list2int_array([(self.vertices[id].background-1) for id in sortedIDs])

        # fill higher order array
        if self.modelType & ghmmwrapper.kHigherOrderEmissions:
            self.cmodel.order = ghmmhelper.list2int_array([self.vertices[id].emission.order for id in sortedIDs])

        # fil label id array
        if self.modelType & ghmmwrapper.kLabeledStates:
            self.cmodel.label_alphabet = self.label_alphabet.WriteCAlphabet()
            self.cmodel.label = ghmmhelper.list2int_array([self.vertices[id].label for id in sortedIDs])

        self.cmodel.model_type = self.modelType

        # create each state
        initial_sum = 0.0
        for i, id in enumerate(sortedIDs):
            self.vertices[id].num = i
            initial_sum += self.vertices[id].initial

        if initial_sum < 1E-14:
            for id in sortedIDs:
                self.vertices[id].initial = 1.0
            initial_sum = float(self.Order())

        for i, id in enumerate(sortedIDs):
            cstate = self.cmodel.getState(i)
            self.vertices[id].initial /= initial_sum
            self.vertices[id].WriteCState(cstate)

        # write to file
        #writeFunction(filename, castFunction(self.cmodel), 1)
        self.cmodel.write_xml(filename)

        self.cmodel = None
