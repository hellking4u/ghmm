# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _ghmmwrapper

def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types



gsl_rng_init = _ghmmwrapper.gsl_rng_init

gsl_rng_timeseed = _ghmmwrapper.gsl_rng_timeseed

time_seed = _ghmmwrapper.time_seed

matrix_d_alloc = _ghmmwrapper.matrix_d_alloc

matrix_d_alloc_copy = _ghmmwrapper.matrix_d_alloc_copy

matrix_d_free = _ghmmwrapper.matrix_d_free

matrix_i_alloc = _ghmmwrapper.matrix_i_alloc

matrix_i_free = _ghmmwrapper.matrix_i_free

matrix_d_print = _ghmmwrapper.matrix_d_print
class sequence_t(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, sequence_t, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, sequence_t, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C sequence_t instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["seq"] = _ghmmwrapper.sequence_t_seq_set
    __swig_getmethods__["seq"] = _ghmmwrapper.sequence_t_seq_get
    if _newclass:seq = property(_ghmmwrapper.sequence_t_seq_get, _ghmmwrapper.sequence_t_seq_set)
    __swig_setmethods__["states"] = _ghmmwrapper.sequence_t_states_set
    __swig_getmethods__["states"] = _ghmmwrapper.sequence_t_states_get
    if _newclass:states = property(_ghmmwrapper.sequence_t_states_get, _ghmmwrapper.sequence_t_states_set)
    __swig_setmethods__["seq_len"] = _ghmmwrapper.sequence_t_seq_len_set
    __swig_getmethods__["seq_len"] = _ghmmwrapper.sequence_t_seq_len_get
    if _newclass:seq_len = property(_ghmmwrapper.sequence_t_seq_len_get, _ghmmwrapper.sequence_t_seq_len_set)
    __swig_setmethods__["seq_label"] = _ghmmwrapper.sequence_t_seq_label_set
    __swig_getmethods__["seq_label"] = _ghmmwrapper.sequence_t_seq_label_get
    if _newclass:seq_label = property(_ghmmwrapper.sequence_t_seq_label_get, _ghmmwrapper.sequence_t_seq_label_set)
    __swig_setmethods__["seq_id"] = _ghmmwrapper.sequence_t_seq_id_set
    __swig_getmethods__["seq_id"] = _ghmmwrapper.sequence_t_seq_id_get
    if _newclass:seq_id = property(_ghmmwrapper.sequence_t_seq_id_get, _ghmmwrapper.sequence_t_seq_id_set)
    __swig_setmethods__["seq_w"] = _ghmmwrapper.sequence_t_seq_w_set
    __swig_getmethods__["seq_w"] = _ghmmwrapper.sequence_t_seq_w_get
    if _newclass:seq_w = property(_ghmmwrapper.sequence_t_seq_w_get, _ghmmwrapper.sequence_t_seq_w_set)
    __swig_setmethods__["seq_number"] = _ghmmwrapper.sequence_t_seq_number_set
    __swig_getmethods__["seq_number"] = _ghmmwrapper.sequence_t_seq_number_get
    if _newclass:seq_number = property(_ghmmwrapper.sequence_t_seq_number_get, _ghmmwrapper.sequence_t_seq_number_set)
    __swig_setmethods__["total_w"] = _ghmmwrapper.sequence_t_total_w_set
    __swig_getmethods__["total_w"] = _ghmmwrapper.sequence_t_total_w_get
    if _newclass:total_w = property(_ghmmwrapper.sequence_t_total_w_get, _ghmmwrapper.sequence_t_total_w_set)
    __swig_setmethods__["state_labels"] = _ghmmwrapper.sequence_t_state_labels_set
    __swig_getmethods__["state_labels"] = _ghmmwrapper.sequence_t_state_labels_get
    if _newclass:state_labels = property(_ghmmwrapper.sequence_t_state_labels_get, _ghmmwrapper.sequence_t_state_labels_set)
    __swig_setmethods__["state_labels_len"] = _ghmmwrapper.sequence_t_state_labels_len_set
    __swig_getmethods__["state_labels_len"] = _ghmmwrapper.sequence_t_state_labels_len_get
    if _newclass:state_labels_len = property(_ghmmwrapper.sequence_t_state_labels_len_get, _ghmmwrapper.sequence_t_state_labels_len_set)

class sequence_tPtr(sequence_t):
    def __init__(self, this):
        _swig_setattr(self, sequence_t, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, sequence_t, 'thisown', 0)
        _swig_setattr(self, sequence_t,self.__class__,sequence_t)
_ghmmwrapper.sequence_t_swigregister(sequence_tPtr)
cvar = _ghmmwrapper.cvar


new_sequence_setPtr = _ghmmwrapper.new_sequence_setPtr

copy_sequence_setPtr = _ghmmwrapper.copy_sequence_setPtr

delete_sequence_setPtr = _ghmmwrapper.delete_sequence_setPtr

sequence_setPtr_assign = _ghmmwrapper.sequence_setPtr_assign

sequence_setPtr_value = _ghmmwrapper.sequence_setPtr_value
class sequence_d_t(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, sequence_d_t, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, sequence_d_t, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C sequence_d_t instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["seq"] = _ghmmwrapper.sequence_d_t_seq_set
    __swig_getmethods__["seq"] = _ghmmwrapper.sequence_d_t_seq_get
    if _newclass:seq = property(_ghmmwrapper.sequence_d_t_seq_get, _ghmmwrapper.sequence_d_t_seq_set)
    __swig_setmethods__["seq_len"] = _ghmmwrapper.sequence_d_t_seq_len_set
    __swig_getmethods__["seq_len"] = _ghmmwrapper.sequence_d_t_seq_len_get
    if _newclass:seq_len = property(_ghmmwrapper.sequence_d_t_seq_len_get, _ghmmwrapper.sequence_d_t_seq_len_set)
    __swig_setmethods__["seq_label"] = _ghmmwrapper.sequence_d_t_seq_label_set
    __swig_getmethods__["seq_label"] = _ghmmwrapper.sequence_d_t_seq_label_get
    if _newclass:seq_label = property(_ghmmwrapper.sequence_d_t_seq_label_get, _ghmmwrapper.sequence_d_t_seq_label_set)
    __swig_setmethods__["seq_id"] = _ghmmwrapper.sequence_d_t_seq_id_set
    __swig_getmethods__["seq_id"] = _ghmmwrapper.sequence_d_t_seq_id_get
    if _newclass:seq_id = property(_ghmmwrapper.sequence_d_t_seq_id_get, _ghmmwrapper.sequence_d_t_seq_id_set)
    __swig_setmethods__["seq_w"] = _ghmmwrapper.sequence_d_t_seq_w_set
    __swig_getmethods__["seq_w"] = _ghmmwrapper.sequence_d_t_seq_w_get
    if _newclass:seq_w = property(_ghmmwrapper.sequence_d_t_seq_w_get, _ghmmwrapper.sequence_d_t_seq_w_set)
    __swig_setmethods__["seq_number"] = _ghmmwrapper.sequence_d_t_seq_number_set
    __swig_getmethods__["seq_number"] = _ghmmwrapper.sequence_d_t_seq_number_get
    if _newclass:seq_number = property(_ghmmwrapper.sequence_d_t_seq_number_get, _ghmmwrapper.sequence_d_t_seq_number_set)
    __swig_setmethods__["total_w"] = _ghmmwrapper.sequence_d_t_total_w_set
    __swig_getmethods__["total_w"] = _ghmmwrapper.sequence_d_t_total_w_get
    if _newclass:total_w = property(_ghmmwrapper.sequence_d_t_total_w_get, _ghmmwrapper.sequence_d_t_total_w_set)

class sequence_d_tPtr(sequence_d_t):
    def __init__(self, this):
        _swig_setattr(self, sequence_d_t, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, sequence_d_t, 'thisown', 0)
        _swig_setattr(self, sequence_d_t,self.__class__,sequence_d_t)
_ghmmwrapper.sequence_d_t_swigregister(sequence_d_tPtr)


new_sequence_d_setPtr = _ghmmwrapper.new_sequence_d_setPtr

copy_sequence_d_setPtr = _ghmmwrapper.copy_sequence_d_setPtr

delete_sequence_d_setPtr = _ghmmwrapper.delete_sequence_d_setPtr

sequence_d_setPtr_assign = _ghmmwrapper.sequence_d_setPtr_assign

sequence_d_setPtr_value = _ghmmwrapper.sequence_d_setPtr_value

sequence_calloc = _ghmmwrapper.sequence_calloc

sequence_d_calloc = _ghmmwrapper.sequence_d_calloc

sequence_clean = _ghmmwrapper.sequence_clean

sequence_d_clean = _ghmmwrapper.sequence_d_clean

sequence_print = _ghmmwrapper.sequence_print

sequence_copy = _ghmmwrapper.sequence_copy

sequence_d_read = _ghmmwrapper.sequence_d_read

sequence_add = _ghmmwrapper.sequence_add

sequence_d_add = _ghmmwrapper.sequence_d_add

sequence_free = _ghmmwrapper.sequence_free

sequence_d_free = _ghmmwrapper.sequence_d_free

sequence_d_print = _ghmmwrapper.sequence_d_print

get_onesequence = _ghmmwrapper.get_onesequence

seq_read = _ghmmwrapper.seq_read

get_seq_ptr = _ghmmwrapper.get_seq_ptr

sequence_d_t_array = _ghmmwrapper.sequence_d_t_array

get_seq_d_ptr = _ghmmwrapper.get_seq_d_ptr

set_seq_d_array = _ghmmwrapper.set_seq_d_array

set_sequence_d_label = _ghmmwrapper.set_sequence_d_label

get_sequence_d_label = _ghmmwrapper.get_sequence_d_label

seq_d_read = _ghmmwrapper.seq_d_read

call_sequence_print = _ghmmwrapper.call_sequence_print

call_sequence_d_print = _ghmmwrapper.call_sequence_d_print
class background_distributions(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, background_distributions, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, background_distributions, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C background_distributions instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["n"] = _ghmmwrapper.background_distributions_n_set
    __swig_getmethods__["n"] = _ghmmwrapper.background_distributions_n_get
    if _newclass:n = property(_ghmmwrapper.background_distributions_n_get, _ghmmwrapper.background_distributions_n_set)
    __swig_setmethods__["m"] = _ghmmwrapper.background_distributions_m_set
    __swig_getmethods__["m"] = _ghmmwrapper.background_distributions_m_get
    if _newclass:m = property(_ghmmwrapper.background_distributions_m_get, _ghmmwrapper.background_distributions_m_set)
    __swig_setmethods__["order"] = _ghmmwrapper.background_distributions_order_set
    __swig_getmethods__["order"] = _ghmmwrapper.background_distributions_order_get
    if _newclass:order = property(_ghmmwrapper.background_distributions_order_get, _ghmmwrapper.background_distributions_order_set)
    __swig_setmethods__["b"] = _ghmmwrapper.background_distributions_b_set
    __swig_getmethods__["b"] = _ghmmwrapper.background_distributions_b_get
    if _newclass:b = property(_ghmmwrapper.background_distributions_b_get, _ghmmwrapper.background_distributions_b_set)

class background_distributionsPtr(background_distributions):
    def __init__(self, this):
        _swig_setattr(self, background_distributions, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, background_distributions, 'thisown', 0)
        _swig_setattr(self, background_distributions,self.__class__,background_distributions)
_ghmmwrapper.background_distributions_swigregister(background_distributionsPtr)

class state(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, state, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, state, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C state instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["pi"] = _ghmmwrapper.state_pi_set
    __swig_getmethods__["pi"] = _ghmmwrapper.state_pi_get
    if _newclass:pi = property(_ghmmwrapper.state_pi_get, _ghmmwrapper.state_pi_set)
    __swig_setmethods__["b"] = _ghmmwrapper.state_b_set
    __swig_getmethods__["b"] = _ghmmwrapper.state_b_get
    if _newclass:b = property(_ghmmwrapper.state_b_get, _ghmmwrapper.state_b_set)
    __swig_setmethods__["order"] = _ghmmwrapper.state_order_set
    __swig_getmethods__["order"] = _ghmmwrapper.state_order_get
    if _newclass:order = property(_ghmmwrapper.state_order_get, _ghmmwrapper.state_order_set)
    __swig_setmethods__["out_id"] = _ghmmwrapper.state_out_id_set
    __swig_getmethods__["out_id"] = _ghmmwrapper.state_out_id_get
    if _newclass:out_id = property(_ghmmwrapper.state_out_id_get, _ghmmwrapper.state_out_id_set)
    __swig_setmethods__["in_id"] = _ghmmwrapper.state_in_id_set
    __swig_getmethods__["in_id"] = _ghmmwrapper.state_in_id_get
    if _newclass:in_id = property(_ghmmwrapper.state_in_id_get, _ghmmwrapper.state_in_id_set)
    __swig_setmethods__["out_a"] = _ghmmwrapper.state_out_a_set
    __swig_getmethods__["out_a"] = _ghmmwrapper.state_out_a_get
    if _newclass:out_a = property(_ghmmwrapper.state_out_a_get, _ghmmwrapper.state_out_a_set)
    __swig_setmethods__["in_a"] = _ghmmwrapper.state_in_a_set
    __swig_getmethods__["in_a"] = _ghmmwrapper.state_in_a_get
    if _newclass:in_a = property(_ghmmwrapper.state_in_a_get, _ghmmwrapper.state_in_a_set)
    __swig_setmethods__["out_states"] = _ghmmwrapper.state_out_states_set
    __swig_getmethods__["out_states"] = _ghmmwrapper.state_out_states_get
    if _newclass:out_states = property(_ghmmwrapper.state_out_states_get, _ghmmwrapper.state_out_states_set)
    __swig_setmethods__["in_states"] = _ghmmwrapper.state_in_states_set
    __swig_getmethods__["in_states"] = _ghmmwrapper.state_in_states_get
    if _newclass:in_states = property(_ghmmwrapper.state_in_states_get, _ghmmwrapper.state_in_states_set)
    __swig_setmethods__["fix"] = _ghmmwrapper.state_fix_set
    __swig_getmethods__["fix"] = _ghmmwrapper.state_fix_get
    if _newclass:fix = property(_ghmmwrapper.state_fix_get, _ghmmwrapper.state_fix_set)
    __swig_setmethods__["label"] = _ghmmwrapper.state_label_set
    __swig_getmethods__["label"] = _ghmmwrapper.state_label_get
    if _newclass:label = property(_ghmmwrapper.state_label_get, _ghmmwrapper.state_label_set)

class statePtr(state):
    def __init__(self, this):
        _swig_setattr(self, state, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, state, 'thisown', 0)
        _swig_setattr(self, state,self.__class__,state)
_ghmmwrapper.state_swigregister(statePtr)

class model(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, model, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, model, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C model instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["N"] = _ghmmwrapper.model_N_set
    __swig_getmethods__["N"] = _ghmmwrapper.model_N_get
    if _newclass:N = property(_ghmmwrapper.model_N_get, _ghmmwrapper.model_N_set)
    __swig_setmethods__["M"] = _ghmmwrapper.model_M_set
    __swig_getmethods__["M"] = _ghmmwrapper.model_M_get
    if _newclass:M = property(_ghmmwrapper.model_M_get, _ghmmwrapper.model_M_set)
    __swig_setmethods__["s"] = _ghmmwrapper.model_s_set
    __swig_getmethods__["s"] = _ghmmwrapper.model_s_get
    if _newclass:s = property(_ghmmwrapper.model_s_get, _ghmmwrapper.model_s_set)
    __swig_setmethods__["prior"] = _ghmmwrapper.model_prior_set
    __swig_getmethods__["prior"] = _ghmmwrapper.model_prior_get
    if _newclass:prior = property(_ghmmwrapper.model_prior_get, _ghmmwrapper.model_prior_set)
    __swig_setmethods__["name"] = _ghmmwrapper.model_name_set
    __swig_getmethods__["name"] = _ghmmwrapper.model_name_get
    if _newclass:name = property(_ghmmwrapper.model_name_get, _ghmmwrapper.model_name_set)
    __swig_setmethods__["model_type"] = _ghmmwrapper.model_model_type_set
    __swig_getmethods__["model_type"] = _ghmmwrapper.model_model_type_get
    if _newclass:model_type = property(_ghmmwrapper.model_model_type_get, _ghmmwrapper.model_model_type_set)
    __swig_setmethods__["silent"] = _ghmmwrapper.model_silent_set
    __swig_getmethods__["silent"] = _ghmmwrapper.model_silent_get
    if _newclass:silent = property(_ghmmwrapper.model_silent_get, _ghmmwrapper.model_silent_set)
    __swig_setmethods__["maxorder"] = _ghmmwrapper.model_maxorder_set
    __swig_getmethods__["maxorder"] = _ghmmwrapper.model_maxorder_get
    if _newclass:maxorder = property(_ghmmwrapper.model_maxorder_get, _ghmmwrapper.model_maxorder_set)
    __swig_setmethods__["emission_history"] = _ghmmwrapper.model_emission_history_set
    __swig_getmethods__["emission_history"] = _ghmmwrapper.model_emission_history_get
    if _newclass:emission_history = property(_ghmmwrapper.model_emission_history_get, _ghmmwrapper.model_emission_history_set)
    __swig_setmethods__["tied_to"] = _ghmmwrapper.model_tied_to_set
    __swig_getmethods__["tied_to"] = _ghmmwrapper.model_tied_to_get
    if _newclass:tied_to = property(_ghmmwrapper.model_tied_to_get, _ghmmwrapper.model_tied_to_set)
    __swig_setmethods__["background_id"] = _ghmmwrapper.model_background_id_set
    __swig_getmethods__["background_id"] = _ghmmwrapper.model_background_id_get
    if _newclass:background_id = property(_ghmmwrapper.model_background_id_get, _ghmmwrapper.model_background_id_set)
    __swig_setmethods__["bp"] = _ghmmwrapper.model_bp_set
    __swig_getmethods__["bp"] = _ghmmwrapper.model_bp_get
    if _newclass:bp = property(_ghmmwrapper.model_bp_get, _ghmmwrapper.model_bp_set)
    __swig_setmethods__["topo_order"] = _ghmmwrapper.model_topo_order_set
    __swig_getmethods__["topo_order"] = _ghmmwrapper.model_topo_order_get
    if _newclass:topo_order = property(_ghmmwrapper.model_topo_order_get, _ghmmwrapper.model_topo_order_set)
    __swig_setmethods__["topo_order_length"] = _ghmmwrapper.model_topo_order_length_set
    __swig_getmethods__["topo_order_length"] = _ghmmwrapper.model_topo_order_length_get
    if _newclass:topo_order_length = property(_ghmmwrapper.model_topo_order_length_get, _ghmmwrapper.model_topo_order_length_set)

class modelPtr(model):
    def __init__(self, this):
        _swig_setattr(self, model, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, model, 'thisown', 0)
        _swig_setattr(self, model,self.__class__,model)
_ghmmwrapper.model_swigregister(modelPtr)


model_free = _ghmmwrapper.model_free

model_read = _ghmmwrapper.model_read

model_print = _ghmmwrapper.model_print

model_from_sequence_ascii = _ghmmwrapper.model_from_sequence_ascii

model_from_sequence = _ghmmwrapper.model_from_sequence

model_copy = _ghmmwrapper.model_copy

model_check = _ghmmwrapper.model_check

model_check_compatibility = _ghmmwrapper.model_check_compatibility

model_generate_from_sequence = _ghmmwrapper.model_generate_from_sequence

model_generate_sequences = _ghmmwrapper.model_generate_sequences

model_likelihood = _ghmmwrapper.model_likelihood

model_prob_distance = _ghmmwrapper.model_prob_distance

reestimate_baum_welch = _ghmmwrapper.reestimate_baum_welch

reestimate_baum_welch_nstep = _ghmmwrapper.reestimate_baum_welch_nstep

reestimate_update_tie_groups = _ghmmwrapper.reestimate_update_tie_groups

reestimate_baum_welch_label = _ghmmwrapper.reestimate_baum_welch_label

reestimate_baum_welch_nstep_label = _ghmmwrapper.reestimate_baum_welch_nstep_label

gradient_descent = _ghmmwrapper.gradient_descent

kbest = _ghmmwrapper.kbest

viterbi = _ghmmwrapper.viterbi

viterbi_logp = _ghmmwrapper.viterbi_logp

discriminative = _ghmmwrapper.discriminative

discrime_compute_performance = _ghmmwrapper.discrime_compute_performance

discrime_modelarray_alloc = _ghmmwrapper.discrime_modelarray_alloc

discrime_modelarray_dealloc = _ghmmwrapper.discrime_modelarray_dealloc

discrime_modelarray_setptr = _ghmmwrapper.discrime_modelarray_setptr

discrime_modelarray_getptr = _ghmmwrapper.discrime_modelarray_getptr

discrime_seqarray_alloc = _ghmmwrapper.discrime_seqarray_alloc

discrime_seqarray_dealloc = _ghmmwrapper.discrime_seqarray_dealloc

discrime_seqarray_setptr = _ghmmwrapper.discrime_seqarray_setptr

discrime_seqarray_getptr = _ghmmwrapper.discrime_seqarray_getptr

foba_forward = _ghmmwrapper.foba_forward

foba_backward = _ghmmwrapper.foba_backward

foba_logp = _ghmmwrapper.foba_logp

model_set_transition = _ghmmwrapper.model_set_transition

foba_forward_lean = _ghmmwrapper.foba_forward_lean

foba_label_forward = _ghmmwrapper.foba_label_forward

foba_label_logp = _ghmmwrapper.foba_label_logp

foba_label_backward = _ghmmwrapper.foba_label_backward

model_alloc_background_distributions = _ghmmwrapper.model_alloc_background_distributions

model_copy_background_distributions = _ghmmwrapper.model_copy_background_distributions

model_free_background_distributions = _ghmmwrapper.model_free_background_distributions

get_emission_index = _ghmmwrapper.get_emission_index

update_emission_history = _ghmmwrapper.update_emission_history

update_emission_history_front = _ghmmwrapper.update_emission_history_front

model_normalize = _ghmmwrapper.model_normalize

model_add_noise = _ghmmwrapper.model_add_noise

model_apply_background = _ghmmwrapper.model_apply_background

new_model = _ghmmwrapper.new_model

arraystate = _ghmmwrapper.arraystate

get_stateptr = _ghmmwrapper.get_stateptr

call_model_print = _ghmmwrapper.call_model_print

call_model_free = _ghmmwrapper.call_model_free

get_model_ptr = _ghmmwrapper.get_model_ptr

cast_model_ptr = _ghmmwrapper.cast_model_ptr

model_label_generate_sequences = _ghmmwrapper.model_label_generate_sequences
class sdstate(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, sdstate, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, sdstate, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C sdstate instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["pi"] = _ghmmwrapper.sdstate_pi_set
    __swig_getmethods__["pi"] = _ghmmwrapper.sdstate_pi_get
    if _newclass:pi = property(_ghmmwrapper.sdstate_pi_get, _ghmmwrapper.sdstate_pi_set)
    __swig_setmethods__["b"] = _ghmmwrapper.sdstate_b_set
    __swig_getmethods__["b"] = _ghmmwrapper.sdstate_b_get
    if _newclass:b = property(_ghmmwrapper.sdstate_b_get, _ghmmwrapper.sdstate_b_set)
    __swig_setmethods__["out_id"] = _ghmmwrapper.sdstate_out_id_set
    __swig_getmethods__["out_id"] = _ghmmwrapper.sdstate_out_id_get
    if _newclass:out_id = property(_ghmmwrapper.sdstate_out_id_get, _ghmmwrapper.sdstate_out_id_set)
    __swig_setmethods__["in_id"] = _ghmmwrapper.sdstate_in_id_set
    __swig_getmethods__["in_id"] = _ghmmwrapper.sdstate_in_id_get
    if _newclass:in_id = property(_ghmmwrapper.sdstate_in_id_get, _ghmmwrapper.sdstate_in_id_set)
    __swig_setmethods__["out_a"] = _ghmmwrapper.sdstate_out_a_set
    __swig_getmethods__["out_a"] = _ghmmwrapper.sdstate_out_a_get
    if _newclass:out_a = property(_ghmmwrapper.sdstate_out_a_get, _ghmmwrapper.sdstate_out_a_set)
    __swig_setmethods__["in_a"] = _ghmmwrapper.sdstate_in_a_set
    __swig_getmethods__["in_a"] = _ghmmwrapper.sdstate_in_a_get
    if _newclass:in_a = property(_ghmmwrapper.sdstate_in_a_get, _ghmmwrapper.sdstate_in_a_set)
    __swig_setmethods__["out_states"] = _ghmmwrapper.sdstate_out_states_set
    __swig_getmethods__["out_states"] = _ghmmwrapper.sdstate_out_states_get
    if _newclass:out_states = property(_ghmmwrapper.sdstate_out_states_get, _ghmmwrapper.sdstate_out_states_set)
    __swig_setmethods__["in_states"] = _ghmmwrapper.sdstate_in_states_set
    __swig_getmethods__["in_states"] = _ghmmwrapper.sdstate_in_states_get
    if _newclass:in_states = property(_ghmmwrapper.sdstate_in_states_get, _ghmmwrapper.sdstate_in_states_set)
    __swig_setmethods__["fix"] = _ghmmwrapper.sdstate_fix_set
    __swig_getmethods__["fix"] = _ghmmwrapper.sdstate_fix_get
    if _newclass:fix = property(_ghmmwrapper.sdstate_fix_get, _ghmmwrapper.sdstate_fix_set)

class sdstatePtr(sdstate):
    def __init__(self, this):
        _swig_setattr(self, sdstate, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, sdstate, 'thisown', 0)
        _swig_setattr(self, sdstate,self.__class__,sdstate)
_ghmmwrapper.sdstate_swigregister(sdstatePtr)

class sdmodel(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, sdmodel, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, sdmodel, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C sdmodel instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["N"] = _ghmmwrapper.sdmodel_N_set
    __swig_getmethods__["N"] = _ghmmwrapper.sdmodel_N_get
    if _newclass:N = property(_ghmmwrapper.sdmodel_N_get, _ghmmwrapper.sdmodel_N_set)
    __swig_setmethods__["M"] = _ghmmwrapper.sdmodel_M_set
    __swig_getmethods__["M"] = _ghmmwrapper.sdmodel_M_get
    if _newclass:M = property(_ghmmwrapper.sdmodel_M_get, _ghmmwrapper.sdmodel_M_set)
    __swig_setmethods__["cos"] = _ghmmwrapper.sdmodel_cos_set
    __swig_getmethods__["cos"] = _ghmmwrapper.sdmodel_cos_get
    if _newclass:cos = property(_ghmmwrapper.sdmodel_cos_get, _ghmmwrapper.sdmodel_cos_set)
    __swig_setmethods__["s"] = _ghmmwrapper.sdmodel_s_set
    __swig_getmethods__["s"] = _ghmmwrapper.sdmodel_s_get
    if _newclass:s = property(_ghmmwrapper.sdmodel_s_get, _ghmmwrapper.sdmodel_s_set)
    __swig_setmethods__["prior"] = _ghmmwrapper.sdmodel_prior_set
    __swig_getmethods__["prior"] = _ghmmwrapper.sdmodel_prior_get
    if _newclass:prior = property(_ghmmwrapper.sdmodel_prior_get, _ghmmwrapper.sdmodel_prior_set)
    __swig_setmethods__["get_class"] = _ghmmwrapper.sdmodel_get_class_set
    __swig_getmethods__["get_class"] = _ghmmwrapper.sdmodel_get_class_get
    if _newclass:get_class = property(_ghmmwrapper.sdmodel_get_class_get, _ghmmwrapper.sdmodel_get_class_set)
    __swig_setmethods__["model_type"] = _ghmmwrapper.sdmodel_model_type_set
    __swig_getmethods__["model_type"] = _ghmmwrapper.sdmodel_model_type_get
    if _newclass:model_type = property(_ghmmwrapper.sdmodel_model_type_get, _ghmmwrapper.sdmodel_model_type_set)
    __swig_setmethods__["silent"] = _ghmmwrapper.sdmodel_silent_set
    __swig_getmethods__["silent"] = _ghmmwrapper.sdmodel_silent_get
    if _newclass:silent = property(_ghmmwrapper.sdmodel_silent_get, _ghmmwrapper.sdmodel_silent_set)

class sdmodelPtr(sdmodel):
    def __init__(self, this):
        _swig_setattr(self, sdmodel, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, sdmodel, 'thisown', 0)
        _swig_setattr(self, sdmodel,self.__class__,sdmodel)
_ghmmwrapper.sdmodel_swigregister(sdmodelPtr)


sdmodel_free = _ghmmwrapper.sdmodel_free

sdmodel_generate_sequences = _ghmmwrapper.sdmodel_generate_sequences

sdmodel_likelihood = _ghmmwrapper.sdmodel_likelihood

smodel_individual_likelihoods = _ghmmwrapper.smodel_individual_likelihoods

sdviterbi = _ghmmwrapper.sdviterbi

sdfoba_forward = _ghmmwrapper.sdfoba_forward

cp_class_change = _ghmmwrapper.cp_class_change

setSwitchingFunction = _ghmmwrapper.setSwitchingFunction

call_sdmodel_free = _ghmmwrapper.call_sdmodel_free

arraysdstate = _ghmmwrapper.arraysdstate

get_sdstateptr = _ghmmwrapper.get_sdstateptr
normal = _ghmmwrapper.normal
normal_pos = _ghmmwrapper.normal_pos
normal_approx = _ghmmwrapper.normal_approx
density_number = _ghmmwrapper.density_number
class sstate(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, sstate, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, sstate, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C sstate instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["pi"] = _ghmmwrapper.sstate_pi_set
    __swig_getmethods__["pi"] = _ghmmwrapper.sstate_pi_get
    if _newclass:pi = property(_ghmmwrapper.sstate_pi_get, _ghmmwrapper.sstate_pi_set)
    __swig_setmethods__["out_id"] = _ghmmwrapper.sstate_out_id_set
    __swig_getmethods__["out_id"] = _ghmmwrapper.sstate_out_id_get
    if _newclass:out_id = property(_ghmmwrapper.sstate_out_id_get, _ghmmwrapper.sstate_out_id_set)
    __swig_setmethods__["in_id"] = _ghmmwrapper.sstate_in_id_set
    __swig_getmethods__["in_id"] = _ghmmwrapper.sstate_in_id_get
    if _newclass:in_id = property(_ghmmwrapper.sstate_in_id_get, _ghmmwrapper.sstate_in_id_set)
    __swig_setmethods__["out_a"] = _ghmmwrapper.sstate_out_a_set
    __swig_getmethods__["out_a"] = _ghmmwrapper.sstate_out_a_get
    if _newclass:out_a = property(_ghmmwrapper.sstate_out_a_get, _ghmmwrapper.sstate_out_a_set)
    __swig_setmethods__["in_a"] = _ghmmwrapper.sstate_in_a_set
    __swig_getmethods__["in_a"] = _ghmmwrapper.sstate_in_a_get
    if _newclass:in_a = property(_ghmmwrapper.sstate_in_a_get, _ghmmwrapper.sstate_in_a_set)
    __swig_setmethods__["out_states"] = _ghmmwrapper.sstate_out_states_set
    __swig_getmethods__["out_states"] = _ghmmwrapper.sstate_out_states_get
    if _newclass:out_states = property(_ghmmwrapper.sstate_out_states_get, _ghmmwrapper.sstate_out_states_set)
    __swig_setmethods__["in_states"] = _ghmmwrapper.sstate_in_states_set
    __swig_getmethods__["in_states"] = _ghmmwrapper.sstate_in_states_get
    if _newclass:in_states = property(_ghmmwrapper.sstate_in_states_get, _ghmmwrapper.sstate_in_states_set)
    __swig_setmethods__["c"] = _ghmmwrapper.sstate_c_set
    __swig_getmethods__["c"] = _ghmmwrapper.sstate_c_get
    if _newclass:c = property(_ghmmwrapper.sstate_c_get, _ghmmwrapper.sstate_c_set)
    __swig_setmethods__["mue"] = _ghmmwrapper.sstate_mue_set
    __swig_getmethods__["mue"] = _ghmmwrapper.sstate_mue_get
    if _newclass:mue = property(_ghmmwrapper.sstate_mue_get, _ghmmwrapper.sstate_mue_set)
    __swig_setmethods__["u"] = _ghmmwrapper.sstate_u_set
    __swig_getmethods__["u"] = _ghmmwrapper.sstate_u_get
    if _newclass:u = property(_ghmmwrapper.sstate_u_get, _ghmmwrapper.sstate_u_set)
    __swig_setmethods__["fix"] = _ghmmwrapper.sstate_fix_set
    __swig_getmethods__["fix"] = _ghmmwrapper.sstate_fix_get
    if _newclass:fix = property(_ghmmwrapper.sstate_fix_get, _ghmmwrapper.sstate_fix_set)
    __swig_setmethods__["mixture_fix"] = _ghmmwrapper.sstate_mixture_fix_set
    __swig_getmethods__["mixture_fix"] = _ghmmwrapper.sstate_mixture_fix_get
    if _newclass:mixture_fix = property(_ghmmwrapper.sstate_mixture_fix_get, _ghmmwrapper.sstate_mixture_fix_set)

class sstatePtr(sstate):
    def __init__(self, this):
        _swig_setattr(self, sstate, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, sstate, 'thisown', 0)
        _swig_setattr(self, sstate,self.__class__,sstate)
_ghmmwrapper.sstate_swigregister(sstatePtr)

class smodel(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, smodel, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, smodel, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C smodel instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["N"] = _ghmmwrapper.smodel_N_set
    __swig_getmethods__["N"] = _ghmmwrapper.smodel_N_get
    if _newclass:N = property(_ghmmwrapper.smodel_N_get, _ghmmwrapper.smodel_N_set)
    __swig_setmethods__["M"] = _ghmmwrapper.smodel_M_set
    __swig_getmethods__["M"] = _ghmmwrapper.smodel_M_get
    if _newclass:M = property(_ghmmwrapper.smodel_M_get, _ghmmwrapper.smodel_M_set)
    __swig_setmethods__["cos"] = _ghmmwrapper.smodel_cos_set
    __swig_getmethods__["cos"] = _ghmmwrapper.smodel_cos_get
    if _newclass:cos = property(_ghmmwrapper.smodel_cos_get, _ghmmwrapper.smodel_cos_set)
    __swig_setmethods__["density"] = _ghmmwrapper.smodel_density_set
    __swig_getmethods__["density"] = _ghmmwrapper.smodel_density_get
    if _newclass:density = property(_ghmmwrapper.smodel_density_get, _ghmmwrapper.smodel_density_set)
    __swig_setmethods__["prior"] = _ghmmwrapper.smodel_prior_set
    __swig_getmethods__["prior"] = _ghmmwrapper.smodel_prior_get
    if _newclass:prior = property(_ghmmwrapper.smodel_prior_get, _ghmmwrapper.smodel_prior_set)
    __swig_setmethods__["s"] = _ghmmwrapper.smodel_s_set
    __swig_getmethods__["s"] = _ghmmwrapper.smodel_s_get
    if _newclass:s = property(_ghmmwrapper.smodel_s_get, _ghmmwrapper.smodel_s_set)

class smodelPtr(smodel):
    def __init__(self, this):
        _swig_setattr(self, smodel, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, smodel, 'thisown', 0)
        _swig_setattr(self, smodel,self.__class__,smodel)
_ghmmwrapper.smodel_swigregister(smodelPtr)


smodel_free = _ghmmwrapper.smodel_free

smodel_read = _ghmmwrapper.smodel_read

smodel_copy = _ghmmwrapper.smodel_copy

smodel_generate_sequences = _ghmmwrapper.smodel_generate_sequences

smodel_prob_distance = _ghmmwrapper.smodel_prob_distance

sfoba_forward = _ghmmwrapper.sfoba_forward

sfoba_backward = _ghmmwrapper.sfoba_backward

sfoba_logp = _ghmmwrapper.sfoba_logp

smodel_alloc_fill = _ghmmwrapper.smodel_alloc_fill

smodel_set_pivector = _ghmmwrapper.smodel_set_pivector

smodel_set_fixvector = _ghmmwrapper.smodel_set_fixvector

smodel_set_transition = _ghmmwrapper.smodel_set_transition

smodel_get_transition = _ghmmwrapper.smodel_get_transition

smodel_set_mean = _ghmmwrapper.smodel_set_mean

smodel_set_variance = _ghmmwrapper.smodel_set_variance

smodel_likelihood = _ghmmwrapper.smodel_likelihood

smodel_sorted_individual_likelihoods = _ghmmwrapper.smodel_sorted_individual_likelihoods

sviterbi = _ghmmwrapper.sviterbi

set_trunc_density = _ghmmwrapper.set_trunc_density

new_smodel = _ghmmwrapper.new_smodel

arraysstate = _ghmmwrapper.arraysstate

get_sstate_ptr = _ghmmwrapper.get_sstate_ptr

call_smodel_free = _ghmmwrapper.call_smodel_free

free_smodel_array = _ghmmwrapper.free_smodel_array

smodel_print_stdout = _ghmmwrapper.smodel_print_stdout

get_sstate = _ghmmwrapper.get_sstate

smodel_array = _ghmmwrapper.smodel_array

get_smodel_ptr = _ghmmwrapper.get_smodel_ptr

set_smodel_ptr = _ghmmwrapper.set_smodel_ptr

cast_smodel_ptr = _ghmmwrapper.cast_smodel_ptr

call_smodel_print = _ghmmwrapper.call_smodel_print
class scluster_t(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, scluster_t, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, scluster_t, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C scluster_t instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["smo"] = _ghmmwrapper.scluster_t_smo_set
    __swig_getmethods__["smo"] = _ghmmwrapper.scluster_t_smo_get
    if _newclass:smo = property(_ghmmwrapper.scluster_t_smo_get, _ghmmwrapper.scluster_t_smo_set)
    __swig_setmethods__["smo_seq"] = _ghmmwrapper.scluster_t_smo_seq_set
    __swig_getmethods__["smo_seq"] = _ghmmwrapper.scluster_t_smo_seq_get
    if _newclass:smo_seq = property(_ghmmwrapper.scluster_t_smo_seq_get, _ghmmwrapper.scluster_t_smo_seq_set)
    __swig_setmethods__["smo_number"] = _ghmmwrapper.scluster_t_smo_number_set
    __swig_getmethods__["smo_number"] = _ghmmwrapper.scluster_t_smo_number_get
    if _newclass:smo_number = property(_ghmmwrapper.scluster_t_smo_number_get, _ghmmwrapper.scluster_t_smo_number_set)
    __swig_setmethods__["seq_counter"] = _ghmmwrapper.scluster_t_seq_counter_set
    __swig_getmethods__["seq_counter"] = _ghmmwrapper.scluster_t_seq_counter_get
    if _newclass:seq_counter = property(_ghmmwrapper.scluster_t_seq_counter_get, _ghmmwrapper.scluster_t_seq_counter_set)
    __swig_setmethods__["smo_Z_MD"] = _ghmmwrapper.scluster_t_smo_Z_MD_set
    __swig_getmethods__["smo_Z_MD"] = _ghmmwrapper.scluster_t_smo_Z_MD_get
    if _newclass:smo_Z_MD = property(_ghmmwrapper.scluster_t_smo_Z_MD_get, _ghmmwrapper.scluster_t_smo_Z_MD_set)
    __swig_setmethods__["smo_Z_MAW"] = _ghmmwrapper.scluster_t_smo_Z_MAW_set
    __swig_getmethods__["smo_Z_MAW"] = _ghmmwrapper.scluster_t_smo_Z_MAW_get
    if _newclass:smo_Z_MAW = property(_ghmmwrapper.scluster_t_smo_Z_MAW_get, _ghmmwrapper.scluster_t_smo_Z_MAW_set)

class scluster_tPtr(scluster_t):
    def __init__(self, this):
        _swig_setattr(self, scluster_t, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, scluster_t, 'thisown', 0)
        _swig_setattr(self, scluster_t,self.__class__,scluster_t)
_ghmmwrapper.scluster_t_swigregister(scluster_tPtr)


scluster_t_free = _ghmmwrapper.scluster_t_free

scluster_random_labels = _ghmmwrapper.scluster_random_labels

scluster_prob = _ghmmwrapper.scluster_prob

scluster_best_model = _ghmmwrapper.scluster_best_model

scluster_update = _ghmmwrapper.scluster_update

scluster_print_likelihood = _ghmmwrapper.scluster_print_likelihood

scluster_out = _ghmmwrapper.scluster_out

scluster_log_aposteriori = _ghmmwrapper.scluster_log_aposteriori

scluster_avoid_empty_smodel = _ghmmwrapper.scluster_avoid_empty_smodel

scluster_hmm = _ghmmwrapper.scluster_hmm

scluster_printl_stdout = _ghmmwrapper.scluster_printl_stdout
class smosqd_t(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, smosqd_t, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, smosqd_t, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C smosqd_t instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["smo"] = _ghmmwrapper.smosqd_t_smo_set
    __swig_getmethods__["smo"] = _ghmmwrapper.smosqd_t_smo_get
    if _newclass:smo = property(_ghmmwrapper.smosqd_t_smo_get, _ghmmwrapper.smosqd_t_smo_set)
    __swig_setmethods__["sqd"] = _ghmmwrapper.smosqd_t_sqd_set
    __swig_getmethods__["sqd"] = _ghmmwrapper.smosqd_t_sqd_get
    if _newclass:sqd = property(_ghmmwrapper.smosqd_t_sqd_get, _ghmmwrapper.smosqd_t_sqd_set)
    __swig_setmethods__["logp"] = _ghmmwrapper.smosqd_t_logp_set
    __swig_getmethods__["logp"] = _ghmmwrapper.smosqd_t_logp_get
    if _newclass:logp = property(_ghmmwrapper.smosqd_t_logp_get, _ghmmwrapper.smosqd_t_logp_set)
    __swig_setmethods__["eps"] = _ghmmwrapper.smosqd_t_eps_set
    __swig_getmethods__["eps"] = _ghmmwrapper.smosqd_t_eps_get
    if _newclass:eps = property(_ghmmwrapper.smosqd_t_eps_get, _ghmmwrapper.smosqd_t_eps_set)
    __swig_setmethods__["max_iter"] = _ghmmwrapper.smosqd_t_max_iter_set
    __swig_getmethods__["max_iter"] = _ghmmwrapper.smosqd_t_max_iter_get
    if _newclass:max_iter = property(_ghmmwrapper.smosqd_t_max_iter_get, _ghmmwrapper.smosqd_t_max_iter_set)

class smosqd_tPtr(smosqd_t):
    def __init__(self, this):
        _swig_setattr(self, smosqd_t, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, smosqd_t, 'thisown', 0)
        _swig_setattr(self, smosqd_t,self.__class__,smosqd_t)
_ghmmwrapper.smosqd_t_swigregister(smosqd_tPtr)


sreestimate_baum_welch = _ghmmwrapper.sreestimate_baum_welch

smosqd_t_array = _ghmmwrapper.smosqd_t_array

set_smosq_t_smo = _ghmmwrapper.set_smosq_t_smo

get_smosqd_t_ptr = _ghmmwrapper.get_smosqd_t_ptr

free_smosqd_t = _ghmmwrapper.free_smosqd_t

int_array = _ghmmwrapper.int_array

set_arrayint = _ghmmwrapper.set_arrayint

get_arrayint = _ghmmwrapper.get_arrayint

free_arrayi = _ghmmwrapper.free_arrayi

double_array = _ghmmwrapper.double_array

set_arrayd = _ghmmwrapper.set_arrayd

get_arrayd = _ghmmwrapper.get_arrayd

free_arrayd = _ghmmwrapper.free_arrayd

long_array = _ghmmwrapper.long_array

set_arrayl = _ghmmwrapper.set_arrayl

get_arrayl = _ghmmwrapper.get_arrayl

free_arrayl = _ghmmwrapper.free_arrayl

char_array = _ghmmwrapper.char_array

set_arraychar = _ghmmwrapper.set_arraychar

get_arraychar = _ghmmwrapper.get_arraychar

double_2d_array = _ghmmwrapper.double_2d_array

double_2d_array_nocols = _ghmmwrapper.double_2d_array_nocols

set_2d_arrayd_col = _ghmmwrapper.set_2d_arrayd_col

set_2d_arrayd = _ghmmwrapper.set_2d_arrayd

get_2d_arrayd = _ghmmwrapper.get_2d_arrayd

get_col_pointer_d = _ghmmwrapper.get_col_pointer_d

double_2d_print = _ghmmwrapper.double_2d_print

cast_ptr_d = _ghmmwrapper.cast_ptr_d

free_2darrayd = _ghmmwrapper.free_2darrayd

int_2d_array_nocols = _ghmmwrapper.int_2d_array_nocols

set_2d_arrayint_col = _ghmmwrapper.set_2d_arrayint_col

get_col_pointer_int = _ghmmwrapper.get_col_pointer_int

set_2d_arrayint = _ghmmwrapper.set_2d_arrayint

get_2d_arrayint = _ghmmwrapper.get_2d_arrayint

cast_ptr_int = _ghmmwrapper.cast_ptr_int

free_2darrayint = _ghmmwrapper.free_2darrayint

freearray = _ghmmwrapper.freearray

