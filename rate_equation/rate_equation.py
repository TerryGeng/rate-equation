import itertools
import numpy as np

from rate_equation.transition_profile import Transition, TransitionProfile
from rate_equation.radiation_field import RadiationField


class RateEquation:
    # d/dt Gn = In - Out
    # In = \sum_j(!=n) Gj \sum_k Rjk \beta_jk \beta_nk
    #    (j: all ground states, k: all excited states)
    # Out = Gn \sum_k Rnk \beta_nk (1 - \beta_nk)
    #    (k: all excited states)

    def __init__(self, trans_profile, radiation_field_profile, detunings):
        self.trans_profile = trans_profile
        self.radiation = radiation_field_profile
        self.detunings = detunings
        self.pump_terms = {es: self.build_pump_terms(es) for es in self.trans_profile.excited_states}

    def build_matrix(self):
        ground_states = self.trans_profile.ground_states
        num_of_gs = len(ground_states)

        mat = np.zeros((num_of_gs, num_of_gs))

        inds = np.arange(num_of_gs)

        for i, j in itertools.product(inds, inds):
            mat[i][j] = self.calculate_matrix_element(ground_states[i], ground_states[j])

        return mat

    def build_pump_terms(self, excited_state):
        #    pump_term_j_k = Gj \sum_k Rjk \beta_jk

        e_to_all_gs = self.trans_profile.get_exc_to_gnd(excited_state)

        terms = {}

        for g_to_e in e_to_all_gs:
            base_freq = self.trans_profile.frequencies[g_to_e.group]
            tot_int = self.radiation.get_effective_intensity(g_to_e,
                                                             base_freq,
                                                             self.detunings,
                                                             self.trans_profile.gamma)
            terms[g_to_e.ground_state] = tot_int * g_to_e.strength

        return terms


    def calculate_matrix_element(self, gs1, gs2):
        # d/dt Gn = In - Out
        # In (j!=n) = \sum_j(!=n) Gj \sum_k Rjk \beta_jk \beta_nk
        #    (j: all ground states, k: all excited states)
        # Out (j==n) = Gn \sum_k Rnk \beta_nk (1 - \beta_nk)
        #    (k: all excited states)

        gs1_to_es = self.trans_profile.get_gnd_to_exc(gs1)

        if gs1 != gs2:
            # `In` term
            _in = 0

            for g_to_e in gs1_to_es:  # sum over k
                pump = self.pump_terms[g_to_e.excited_state].get(gs2, 0)
                _in += pump * g_to_e.strength

            return _in
        else:
            # `Out term`
            _out = 0

            for g_to_e in gs1_to_es:  # sum over k
                pump = self.pump_terms[g_to_e.excited_state].get(gs2, 0)
                _out += pump * (1 - g_to_e.strength)

            return -1 * _out
