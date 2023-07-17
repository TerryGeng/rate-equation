import numpy as np

from rate_equation.transition_profile import State


class Detuning:
    def get_detuning(self, field_freq, ground_state, excited_state):  # get detuning normalized by gamma
        raise NotImplementedError


class ZeemanDetuning(Detuning):
    def __init__(self, g_factors, b_field):  # B in Tesla, gamma in Hz (w/o 2\pi)
        self.g_factors = g_factors  # dict: hyperfine label -> g factor
        self.b_field = b_field

    def get_detuning(self, field_freq, ground_state, excited_state):
        from scipy.constants import physical_constants, h
        mu_B = physical_constants['Bohr magneton'][0]  # J / T

        gs_det = mu_B * self.g_factors[ground_state.hyperfine] \
                * ground_state.m * self.b_field / h

        es_det = mu_B * self.g_factors[excited_state.hyperfine] \
                * excited_state.m * self.b_field / h

        return es_det - gs_det


class DopplerDetuning(Detuning):
    def __init__(self, velocity):  # velocity in m/s
        self.velocity = velocity

    def get_detuning(self, field_freq, ground_state, excited_state):
        # opposite direction light creates positive detuning
        from scipy.constants import c
        return -1 * field_freq / c * self.velocity
