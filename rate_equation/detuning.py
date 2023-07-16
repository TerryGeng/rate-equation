import numpy as np

from rate_equation.transition_profile import State


class Detuning:
    def get_detuning(self, ground_state, excited_state):  # get detuning normalized by gamma
        raise NotImplementedError


class ZeemanDetuning(Detuning):
    def __init__(self, g_factors, b_field):  # B in Tesla, gamma in Hz (w/o 2\pi)
        self.g_factors = g_factors  # dict: hyperfine label -> g factor
        self.b_field = b_field

    def get_detuning(self, ground_state, excited_state):
        from scipy.constants import physical_constants, h
        mu_B = physical_constants['Bohr magneton'][0]  # J / T

        gs_det = mu_B * self.g_factors[ground_state.hyperfine] \
                * ground_state.m * self.b_field / h

        es_det = mu_B * self.g_factors[excited_state.hyperfine] \
                * excited_state.m * self.b_field / h

        return es_det - gs_det


class DopplerDetuning(Detuning):
    def __init__(self, wavelength, velocity):  # wavelength in m, gamma in Hz (w/o 2\pi)
        self.velocity = velocity
        self.k = 1 / wavelength

    def get_detuning(self, ground_state, excited_state):
        # opposite direction light creates positive detuning
        return -1 * self.k * self.velocity
