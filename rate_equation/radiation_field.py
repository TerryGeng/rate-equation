import numpy as np
from collections import namedtuple

RadiationField = namedtuple("RadiationField", ["frequency", "delta_m", "normalized_intensity"])


# some utilities

def frequency(base_wavelength):
    from scipy.constants import c

    return c / base_wavelength


class RadiationFieldProfile:
    def __init__(self, fields):
        self.fields = fields

    def get_effective_scattering_rate(self, transition, base_frequency, detunings, gamma):
        from scipy.constants import c

        group = transition.group
        delta_m = transition.delta_m

        tot_Gamma_p = 0

        for field in self.fields:
            if field.delta_m != delta_m:
                continue

            field_freq = field.frequency

            det = base_frequency - field_freq

            det += sum([ det.get_detuning(field_freq, transition.ground_state, transition.excited_state) \
                    for det in detunings ])

            i_sat_ratio = field.normalized_intensity

            tot_Gamma_p += np.pi * gamma * i_sat_ratio / (1 + i_sat_ratio + (2*det / gamma)**2)

        return tot_Gamma_p

