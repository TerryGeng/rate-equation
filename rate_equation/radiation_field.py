import numpy as np
from collections import namedtuple

RadiationField = namedtuple("RadiationField", ["wavelength", "detuning", "delta_m", "normalized_intensity"])


# some utilities

def wavelength_detuned_by_hz(base_wavelength, *hz_args):
    from scipy.constants import c

    f0 = c / base_wavelength

    f = f0 + sum([c / hz for hz in hz_args])

    return c / f


class RadiationFieldProfile:
    def __init__(self, fields):
        self.fields = fields

    def get_effective_intensity(self, transition, base_frequency, detunings, gamma):
        from scipy.constants import c

        group = transition.group
        delta_m = transition.delta_m

        tot_int = 0

        for field in self.fields:
            if field.delta_m != delta_m:
                continue

            field_freq = c / field.wavelength + field.detuning

            det = base_frequency - field_freq

            det += sum([ det.get_detuning(transition.ground_state, transition.excited_state) \
                    for det in detunings ])

            i_sat_ratio = field.normalized_intensity

            tot_int += np.pi * gamma * i_sat_ratio / (1 + i_sat_ratio + (2*det / gamma)**2)

        return tot_int

