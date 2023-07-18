import numpy as np

from rate_equation.transition_profile import state
from rate_equation.detuning import ZeemanDetuning, DopplerDetuning


class TestDetuning:
    def test_zeeman(self):
        # Sodium D2 line
        det = ZeemanDetuning(
                g_factors={"G": 1/2, "E": 2/3},
                b_field=0.06
                )

        # This particular number is grabbed from PhysRevLett.48.596 p.598
        assert round(det.get_detuning(None, state("G2"), state("E3")) / 1e6) == 840

    def test_doppler(self):
        det = DopplerDetuning(
                velocity=-1000
                )

        # This particular number is grabbed from PhysRevLett.48.596 p.597
        assert round(det.get_detuning(508.3331958e12, state("G2"), state("E3")) / 1e6) == -1696


