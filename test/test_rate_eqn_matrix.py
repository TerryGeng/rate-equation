import numpy as np

from rate_equation.transition_profile import (Transition, TransitionProfile, State,
                                              state, transition, group)
from rate_equation.rate_equation import RateEquation
from rate_equation.radiation_field import RadiationFieldProfile, RadiationField

class TestTransitions:
    def _create_87Rb_trans(self):
        # 87Rb D2-line, Fg=2 -> Fe=3
        trans = TransitionProfile(
                ground_states=[state(s) for s in ["G2", "G1", "G0", "G-1", "G-2"]],
                excited_states=[state(s) for s in ["E3", "E2", "E1", "E0", "E-1", "E-2", "E-3"]],
                transitions=[
                    # Metcalf, Appendix D
                    transition("G-2", "E-3", 60),
                    transition("G-2", "E-2", 20),
                    transition("G-2", "E-1", 4),
                    transition("G-1", "E-2", 40),
                    transition("G-1", "E-1", 32),
                    transition("G-1", "E0",  12),
                    transition("G0", "E-1",  24),
                    transition("G0", "E0",   36),
                    transition("G0", "E1",   24),
                    transition("G1", "E0",   12),
                    transition("G1", "E1",   32),
                    transition("G1", "E2",   40),
                    transition("G2", "E1",   4),
                    transition("G2", "E2",   20),
                    transition("G2", "E3",   60),
                    ],
                frequencies={
                    group("G->E"): 384.2304844685e12

                    },
                gamma=6.0666e6
                )

        return trans

    def _create_rate_eqn(self):
        # pumping Fg=2 with sigma-plus light to Fe=3
        # Atonche 2017, example 1

        trans = self._create_87Rb_trans()
        fields = RadiationFieldProfile([ RadiationField(
            frequency=384.2304844685e12,
            delta_m=+1,
            normalized_intensity=1
            ) ])
        detunings=[]

        return RateEquation(trans, fields, detunings)

    def test_matrix_element(self):
        rate_eqn = self._create_rate_eqn()

        gamma = 6.0666e6

        assert np.isclose(rate_eqn.calculate_matrix_element(state("G2"), state("G2")),   np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G2"), state("G1")),   np.pi * gamma * 0.5 * 50/225)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G2"), state("G0")),   np.pi * gamma * 0.5 * 6/225)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G2"), state("G-1")),  np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G2"), state("G-2")),  np.pi * gamma * 0.5 * 0)

        assert np.isclose(rate_eqn.calculate_matrix_element(state("G1"), state("G2")),   np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G1"), state("G1")),   np.pi * gamma * 0.5 * -50/225)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G1"), state("G0")),   np.pi * gamma * 0.5 * 48/225)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G1"), state("G-1")),  np.pi * gamma * 0.5 * 9/225)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G1"), state("G-2")),  np.pi * gamma * 0.5 * 0)

        assert np.isclose(rate_eqn.calculate_matrix_element(state("G0"), state("G2")),   np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G0"), state("G1")),   np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G0"), state("G0")),   np.pi * gamma * 0.5 * -54/225)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G0"), state("G-1")),  np.pi * gamma * 0.5 * 27/225)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G0"), state("G-2")),  np.pi * gamma * 0.5 * 6/225)

        assert np.isclose(rate_eqn.calculate_matrix_element(state("G-1"), state("G2")),  np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G-1"), state("G1")),  np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G-1"), state("G0")),  np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G-1"), state("G-1")), np.pi * gamma * 0.5 * -36/225)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G-1"), state("G-2")), np.pi * gamma * 0.5 * 8/225)

        assert np.isclose(rate_eqn.calculate_matrix_element(state("G-2"), state("G2")),  np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G-2"), state("G1")),  np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G-2"), state("G0")),  np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G-2"), state("G-1")), np.pi * gamma * 0.5 * 0)
        assert np.isclose(rate_eqn.calculate_matrix_element(state("G-2"), state("G-2")), np.pi * gamma * 0.5 * -14/225)

    def test_matrix(self):
        rate_eqn = self._create_rate_eqn()
        mat = rate_eqn.build_matrix()

        gamma = 6.0666e6

        truth = np.array(
                [
                    [   0,  50,   6,   0,   0 ],
                    [   0, -50,  48,   9,   0 ],
                    [   0,   0, -54,  27,   6 ],
                    [   0,   0,   0, -36,   8 ],
                    [   0,   0,   0,   0, -14 ],
                ]) / 225

        assert np.allclose(mat, truth*np.pi*gamma*0.5)



