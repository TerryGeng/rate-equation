from rate_equation.transition_profile import (Transition, TransitionProfile, State,
                                              state, transition, group)

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
                wavelengths={
                    group("G->E"): 780.241209686e-9
                    },
                gamma=6.0666e6
                )

        return trans

    def _has_transition(self, trans_list, gs, es):
        return any(filter(lambda t: t.ground_state == gs and t.excited_state == es, trans_list))

    def _get_transition(self, g_to_e, gs, es):
        ets = g_to_e[gs]
        it = filter(lambda t: t.excited_state == es, ets)

        try:
            return next(it)
        except StopIteration:
            return None

    def test_parse_state_label(self):
        assert state("G1") == State("G", +1)
        assert state("G-2") == State("G", -2)
        assert state("E1,3") == State("E1", +3)
        assert state("E1,-2") == State("E1", -2)


    def test_exc_to_gnd_trans_map(self):
        trans_map = self._create_87Rb_trans()

        e_to_g = trans_map._exc_to_gnd_trans

        assert len(e_to_g[state("E2")]) == 2
        assert len(e_to_g[state("E-2")]) == 2
        assert self._has_transition(e_to_g[state("E2")], state("G1"), state("E2"))
        assert self._has_transition(e_to_g[state("E2")], state("G2"), state("E2"))

        assert len(e_to_g[state("E0")]) == 3
        assert self._has_transition(e_to_g[state("E0")], state("G-1"), state("E0"))
        assert self._has_transition(e_to_g[state("E0")], state("G0"), state("E0"))
        assert self._has_transition(e_to_g[state("E0")], state("G1"), state("E0"))

    def test_gnd_to_trans_map(self):
        trans_map = self._create_87Rb_trans()

        g_to_e = trans_map._gnd_to_exc_trans

        assert len(g_to_e[state("G0")]) == 3
        assert self._has_transition(g_to_e[state("G0")], state("G0"), state("E-1"))
        assert self._has_transition(g_to_e[state("G0")], state("G0"), state("E0"))
        assert self._has_transition(g_to_e[state("G0")], state("G0"), state("E1"))

        assert len(g_to_e[state("G-2")]) == 3
        assert len(g_to_e[state("G2")]) == 3
        assert self._has_transition(g_to_e[state("G-2")], state("G-2"), state("E-3"))
        assert self._has_transition(g_to_e[state("G-2")], state("G-2"), state("E-2"))
        assert self._has_transition(g_to_e[state("G-2")], state("G-2"), state("E-1"))

    def test_normalization_to_1(self):
        trans_map = self._create_87Rb_trans()

        e_to_g = trans_map._exc_to_gnd_trans

        for trans in e_to_g.values():
            assert sum([t.strength for t in trans]) == 1

    def test_normalization_numbers(self):
        trans_map = self._create_87Rb_trans()

        g_to_e = trans_map._gnd_to_exc_trans

        # Transition(state("G-2"), state("E-3"), -1, 1),
        # Transition(state("G-2"), state("E-2"), 0, 1/3),
        # Transition(state("G-2"), state("E-1"), 1, 1/15),
        assert self._get_transition(g_to_e, state("G-2"), state("E-3")).strength == 1
        assert self._get_transition(g_to_e, state("G-2"), state("E-2")).strength == 1/3
        assert self._get_transition(g_to_e, state("G-2"), state("E-1")).strength == 1/15

        # Transition(state("G-1"), state("E-2"), -1, 2/3),
        # Transition(state("G-1"), state("E-1"), 0, 8/15),
        # Transition(state("G-1"), state("E0"), 1, 1/5),
        assert self._get_transition(g_to_e, state("G-1"), state("E-2")).strength == 2/3
        assert self._get_transition(g_to_e, state("G-1"), state("E-1")).strength == 8/15
        assert self._get_transition(g_to_e, state("G-1"), state("E0")).strength == 1/5

        # Transition(state("G0"), state("E-1"), -1, 2/5),
        # Transition(state("G0"), state("E0"), 0, 3/5),
        # Transition(state("G0"), state("E1"), 1, 2/5),
        assert self._get_transition(g_to_e, state("G0"), state("E-1")).strength == 2/5
        assert self._get_transition(g_to_e, state("G0"), state("E0")).strength == 3/5
        assert self._get_transition(g_to_e, state("G0"), state("E1")).strength == 2/5

        # Transition(state("G1"), state("E0"), -1, 1/5),
        # Transition(state("G1"), state("E1"), 0, 8/15),
        # Transition(state("G1"), state("E2"), 1, 2/3),
        assert self._get_transition(g_to_e, state("G1"), state("E0")).strength == 1/5
        assert self._get_transition(g_to_e, state("G1"), state("E1")).strength == 8/15
        assert self._get_transition(g_to_e, state("G1"), state("E2")).strength == 2/3

        # Transition(state("G2"), state("E1"), -1, 1/15),
        # Transition(state("G2"), state("E2"), 0, 1/3),
        # Transition(state("G2"), state("E3"), 1, 1),
        assert self._get_transition(g_to_e, state("G2"), state("E1")).strength == 1/15
        assert self._get_transition(g_to_e, state("G2"), state("E2")).strength == 1/3
        assert self._get_transition(g_to_e, state("G2"), state("E3")).strength == 1
