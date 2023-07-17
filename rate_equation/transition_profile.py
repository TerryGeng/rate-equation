import re
from collections import namedtuple
import itertools

# delta_m - in [-1, 0, 1], assuming dipole transition
# strength - can be unnormalized
# ex.
#   trans_g0_e1 = Transition("G0", "E1", "G->E", 1, 25)
State = namedtuple("State", ["hyperfine", "m"])
TransitionGroupLabel = namedtuple("TransitionGroupLabel", ["ground_state_hyperfine", "excited_state_hyperfine"])
Transition = namedtuple("Transition", ["ground_state", "excited_state", "group", "delta_m", "strength"])


def state(label):
    mg = re.match(r"([^-,]+),{0,1}(-{0,1}[\d+])", label)
    assert mg, f"Cannot parse state label {label}"

    hyperfine, m = mg[1], int(mg[2])

    return State(hyperfine, m)


def transition(ground_state_label, excited_state_label, strength):
    gs = state(ground_state_label)
    es = state(excited_state_label)

    return Transition(
            ground_state=gs,
            excited_state=es,
            group=TransitionGroupLabel(gs.hyperfine, es.hyperfine),
            delta_m=es.m-gs.m,
            strength=strength
            )

def group(label):
    mg = re.match(r"(.*)->(.*)", label)
    assert mg, f"Cannot parse transition group {label}"

    return TransitionGroupLabel(mg[1], mg[2])


class TransitionProfile:
    def __init__(self, ground_states, excited_states, transitions, frequencies, gamma):
        from scipy.constants import c

        self.ground_states = ground_states  # ground state labels
        self.excited_states = excited_states  # excited state labels
        self.gamma = gamma  # decay rate, w/o 2 \pi
        self.frequencies = frequencies # Hz, w/o 2 \pi

        # All transition strengths are normalized to make those a
        # excited states sums to 1.
        # Consistently, for cycling transition, the normalized transition
        # strength will be 1.
        # The directly result of this convention is that the pumping Rabi
        # frequency to provide is naturally the Rabi rate of the cycling
        # transition.
        # There are also symmetry (Steck, sodium number eqn. 40) that guarantees
        # the sum transition strengths of any ground state to all possible excited
        # state (vice versa) sum to a fixed number. So this can be a consistency check.
        self._exc_to_gnd_trans = self._normalize_trans_strength(
                self._build_exc_to_gnd_map(transitions))  # indexed by excited state label

        self._gnd_to_exc_trans = self._build_gnd_to_exc_map(
                itertools.chain.from_iterable(self._exc_to_gnd_trans.values()))

    @staticmethod
    def _build_exc_to_gnd_map(transitions):
        _exc_to_gnd_trans = {}
        for trans in transitions:
            if trans.excited_state not in _exc_to_gnd_trans:
                _exc_to_gnd_trans[trans.excited_state] = []

            _exc_to_gnd_trans[trans.excited_state].append(trans)

        return _exc_to_gnd_trans

    @staticmethod
    def _normalize_trans_strength(exc_to_gnd_trans):
        # normalize transition strength, so that for a given excited state,
        # all transitions, who link it to the ground states, sum up to 1
        strength_sum = None
        for exc_state, trans in exc_to_gnd_trans.items():
            _strength_sum = sum([t.strength for t in trans])
            assert strength_sum is None or strength_sum == _strength_sum, \
                    f"Inconsistent transition strength related to {exc_state}."  # see above comment
            strength_sum = _strength_sum

            new_trans = []
            for t in trans:
                new_trans.append(t._replace(strength=t.strength / strength_sum))

            exc_to_gnd_trans[exc_state] = new_trans

        return exc_to_gnd_trans

    @staticmethod
    def _build_gnd_to_exc_map(transitions):
        _gnd_to_exc_trans = {}
        for trans in transitions:
            if trans.ground_state not in _gnd_to_exc_trans:
                _gnd_to_exc_trans[trans.ground_state] = []

            _gnd_to_exc_trans[trans.ground_state].append(trans)

        return _gnd_to_exc_trans

    def get_gnd_to_exc(self, gs):
        assert gs in self.ground_states

        return self._gnd_to_exc_trans[gs] if gs in self._gnd_to_exc_trans else []

    def get_exc_to_gnd(self, es):
        assert es in self.excited_states

        return self._exc_to_gnd_trans[es] if es in self._exc_to_gnd_trans else []
