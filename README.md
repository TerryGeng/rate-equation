# `rate_equation`: Rate equation model implemented in python

![Figure: Ground state population vs. laser detuning, under given magnetic
field](./screenshot.png)

_Figure: Ground state population vs. laser detuning, under given magnetic
field_

This package implements rate equation model in python. Rate equation model is
commonly used to simulate the evolution of atomic populations in different
ground states, under given external radiation (e.g. lasers) and magnetic
fields.

Rate equation model assumes the time atoms spend in excited states is
negligible. That is, once excited, they immediately decay back to ground
states. Therefore, all coherence terms are discarded, leaving a set of
equations describing the population transfer from one ground state to others.
This assumption is equivalent to assuming the timescale of interest is much
larger than the decay rate of studied energy levels, while the laser intensity
is way below saturation intensity.

A very good reference for this approach can be found at:
[F. Atoneche and A. Kastberg, Simplified Approach for Quantitative Calculations
of Optical Pumping, Eur. J. Phys. 38, 045703 (2017)](
https://doi.org/10.1088/1361-6404/aa6e6f).

## Features

This package

- Takes hyperfine structure and transition strength as input, normalizes
  transition strength properly and produces the rate equation matrix $G$.

- Supports configurations that include multiple radiation fields with different
  frequencies and polarizations.

- Provides flexible ways of defining detuning terms (Zeeman shift and Doppler
  shift).

- Each part is [individually tested](./test/) against published results to ensure
  correctness.

## Limitations

- This package is designed to work with hyperfine states $I, J, F, m_F$ as the
  eigenstates of the unperturbed Hamiltonian. However, under strong magnetic
  field, hyperfine states are no longer eigenstates and the system is better
  described by $I, J, m_J, m_I$ (Paschen-Back effect). The only way to
  correctly deal with this is to diagonalize the Hamiltonian under all magnetic
  field configurations and use the acquired eigenstates to perform calculation.
  This package is not well geared to achieve this.

## Usage And Examples

In the [example.ipynb](./example.ipynb) IPython notebook, I reproduced all relevant
figures in [Atoneche (2017)](https://doi.org/10.1088/1361-6404/aa6e6f).

