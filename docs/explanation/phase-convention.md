# Eigenvector phase convention

`brille` finds eigenvectors at any $\mathbf{Q}$ by rotating and translating those stored in the irreducible Brillouin zone. For phonon eigenvectors ([`RotatesLike.Gamma`][brille._brille.RotatesLike]) this is correct only when the eigenvectors follow the *cell* phase convention. Eigenvectors in the other common convention, the *atom* convention, give wrong eigenvectors, and so wrong structure factors, without any error or warning.

## The two conventions

With force constants $\Phi(0k; l k')$ between atom $k$ in the origin cell and atom $k'$ in cell $l$, lattice vectors $\mathbf{R}_l$ and atom positions $\mathbf{x}_k$, the dynamical matrix is built with one of

$$
\begin{aligned}
D^\text{cell}_{kk'}(\mathbf{q}) &= \sum_l \Phi(0k; lk')\, e^{i\mathbf{q}\cdot\mathbf{R}_l} \\
D^\text{atom}_{kk'}(\mathbf{q}) &= \sum_l \Phi(0k; lk')\, e^{i\mathbf{q}\cdot(\mathbf{R}_l + \mathbf{x}_{k'} - \mathbf{x}_k)}
\end{aligned}
$$

Both give the same frequencies. Their eigenvectors differ by a phase for each atom,

$$
\mathbf{e}^\text{cell}_k(\mathbf{q}) = e^{i\mathbf{q}\cdot\mathbf{x}_k}\, \mathbf{e}^\text{atom}_k(\mathbf{q}),
$$

so the cell-convention eigenvectors are periodic in reciprocal space, $\mathbf{e}^\text{cell}(\mathbf{q}+\mathbf{G}) = \mathbf{e}^\text{cell}(\mathbf{q})$, while the atom-convention ones pick up $e^{-i\mathbf{G}\cdot\mathbf{x}_k}$. `brille` relies on that periodicity.

## Which convention your data use

- [Euphonic](https://github.com/pace-neutrons/Euphonic) uses the cell convention, so its eigenvectors can be given to `brille` directly.
- [phonopy](https://phonopy.github.io/phonopy/) builds its dynamical matrix with the atom convention (checked for phonopy 4.6), so its eigenvectors must be converted first.

For other codes, check how the phase in the dynamical matrix is defined: a phase from lattice vectors alone is the cell convention, and one from vectors between atoms is the atom convention.

## Converting atom-convention eigenvectors

Multiply each atom's part of every eigenvector by $e^{2\pi i\,\mathbf{q}\cdot\mathbf{x}_k}$, with $\mathbf{q}$ in reciprocal lattice units and $\mathbf{x}_k$ in fractional coordinates, before filling the grid:

``` python
import numpy as np

# q: (n_q, 3) grid points in r.l.u.; x: (n_atoms, 3) fractional positions
# atom_vecs: (n_q, n_modes, n_atoms, 3) eigenvectors in the atom convention
phase = np.exp(2j * np.pi * q @ x.T)                # (n_q, n_atoms)
cell_vecs = atom_vecs * phase[:, None, :, None]
```

The frequencies need no change.
