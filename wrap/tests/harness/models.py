"""Reference models that can be evaluated exactly at any q.

:py:class:`ReferenceModel` is the interface every model implements (requirement
R4 in ``DEFICIENCIES.md``). A model returns :py:class:`Modes`, which carry the
inner-product metric their eigenvectors are orthonormal in (requirement R1):
``None`` means the ordinary inner product, as for phonons; a Bogoliubov model
would return ``η = diag(1,…,1,−1,…,−1)``.
"""
from typing import NamedTuple, Optional, Protocol

import numpy as np
from scipy import sparse

from .crystals import Crystal


class Modes(NamedTuple):
    """Eigen-solutions at a set of q points.

    Attributes
    ----------
    values : (Q, B) real array
        Mode energies, ascending along B at every q.
    vectors : (Q, B, N, 3) complex array
        Eigenvectors, one 3-vector per atom, Cartesian (angstrom frame of the
        crystal lattice), orthonormal under ``metric``.
    metric : (B,) array or None
        Diagonal of the inner-product metric; ``None`` for the identity.
    """

    values: np.ndarray
    vectors: np.ndarray
    metric: Optional[np.ndarray] = None


class ReferenceModel(Protocol):
    """Anything brille's interpolation can be checked against.

    Implementations state how their eigenvectors transform (``rotates_like``,
    one of the names in :py:data:`harness.adapter.ROTATES_LIKE`) and whether
    they obey time-reversal symmetry.
    """

    crystal: Crystal
    rotates_like: str
    time_reversal_invariant: bool

    def modes(self, q_rlu: np.ndarray) -> Modes:
        """Exact modes at ``q_rlu``, shape (Q, 3), in reciprocal lattice units."""
        ...


class BornVonKarman:
    """Phonons of a central- plus transverse-spring model.

    Every pair of atoms closer than ``cutoff`` is joined by a spring whose
    longitudinal and transverse stiffnesses depend only on the distance and on
    the two species. Such a model has at least the symmetry of the crystal,
    and is time-reversal invariant because its force constants are real.

    Parameters
    ----------
    crystal : Crystal
    convention : {"cell", "atom"}
        Phase convention of the dynamical matrix:

        * ``"cell"``: ``D_kk'(q) = Σ_l Φ(0k; l k') exp(2πi q·l)``. Eigenvectors
          are periodic in q.
        * ``"atom"``: an extra ``exp(2πi q·(x_k' − x_k))``. Eigenvectors pick up
          ``exp(2πi G·x_k)`` under ``q → q + G``.
    cutoff : float
        Spring range, in units of the mean interatomic spacing.
    """

    rotates_like = "gamma"
    time_reversal_invariant = True

    def __init__(self, crystal: Crystal, convention="cell", cutoff=1.8, masses=None, stiffness=None):
        if convention not in ("cell", "atom"):
            raise ValueError(f"Unknown convention {convention!r}")
        self.crystal = crystal
        self.convention = convention
        nspecies = crystal.types.max() + 1
        self.masses = np.asarray(masses if masses is not None else 1.0 + 1.3 * np.arange(nspecies))
        c = np.asarray(stiffness if stiffness is not None else 1.0 + 0.3 * np.add.outer(np.arange(nspecies), np.arange(nspecies)))
        self._build(c, cutoff)

    def _build(self, stiffness, cutoff):
        cr = self.crystal
        n = cr.natoms
        a = cr.lattice
        d0 = (abs(np.linalg.det(a)) / n) ** (1 / 3)
        rc = cutoff * d0
        reach = np.ceil(rc * np.linalg.norm(np.linalg.inv(a), axis=0)).astype(int) + 1
        cells = np.array(np.meshgrid(*[np.arange(-m, m + 1) for m in reach], indexing="ij")).reshape(3, -1).T

        k, kp, ls, blocks = [], [], [], []
        eye = np.eye(3)
        for i in range(n):
            dx = cr.positions[None, :, None, :] - cr.positions[i] + cells[None, None, :, :]
            dx = dx[0]  # (N, C, 3)
            r = dx @ a
            d = np.linalg.norm(r, axis=-1)
            j, c = np.nonzero((d > 1e-9) & (d < rc))
            rhat = r[j, c] / d[j, c, None]
            rr = np.einsum("pa,pb->pab", rhat, rhat)
            s = stiffness[cr.types[i], cr.types[j]]
            k_long = s * np.exp(-d[j, c] / d0)
            k_tran = 0.35 * s * np.exp(-1.5 * d[j, c] / d0)
            phi = -(k_long[:, None, None] * rr + k_tran[:, None, None] * (eye - rr))
            k.append(np.full(len(j), i))
            kp.append(j)
            ls.append(cells[c])
            blocks.append(phi)
            # acoustic sum rule: the self term balances every spring on atom i
            k.append([i])
            kp.append([i])
            ls.append(np.zeros((1, 3), dtype=int))
            blocks.append(-phi.sum(axis=0, keepdims=True))
        self._k = np.concatenate(k)
        self._kp = np.concatenate(kp)
        self._l = np.concatenate(ls)
        self._phi = np.concatenate(blocks).reshape(-1, 9)
        inv_sqrt_m = 1 / np.sqrt(self.masses[cr.types])
        self._phi = self._phi * (inv_sqrt_m[self._k] * inv_sqrt_m[self._kp])[:, None]
        shift = self._l.astype(float)
        if self.convention == "atom":
            shift = shift + cr.positions[self._kp] - cr.positions[self._k]
        self._shift = shift
        npairs = len(self._k)
        self._scatter = sparse.csr_matrix(
            (np.ones(npairs), (self._k * n + self._kp, np.arange(npairs))), shape=(n * n, npairs)
        )

    def dynamical_matrix(self, q_rlu):
        q = np.atleast_2d(q_rlu)
        n = self.crystal.natoms
        phase = np.exp(2j * np.pi * (q @ self._shift.T))  # (Q, P)
        d = np.empty((len(q), n * n, 9), dtype=complex)
        for c in range(9):
            d[:, :, c] = (self._scatter @ (phase * self._phi[:, c]).T).T
        d = d.reshape(len(q), n, n, 3, 3).transpose(0, 1, 3, 2, 4).reshape(len(q), 3 * n, 3 * n)
        return 0.5 * (d + d.conj().transpose(0, 2, 1))

    def modes(self, q_rlu):
        w2, u = np.linalg.eigh(self.dynamical_matrix(q_rlu))
        n = self.crystal.natoms
        vectors = u.transpose(0, 2, 1).reshape(len(w2), 3 * n, n, 3)
        return Modes(values=w2, vectors=vectors)
