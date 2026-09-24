"""Degeneracy-aware comparison of modes.

Eigenvectors are only defined up to a phase, and inside a degenerate group
only up to a unitary mixing. Comparing projectors onto each group removes
both ambiguities. The metric parameter (requirement R1 in ``DEFICIENCIES.md``)
is the identity for phonons; for Bogoliubov modes the projector becomes
``Σ v v† η``.
"""
from typing import NamedTuple

import numpy as np


class ModeComparison(NamedTuple):
    """Per-q discrepancies between reference and test modes.

    Attributes
    ----------
    value_error : (Q,) array
        Largest absolute eigenvalue difference.
    projector_error : (Q,) array
        Largest Frobenius-norm difference of any degenerate-group projector.
    norm_error : (Q,) array
        Largest deviation of a test eigenvector's norm from one.
    """

    value_error: np.ndarray
    projector_error: np.ndarray
    norm_error: np.ndarray

    def worst(self):
        return {name: float(np.max(getattr(self, name))) for name in self._fields}


def degenerate_groups(values, rtol=1e-6, atol=1e-8):
    """Split ascending ``values`` into runs of (near-)equal values; returns slices.

    The tolerance is relative to the whole spectrum, not to each value: rounding
    in the eigensolver scales with the largest eigenvalue, so two small values
    can be degenerate to within rounding while differing by more than ``rtol``
    of themselves.
    """
    tol = atol + rtol * np.max(np.abs(values))
    groups, start = [], 0
    for i in range(1, len(values) + 1):
        if i == len(values) or values[i] - values[i - 1] > tol:
            groups.append(slice(start, i))
            start = i
    return groups


def _projector(vectors, metric):
    v = vectors.reshape(len(vectors), -1)
    p = v.T @ v.conj()
    return p if metric is None else p * metric[None, :]


def compare_modes(reference, values, vectors, metric=None, rtol=1e-6, atol=1e-8):
    """Compare test modes against reference modes, q by q.

    Parameters
    ----------
    reference : Modes
        Exact modes; their ascending ``values`` define the degenerate groups.
    values : (Q, B) array
    vectors : (Q, B, N, 3) array
        Test modes, branch-aligned with ``reference``, in the same frame.
    metric : (3N,) array or None
        Inner-product metric; defaults to ``reference.metric``.
    """
    metric = reference.metric if metric is None else metric
    nq = len(reference.values)
    value_error = np.abs(np.asarray(values).reshape(nq, -1) - reference.values).max(axis=1)
    norms = np.linalg.norm(vectors.reshape(nq, vectors.shape[1], -1), axis=-1)
    norm_error = np.abs(norms - 1).max(axis=1)
    projector_error = np.zeros(nq)
    for iq in range(nq):
        for g in degenerate_groups(reference.values[iq], rtol, atol):
            diff = _projector(reference.vectors[iq, g], metric) - _projector(vectors[iq, g], metric)
            projector_error[iq] = max(projector_error[iq], np.linalg.norm(diff))
    return ModeComparison(value_error, projector_error, norm_error)
