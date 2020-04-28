=====================================
Interpolating Sodium Chloride phonons
=====================================

Sodium Chloride forms face-centred cubic crystals in Space group
:math:`F m \bar{3} m` (Hall symbol :math:`-F 4 2 3`) with lattice
parameter :math:`a \approx 5.69` Å.

Elemental niobium forms body-centred cubic single crystals in Space group
:math:`I m \bar{3} m` with lattice parameter :math:`a = 3.3004` Å.

:py:mod:`Euphonic` is a project which can take force constant information at
gridded positions in :math:`\mathbf{Q}` and use Fourier interpolation to
approximate the grand dynamical matrix at arbitrary :math:`\mathbf{Q}`.
It can then solve the eigenvalue problem to determine the eigenvalues,
:math:`\omega_i^2(\mathbf{Q})`, and eigenvectors,
:math:`\hat{\epsilon}_{ij}(\mathbf{Q})`, which are the squared energy
and atom-resolved displacement vectors for the
:math:`{i=1,…,3 N_\text{atom}}` phonon branches,
where :math:`{j=1,\ldots,N_\text{atom}}`.

In the specific case of NaCl, the grand dynamical matrix is a
:math:`6\times6` matrix and the eigenvalue problem can be solved quickly
on a modern computer.
Still, it may be advantageous to minimise the work performed by
:py:mod:`Euphonic` and let :py:mod:`brille` do the heavy lifting.

The following is :download:`available in script form<tutorial_00>.py` as well.

The interface between :py:mod:`brille` and :py:mod:`Euphonic` is implemented as
a class, :py:class:`brille.euphonic.BrEu`.



The first step in creating an appropriate interpolation grid object is to load
the force constant data from, e.g., a CASTEP binary file, an appropriate file
for NaCl :download:`obtained here<NaCl.castep_bin>`

.. code:: python

    from pathlib import Path
    import numpy as np
    from euphonic.data.interpolation import InterpolationData
    from brille.euphonic as breu

    SEED = Path('path','to','nb','castep','file')
    ID = InterpolationData.from_castep(str(SEED))
    be = breu(ID, trellis=True, max_volume=0.0001, parallel=False)
