"""Give a lattice its symmetry: by name, by generators, or by explicit operations.

The sections between ``--8<-- [start:name]`` and ``[end:name]`` markers are shown
in docs/how-to/symmetry.md, and the tests run this script.
"""
import numpy as np

# --8<-- [start:name]
from brille import BrillouinZone, Lattice

a = 5.69   # angstrom, about NaCl's lattice constant
by_hall = Lattice(([a, a, a], [90, 90, 90]), spacegroup="-F 4 2 3")
by_hermann_mauguin = Lattice(([a, a, a], [90, 90, 90]), spacegroup="Fm-3m")
# --8<-- [end:name]

# --8<-- [start:generators]
from brille import HallSymbol

generators = HallSymbol("-F 4 2 3").generators
by_generators = Lattice(([a, a, a], [90, 90, 90]), symmetry=generators)
print(f"{generators.size} generators make {generators.generate().size} operations")
# --8<-- [end:generators]

# --8<-- [start:explicit]
operations = generators.generate()
rotations, translations = operations.W, operations.w   # (N, 3, 3) integers, (N, 3)
by_operations = Lattice(([a, a, a], [90, 90, 90]), symmetry=(rotations, translations))
# --8<-- [end:explicit]

# --8<-- [start:cif]
from brille import Symmetry

# one generator of Pnma, as it might appear in a CIF file's _symmetry_equiv_pos_as_xyz
mirror = Symmetry("x, 1/2-y, 1/2+z")
same = Symmetry([[[1, 0, 0], [0, -1, 0], [0, 0, 1]]], [[0, 1 / 2, 1 / 2]])
print(f"xyz and matrix forms are the same: {mirror == same}")
orthorhombic = Lattice(([4.0, 5.0, 6.0], [90, 90, 90]), symmetry="x,y,z;-x,-y,-z;x,1/2-y,1/2+z")
# --8<-- [end:cif]

# --8<-- [start:check]
for name, lattice in (("Hall symbol", by_hall), ("Hermann-Mauguin", by_hermann_mauguin),
                      ("generators", by_generators), ("operations", by_operations)):
    zone = BrillouinZone(lattice)
    print(f"{name:16s} {len(lattice.pointgroup.W):3d} point operations, "
          f"irreducible zone {zone.ir_polyhedron.volume:.6f} Å⁻³")
# --8<-- [end:check]

volumes = [BrillouinZone(x).ir_polyhedron.volume for x in (by_hall, by_hermann_mauguin, by_generators, by_operations)]
assert np.allclose(volumes, volumes[0])
assert mirror == same
assert BrillouinZone(orthorhombic).ir_polyhedron.volume > 0
