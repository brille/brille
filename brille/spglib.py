import warnings
import numpy as np
import spglib

import brille as b

class BrCell:
    def __init__(self, lattice_vectors, atom_positions, atom_types):
        self.lat = lattice_vectors
        self.pos = atom_positions
        self.typ = atom_types
    def cell(self):
        return (self.lat, self.pos, self.typ)
    def isapprox(self, other):
        return b.Direct(self.lat).isapprox(b.Direct(other.lat))
    def isorthogonal(self):
        return np.allclose(np.diag(np.diag(self.lat)), self.lat)

class BrSpgl:
    def __init__(self, lattice_vectors, atom_positions, atom_types, hall=None):
        if not issubclass(atom_types.dtype.type, np.integer):
            _, atom_types = np.unique(atom_types, return_inverse=True)
        self.input = BrCell(lattice_vectors, atom_positions, atom_types)
        self.primitive = BrCell(*spglib.find_primitive(self.input.cell()))
        self.input_is_primitive = self.input.isapprox(self.primitive)
        if hall is None:
            # Relying on spglib to determine the symmetry of this crystal has
            # the disadvantage that it might introduce an axis rotation, so
            # that even if the input and primitive cells are equivalent the
            # *wrong* Hall symbol will be passed to brille. If the input is not
            # the primitive cell there might *still* be a rotation of axes, but
            # this needs to be investigated.
            self.spglib_dict = spglib.get_symmetry_dataset(self.input.cell())
            self.hall = self.spglib_dict['hall']
            self.spglib_used = True
        else:
            self.hall = hall # no error checking
            self.spglib_used = False

    # pylint: disable=no-member
    def get_primitive_transform(self):
        return b.PrimitiveTransform(self.spglib_dict['hall_number'])

    def transform_needed(self):
        if not self.spglib_used:
            return False
        pt = self.get_primitive_transform()
        return pt.does_anything and self.input_is_primitive

    def spglib_introduced_rotation(self):
        if not self.spglib_used:
            return False
        rot = self.spglib_dict['std_rotation_matrix']
        return not np.allclose(np.diag(np.diag(rot)), rot)

    def spglib_introduced_transformation(self):
        if not self.spglib_used:
            return False
        trn = self.spglib_dict['transformation_matrix']
        return not np.allclose(np.diag(np.diag(trn)), trn)

    def get_input_basis(self):
        return self.input.lat
    def get_spglib_primitive_basis(self):
        return self.primtive.lat
    def get_conventional_basis(self):
        if self.transform_needed():
            # the transformation matrix relies on the spglib-derived
            # symmetry information, so we must use the spglib lattice
            pt = self.get_primitive_transform()
            basis = np.matmul(pt.invPt, self.primitive.lat)
            return np.round(basis, 10) # this might not be necessary for spglib results
        else:
            # Rotating the (x,y,z) representation of the lattice vectors should
            # not make any difference to brille since it only cares about
            # basis vector lengths and mutual angles. A lattice 'transformation'
            # might cause a problem, however.
            if self.spglib_introduced_transformation():
                warnings.warn('Spglib introduced a transformation. Verify before proceeding')
            return self.input.lat

    def get_inverse_input_basis(self):
        return np.linalg.inv(self.get_input_basis())
    def get_inverse_spglib_primitive_basis(self):
        return np.linalg.inv(self.get_spglib_primitive_basis())
    def get_inverse_conventional_basis(self):
        return np.linalg.inv(self.get_conventional_basis())

    def get_input_atom_positions(self):
        return self.input.pos
    def get_conventional_atom_positions(self):
        if self.transform_needed():
            pt = self.get_primitive_transform()
            atom_pos = np.einsum('ai,ij->aj', self.primitive.pos, pt.invPt)
            if pt.does_anything and atom_pos.shape[0]>1:
                warnings.warn('The transformation of atom positions has not been verified')
            return atom_pos
        else:
            if self.spglib_introduced_transformation():
                warnings.warn('Spglib introduced a transformation. Verify before proceeding')
            # Note that a 'standard' rotation doesn't effect the atom positions
            # since they're expressed in units *of the lattice* (i.e., dlu)
            return self.input.pos

    def get_input_atom_index(self):
        return self.input.typ
    def get_conventional_atom_index(self):
        if self.input_is_primitive:
            if not np.allclose(self.input.typ, self.primitive.typ):
                warnings.warn('Spglib introduced an atom permutation. Verify before proceeding')
            return self.primitive.typ
        else:
            return self.input.typ

    def get_input_cell(self):
        return (self.get_input_basis(), self.get_input_atom_positions(), self.get_input_atom_index())
    def get_conventional_cell(self):
        return (self.get_conventional_basis(), self.get_conventional_atom_positions(), self.get_conventional_atom_index())

    def get_input_Direct(self):
        return b.Direct(self.get_input_basis(), self.get_input_atom_positions(), self.get_input_atom_index(), self.hall)
    def get_conventional_Direct(self):
        return b.Direct(self.get_conventional_basis(), self.get_conventional_atom_positions(), self.get_conventional_atom_index(), self.hall)

    def get_input_BrillouinZone(self):
        return b.BrillouinZone(self.get_input_Direct().star)
    def get_conventional_BrillouinZone(self):
        return b.BrillouinZone(self.get_conventional_Direct().star)

    def conventional_to_input_Q(self, qpts):
        if self.transform_needed():
            pt = self.get_primitive_transform()
            qpts = np.einsum('ij,kj->ki', pt.P, qpts)
        # should this be before or after the Conventional to Primitive transformation?
        if self.spglib_introduced_transformation():
            invtrn = np.linalg.inv(self.spglib_dict['transformation_matrix'])
            qpts = np.einsum('ij,ki->kj', invtrn, qpts) # q' = (R⁻¹)ᵀ q
            warnings.warn('Spglib introduced a transformation. Verify before proceeding.')
        return qpts

    def orthogonal_to_conventional_eigenvectors(self, vecs):
        # some quantities are expressed in an orthogonal coordinate system, e.g.
        # eigenvectors, but brille needs them expressed in the units of the
        # conventional lattice which it uses.
        #
        # we need to convert the eigenvector components from units of
        # (x,y,z) to (a,b,c) via the inverse of the basis vectors:
        # For column vectors ⃗x and ⃗a and A ≡ self._basis()
        #        ⃗xᵀ = ⃗aᵀ A
        # which can be inverted to find ⃗a from ⃗x
        #       ⃗a = (A⁻¹)ᵀ ⃗x
        # A⁻¹ is (3,3) and the eigenvectors are (n_pt, n_br, n_io, 3)
        # we want to perform the matrix product of A⁻¹ with vecs[i,j,k,:]
        # which can be most easily accomplished using numpy's einsum
        # which even lets us skip the explicit transpose of A⁻¹
        return np.einsum('ba,ijkb->ijka', self.get_inverse_conventional_basis(), vecs)

    def conventional_to_orthogonal_eigenvectors(self, vecs):
        return np.einsum('ba,ijkb->ijka', self.get_conventional_basis(), vecs)
