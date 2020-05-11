# Copyright 2020 Greg Tucker
#
# This file is part of brille.
#
# brille is free software: you can redistribute it and/or modify it under the
# terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# brille is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.
# You should have received a copy of the GNU Affero General Public License
# along with brille. If not, see <https://www.gnu.org/licenses/>.

"""
Read symmetry information from ``CASTEP`` binary files
------------------------------------------------------

.. currentmodule:: brille.castep

.. autosummary::
    :toctree: _generate
"""
# This file borrows heavily from Euphonic's CASTEP binary reader
# https://github.com/pace-neutrons/Euphonic/blob/v0.2.2/euphonic/_readers/_castep.py

import re
import os
import struct
import numpy as np

def read_castep_bin_symmetry(seedname, path):
    """
    Reads symmetry data from a .castep_bin or .check file and returns it in a dictionary

    Parameters
    ----------
    seedname : str
        Seedname of file(s) to read
    path : str
        Path to dir containing the file(s), if in another directory

    Returns
    -------
    data_dict : dict
        A dict with the following keys: 'symmetry_operations', 'symmetry_'disps'
    """
    file = os.path.join(path, seedname + '.castep_bin')
    if not os.path.isfile(file):
        print(
            '{:s}.castep_bin file not found, trying to read {:s}.check'
            .format(seedname, seedname))
        file = os.path.join(path, seedname + '.check')

    with open(file, 'rb') as f:
        int_type = '>i4'
        float_type = '>f8'
        header = ''
        first_sym_read = True
        while header.strip() != b'END':
            header = _read_entry(f)
            if header.find(b'NUM_CRYSTAL_SYMMETRY_OPERATIONS') is 0:
                # There are two sets of symmetry information for some reason
                # hopefully the first one is what we want
                if first_sym_read:
                    n_ops = _read_entry(f, int_type)
                    ops, disps = _read_cell_symmetry(n_ops, f, int_type, float_type)
                    first_sym_read = False

    data_dict = {}
    data_dict['symmetry_operations'] = ops
    data_dict['symmetry_disps'] = disps

    return data_dict


def _read_cell_symmetry(n_ops, file_obj, int_type, float_type):
    """
    Read cell symmetry data from a .castep_bin or .check file

    Parameters
    ----------
    f : file object
        File object in read mode for the .castep_bin or .check file
    int_type : str
        Python struct format string describing the size and endian-ness of
        ints in the file
    float_type : str
        Python struct format string describing the size and endian-ness of
        floats in the file

        Returns
    -------
    n_ops : int
        Number of symmetry operations
    ops : (n_ops, 3, 3) integer ndarray
        The rotation part of each symmetry operation
    disps : (n_ops, 3) float ndarray
        The translation part of each symmetry operation if the cell is symmorphic
        otherwise all entries are zero(?).
    """
    header = b''
    while header.find(b'SYMMETRY_TOL') is -1:
        header = _read_entry(file_obj)
        if header.find(b'CRYSTAL_SYMMETRY_OPERATIONS') is 0:
            ops = np.moveaxis(np.reshape(_read_entry(file_obj, float_type), (n_ops, 3, 3)), 2, 1)
        elif header.find(b'CELL_SYMMORPHIC') is 0:
            is_symmorphic = _read_entry(file_obj, int_type)
        elif header.find(b'CRYSTAL_SYMMETRY_DISPS') is 0:
            divs = np.reshape(_read_entry(file_obj, float_type), (n_ops, 3))
    return ops, divs


def _read_entry(file_obj, dtype=''):
    """
    Read a record from a Fortran binary file, including the beginning
    and end record markers and return the data inbetween

    Parameters
    ----------
    f : file object
        File object in read mode for the Fortran binary file
    dtype : str, optional, default ''
        String determining what order and type to unpack the bytes as. See
        'Format Strings' in Python struct documentation

    Returns
    -------
    data : str, int, float or ndarray
        Data type returned depends on dtype specified. If dtype is not
        specified, return type is a string. If there is more than one
        element in the record, it is returned as an ndarray of floats or
        integers

    """
    def record_mark_read(file_obj):
        # Read 4 byte Fortran record marker
        rawdata = file_obj.read(4)
        if rawdata == b'':
            raise EOFError(
                'Problem reading binary file: unexpected EOF reached')
        return struct.unpack('>i', rawdata)[0]

    begin = record_mark_read(file_obj)
    if dtype:
        n_bytes = int(dtype[-1])
        n_elems = int(begin/n_bytes)
        if n_elems > 1:
            data = np.fromfile(file_obj, dtype=dtype, count=n_elems)
            if 'i' in dtype:
                data = data.astype(np.int32)
            elif 'f' in dtype:
                data = data.astype(np.float64)
            elif '?' in dtype:
                data = data.astype(np.bool)
        else:
            if 'i' in dtype:
                data = struct.unpack('>i', file_obj.read(begin))[0]
            elif 'f' in dtype:
                data = struct.unpack('>d', file_obj.read(begin))[0]
            elif '?' in dtype:
                data = struct.unpack('>?', file_obj.read(begin))[0]
            else:
                data = file_obj.read(begin)
    else:
        data = file_obj.read(begin)
    end = record_mark_read(file_obj)
    if begin != end:
        raise IOError("""Problem reading binary file: beginning and end
                         record markers do not match""")

    return data
