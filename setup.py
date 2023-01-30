import os
import re
import sys
import pkgutil
from sysconfig import get_platform
from subprocess import check_output, check_call
from pathlib import Path

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

if pkgutil.find_loader('packaging') is None:
    from distutils.version import LooseVersion as Version
else:
    from packaging.version import Version

# We can use cmake provided from pip which (normally) gets installed at /bin
# Except that in the manylinux builds it's placed at /opt/python/[version]/bin/
# (as a symlink at least) which is *not* on the path.
# If cmake is a known module, import it and use it to tell us its binary directory
if pkgutil.find_loader('cmake') is not None:
    import cmake

    CMAKE_BIN = cmake.CMAKE_BIN_DIR + os.path.sep + 'cmake'
else:
    CMAKE_BIN = 'cmake'


def get_cmake():
    return CMAKE_BIN


# We want users to be able to specify to *not* use HDF5 for object IO.
# Disable HDF5 IO by passing `--no-hdf` when calling python setup.py.
USE_HDF5 = True


def is_vsc():
    platform = get_platform()
    return platform.startswith("win")


def is_mingw():
    platform = get_platform()
    return platform.startswith("mingw")


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

PACKAGE_ROOT = Path(__file__).absolute().parent

if PACKAGE_ROOT != Path(os.getcwd()):
    raise RuntimeError(f"{PACKAGE_ROOT} != {os.getcwd()}")

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = check_output([get_cmake(), '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build" +
                               " the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        rex = r'version\s*([\d.]+)'
        cmake_version = Version(re.search(rex, out.decode()).group(1))
        if cmake_version < Version('3.18.2'):
            raise RuntimeError("CMake >= 3.18.2 is required")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.dirname(self.get_ext_fullpath(ext.name))
        extdir = os.path.abspath(extdir)
        cmake_args = []
        if is_vsc():
            if sys.maxsize > 2 ** 32:
                cmake_args += ['-A', 'x64']
            else:
                cmake_args += ['-A', 'Win32']

        if is_mingw():
            cmake_args += ['-G', 'Unix Makefiles']  # Must be two entries to work

        cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                       '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        # cfg = 'Debug' if self.debug else 'RelWithDebInfo'
        build_args = ['--config', cfg, '--target', '_brille']

        # make sure all library files end up in one place
        cmake_args += ["-DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE"]
        cmake_args += ["-DCMAKE_INSTALL_RPATH={}".format("$ORIGIN")]

        if not USE_HDF5:
            cmake_args += ["-DBRILLE_HDF5=FALSE"]

        if is_vsc():
            cmake_lib_out_dir = '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'
            cmake_args += [cmake_lib_out_dir.format(cfg.upper(), extdir)]
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '/m:4']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j']

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        check_call([get_cmake(), str(PACKAGE_ROOT)] + cmake_args, cwd=self.build_temp)
        check_call([get_cmake(), '--build', '.', '--target', "_brille"] + build_args, cwd=self.build_temp)


with open(PACKAGE_ROOT.joinpath('README.md'), 'r') as fh:
    LONG_DESCRIPTION = fh.read()

if "--use-hdf5" in sys.argv:
    USE_HDF5 = True
    sys.argv.remove("--use-hdf5")
if "--no-hdf5" in sys.argv:
    USE_HDF5 = False
    sys.argv.remove("--no-hdf5")

setup(
    name='brille',
    author='Greg Tucker',
    author_email='gregory.tucker@ess.eu',
    description='Irreducible Brillouin zone symmetry and interpolation.',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    ext_modules=[CMakeExtension('brille._brille')],
    packages=find_packages(str(PACKAGE_ROOT)),
    extras_require={'plotting': ['matplotlib>=2.2.0', ], 'vis': ['pyglet>=1.5.27', 'vispy>=0.12.1', ]},
    cmdclass=dict(build_ext=CMakeBuild),
    url="https://github.com/brille/brille",
    zip_safe=False,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: Microsoft :: Windows :: Windows 10",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
    ]
)
