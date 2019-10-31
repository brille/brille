import os
import re
import sys
import subprocess
from sysconfig import get_platform
from subprocess import CalledProcessError, check_output, check_call
from distutils.version import LooseVersion
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from write_version_info import get_version_info

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


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build" +
                               " the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if is_vsc():
            rex = r'version\s*([\d.]+)'
            cmake_version = LooseVersion(re.search(rex, out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.dirname(self.get_ext_fullpath(ext.name))
        extdir = os.path.abspath(extdir)
        cmake_args = []
        if is_vsc():
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            else:
                cmake_args += ['-A', 'Win32']
            # Try to be clever about finding pybind11 on Windows:
            if 'VCPKG' in os.environ:
                vcpkg = os.environ['VCPKG']
                vcpkg += '\\scripts\\buildsystems\\vcpkg.cmake'
                cmake_args += ['-DCMAKE_TOOLCHAIN_FILE='+vcpkg]

        if is_mingw():
            cmake_args += ['-G','Unix Makefiles'] # Must be two entries to work

        cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                       '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        cfg = 'Debug'
        build_args = ['--config', cfg]

        # make sure all library files end up in one place
        cmake_args += ["-DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE"]
        cmake_args += ["-DCMAKE_INSTALL_RPATH={}".format("$ORIGIN")]

        if is_vsc():
            cmake_lib_out_dir = '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'
            cmake_args += [cmake_lib_out_dir.format(cfg.upper(), extdir)]
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '/m:4']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j']

        env = os.environ.copy()
        cxxflags = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''), self.distribution.get_version())
        env['CXXFLAGS'] = cxxflags
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        check_call(
            ['cmake', ext.sourcedir] + cmake_args,
            cwd=self.build_temp, env=env)
        check_call(
            ['cmake', '--build', '.', '--target', "_symbz"] + build_args,
            cwd=self.build_temp)


with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

KEYWORDARGS = dict(
    name='symbz',
    version=get_version_info()[3],
    author='Greg Tucker',
    author_email='greg.tucker@stfc.ac.uk',
    description='First Brillouin zone symmetry and interpolation.',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    ext_modules=[CMakeExtension('symbz._symbz')],
    packages=find_packages(),
    cmdclass=dict(build_ext=CMakeBuild),
    url="https://github.com/g5t/symbz",
    zip_safe=False,
)

try:
    setup(**KEYWORDARGS)
except CalledProcessError:
    print("Failed to build the extension!")
