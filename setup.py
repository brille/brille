from skbuild import setup
from setuptools import find_packages 

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name='brille',
    author='Greg Tucker',
    author_email='gregory.tucker@ess.eu',
    description='Irreducible Brillouin zone symmetry and interpolation.',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    extras_require={'interactive': ['matplotlib>=2.2.0', ], },
    url="https://github.com/brille/brille",
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
