import os
import sys
from pathlib import Path

def found(chain):
    from importlib.util import find_spec as find
    for x in chain:
        if find(x) is None:
            return True
    return False

def load(chains, prefer_installed=False, search=None, cmake_mode=True):
    from importlib import import_module
    if prefer_installed:
        for chain in chains:
            if found(chain):
                return import_module(chain[0])

    # We don't prefer installed modules, or we did not find the installed module.
    # Now look locally
    if search is None:
        search = [Path(), ]
        # search = [Path(), Path('..')]
    if cmake_mode and os.environ.get('CMAKE_CONFIG_TYPE'):
        config = os.environ.get('CMAKE_CONFIG_TYPE')
        for locale in search:
            if Path(locale, config).exists():
                search.append(Path(locale, config))

    sys.path[:0] = [str(locale.absolute()) for locale in search]

    for chain in chains:
        if found(chain):
            return import_module(chain[0])

    msg = "None of \n\t" + '\n\t'.join(chains)
    msg += " installed or found in\n\t" if prefer_installed else " found in\n\t"
    msg += "\n\t".join([str(locale.absolute()) for locale in search])
    raise ModuleNotFoundError(msg)


# use as, e.g.,
#  brille = load((['_brille',], ['brille', 'brille._brille']], search=[Path(), Path('..')], prefer_installed=True)
# to look for:
#       installed _brille  -- never true
#       installed brille with submodule brille._brille
#       local _brille
#       local brille with submodule brille._brille