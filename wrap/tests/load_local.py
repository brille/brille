import os
import sys
from pathlib import Path

def found(chain):
    from importlib.util import find_spec as find
    for i in range(len(chain)):
        part = '.'.join(chain[:i+1])
        if find(part) is None:
            return False
    return True

def _load_first_found(modules_with_required_submodules):
    from importlib import import_module
    for module_with_required_submodules in modules_with_required_submodules:
        chain = module_with_required_submodules.split('.')
        if found(chain):
            return import_module(chain[0])
    return None

def load(modules_with_required_submodules, prefer_installed=False, search=None, cmake_mode=True):
    mod = _load_first_found(modules_with_required_submodules) if prefer_installed else None

    if mod is None:
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

        mod = _load_first_found(modules_with_required_submodules)

    if mod is None:
        msg = "None of \n\t" + '\n\t'.join(modules_with_required_submodules)
        msg += " installed or found in\n\t" if prefer_installed else " found in\n\t"
        msg += "\n\t".join([str(locale.absolute()) for locale in search])
        raise ModuleNotFoundError(msg)

    return mod


# use as, e.g.,
#  brille = load((['_brille',], ['brille', 'brille._brille']], search=[Path(), Path('..')], prefer_installed=True)
# to look for:
#       installed _brille  -- never true
#       installed brille with submodule brille._brille
#       local _brille
#       local brille with submodule brille._brille