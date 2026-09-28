"""Griffe extension for brille's documentation.

brille's compiled module is documented from type stubs made by pybind11-stubgen.
Two things in them need adjusting before mkdocstrings renders them:

- An overloaded method or function appears in a stub only as its
  ``@overload`` signatures, which griffe keeps apart from the members, so it
  would not be documented at all. Each becomes a member, with its overloads.
- The docstrings were written for Sphinx, so they use roles such as
  ``:py:meth:`~brille._brille.BZMeshQdc.refine```. These become Markdown
  cross-references, and ``:math:`x``` becomes ``$x$``.
"""
import re

import griffe

ROLE = re.compile(r":py:(?:meth|class|attr|func|mod|obj|data|const|exc):`(~?)([^`]+)`")
MATH = re.compile(r":math:`([^`]+)`")
REF = re.compile(r":ref:`([^`<]+?)(?:\s*<[^>]+>)?`")


def _objects(obj):
    yield obj
    for member in obj.members.values():
        if not member.is_alias:
            yield from _objects(member)


def _crossref(match, obj):
    short, written = match.group(1), match.group(2)
    target = written.rstrip("()")
    module = obj.module
    if "." not in target:
        # as Sphinx resolves it: a member of the enclosing class, else of the module
        cls = obj if obj.is_class else obj.parent
        while cls is not None and not cls.is_class:
            cls = cls.parent
        if cls is not None and target in cls.all_members:
            target = f"{cls.path}.{target}"
        elif target in module.members:
            target = f"{module.path}.{target}"
        else:
            return f"`{written}`"      # not documented here
    elif target.split(".", 1)[0] in module.members and not module.members[target.split(".", 1)[0]].is_alias:
        target = f"{module.path}.{target}"   # relative to the module, e.g. Lattice.pointgroup
    # the dd and cc grid types are shown without members, which are the dc type's
    target = re.sub(r"(BZ(?:Mesh|Trellis|Nest)Q)(?:dd|cc)\.", r"\1dc.", target)
    if target.startswith(obj.package.path + ".") and not _exists(obj.package, target):
        return f"`{written.lstrip('~')}`"    # names something brille no longer has
    shown = target.rsplit(".", 1)[-1] if short or "." not in written else target
    return f"[`{shown}`][{target}]"


def _exists(package, path):
    try:
        package[path[len(package.path) + 1:]]
    except (KeyError, ValueError, griffe.AliasResolutionError, griffe.CyclicAliasError):
        return False
    return True


def _convert(obj):
    if obj.docstring is not None and ":" in obj.docstring.value:
        text = ROLE.sub(lambda m: _crossref(m, obj), obj.docstring.value)
        text = MATH.sub(r"$\1$", text)
        obj.docstring.value = REF.sub(r"\1", text)
        obj.docstring.__dict__.pop("parsed", None)   # parsed before this change, perhaps


class BrilleDocs(griffe.Extension):
    def on_package(self, *, pkg, **kwargs):
        for obj in _objects(pkg):
            if obj.is_class or obj.is_module:
                for name, overloads in list(obj.overloads.items()):
                    if name not in obj.members and overloads:
                        first = overloads[0]
                        function = griffe.Function(
                            name, parameters=first.parameters, returns=first.returns,
                            docstring=first.docstring, lineno=first.lineno, endlineno=first.endlineno,
                        )
                        function.overloads = overloads
                        obj.set_member(name, function)
        # after the members exist, so that cross-references to them resolve
        for obj in _objects(pkg):
            _convert(obj)
            if obj.is_function:
                for overload in obj.overloads or []:
                    if overload.parent is None:
                        overload.parent = obj.parent
                    _convert(overload)
