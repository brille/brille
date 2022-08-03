from sphinx.application import Sphinx
from typing import Any, Dict, Optional
from docutils.nodes import document


def _link(text, url):
    return {'type': 'link', 'text': text, 'url': url}


def _button(content, text, icon):
    return {'type': 'group', 'buttons': content, 'icon': icon, 'text': text}


def make_add_launch_buttons(version, version_info):
    
    def add_launch_buttons(app: Sphinx, pagename: str, templatename: str, context: Dict[str, Any], doctree: Optional[document]):
        # grab the current header buttons
        header_buttons = context["header_buttons"]

        # related project links:
        # -- those in the brille organization --
        base = "https://brille.github.io"
        l1 = [_link(name, f"{base}{ext}") for name, ext in (("brille",""), ("brillem","/brillem"), ("brilleu", "/brilleu"))]
        # -- PACE members --
        l1.append(_link('PACE', "https://github.com/pace-neutrons"))
        l1.append(_link('Euphonic', "https://euphonic.readthedocs.io"))
        # add the related project links button
        header_buttons.append(_button(l1, 'Related', 'fa fa-link'))
        
        # (earlier) version links:
        # -- decide if this version is the latest release, a deveopment version, or outdated --
        prerelease = version_info.is_prerelease(version)
        stable = version_info.is_released(version) and version_info.is_latest(version)

        # -- build the list of links to released versions -- 
        variants = [] if stable else [_link('stable', f"{base}/stable")]
        
        # -- get the list of all Major.Minor numbers (after 0.2) --
        releases = version_info.releases(first='0.3')
        # -- add each to the links, making the current version bold (if it is released) --
        for release in releases:
            if release == version:
                variants.append(_link(f"<strong>{release}</strong>", f"{base}/{release}"))
            else:
                variants.append(_link(str(release), f"{base}/{release}"))

        extra = 'latest' if prerelease else 'outdated'
        current = f'{version}' if stable else f'{version} ({extra})'

        # add the links button
        header_buttons.append(_button(variants, current, 'fa fa-history'))

    return add_launch_buttons

