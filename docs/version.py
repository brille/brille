import sys
from typing import List
from packaging.version import parse, Version, LegacyVersion
import requests
import argparse


class V:
    def __init__(self, v: Version):
        self.major = v.major
        self.minor = v.minor

    def __str__(self):
        return f"v{self.major}.{self.minor}"

    def __repr__(self):
        return f"V({self.major}, {self.minor})"

    def _to_V(self, v):
        if isinstance(v, str):
            v = parse(v)
        return V(v) if isinstance(v, Version) else v

    def __eq__(self, other):
        other = self._to_V(other)
        return (self.major, self.minor) == (other.major, other.minor)

    def __lt__(self, other):
        other = self._to_V(other)
        return (self.major, self.minor) < (other.major, other.minor)
    
    def __le__(self, other):
        other = self._to_V(other)
        return (self.major, self.minor) <= (other.major, other.minor)
    
    def __hash__(self):
        return hash((self.major, self.minor))


def _get_releases(repo: str, organization: str = 'brille') -> List[Version]:
    """Return reversed sorted list of release tag names."""
    r = requests.get(f'https://api.github.com/repos/{organization}/{repo}/releases')
    if r.status_code != 200:
        return []
    data = r.json()
    return sorted([parse(e['tag_name']) for e in data if not e['draft']], reverse=True)


class VersionInfo:

    def __init__(self, repo: str, organization: str = 'brille'):
        self._releases = _get_releases(repo=repo, organization=organization)

    def _to_version(self, version) -> Version:
        if isinstance(version, str):
            version = parse(version)
            # When not building for a tagged release we may get, e.g., 'main'.
            # Pretend this means the current latest release.
            if isinstance(version, LegacyVersion):
                return self._releases[0]
        return version

    def minor_releases(self, first: str = '0.1') -> List[str]:
        """Return set of minor releases in the form '1.2'.

        `first` gives the first release to be included. By default '0.0' releases are
        not included.
        """
        first = parse(first)
        releases = [r for r in self._releases if r >= first]
        releases = sorted(set((r.major, r.minor) for r in releases), reverse=True)
        return [f'{major}.{minor}' for major, minor in releases]

    def releases(self, first=None) -> List[V]:
        if first is None:
            first = parse('0.0.1')
        if not isinstance(first, Version):
            first = parse(first)
        rs = [r for r in self._releases if r >= first]
        return sorted(set(V(r) for r in rs), reverse=True)

    def is_latest(self, version: str) -> bool:
        """Return True if `version` has the same or larger major and minor as the
        latest release.
        """
        version = self._to_version(version)
        latest = self._releases[0]
        return (latest.major, latest.minor) <= (version.major, version.minor)

    def is_new(self, version: str) -> bool:
        """Return True if `version` is a new major or minor release."""
        version = self._to_version(version)
        releases = [r for r in self._releases if r != version]
        latest = releases[0]
        return (latest.major, latest.minor) < (version.major, version.minor)

    def is_released(self, version: str) -> bool:
        """Return True if `version` is in the known releases"""
        version = self._to_version(version)
        return any([r == version for r in self._releases])

    def is_prerelease(self, version: str) -> bool:
        """Return True if `version` is not in the known releases and it is greater than the latest"""
        return self._to_version(version) > self._releases[0] and not self.is_released(version)

    def target(self, version: str) -> str:
        version = self._to_version(version)
        if self.is_latest(version):
            return ''
        else:
            return f'release/{version.major}.{version.minor}'

    def replaced(self, version: str) -> str:
        version = self._to_version(version)
        # If we release 1.2 we want to find 1.1
        for release in self._releases:
            if (release.major, release.minor) < (version.major, version.minor):
                return release


def main(repo: str, action: str, version: str) -> int:
    info = VersionInfo(repo=repo)
    if action == 'is-latest':
        print(info.is_latest(version))
    elif action == 'is-new':
        print(info.is_new(version))
    elif action == 'get-replaced':
        print(info.replaced(version))
    elif action == 'get-target':
        print(info.target(version))
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--repo', dest='repo', required=True, help='Repository name')
    parser.add_argument('--action',
                        choices=['is-latest', 'is-new', 'get-replaced', 'get-target'],
                        required=True,
                        help='Action to perform: Check whether this major or minor '
                        'release exists or is new (is-latest), check whether this is a '
                        'new major or minor release (is-new), get the version this is '
                        'replacing (get-replaced), get the target folder for '
                        'publishing the docs (get-target). In all cases the '
                        'patch/micro version is ignored.')
    parser.add_argument('--version',
                        dest='version',
                        required=True,
                        help='Version the action refers to')

    args = parser.parse_args()
    sys.exit(main(**vars(args)))
