# Build the documentation pages in a reproducible environment

ENTRYPOINT="""#!/bin/sh
set -eu
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ! This file is auto-generated and will be overwritten !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function echo_run {{
    echo "-------------------------------------------------"
	echo "$" "$@"
    echo "-------------------------------------------------"
	"$@"
}}

export PATH={path}:$PATH
export CONAN_USER_HOME=/conan

# Configure the project, specifying which python interpreter and library to use
echo_run echo $CONAN_USER_HOME
echo_run ls -la $CONAN_USER_HOME
echo_run {python} -m pip wheel -w /build --no-deps /source

echo_run find /build -type f -iname "*-linux*.whl" -exec sh -c "auditwheel repair '{{}}' -w \$(dirname '{{}}')" \;

"""

PAGES_ENTRYPOINT="""#!/bin/sh
set -eu
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ! This file is auto-generated and will be overwritten !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function echo_run {{
    echo "-------------------------------------------------"
	echo "$" "$@"
    echo "-------------------------------------------------"
	"$@"
}}

sphinx_doctree=/conan/.doctree
echo_run mkdir -p $sphinx_doctree

html_dir=`mktemp -d`

echo_run git config --global --add safe.directory /source
echo_run cd /source

echo_run ls /wheelhouse
echo_run python3 -m pip install /wheelhouse/{wheel}

echo_run mkdir -p /source/docs/_build/doxygenxml

echo_run sphinx-build -b html /source/docs $html_dir -E -d $sphinx_doctree

echo_run mkdir -p /build/pages
echo_run rsync -a --delete "${{html_dir}}/" /build/pages

"""

def find_wheel(p):
    return [x.name for x in p.iterdir() if x.is_file() and 'cp39-musllinux' in str(x)]


def main(build_wheel, wheelhouse, output):
    from config import get_client, get_image, get_volumes, get_folder, write_entrypoint
    from python_on_whales import exceptions

    client = get_client()
    folder = get_folder(wheelhouse)

    if not build_wheel:
        wheel = find_wheel(folder)
    else:
        wheel = []

    if len(wheel) != 1:
        # Build the musllinux wheel
        image = get_image(client, 'quay.io/pypa/musllinux_1_1_x86_64')
        volumes = get_volumes(client, {'build': folder, 'conan': None})
        print(volumes)
        write_entrypoint(ENTRYPOINT, folder, 'cp39-cp39')
        try:
            result = client.run(image, ['sh', '/build/entrypoint.sh'], volumes=volumes, tty=True)
            print(result)
        except exceptions.DockerException as ex:
            raise RuntimeError()
        wheel = find_wheel(folder)

    # Use the built wheel to build the documentation
    if len(wheel) != 1:
        raise RuntimeError(f"Failed to find single musllinux wheel; found {wheel} instead")
    image = get_image(client, 'docker.io/brille/docact:v5.0.0')
    volumes = get_volumes(client, {'wheelhouse': folder, 'build': get_folder(output), 'conan': None})
    write_entrypoint(PAGES_ENTRYPOINT, folder, sub={'wheel': wheel[0]})
    try:
        result = client.run(image, ['sh', '/wheelhouse/entrypoint.sh'], volumes=volumes, tty=True)
    except exceptions.DockerException:
        raise RuntimeError()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Build image-based Sphinx documentation locally for brille')
    parser.add_argument('-s', '--skip-wheel', action='store_true', help='Skip building the musl linux wheel if it already exists in WHEELHOUSE')
    parser.add_argument('-w', '--wheelhouse', type=str, default='wheelhouse', help='Wheelhouse directory, default="wheelhouse"')
    parser.add_argument('-o', '--output', type=str, default='pages_build', help='Output HTML paged directory, default="pages_build"')
    args = parser.parse_args()

    main(not args.skip_wheel, args.wheelhouse, args.output)

