# Build and test the CMake project in a manylinux container using a persistent mounted folder


def get_client_name():
    from shutil import which
    # python_on_whales supported clients, sorted by our preference
    commands = ('podman', 'docker', 'nerdctl')
    for command in commands:
        if which(command) is not None:
            return command
    raise RuntimeError("One of [docker, podman, nerdctl] required to use virtual-image based builder")


def get_client():
    from python_on_whales import DockerClient
    return DockerClient(client_call=[get_client_name()])


def get_image(client):
    IMAGE = 'quay.io/pypa/manylinux2014_x86_64'
    quiet = client.image.exists(IMAGE)
    img = client.image.pull(IMAGE, quiet=quiet)
    return img


def get_build_folder():
    from pathlib import Path
    build_folder = Path(__file__).parent.joinpath('cmake_build')
    if not build_folder.exists():
        build_folder.mkdir()
    return build_folder


def get_volumes(client):
    from pathlib import Path
    # ensure the cmake directory exists, relative to this file path
    source_folder = Path(__file__).parent.parent
    volumes = [[get_build_folder(), '/build'], [source_folder, '/source']]

    name = client.client_config.client_call[0]
    if 'podman' in name:
        for volume in volumes:
            volume.append('Z')
    # the python_on_whales volumes keyword requires list(tuple(...)):    
    return [tuple(volume) for volume in volumes]


ENTRYPOINT="""#!/bin/bash
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

# Install HDF5 development files
echo_run yum install -y hdf5-devel

# Ensure the python interpreter has setuptools_scm for versioning 
echo_run {python} -m pip install setuptools_scm

# Configure the project, specifying which python interpreter and library to use
echo_run cmake -S /source -B /build -DPYTHON_EXECUTABLE={python}

# Build the testing project
echo_run cmake --build /build --target tester --parallel {count}

# Run the tests
echo_run ctest --test-dir /build --parallel {count} --rerun-failed --output-on-failure
"""


def write_entrypoint():
    from os import cpu_count
    path = "/opt/python/cp310-cp310/bin"
    python = "/opt/python/cp310-cp310/bin/python" 
    count = cpu_count()
    filepath = get_build_folder().joinpath('entrypoint.sh')
    with open(filepath, 'w') as file:
        file.writelines(ENTRYPOINT.format(path=path, python=python, count=1 if count is None else count >> 1))


def main():
    from python_on_whales import exceptions
    write_entrypoint()

    client = get_client()
    image = get_image(client)
    volumes = get_volumes(client)
    
    try:
        result = client.run(image, ['sh', '/build/entrypoint.sh'], volumes=volumes, tty=True)
        print(result)
    except exceptions.DockerException as ex:
        raise RuntimeError()


if __name__ == '__main__':
    main()

