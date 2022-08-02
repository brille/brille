# Build and test the CMake project in a manylinux container using a persistent mounted folder

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
# echo_run yum install -y hdf5-devel

# Ensure the python interpreter has setuptools_scm for versioning 
echo_run {python} -m pip install setuptools_scm conan

# Configure the project, specifying which python interpreter and library to use
echo_run export CONAN_USER_HOME=/conan
echo_run cmake -S /source -B /build -DPYTHON_EXECUTABLE={python}

# Build the testing project
echo_run cmake --build /build --target tester --parallel {count}

# Run the tests
echo_run ctest --test-dir /build --parallel {count} --rerun-failed --output-on-failure
"""


def main():
    from config import get_client, get_image, get_volumes, get_folder, write_entrypoint
    from python_on_whales import exceptions

    client = get_client()
    image = get_image(client)
    folder = get_folder('cmake_build')
    volumes = get_volumes(client, {'build': folder, 'conan': None})
    
    write_entrypoint(ENTRYPOINT, folder)
    try:
        result = client.run(image, ['sh', '/build/entrypoint.sh'], volumes=volumes, tty=True)
        print(result)
    except exceptions.DockerException as ex:
        raise RuntimeError()


if __name__ == '__main__':
    main()

