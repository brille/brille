# Build and test the python module using CIBuildWheel

ENTRYPOINT = r"""#!/bin/bash
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
echo_run {python} -m pip wheel -w /build --no-deps /source

echo_run find /build -type f -iname "*-linux*.whl" -exec sh -c "auditwheel repair '{{}}' -w \$(dirname '{{}}')" \;

"""


def main(variant=None):
    from config import get_client, get_image, get_volumes, get_folder, write_entrypoint
    from python_on_whales import exceptions

    client = get_client()
    image = get_image(client)
    folder = get_folder("wheelhouse")
    volumes = get_volumes(client, {"build": folder, "conan": None})

    write_entrypoint(ENTRYPOINT, folder, variant=variant)
    try:
        result = client.run(
            image, ["sh", "/build/entrypoint.sh"], volumes=volumes, tty=True
        )
        print(result)
    except exceptions.DockerException as ex:
        raise RuntimeError(ex)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-p", "--python", type=str, default="cp312-cp312")
    args = parser.parse_args()
    main(variant=args.python)
