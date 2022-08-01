

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


def get_folder(name):
    from pathlib import Path
    build_folder = Path(__file__).parent.joinpath(name)
    if not build_folder.exists():
        build_folder.mkdir()
    return build_folder


def get_volumes(client, build_folder):
    from pathlib import Path
    # ensure the cmake directory exists, relative to this file path
    source_folder = Path(__file__).parent.parent
    conan_folder = get_folder('conan')
    volumes = [[build_folder, '/build'], [source_folder, '/source'], [conan_folder, '/conan']]

    name = client.client_config.client_call[0]
    if 'podman' in name:
        for volume in volumes:
            volume.append('Z')
    # the python_on_whales volumes keyword requires list(tuple(...)):    
    return [tuple(volume) for volume in volumes]


def write_entrypoint(entrypoint_string, folder):
    from os import cpu_count
    path = "/opt/python/cp310-cp310/bin"
    python = "/opt/python/cp310-cp310/bin/python" 
    count = cpu_count()
    filepath = folder.joinpath('entrypoint.sh')
    with open(filepath, 'w') as file:
        file.writelines(entrypoint_string.format(path=path, python=python, count=1 if count is None else count >> 1))


