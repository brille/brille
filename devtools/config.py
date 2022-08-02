

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


def get_image(client, image=None):
    if image is None:
        image = 'quay.io/pypa/manylinux2014_x86_64'
    quiet = client.image.exists(image)
    return client.image.pull(image, quiet=quiet)


def get_folder(name):
    from pathlib import Path
    folder = Path(__file__).parent.joinpath(name)
    if not folder.exists():
        folder.mkdir()
    return folder


def get_volumes(client, folders: dict):
    from pathlib import Path
    # ensure the cmake directory exists, relative to this file path
    folders['source'] = Path(__file__).parent.parent
    for key in [k for k, v in folders.items() if v is None]:
        folders[key] = get_folder(key)
    volumes = [[v, f"/{k}"] for k, v in folders.items()]
    name = client.client_config.client_call[0]
    if 'podman' in name:
        for volume in volumes:
            volume.append('Z')
    # the python_on_whales volumes keyword requires list(tuple(...)):    
    return [tuple(volume) for volume in volumes]


def write_entrypoint(entrypoint_string, folder, variant=None, sub=None):
    from os import cpu_count
    if variant is None:
        variant = 'cp310-cp310'
    if sub is None:
        sub = {}
    path = f"/opt/python/{variant}/bin"
    python = f"{path}/python" 
    count = cpu_count()
    filepath = folder.joinpath('entrypoint.sh')
    with open(filepath, 'w') as file:
        file.writelines(entrypoint_string.format(path=path, python=python, count=1 if count is None else count >> 1, **sub))


