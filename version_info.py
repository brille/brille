import subprocess
import sys
import platform
import git
from datetime import datetime

def is_git_repo():
    try:
        _ = git.Repo().git_dir
        return True
    except git.exc.InvalidGitRepositoryError:
        return False

def version_info():
    if is_git_repo():
        repo = git.Repo()
        git_branch = repo.active_branch.name
        git_revision = repo.head.commit.hexsha
    else:
        git_branch = "non-git"
        git_revision = ""
    	
    def read_version():
        with open("VERSION") as f:
            return f.readline().strip()

    build_datetime = datetime.now().isoformat(timespec='minutes')
    version_number = read_version()
    hostname = platform.node()
    return git_revision, git_branch, build_datetime, version_number, hostname

def version_number():
    sys.stdout.write(version_info()[3])

if __name__ =="__main__":

    output_file = sys.argv[1]
    with open(output_file, "w") as fout:
        fout.write("""#pragma once
namespace brille{{ namespace version{{
    auto constexpr git_revision = u8"{0}";
    auto constexpr git_branch = u8"{1}";
    auto constexpr build_datetime = u8"{2}";
    auto constexpr version_number = u8"{3}";
    auto constexpr build_hostname = u8"{4}";
}} }}

""".format(*version_info()))
