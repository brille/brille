import subprocess
import sys
import platform
from datetime import datetime

def version_info():
    try:
        git_revision = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("utf-8") .split("\n")[0]
        git_branch = subprocess.check_output(["git", "rev-parse","--abbrev-ref", "HEAD"]).decode("utf-8").split("\n")[0]
    except (subprocess.CalledProcessError, OSError):
        git_revision = ""
        git_branch = "non-git"

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
        fout.write("""#ifndef BRILLE_VERSION_HPP_
#define BRILLE_VERSION_HPP_
//! \\file
namespace brille::version{{
    //! `brille` git repository revision information at build time
    auto constexpr git_revision = u8"{0}";
    //! `brille` git repository branch at build time
    auto constexpr git_branch = u8"{1}";
    //! build date and time in YYYY-MM-DDThh:mm format
    auto constexpr build_datetime = u8"{2}";
    //! `brille` version
    auto constexpr version_number = u8"{3}";
    //! hostname of the build machine
    auto constexpr build_hostname = u8"{4}";
}}
#endif
""".format(*version_info()))
