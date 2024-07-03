import subprocess
import sys
import platform
from datetime import datetime


def git_revision():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("utf-8") .split("\n")[0]
    except (subprocess.CalledProcessError, OSError):
        return ""


def git_branch():
    try:
        return subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]).decode("utf-8").split("\n")[0]
    except (subprocess.CalledProcessError, OSError):
        return "non-git"


def build_datetime():
    return datetime.now().isoformat(timespec='minutes')


def version_number():
    with open("VERSION") as file:
        return file.readline().strip()


def hostname():
    return platform.node()


def version_info():
    return git_revision(), git_branch(), build_datetime(), version_number(), hostname()


if __name__ == "__main__":

    output_file = sys.argv[1]
    with open(output_file, "w") as file_out:
        file_out.write(f"""#ifndef POLYSTAR_VERSION_HPP_
#define POLYSTAR_VERSION_HPP_
//! \\file
namespace polystar::version{{
    //! `polystar` git repository revision information at build time
    auto constexpr git_revision = u8"{git_revision()}";
    //! `polystar` git repository branch at build time
    auto constexpr git_branch = u8"{git_branch()}";
    //! build date and time in YYYY-MM-DDThh:mm format
    auto constexpr build_datetime = u8"{build_datetime()}";
    //! `polystar` version
    auto constexpr version_number = u8"{version_number()}";
    //! hostname of the build machine
    auto constexpr build_hostname = u8"{hostname()}";
}}
#endif
""")
