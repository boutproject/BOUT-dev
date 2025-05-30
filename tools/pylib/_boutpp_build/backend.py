#!/bin/python3

import os  # corelib
import glob  # corelib
import hashlib  # corelib
import base64  # corelib
import tempfile  # corelib
import subprocess  # corelib
import re  # corelib
import pathlib  # corelib
import contextlib  # corelib

try:
    import packaging.tags  # packaging
except:
    packaging = None


useLocalVersion = True
pkgname = "boutpp"
version = None


def getversion():
    """
    Get the current BOUT++ version.

    Use a version similar to setuptools_scm but fancier. Similar to cmake's
    versioning.
    """
    global version
    if version is None:
        with contextlib.suppress(KeyError):
            # 0. Check whether version is set via environment variable
            version = os.environ["BOUT_PRETEND_VERSION"]
            return version

        _bout_previous_version = "v5.1.1"
        _bout_next_version = "v5.2.0"

        try:
            try:
                # 1. Check whether we are at a tag
                version = run2("git describe --exact-match --tags HEAD").strip()
            except subprocess.CalledProcessError:
                # 2. default mode, try to derive version from previous tag
                tmp = run2(
                    f"git describe --tags --match={_bout_previous_version}"
                ).strip()
                tmp = re.sub(
                    f"{_bout_previous_version}-", f"{_bout_next_version}.dev", tmp
                )
                if useLocalVersion:
                    tmp = re.sub("-", "+", tmp)
                else:
                    tmp = re.sub("-.*", "", tmp)
                version = tmp
            with open("_version.txt", "w") as f:
                f.write(version + "\n")
        except subprocess.CalledProcessError:
            try:
                # 3. Check whether there is a _version - e.g. we have a tarball
                with open("_version.txt") as f:
                    version = f.read().strip()
            except FileNotFoundError:
                # 4. Maybe not released yet, but version already bumped?
                #    Things are messy here, so always assume useLocalVersion
                try:
                    # 4.1 us proper hash
                    hash = "g" + run2('git log -n 1 --pretty=format:"%h"')
                except subprocess.CalledProcessError:
                    # 4.2 fallback
                    hash = "unknown"
                version = _bout_previous_version + ".rc+" + hash
                with open("_version.txt", "w") as f:
                    f.write(version + "\n")
    return version


def run(cmd):
    """
    Run a command without capturing output, so it gets printed.
    """
    print(f"running `{cmd}`")
    ret = os.system(cmd)
    if ret != 0:
        raise subprocess.CalledProcessError(ret, cmd)


def run2(cmd):
    """
    Run a command and return standard-out
    """
    return subprocess.run(
        cmd, capture_output=True, shell=True, check=True
    ).stdout.decode()


def hash(fn):
    """
    Calculate the hash of a file
    """
    sha256_hash = hashlib.sha256()
    with open(fn, "rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return f"sha256={base64.urlsafe_b64encode(sha256_hash.digest()).decode()[:-1]}"


def size(fn):
    """
    Return the size of a file in bytes
    """
    return os.path.getsize(fn)


def gettag():
    """
    Get the platform tag
    """
    thisos = list(packaging.tags.platform_tags())[-1]
    tag = "-".join(str(next(packaging.tags.sys_tags())).split("-")[:2] + [thisos])
    return tag


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    """
    Build the wheel.

    Calls cmake internaly.
    """
    if metadata_directory:
        parse(f"{metadata_directory}/METADATA")
    opts = ""
    if config_settings is not None:
        global useLocalVersion, pkgname
        for k, v in config_settings.items():
            if k == "sdist":
                continue
            if k == "useLocalVersion":
                useLocalVersion = False
                continue
            if k == "nightly":
                useLocalVersion = False
                pkgname = "boutpp-nightly"
                continue
            if v:
                opts += f" {k}={v}"
            else:
                opts += f" {k}=ON"
    tag = gettag()
    whlname = f"{pkgname.replace('-', '_')}-{getversion()}-{tag}.whl"
    trueprefix = f"{os.getcwd()}/_wheel_install/"
    prefix = f"{trueprefix}/boutpp/"
    run(
        "cmake -S . -B _wheel_build/ -DBOUT_ENABLE_PYTHON=ON"
        + f" -DCMAKE_INSTALL_PREFIX={prefix} -DCMAKE_INSTALL_LIBDIR={prefix}"
        + f" -DCMAKE_INSTALL_PYTHON_SITEARCH={trueprefix} -DCMAKE_INSTALL_RPATH=$ORIGIN"
        + opts
    )
    run(f"cmake --build  _wheel_build/ -j {os.cpu_count()}")
    run("cmake --install _wheel_build/")
    distinfo = f"_wheel_install"
    prepare_metadata_for_build_wheel("_wheel_install", record=True)

    # Do not add --symlink as python's does not extract that as symlinks
    run(f"cd {trueprefix} ; zip  {wheel_directory}/{whlname} . -rq")
    return whlname


def build_sdist(sdist_directory, config_settings=None):
    """
    Create an archive of the code including some metadata files
    """
    print(config_settings, sdist_directory)
    enable_gz = True
    enable_xz = False
    external = {"fmt", "mpark.variant"}
    if config_settings is not None:
        global useLocalVersion, pkgname
        for k, v in config_settings.items():
            if k == "sdist":
                if v == "onlygz":
                    enable_gz = True
                    enable_xz = False
                elif v == "onlyxz":
                    enable_xz = True
                    enable_gz = False
                elif v == "both":
                    enable_xz = True
                else:
                    raise ValueError(f"unknown option {v} for {k}")
            if k == "dist":
                enable_xz = True
                pkgname = "BOUT++"
                external.add("googletest")
            if k == "useLocalVersion":
                useLocalVersion = False
            if k == "nightly":
                useLocalVersion = False
                pkgname = "boutpp-nightly"
    prefix = f"{pkgname.replace('-', '_')}-{getversion()}"
    fname = f"{prefix}.tar"
    run(f"git archive HEAD --prefix {prefix}/ -o {sdist_directory}/{fname}")
    _, tmp = tempfile.mkstemp(suffix=".tar")
    for ext in sorted(external):
        run(
            f"git archive --remote=externalpackages/{ext} HEAD --prefix  {prefix}/externalpackages/{ext}/ --format=tar > {tmp}"
        )
        run(f"tar -Af {sdist_directory}/{fname} {tmp}")
        run(f"rm {tmp}")

    with open(tmp, "w") as f:
        f.write(
            f"""Metadata-Version: 2.1
Name: {pkgname}
Version: {getversion()}
License-File: COPYING
"""
        )
    run(
        f"tar --append -f {sdist_directory}/{fname} _version.txt --xform='s\\_version.txt\\{prefix}/_version.txt\\'"
    )
    run(
        f"tar --append -f {sdist_directory}/{fname} {tmp} --xform='s\\{tmp[1:]}\\{prefix}/PKG-INFO\\'"
    )

    if enable_gz:
        run(f"gzip --force --best --keep {sdist_directory}/{fname}")
        if not enable_xz:
            fname += ".gz"
    if enable_xz:
        run(f"rm {sdist_directory}/{fname}.xz -f")
        run(f"xz --best {sdist_directory}/{fname}")
        fname += ".xz"
    return fname


def get_requires_for_build_sdist(config_settings=None):
    return []


def get_requires_for_build_wheel(config_settings=None):
    return ["packaging", "cython", "jinja2", "numpy"]


def mkdir_p(path):
    return pathlib.Path(path).mkdir(parents=True, exist_ok=True)


def prepare_metadata_for_build_wheel(
    metadata_directory, config_settings=None, record=False
):
    """
    Create dist-info directory with files
    """
    if not record and os.path.isfile("PKG-INFO"):
        parse("PKG-INFO")

    thisdir = f"{pkgname.replace('-', '_')}-{getversion()}.dist-info"
    distinfo = f"{metadata_directory}/{thisdir}"
    mkdir_p(distinfo)
    with open(f"{distinfo}/METADATA", "w") as f:
        f.write(
            f"""Metadata-Version: 2.1
Name: {pkgname}
Version: {getversion()}
License-File: COPYING
"""
        )
    run(f"cp LICENSE {distinfo}/COPYING")
    run(f"cp LICENSE.GPL {distinfo}/COPYING.GPL")
    with open(f"{distinfo}/WHEEL", "w") as f:
        f.write(
            f"""Wheel-Version: 1.0
Generator: boutpp_custom_build_wheel ({getversion()})
Root-Is-Purelib: false
Tag: {gettag()}
"""
        )

    if record:
        with open(f"{distinfo}/RECORD", "w") as f:
            for fn in glob.iglob("_wheel_install/**", recursive=True):
                if not os.path.isfile(fn):
                    continue
                fn0 = fn.removeprefix("_wheel_install/")
                if fn0 != f"{distinfo}/RECORD":
                    f.write(f"{fn0},{hash(fn)},{size(fn)}\n")
                else:
                    f.write(f"{fn0},,\n")
    return thisdir


def parse(fn):
    with open(fn) as f:
        global pkgname, version
        for line in f:
            if line.startswith("Name:"):
                pkgname = line[5:].strip()
            if line.startswith("Version:"):
                version = line[8:].strip()


def nightly():
    """
    Build the python sdist for upload to boutpp-nightly.
    """
    return build_sdist(os.getcwd() + "/dist/", dict(nightly=True))


def sdist():
    """
    Build the python sdist
    """
    return build_sdist(os.getcwd() + "/dist/")


def wheel():
    """
    Build the python binary wheel
    """
    return build_wheel(os.getcwd() + "/dist/")


def dist():
    """
    Build an archive for BOUT++ release
    """
    return build_sdist(os.getcwd(), config_settings=dict(dist=True))


def help():
    """
    Print this help
    """
    table = []
    for k, v in todos.items():
        try:
            doc = v.__doc__.strip()
            doc = " : " + doc
        except:
            doc = ""
        table.append((k, doc))
    maxkey = max([len(k) for k, _ in table])
    fmt = f"   %-{maxkey}s%s"
    print(f"{sys.argv[0]} [command] [command]")
    for row in table:
        print(fmt % row)


def printVersion():
    """
    print the version
    """
    print(getversion())


todos = dict(
    nightly=nightly,
    sdist=sdist,
    wheel=wheel,
    dist=dist,
    version=printVersion,
    help=help,
)
todos.update(
    {
        "--help": help,
        "-?": help,
        "-h": help,
    }
)

if __name__ == "__main__":
    import sys

    for todo in sys.argv[1:]:
        if todo not in todos:
            help()
            sys.exit(1)
        todos[todo]()
