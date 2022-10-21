import os  # corelib
import glob  # corelib
import hashlib  # corelib
import base64  # corelib
import tempfile  # corelib
import subprocess  # corelib

try:
    import packaging.tags  # packaging
except:
    packaging = None


def run(cmd):
    print(f"running `{cmd}`")
    ret = os.system(cmd)
    assert ret == 0, f"{cmd} failed with {ret}"


useLocalVersion = True
pkgname = "boutpp"


def getversion(_cache={}):
    if "r" not in _cache:
        try:
            tmp = run2("git describe --tags --match=v4.0.0")
            if useLocalVersion:
                tmp = run2(
                    f"echo {tmp.strip()} | sed -e s/v4.0.0-/v5.0.0.dev/ -e s/-/+/"
                )
            else:
                tmp = run2(
                    f"echo {tmp.strip()} | sed -e s/v4.0.0-/v5.0.0.dev/ -e s/-.*//"
                )
            _cache["r"] = tmp
            with open("_version.txt", "w") as f:
                f.write(_cache["r"])
        except AssertionError:
            _cache["r"] = run2("cat _version.txt")
    return _cache["r"].strip()


def run2(cmd):
    child = subprocess.Popen(
        cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True
    )
    output = child.stdout.read().decode("utf-8", "ignore")
    child.communicate()
    assert child.returncode == 0, f"{cmd} failed with {child.returncode}"
    return output


def hash(fn):
    sha256_hash = hashlib.sha256()
    with open(fn, "rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return f"sha256={base64.urlsafe_b64encode(sha256_hash.digest()).decode()[:-1]}"


def size(fn):
    return os.path.getsize(fn)


def gettag():
    thisos = list(packaging.tags.platform_tags())[-1]
    tag = "-".join(str(next(packaging.tags.sys_tags())).split("-")[:2] + [thisos])
    return tag


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    print(config_settings, metadata_directory)
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
    whlname = f"{pkgname}-{getversion()}-{tag}.whl"
    trueprefix = f"{os.getcwd()}/_wheel_install/"
    prefix = f"{trueprefix}/boutpp/"
    run(
        "cmake -S . -B _wheel_build/ -DBOUT_ENABLE_PYTHON=ON "
        + f" -DCMAKE_INSTALL_PREFIX={prefix} -DCMAKE_INSTALL_LIBDIR={prefix} -DCMAKE_INSTALL_PYTHON_SITEARCH={trueprefix}"
        + opts
    )
    run(f"cmake --build  _wheel_build/ -j {os.cpu_count()}")
    run("cmake --install _wheel_build/")
    distinfo = f"_wheel_install"
    prepare_metadata_for_build_wheel("_wheel_install", record=True)

    # Do not add --symlink as python's does not extract that as symlinks
    run(f"cd {trueprefix} ; zip  {wheel_directory}/{whlname} . -rq")
    # cmd = f"git archive HEAD -o {wheel_directory}/{whlname}"
    # run(cmd)
    return whlname


def build_sdist(sdist_directory, config_settings=None):
    print(config_settings, sdist_directory)
    enable_gz = True
    enable_xz = False
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
            if k == "useLocalVersion":
                useLocalVersion = False
            if k == "nightly":
                useLocalVersion = False
                pkgname = "boutpp-nightly"
    prefix = f"{pkgname}-{getversion()}"
    fname = f"{prefix}.tar"
    run(f"git archive HEAD --prefix {prefix}/ -o {sdist_directory}/{fname}")
    _, tmp = tempfile.mkstemp(suffix=".tar")
    for ext in "fmt", "mpark.variant":
        run(
            f"git archive --remote=externalpackages/{ext} HEAD --prefix  {prefix}/externalpackages/{ext}/ --format=tar > {tmp}"
        )
        run(f"tar -Af {sdist_directory}/{fname} {tmp}")
        run(f"rm {tmp}")

    # _, tmpd = tempfile.mkdtemp()
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
    full = ""
    for k in path.split("/"):
        full += k + "/"
        try:
            os.mkdir(full)
        except FileExistsError:
            pass


def prepare_metadata_for_build_wheel(
    metadata_directory, config_settings=None, record=False
):
    thisdir = f"{pkgname}-{getversion()}.dist-info"
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


def nightly():
    return build_sdist(os.getcwd() + "/dist/", dict(nightly=True))


def sdist():
    return build_sdist(os.getcwd() + "/dist/")


def wheel():
    return build_wheel(os.getcwd() + "/dist/")


if __name__ == "__main__":
    import sys

    todos = dict(nightly=nightly, sdist=sdist, wheel=wheel)
    for todo in sys.argv[1:]:
        todos[todo]()
