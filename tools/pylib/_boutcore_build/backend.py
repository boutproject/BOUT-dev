import os  # corelib

try:
    import packaging.tags  # packaging
except:
    packaging = None
import glob  # corelib
import hashlib  # corelib
import base64  # corelib


def run(cmd):
    print(f"running `{cmd}`")
    ret = os.system(cmd)
    assert ret == 0


def getversion():
    return "v5.0.0.alpha.dev10000"


def hash(fn):
    sha256_hash = hashlib.sha256()
    with open(fn, "rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return f"sha256={base64.urlsafe_b64encode(sha256_hash.digest()).decode().removesuffix('=')}"


def size(fn):
    return os.path.getsize(fn)


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    print(config_settings, metadata_directory)
    opts = ""
    for k, v in config_settings.items():
        if v:
            opts += f" {k}={v}"
        else:
            opts += f" {k}=ON"
    print(wheel_directory)
    thisos = list(packaging.tags.platform_tags())[-1]
    tag = "-".join(str(next(packaging.tags.sys_tags())).split("-")[:2] + [thisos])
    whlname = f"boutcore-{getversion()}-{tag}.whl"
    trueprefix = f"{os.getcwd()}/_wheel_install/"
    prefix = f"{trueprefix}/boutcore/"
    run(
        "cmake -S . -B _wheel_build/ -DBOUT_ENABLE_PYTHON=ON "
        + f" -DCMAKE_INSTALL_PREFIX={prefix} -DCMAKE_INSTALL_LIBDIR={prefix} -DCMAKE_INSTALL_PYTHON_SITEARCH={trueprefix}"
        + opts
    )
    run("cmake --build  _wheel_build/ -j 4")
    run("cmake --install _wheel_build/")
    distinfo = f"_wheel_install/boutcore-{getversion()}.dist-info"
    try:
        os.mkdir(distinfo)
    except FileExistsError:
        pass
    with open(f"{distinfo}/METADATA", "w") as f:
        f.write(
            f"""Metadata-Version: 2.1
Name: boutcore
Version: {getversion()}
License-File: COPYING
"""
        )
    run(f"cp LICENSE {distinfo}/COPYING")
    run(f"cp LICENSE.GPL {distinfo}/COPYING.GPL")
    with open(f"{distinfo}/WHEEL", "w") as f:
        f.write(
            f"""Wheel-Version: 1.0
Generator: boutcore_custom_build_wheel ({getversion()})
Root-Is-Purelib: false
Tag: {tag}
"""
        )

    with open(f"{distinfo}/RECORD", "w") as f:
        for fn in glob.iglob("_wheel_install/**", recursive=True):
            if not os.path.isfile(fn):
                continue
            fn0 = fn.removeprefix("_wheel_install/")
            if fn0 != f"{distinfo}/RECORD":
                f.write(f"{fn0},{hash(fn)},{size(fn)}\n")
            else:
                f.write(f"{fn0},,\n")

    run(f"cd {trueprefix} ; zip  {wheel_directory}/{whlname} . -rq --symlinks")
    # cmd = f"git archive HEAD -o {wheel_directory}/{whlname}"
    # run(cmd)
    return whlname


def build_sdist(sdist_directory, config_settings=None):
    print(config_settings)
    print(sdist_directory)
    prefix = f"boutcore-{getversion()}"
    prefix = "boutcore-v5.0.0"
    name = f"{prefix}.tar.gz"
    run(f"git archive HEAD --prefix {prefix}/ -o {sdist_directory}/{name}")
    # run(f"tar --append -f {sdist_directory}/{name} _bout_version.txt")
    return name


def get_requires_for_build_sdist(config_settings=None):
    return []


def get_requires_for_build_wheel(config_settings=None):
    return ["packaging"]
