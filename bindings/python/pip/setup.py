"""Build of the mgis pip package.

MGIS is configured, built and installed with cmake, then the installation
is laid out as follows:

- the python modules (the mgis package) are placed at the root of the
  package,
- the libraries, headers and data are placed in the bin, lib, include and
  share directories of the mgis module, which is thus the MGISHOME
  directory, as a regular installation of MGIS.

MGIS is built against the tfel-mfront wheel declared in pyproject.toml,
which pip installs in the isolated build environment. This script locates
that installation, mandates it to cmake (require-tfel and
MGIS_REQUIRED_TFEL_VERSION), makes the binaries of the wheel look for the
TFEL libraries in the sibling tfel wheel at run time, and records the
version of tfel-mfront in the package, where importing it checks it.

The TFEL libraries are never copied into this wheel: a behaviour loaded in
python would then load a second copy of TFEL (cf. pyproject.toml).
"""

import importlib.metadata
import importlib.util
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.dist import Distribution

pip_directory = Path(__file__).resolve().parent
mgis_sources = pip_directory.parents[2]

# The tfel-mfront wheel this package is built against and depends on.
# pyproject.toml cannot compute the pin, so it is written twice there, and
# check_tfel_dependency verifies that the two pins and this constant agree.
# Bump the three values together.
TFEL_DEPENDENCY = "tfel-mfront==5.2.0.dev0"

# python code prepended to the __init__.py script of the mgis package:
# MGISHOME, the DLL directories on windows, and the check of the version of
# tfel-mfront. The version of TFEL of the installation the wheel was built
# against is recorded at build time (@TFEL_MFRONT_VERSION@).
pip_setup = '''\
# added by the mgis pip package
import os as _os

_os.environ["MGISHOME"] = _os.path.dirname(_os.path.abspath(__file__))
if _os.name == "nt":
    # since Python 3.8, the directories listed in PATH are not searched for
    # the DLLs needed by the extension modules: the bin directory of this
    # package, where MFrontGenericInterface.dll is installed, and the bin
    # directory of the tfel-mfront wheel, where the TFEL DLLs are, are
    # registered. importing mgis does not necessarily run tfel/__init__.py.
    import importlib.util as _importlib_util

    for _name in ("mgis", "tfel"):
        _spec = _importlib_util.find_spec(_name)
        for _directory in (_spec.submodule_search_locations if _spec else ()):
            _directory = _os.path.join(_directory, "bin")
            if _os.path.isdir(_directory):
                _os.add_dll_directory(_directory)

_TFEL_MFRONT_VERSION = "@TFEL_MFRONT_VERSION@"

# The shared libraries of TFEL have a versioned soname, such as
# libTFELMath.so.5.2.0-dev, so importing a module of this wheel fails with
# a obscure 'libTFELXXX.so.5.2.0-dev: cannot open shared object file' error
# when the installed tfel-mfront does not have the version this wheel was
# built against. Check the version of the installed tfel-mfront package
# instead, which also covers the cases where the libraries are not the ones
# of a wheel (TFELHOME set externally), and which is what cmake verified at
# configuration time (require-tfel and MGIS_REQUIRED_TFEL_VERSION).
def _check_tfel_mfront():
    from importlib import metadata

    try:
        _version = metadata.version("tfel-mfront")
    except metadata.PackageNotFoundError:
        _version = None
    if _version != _TFEL_MFRONT_VERSION:
        raise RuntimeError(
            "this mgis wheel was built against tfel-mfront "
            + _TFEL_MFRONT_VERSION + ", but version " + str(_version)
            + " is installed")


_check_tfel_mfront()
'''


def get_mgis_version():
    """return the MGIS version, as defined in the main CMakeLists.txt"""
    cmakelists = (mgis_sources / "CMakeLists.txt").read_text()
    version = re.search(
        r'project\("mfront-generic-interface"\s+VERSION\s+(\d+)\.(\d+)\.(\d+)',
        cmakelists)
    if version is None:
        raise RuntimeError("MGIS version not found in CMakeLists.txt")
    version = ".".join(version.groups())
    if re.search(
            r"set\(MFrontGenericInterfaceSupportDevelopmentVersion\s+ON\)",
            cmakelists):
        return version + "-dev"
    return version


def check_tfel_dependency():
    """check that the pins of pyproject.toml agree with TFEL_DEPENDENCY"""
    pins = set(re.findall(r'"(tfel-mfront==[^"]+)"',
                          (pip_directory / "pyproject.toml").read_text()))
    if pins != {TFEL_DEPENDENCY}:
        raise RuntimeError(
            "inconsistent TFEL dependency: pyproject.toml pins "
            + str(sorted(pins)) + " while setup.py uses '"
            + TFEL_DEPENDENCY + "'")


def get_tfelhome():
    """path of the tfel module installed in the build environment. The
    package is not imported: find_spec does not load any library"""
    spec = importlib.util.find_spec("tfel")
    if spec is None or not spec.submodule_search_locations:
        raise RuntimeError(
            "the tfel-mfront package, declared in [build-system].requires "
            "of pyproject.toml, is not installed in the build environment")
    return Path(next(iter(spec.submodule_search_locations)))


def get_tfel_version(tfelhome):
    """return the version of TFEL as reported by the tfel-config program of
    the tfel-mfront installation, in the same way as cmake's FindTFEL"""
    tfel_config = tfelhome / "bin" / (
        "tfel-config.exe" if os.name == "nt" else "tfel-config")
    output = subprocess.run([str(tfel_config), "--version"], check=True,
                            stdout=subprocess.PIPE, text=True).stdout
    return output.split()[1]


class BinaryDistribution(Distribution):
    # compiled modules inside: platform-specific wheel, not py3-none-any
    def has_ext_modules(self):
        return True


def build_and_install_mgis(build_directory, install_directory, tfelhome,
                           tfel_version):
    import pybind11

    # copy of pybind11 at a fixed path: the isolated builds install it in a
    # directory whose name changes at each build, which prevents ccache from
    # reusing the compilations of the bindings
    pybind11_directory = Path(build_directory).parent / "pybind11"
    shutil.rmtree(pybind11_directory, ignore_errors=True)
    shutil.copytree(Path(pybind11.__file__).parent, pybind11_directory)

    # every option is set explicitly: the content of the wheel shall not
    # depend on what the build environment happens to provide (HDF5...)
    cmake_options = [
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_INSTALL_PREFIX=" + str(install_directory),
        "-DMGIS_PYTHON_MODULES_INSTALL_DIRECTORY=site-packages/mgis",
        "-DPython_EXECUTABLE=" + sys.executable,
        "-Dpybind11_DIR=" + str(pybind11_directory / "share" / "cmake" /
                                "pybind11"),
        "-Denable-python-bindings=ON",
        "-Denable-pybind11-stubgen=ON",
        "-Denable-portable-build=ON",
        "-Denable-testing=OFF",
        "-Denable-website=OFF",
        "-Denable-doxygen-doc=OFF",
        "-Denable-c-bindings=OFF",
        "-Denable-fortran-bindings=OFF",
        "-Denable-fenics-bindings=OFF",
        "-Denable-julia-bindings=OFF",
        "-Denable-cuda-algorithms=OFF",
        "-Denable-openmp=OFF",
        "-Denable-exceptions=OFF",
        "-Denable-hdf5-support=OFF",
        "-Denable-hdf5-automatic-support=OFF",
        "-Ddisable-tfel=OFF",
        # the wheel can never be built against another TFEL than the one it
        # depends on: cmake checks TFEL_FOUND, MGIS_HAVE_MFRONT_SUPPORT,
        # TFEL_PYTHON_BINDINGS and MGIS_REQUIRED_TFEL_VERSION
        "-Drequire-tfel=ON",
        "-DMGIS_REQUIRED_TFEL_VERSION=" + tfel_version,
        "-DTFEL_DIR=" + str(tfelhome / "share" / "tfel" / "cmake"),
    ]
    if os.name == "nt":
        cmake_options += ["-A", "x64"]
    jobs = os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL", str(os.cpu_count() or 1))
    environment = dict(os.environ, TFELHOME=str(tfelhome))
    subprocess.run(["cmake", "-S", str(mgis_sources), "-B", str(build_directory)]
                   + cmake_options, check=True, env=environment)
    subprocess.run(["cmake", "--build", str(build_directory), "--config",
                    "Release", "--parallel", jobs], check=True, env=environment)
    subprocess.run(["cmake", "--install", str(build_directory), "--config",
                    "Release"], check=True, env=environment)


def is_binary(path, magic_numbers):
    with open(path, "rb") as f:
        return f.read(4) in magic_numbers


def materialize_symlinks(package_directory):
    """replace the symbolic links of the wheel by copies of the files they
    point to: wheels do not store symbolic links (MGIS installs its library
    unversioned, so there is none in practice)"""
    for path in package_directory.rglob("*"):
        if path.is_symlink():
            target = path.resolve()
            path.unlink()
            shutil.copy2(target, path)


def set_relative_runtime_paths(package_directory):
    """make all the binaries of the wheel look for the libraries in
    mgis/lib and in the lib directory of the tfel-mfront wheel.

    The tfel directory is not part of this wheel: it is the python package
    of the tfel-mfront wheel, installed in site-packages next to this one.
    The relative path of a binary to that directory is hence computed as if
    the root of the wheel was a site-packages directory, which it is once
    installed, whereas the build directory of the package is not.

    auditwheel and delocate keep the libraries already in the wheel, and
    the TFEL ones are excluded from the repair, so that the runtime paths
    below are used to find them: see pyproject.toml.

    On linux, DT_RPATH is used rather than DT_RUNPATH, since it is
    searched before LD_LIBRARY_PATH: the binaries of the package use the
    libraries of the package, even if LD_LIBRARY_PATH points to another
    installation, whose libraries may have the same names."""
    fake_site_packages = package_directory.parent / "site-packages"
    if sys.platform.startswith("linux"):
        magic_numbers = (b"\x7fELF",)
    elif sys.platform == "darwin":
        magic_numbers = (b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe")
    else:
        return
    for path in sorted(package_directory.rglob("*")):
        if path.is_symlink() or not path.is_file():
            continue
        if not is_binary(path, magic_numbers):
            continue
        # where the directory of the binary will be once the wheel is
        # installed, simulated in a site-packages directory
        directory = fake_site_packages / path.parent.relative_to(
            package_directory)
        rpaths = []
        for library_directory in (fake_site_packages / "mgis" / "lib",
                                  fake_site_packages / "tfel" / "lib"):
            relative_path = os.path.relpath(library_directory, directory)
            if sys.platform.startswith("linux"):
                rpaths.append("$ORIGIN" if relative_path == "."
                              else "$ORIGIN/" + relative_path)
            elif relative_path != ".":
                rpaths.append("@loader_path/" + relative_path)
        if sys.platform.startswith("linux"):
            subprocess.run(["patchelf", "--force-rpath", "--set-rpath",
                            ":".join(dict.fromkeys(rpaths)), str(path)],
                           check=True)
        else:
            load_commands = subprocess.run(["otool", "-l", str(path)],
                                           check=True, capture_output=True,
                                           text=True).stdout.splitlines()
            current_rpaths = [load_commands[i + 2].split()[1]
                              for i, line in enumerate(load_commands)
                              if line.strip() == "cmd LC_RPATH"]
            options = []
            for rpath in current_rpaths:
                options += ["-delete_rpath", rpath]
            for rpath in rpaths:
                options += ["-add_rpath", rpath]
            subprocess.run(["install_name_tool"] + options + [str(path)],
                           check=True)
            # install_name_tool invalidates the signature of the binary
            subprocess.run(["codesign", "--force", "--sign", "-", str(path)],
                           check=True)


def layout_package(install_directory, package_directory,
                   tfel_mfront_version):
    """lay out the MGIS installation as a python package"""
    shutil.copytree(install_directory / "site-packages", package_directory,
                    dirs_exist_ok=True)
    mgis_module = package_directory / "mgis"
    for directory in ("bin", "lib", "include", "share"):
        if (install_directory / directory).is_dir():
            shutil.copytree(install_directory / directory,
                            mgis_module / directory, dirs_exist_ok=True)
    init_file = mgis_module / "__init__.py"
    init_file.write_text(
        pip_setup.replace("@TFEL_MFRONT_VERSION@", tfel_mfront_version)
        + init_file.read_text())
    materialize_symlinks(package_directory)
    set_relative_runtime_paths(package_directory)


class build_mgis(build_ext):
    """build MGIS with cmake and fill the build directory of the package.

    build_ext is used since it is the only build command run for a
    distribution without python packages declared to setuptools"""

    def run(self):
        check_tfel_dependency()
        tfelhome = get_tfelhome()
        tfel_version = get_tfel_version(tfelhome)
        tfel_mfront_version = importlib.metadata.version("tfel-mfront")
        # with python 3.9 on Windows, resolve keeps relative the paths which
        # do not exist yet
        build_temp = Path(self.build_temp).absolute().resolve()
        install_directory = build_temp / "install"
        build_and_install_mgis(build_temp / "build", install_directory,
                               tfelhome, tfel_version)
        layout_package(install_directory, Path(self.build_lib),
                       tfel_mfront_version)


mgis_version = get_mgis_version()

setup(
    # 3.1.0-dev is not a valid PEP 440 version, 3.1.0.dev0 is
    version=mgis_version.replace("-dev", ".dev0"),
    distclass=BinaryDistribution,
    cmdclass={"build_ext": build_mgis},
)
