"""Test of an installed mgis wheel, run by cibuildwheel.

On Linux and macOS, a behaviour is compiled with the mfront program of the
tfel-mfront wheel, integrated from python through mgis.behaviour and loaded
through a project using find_package(MFrontGenericInterface) from the
wheel, which checks that the executables, the headers, the cmake package
and the libraries of the two wheels are found and are consistent. On
Windows, building behaviours requires Visual Studio and is not tested.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
from importlib import metadata

import mgis
import mgis.behaviour as mgis_bv
import mgis.function as mgis_fct
import mgis.material_property as mgis_mp
import mgis.model as mgis_model
import tfel  # sets TFELHOME_<version>, used by MFrontGenericInterfaceConfig

# the import of mgis checked that the tfel-mfront wheel this wheel was
# built against is the one installed; check it again explicitly
assert mgis._TFEL_MFRONT_VERSION == metadata.version("tfel-mfront")

mgishome = os.path.dirname(os.path.abspath(mgis.__file__))

# the TFEL libraries must not have been vendored into the mgis wheel:
# they are provided by the tfel-mfront wheel, and a copy under a mangled
# name would make a behaviour load two copies of TFEL
vendored = re.compile(
    r"(lib)?(TFEL|MTest|MFM|MFront(?!GenericInterface))"
    r".*\.(so|dylib|dll)(\.|$)")
for directory, _, files in os.walk(mgishome):
    for name in files:
        if vendored.search(name):
            sys.exit("the mgis wheel must not contain the library " + name)

# the stubs of the python modules are installed next to them, as well as
# the py.typed marker
for stub in ("_mgis.pyi", "behaviour.pyi", "material_property.pyi",
             "model.pyi", "function.pyi", "py.typed"):
    if not os.path.isfile(os.path.join(mgishome, stub)):
        sys.exit("missing " + stub)

if os.name == "nt":
    sys.exit(0)

behaviour = """\
@DSL Implicit;
@Behaviour Elasticity;
@Brick StandardElasticity;
@ElasticMaterialProperties {150e9, 0.3};
"""

library = "libBehaviour.dylib" if sys.platform == "darwin" else "libBehaviour.so"

tfel_version = subprocess.run([shutil.which("tfel-config"), "--version"],
                              check=True, stdout=subprocess.PIPE,
                              text=True).stdout.split()[1]

with tempfile.TemporaryDirectory() as directory:
    os.chdir(directory)
    with open("Elasticity.mfront", "w") as f:
        f.write(behaviour)
    subprocess.run([shutil.which("mfront"), "--obuild", "--interface=generic",
                    "Elasticity.mfront"], check=True)

    # integrate the behaviour from python: with the full strain tensor
    # imposed ({EXX: 1e-3, the other components zero}), the stress is
    # (lambda + 2 mu) EXX for E = 150 GPa and nu = 0.3
    b = mgis_bv.load(os.path.join("src", library), "Elasticity",
                     mgis_bv.Hypothesis.Tridimensional)
    assert b.tfel_version == tfel_version, \
        "the behaviour was built with TFEL " + b.tfel_version + \
        ", not " + tfel_version
    d = mgis_bv.BehaviourData(b)
    d.dt = 1
    mgis_bv.setExternalStateVariable(d.s1, "Temperature", 293.15)
    o = mgis_bv.getVariableOffset(b.gradients, "Strain", b.hypothesis)
    mgis_bv.update(d)
    d.s1.gradients[o] = 1e-3
    mgis_bv.integrate(d, b)
    sxx = float(d.s1.thermodynamic_forces[o])
    e, nu = 150e9, 0.3
    expected = e * (1 - nu) / ((1 + nu) * (1 - 2 * nu)) * 1e-3
    print("SXX:", sxx, "expected:", expected)
    if abs(sxx - expected) > 1e-3 * abs(expected):
        sys.exit("unexpected stress value")

    # the cmake package exported by the wheel must be usable: MFrontGeneric
    # InterfaceConfig.cmake finds TFEL through the TFELHOME_<version>
    # variable set by the import of tfel above
    with open("check.cxx", "w") as f:
        f.write("""\
#include <iostream>
#include <MGIS/Behaviour/Behaviour.hxx>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    return 1;
  }
  const auto b = mgis::behaviour::load(argv[1], argv[2],
                                       mgis::behaviour::Hypothesis::\
TRIDIMENSIONAL);
  std::cout << b.behaviour << std::endl;
  return 0;
}
""")
    with open("CMakeLists.txt", "w") as f:
        f.write("""\
cmake_minimum_required(VERSION 3.20)
project(mgis_wheel_check CXX)
find_package(MFrontGenericInterface REQUIRED)
add_executable(check check.cxx)
target_link_libraries(check PRIVATE mgis::MFrontGenericInterface)
""")
    subprocess.run([shutil.which("cmake"), "-S", ".", "-B", "build",
                    "-DCMAKE_BUILD_TYPE=Release",
                    "-DCMAKE_PREFIX_PATH=" + mgishome], check=True)
    subprocess.run([shutil.which("cmake"), "--build", "build"], check=True)
    check = os.path.join("build", "check" + (".exe" if os.name == "nt" else ""))
    out = subprocess.run([check, os.path.join(directory, "src", library),
                          "Elasticity"], check=True, stdout=subprocess.PIPE,
                         text=True).stdout.strip()
    if out != "Elasticity":
        sys.exit("unexpected output of the cmake check: " + out)

print("the mgis wheel is consistent")
