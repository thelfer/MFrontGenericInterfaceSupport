# Julia version of IntegrateTest.py, followed by a check that behaviour
# data remain valid when the Julia object of their behaviour is collected

module mgis
module behaviour
using CxxWrap
@wrapmodule(() -> ENV["MGIS_JULIA_LIBRARY"], :define_mgis_behaviour_module)
function __init__()
  @initcxx
end
end
end

using Test
const mgis_bv = mgis.behaviour

# path to the test library
lib = ENV["MGIS_TEST_BEHAVIOURS_LIBRARY"]
h = mgis_bv.Tridimensional
# reference values
pref = [0, 1.3523277308229e-11, 1.0955374667213e-07, 5.5890770166084e-06,
        3.2392193670428e-05, 6.645865307584e-05, 9.9676622883138e-05,
        0.00013302758358953, 0.00016635821069889, 0.00019969195920296,
        0.00023302522883648, 0.00026635857194317, 0.000299691903777,
        0.0003330252373404, 0.00036635857063843, 0.00039969190397718,
        0.00043302523730968, 0.00046635857064314, 0.00049969190397646,
        0.00053302523730979, 0.00056635857064313]
# comparison criterion
eps = 1.e-12

b = mgis_bv.load(lib, "Norton", h)
d = mgis_bv.BehaviourData(b)
o = mgis_bv.get_variable_offset(mgis_bv.get_internal_state_variables(b),
                                "EquivalentViscoplasticStrain", h)
# strain increment per time step
de = 5.e-5
# time step
mgis_bv.set_time_increment!(d, 180)
# setting the temperature
mgis_bv.set_external_state_variable!(mgis_bv.get_s1(d), "Temperature", 293.15)
# copy d.s1 in d.s0
mgis_bv.update(d)
mgis_bv.get_gradients(mgis_bv.get_s1(d))[1] = de
# equivalent plastic strain
p = [mgis_bv.get_internal_state_variables(mgis_bv.get_s0(d))[o]]
# integrate the behaviour
for i in 1:20
  mgis_bv.integrate(d, b)
  mgis_bv.update(d)
  mgis_bv.get_gradients(mgis_bv.get_s1(d))[1] += de
  push!(p, mgis_bv.get_internal_state_variables(mgis_bv.get_s1(d))[o])
end
# check-in results
@test all(abs.(p - pref) .< eps)

# behaviour data whose behaviour is collected
function make_data()
  b = mgis_bv.load(lib, "Norton", h)
  return mgis_bv.BehaviourData(b)
end
d2 = make_data()
GC.gc(true)
mgis_bv.set_external_state_variable!(mgis_bv.get_s1(d2), "Temperature", 293.15)
@test mgis_bv.get_external_state_variables(mgis_bv.get_s1(d2))[1] == 293.15
