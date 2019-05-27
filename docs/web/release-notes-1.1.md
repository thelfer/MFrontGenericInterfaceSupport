% MFrontGenericInterfaceSupport Version 1.1 
% Thomas Helfer

# New functionalities

## Finite strain behaviour options

By default, finite strain behaviours use:

- the Cauchy stress measure.
- the derivative of the Cauchy stress with respect to the deformation
  gradient as tangent operator.

The `FiniteStrainBehaviourOptions` data structure allows to specify
which stress measure and which tangent operator is expected.

Here is an example of this structure may be used:

~~~~{.cxx}
auto o = FiniteStrainBehaviourOptions{};
o.stress_measure = FiniteStrainBehaviourOptions::PK1;
const auto b =
    load(o, "src/libBehaviour.so", "SaintVenantKirchhoffElasticity",
         Hypothesis::TRIDIMENSIONAL);
~~~~

## `Julia` bindings

~~~~{.julia}
h = mgis.behaviour.Axisymmetrical
b = mgis.behaviour.load("src/libBehaviour.so","Norton2",h)

# meta data
@test mgis.behaviour.get_source(b)=="Norton2.mfront"
@test mgis.behaviour.get_hypothesis(b)==mgis.behaviour.Axisymmetrical
@test mgis.behaviour.get_symmetry(b)==mgis.behaviour.Isotropic
@test mgis.behaviour.get_kinematic(b)==mgis.behaviour.SmallStrainKinematic
@test mgis.behaviour.get_kinematic(b)==mgis.behaviour.SmallStrainKinematic
@test mgis.behaviour.get_behaviour_type(b)== 
    mgis.behaviour.StandardStrainBasedBehaviour
~~~~
