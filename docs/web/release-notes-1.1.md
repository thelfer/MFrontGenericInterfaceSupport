% MFrontGenericInterfaceSupport Version 1.1 
% Thomas Helfer
% 17/12/2019

Version 1.1 of `MFrontGenericInterfaceSupport` is compatible with the
Version 3.3 of `TFEL/MFront`.

This version has been used in many projects:

- `Mefisto`. See it in action [here](https://www.youtube.com/watch?v=7YvqhzGagt0).
- [Modelling non-linear constitutive material laws in FEniCS with
  MFront](https://github.com/thelfer/tfel-doc/tree/master/MFrontUserDays/FifthUserDays/FEniCS-Bleyer-MFrontUserDays-2019.pdf),
  by J. Bleyer (Laboratoire Navier, UMR 8205, École des Ponts
  ParisTech-IFSTTAR-CNRS). See also
  [@bleyer_fenics_2019;@bleyer_elasto-plastic_2019]. See it in action
  [here](https://www.youtube.com/watch?v=LyRq6_Q15D0).
- [MFront and OpenGeoSys Connecting two open-source initiatives for
  simulations in environmental geosciences and energy
  geotechnics](https://github.com/thelfer/tfel-doc/tree/master/MFrontUserDays/FifthUserDays/OpenGeoSys-Nagel-MFrontUserDays-2019.pdf),
  by T. Nagel, F. Parisio, D. Naumov, C. Lehmann, and O. Kolditz
  (Technische Universität Bergakademie Freiberg/Helmholtz Zentrum für
  Umweltforschung GmbH/Technische Universität Dresden/Competence Centre
  for Environmental Geosciences). See it action [here](https://www.youtube.com/watch?v=juWMIkJ64iE).
- [Xper: a software dedicated to the fracture of nonlinear heterogeneous
  materials. Coupling with MFront using
  MGIS](https://github.com/thelfer/tfel-doc/tree/master/MFrontUserDays/FifthUserDays/Xper-Perales-MFrontUserDays-2019.pdf),
  by F. Péralès (IRNS).
- `JuliaFEM`. See @frondelius_mfrontinterface_2019 for details.
- `Kratos Multiphysics`. See it in action [here](https://www.youtube.com/watch?v=402rqrygT4k).

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

## Add support for generalised behaviours

Generalised behaviours in `MFront` can be based on multiple tangent
operator blocks.

## Markdown output

The `print_markown` function allows printing using a markdown format:

- a detailled (verbose) description of a behaviour
- a detailled (verbose) description of the data associated with an integration point.
- a detailled (verbose) description of the state of an integration point.

## Specify the type of integration to be performed in the `integrate` function

~~~~{.cxx}
enum struct IntegrationType {
  PREDICTION_TANGENT_OPERATOR = -3,
  PREDICTION_SECANT_OPERATOR = -2,
  PREDICTION_ELASTIC_OPERATOR = -1,
  INTEGRATION_NO_TANGENT_OPERATOR = 0,
  INTEGRATION_ELASTIC_OPERATOR = 1,
  INTEGRATION_SECANT_OPERATOR = 2,
  INTEGRATION_TANGENT_OPERATOR = 3,
  INTEGRATION_CONSISTENT_TANGENT_OPERATOR = 4
};  // end of enum IntegrationType
~~~~

## The `MaterialDataManagerInitializer` and `MaterialStateManagerInitializer` data structures

The `MaterialDataManagerInitializer` data structure is in charge of
holding information on how a material data manager shall be initialized.
It may contain pointers to externally allocated data, that won't be
handled by the final data manager. If a pointer is not initialized, the
material data manager will allocate and handle memory internally.

The `MaterialStateManagerInitializer` data structure is in charge of
holding information on how a material state manager shall be
initialized. It may contain pointers to externally allocated data, that
won't be handled by the final state manager. If a pointer is not
initialized, the material state manager will allocate and handle memory
internally.

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

# References