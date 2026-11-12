from dolfin import *
import ufl
import mgis.fenics as mf
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from gradam import *

behaviours = [
              "FCCMericCailletaudSingleCrystalViscoPlasticity",
              "FCCMericCailletaudFiniteStrainSingleCrystalViscoPlasticity",
             ]

## material parameters
# elasticity constants
E  = 1.
nu = 0.3
G  = E/(1+nu)/2
# initial critical resolved shear stress
tau0 = E/1000.
# viscosity parameters
n_visc = 5.
K_visc = tau0/10.
# isotropic hardening
Q = 1.*tau0
b = 10.
# kinematic hardening
C = 0.
d = 0.

# displacement rate v, final time T and number of loading increments Nincr
v = 1.
T = [.25,.25]
Nincr = 100

# geometry - mesh
xdmf_mesh_file = XDMFFile(MPI.comm_world,"meshes/n20-id1.xdmf")
mesh = Mesh()
xdmf_mesh_file.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 3) 
xdmf_mesh_file.read(mvc, "grains")

grains = cpp.mesh.MeshFunctionSizet(mesh, mvc)
facets = MeshFunction("size_t", mesh, dim=2)
dx = Measure('dx', subdomain_data=grains)

# create boundaries where to impose boundary conditions
class X0(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.) and on_boundary
class X1(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 1.) and on_boundary

X0().mark(facets, 1)
X1().mark(facets, 2)

# create function space
Vu  = VectorFunctionSpace(mesh, "CG", degree=1, dim=3)
Vmf = FunctionSpace(mesh, "DG", 0)

# set the rotation matrix in each grain based on Bunge-Euler angles (ZXZ convention) or rotation matrix columns (1st and 2nd)
def make_rotation_matrix(dim,
                         ori_file,
                         V,
                         mf,
                         convention='from_euler'): 
    '''
    Create the rotation matrix R for each grain. By convention, the rotation is done from sample frame
    to crystal frame i.e. V_crystal = R*V_sample).
    
    Inputs:  * dim: int equal to 2 or 3. If 2, only use phi1 (first angle) to compute the rotation matrix,
                                         while Phi and phi2 are set to 0.
                                         If 3, use phi1, Phi and phi2 to compute the rotation matrix.
             * ori_file: string containing the path to the text file defining Euler angles (in degrees) or
                         1st and 2nd column vectors of the rotation matrix (the 3rd column is deduced).
                         If Euler angles are used, each line contains 3 angles phi1, Phi and phi2 separated
                         by spaces. If vectors are used, the first 3 values define the first column of R, the
                         next threee values define the second column of R. Vectors need not be normalized. The 
                         orientation defined on line 1 defines the rotation matrix of grain with tag 1, ... 
                         orientation defined on line n defines the rotation matrix of grain with tag n. If the
                         file does not exist, it is assumed that the whole mesh has the tag 1 and the default 
                         orientation will be phi1=0, Phi=0, phi2=0.
             * V: FunctionSpace on which the rotation matrix is defined.
             * mf: MeshFunction containing the grain tags.
             * convention: string equal to 'from_euler' or 'from_vector' defining the orientation convention used.

    Outputs: * R: Function in the FunctionSpace V defining the local rotation matrix.
             * phi1: Function in the FunctionSpace V defining the local 1st Bunge-Euler angle.
             * Phi:  Function in the FunctionSpace V defining the local 2nd Bunge-Euler angle.
             * phi2: Function in the FunctionSpace V defining the local 3rd Bunge-Euler angle.
    '''
    # get orientations from file
    if os.path.isfile(ori_file):
        print("Orientation file found!")
        ori = np.loadtxt(ori_file)
    else:
        print("Orientation file not found! using phi1=0, Phi=0, phi2=0")
        ori = np.array([0., 0., 0.])
    if convention=='from_euler':
        if (ori.size == 3): # degenerate case with a single grain
            ori = [ori, ori] # duplicate single crystal orientation in order to get a 2D list as well
    elif convention=='from_vector': 
        if (ori.size == 6): # degenerate case with a single grain
            ori = [ori, ori] # duplicate single crystal orientation in order to get a 2D list as well
        angles = []
        for o in ori:
            V1 = np.array(o[0:2])
            V2 = np.array(o[3:5])
            V1 = V1/np.linalg.norm(V1)
            V2 = V2/np.linalg.norm(V2)
            V3 = np.cross(V1,V2)
            R_ = np.hstack((V1,V2,V3))
            r = Rotation.from_matrix(np.transpose(R_))
            phi = r.as_euler('ZXZ', degrees=True)
            angles.append(phi)
        ori = angles
    else:
        raise Exception("The '%s' convention is not defined." % convention)
    class Angle(UserExpression):
        def __init__(self, mf, eulerAngleIndex, **kwargs):
            super().__init__(**kwargs)
            self.mf = mf
            self.index = eulerAngleIndex
        def eval_cell(self, values, x, cell):
            k = self.mf[cell.index]-1
            values[0] = pi*ori[k][self.index]/180.
        def value_shape(self):
            return () 

    phi1_, Phi_, phi2_ = Angle(mf,0), Angle(mf,1), Angle(mf,2)
    phi1, Phi, phi2 = Function(V, name='phi1'), Function(V, name='Phi'), Function(V, name='phi2')
    phi1.interpolate(phi1_)
    Phi.interpolate(Phi_)
    phi2.interpolate(phi2_)
    if (dim==2):
        R = as_matrix([[ cos(phi1),  sin(phi1), 0.],
                       [-sin(phi1),  cos(phi1), 0.], 
                       [        0.,         0., 1.]])
    else:
        R = as_matrix([[cos(phi1)*cos(phi2)-sin(phi1)*sin(phi2)*cos(Phi), 
                        sin(phi1)*cos(phi2)+cos(phi1)*sin(phi2)*cos(Phi),
                        sin(phi2)*sin(Phi)],
                       [-cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(Phi),
                        -sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(Phi),
                        cos(phi2)*sin(Phi)],
                       [sin(phi1)*sin(Phi),
                        -cos(phi1)*sin(Phi),
                        cos(Phi)]])
    return R, phi1, Phi, phi2   

R, phi1, Phi, phi2 = make_rotation_matrix(3,'meshes/n20-id1.ori',Vmf,grains,convention='from_euler')

# set boundary conditions
loading = Expression("q*t", q=v, t=0, degree=2)
bcs = [ 
       DirichletBC(Vu.sub(0), 0, facets, 1),
       DirichletBC(Vu.sub(0), loading, facets, 2),
       DirichletBC(Vu, Constant((0,0,0)), 
         CompiledSubDomain('near(x[0], 0.)*near(x[1], 0.)*near(x[2], 0.)'), 
         method='pointwise'),
       DirichletBC(Vu.sub(1), Constant(0), 
         CompiledSubDomain('near(x[0], 1.)*near(x[1], 0.)*near(x[2], 0.)'), 
         method='pointwise'),
       DirichletBC(Vu.sub(2), Constant(0), 
         CompiledSubDomain('near(x[0], 1.)*near(x[1], 0.)*near(x[2], 0.)'), 
         method='pointwise'),
       DirichletBC(Vu.sub(1), Constant(0), 
         CompiledSubDomain('near(x[0], 1.)*near(x[1], 0.)*near(x[2], 1.)'), 
         method='pointwise'),
      ]

hypothesis = '3d'
results = np.zeros((2,Nincr+1+11,2))
# solve 3D problem at small and finite strains
for (b,behaviour) in enumerate(behaviours):
    # setup DOFs functions
    u = Function(Vu, name="Displacement")

    # setup problem
    mat_prop = {
                "E": E,
                "nu": nu,
                "G": G,
                "n": n_visc,
                "K": K_visc,
                "Q": Q,
                "b": b,
                "tau0": tau0,
                "d": d,
                "C": C,
               }
    material = mf.MFrontNonlinearMaterial("./src/libBehaviour.so",
                                          behaviour,
                                          hypothesis=hypothesis,
                                          material_properties=mat_prop,
                                          rotation_matrix=R.T)
    problem = mf.MFrontNonlinearProblem(u, material, quadrature_degree=1, bcs=bcs)
    '''
    prm = problem.solver.parameters
    prm['linear_solver'] = 'cg' #'mumps' #'cg'
    tol = 1e-5
    prm["relative_tolerance"] = tol
    prm['absolute_tolerance'] = tol
    '''
    
    ''''''
    problem.solver = PETScSNESSolver('newtonls')
    prm = problem.solver.parameters
    prm['line_search'] =  'nleqerr' 
    prm['linear_solver'] = 'cg' #'mumps'
    prm['preconditioner'] = 'amg'
    prm['krylov_solver']['nonzero_initial_guess'] = False # True
    prm['maximum_iterations'] = 50
    #prm['report'] = False
    tol = 1e-6
    prm["solution_tolerance"] = tol
    prm["relative_tolerance"] = tol
    prm['absolute_tolerance'] = tol
    ''''''

    if behaviour=="FCCMericCailletaudSingleCrystalViscoPlasticity":
        epsel = problem.get_state_variable("ElasticStrain")
        assert (ufl.shape(epsel)==(6, ))
        assert (epsel.ufl_element().family() == "Quadrature")
    elif behaviour=="FCCMericCailletaudFiniteStrainSingleCrystalViscoPlasticity":
        Fp = problem.get_state_variable("PlasticPartOfTheDeformationGradient")
        assert (ufl.shape(Fp) == (9,))
        assert (Fp.ufl_element().family() == "Quadrature")

    file_results = XDMFFile("results/{}_{}_polycrystal_results.xdmf".format(behaviour,hypothesis))
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    
    load_steps = np.hstack( (
                             np.linspace(0., 0.005, 11),
                             np.linspace(0.006, T[b]+0.006, Nincr+1),
                            )
                          )
    
    for (i, t) in enumerate(load_steps):
        if i==0:
            dt = 0
        else:
            dt = load_steps[i] - load_steps[i-1]
        print("Increment: {}, Time = {:6.5f}, dt = {:6.5f}".format(i,t,dt))
        loading.t = t
        problem.dt = dt
        problem.solve(u.vector())

        # save fields for post-processing
        file_results.write(u, t)
        file_results.write(phi1, t)
        file_results.write(Phi, t)
        file_results.write(phi2, t)
        if behaviour=="FCCMericCailletaudSingleCrystalViscoPlasticity":
            epsel  = problem.get_state_variable("ElasticStrain", project_on=("DG", 0), as_tensor=True)
            epspl  = problem.get_state_variable("PlasticStrain", project_on=("DG", 0), as_tensor=True)
            strain = problem.get_gradient("Strain", project_on=("DG", 0), as_tensor=True)
            cauchy = problem.get_flux("Stress", project_on=("DG", 0), as_tensor=True)
            file_results.write(epsel, t)
            file_results.write(epspl, t)
            file_results.write(strain, t)
            file_results.write(cauchy, t)
            eps_11 = assemble(strain[0,0]*dx)
            sig_11 = assemble(cauchy[0,0]*dx)
            results[b,i,:] = (eps_11,sig_11)
        elif behaviour=="FCCMericCailletaudFiniteStrainSingleCrystalViscoPlasticity":
            Fp  = problem.get_state_variable("PlasticPartOfTheDeformationGradient", project_on=("DG", 0), as_tensor=True)
            Fe  = problem.get_state_variable("ShiftedElasticPartOfTheDeformationGradient", project_on=("DG", 0), as_tensor=True)
            F   = problem.get_gradient("DeformationGradient", project_on=("DG", 0), as_tensor=True)
            PK1 = problem.get_flux("FirstPiolaKirchhoffStress", project_on=("DG", 0), as_tensor=True)
            file_results.write(Fp, t)
            file_results.write(Fe, t)
            file_results.write(F, t)
            file_results.write(PK1, t)
            F_11   = assemble(F[0,0]*dx) - 1
            PK1_11 = assemble(PK1[0,0]*dx)
            results[b,i,:] = (F_11,PK1_11)
        for s in range(12):
            PlasticSlip = problem.get_state_variable("PlasticSlip[%s]"%s, project_on=("DG", 0))
            EquivalentViscoplasticSlip = problem.get_state_variable("EquivalentViscoplasticSlip[%s]"%s, project_on=("DG", 0))
            BackStrain = problem.get_state_variable("BackStrain[%s]"%s, project_on=("DG", 0))
            ResolvedShearStress = problem.get_state_variable("ResolvedShearStress[%s]"%s, project_on=("DG", 0))
            file_results.write(PlasticSlip, t)
            file_results.write(EquivalentViscoplasticSlip, t)
            file_results.write(BackStrain, t)
            file_results.write(ResolvedShearStress, t)

        print("\n")

xlabel = r"$\bar{\varepsilon}_{11}$ or $\bar{F}_{11}-1$"
ylabel = r"$\bar{\sigma}_{11}/\tau_0$ or $\bar{S}_{11}/\tau_0$"
plt.figure()
plt.title(r"Polycrystal simple tension") 
plt.plot(results[0,:, 0], results[0,:, 1]/tau0, "-o", markevery=3, label="Small strains")
plt.plot(results[1,:, 0], results[1,:, 1]/tau0, "-s", markevery=3, label="Finite strains")
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.xlim([0.,0.25])
plt.ylim([0.,4.0])
plt.legend()
plt.show()