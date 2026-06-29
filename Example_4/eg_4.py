from mpi4py import MPI
from petsc4py import PETSc

import dolfinx
import basix
import numpy as np
import tqdm.autonotebook
from dolfinx import fem
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.mesh import create_mesh, locate_entities_boundary, meshtags
import ufl

from dolfinx.io import VTXWriter

# define mesh
indices = np.linspace(0, 1, num=100)
gdim, shape, degree = 1, "interval", 1
domain = ufl.Mesh(basix.ufl.element("Lagrange", shape, degree, shape=(gdim,)))
mesh_points = np.reshape(indices, (len(indices), 1))
indexes = np.arange(mesh_points.shape[0])
cells = np.stack((indexes[:-1], indexes[1:]), axis=-1)
my_mesh = create_mesh(comm=MPI.COMM_WORLD, cells=cells, x=mesh_points, e=domain)
fdim = my_mesh.topology.dim - 1

# Define function space and functions
element_CG = basix.ufl.element(
    basix.ElementFamily.P,
    my_mesh.basix_cell(),
    1,
    basix.LagrangeVariant.equispaced,
)
element_DG = basix.ufl.element(
    "DG",
    my_mesh.basix_cell(),
    1,
    basix.LagrangeVariant.equispaced,
)
elements = basix.ufl.mixed_element([element_CG, element_DG])

V = fem.functionspace(my_mesh, elements)
u = fem.Function(V)
u_n = fem.Function(V)

cm, ct = ufl.split(u)
cm_n, ct_n = ufl.split(u_n)
V1, V2 = V.sub(0), V.sub(1)
v_cm, v_ct = ufl.TestFunction(V)
cm_pp, ct_pp = u.sub(0).collapse(), u.sub(1).collapse()
_, map_cm_to_u = V.sub(0).collapse()
_, map_ct_to_u = V.sub(1).collapse()


# define boundary conditions
num_facets = my_mesh.topology.index_map(fdim).size_local
mesh_facet_indices = np.arange(num_facets, dtype=np.int32)
tags_facets = np.full(num_facets, 0, dtype=np.int32)

entities_left = locate_entities_boundary(my_mesh, fdim, lambda x: np.isclose(x[0], 0))
entities_right = locate_entities_boundary(my_mesh, fdim, lambda x: np.isclose(x[0], 1))
tags_facets[entities_left] = 1
tags_facets[entities_right] = 2

facet_meshtags = meshtags(my_mesh, fdim, mesh_facet_indices, tags_facets)
my_mesh.topology.create_connectivity(fdim, my_mesh.topology.dim)

left_facets = facet_meshtags.find(1)
left_dofs_c1 = fem.locate_dofs_topological(V.sub(0), fdim, left_facets)
right_facets = facet_meshtags.find(2)
right_dofs_c1 = fem.locate_dofs_topological(V.sub(0), fdim, right_facets)

bc_left_cm = fem.dirichletbc(
    fem.Constant(my_mesh, PETSc.ScalarType(100)), left_dofs_c1, V1
)
bc_right_cm = fem.dirichletbc(fem.Constant(my_mesh, PETSc.ScalarType(0)), right_dofs_c1, V1)
bcs = [bc_left_cm, bc_right_cm]

# Define variational problem
D = 0.1
dt = 0.05
k = 0.0001
p = 0.01
n = 1

reaction = -k * cm * (n - ct) + p * ct
F = ufl.dot(D * ufl.grad(cm), ufl.grad(v_cm)) * ufl.dx - reaction * v_cm * ufl.dx
F += ((cm - cm_n) / dt) * v_cm * ufl.dx

F += reaction * v_ct * ufl.dx
F += ((ct - ct_n) / dt) * v_ct * ufl.dx


# define solver
petsc_options = {
    "snes_type": "newtonls",
    "snes_linesearch_type": "none",
    "snes_stol": np.sqrt(np.finfo(dolfinx.default_real_type).eps)
    * 1e-2,
    "snes_atol": 1e-10,
    "snes_rtol": 1e-10,
    "snes_max_it": 30,
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
}
solver = NonlinearProblem(
    F,
    u,
    bcs=bcs,
    petsc_options=petsc_options,
    petsc_options_prefix="festim_solver",
)
snes = solver.solver
prefix = snes.getOptionsPrefix()
opts = PETSc.Options()
for k in petsc_options.keys():
    del opts[f"{prefix}{k}"]

writer1 = VTXWriter(MPI.COMM_WORLD, "trapping_cm.bp", cm_pp, "BP5")
writer2 = VTXWriter(MPI.COMM_WORLD, "trapping_ct.bp", ct_pp, "BP5")

final_time = 10
t = 0
progress = tqdm.autonotebook.tqdm(
        desc="Solving H transport problem", total=final_time, unit_scale=True
    )
while t < final_time:
    solver.solve()

    u_n.x.array[:] = u.x.array

    cm_pp.x.array[:] = u.x.array[map_cm_to_u]
    ct_pp.x.array[:] = u.x.array[map_ct_to_u]
    
    writer1.write(t)
    writer2.write(t)

    t += dt

    progress.update(dt)

writer1.close()
writer2.close()
