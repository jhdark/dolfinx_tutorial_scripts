from mpi4py import MPI
from petsc4py import PETSc

import dolfinx
import basix
import numpy as np
from dolfinx import fem
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.mesh import create_mesh
from dolfinx import plot
import ufl

import pyvista


def plot_function(function):

    topology, cell_types, geometry = plot.vtk_mesh(function.function_space)
    u_grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    u_grid.point_data["c"] = function.x.array.real
    u_grid.set_active_scalars("c")
    u_plotter = pyvista.Plotter()

    u_plotter.add_mesh(u_grid, cmap="viridis", show_edges=False)

    u_plotter.show()


### write your script here ###
