# -*- coding: utf-8 -*-

#-----------------------------------------------------------------
#
# name: thomas brunner
# email: brunner.th@hotmail.com
# nr: 12018550
#
#-----------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from meshpy.tet import MeshInfo, build
from mpl_toolkits.mplot3d import Axes3D
import meshpy.triangle as triangle
from meshing import UnitSquareMesh, CSVToMesh
from geometry import *
from local_system_assembling import *
from global_assembler import *
import matplotlib.tri as mtri
import scipy.interpolate as ip 
import boundary_condition_classes as BC
plt.rcParams["figure.dpi"] = 300
plt.rcParams["figure.figsize"] = (10, 7)

############################## mesh generation / import #######################

mesh_path = r"C:\Users\brunn\Documents\GitHub\FEMPy\Mesh_files/UnitSquareVeryFine.csv"

p_unitsquare = CSVToMesh(mesh_path)

points = p_unitsquare[0]
simplices = p_unitsquare[1]
hull = p_unitsquare[2]

############################## geometry handling ##############################

Anode = BoundaryPointsInRectangle(0.3,0.35,0.3,0.7, points)
Kathode = BoundaryPointsInRectangle(0.65,0.7,0.3,0.7, points)
Dielectric = BoundaryPointsInRectangle(0.35,0.65,0.3,0.7, points)

############################## parameters #####################################

k1_list = np.zeros((len(points)))
k2_list = np.zeros((len(points)))
rho_list = np.zeros((len(points)))
f_list = np.zeros((len(points)))

k1_list.fill(1)
k2_list.fill(1)

################################ boundary conditions ##########################

Dirichlet_BC = BC.DirichletBoundaryCondition(len(points))
Cauchy_BC = BC.CauchyBoundaryCondition(len(points))

Dirichlet_BC.add_dirichlet_points(Anode, 5)
Dirichlet_BC.add_dirichlet_points(Kathode, 0)


############################### functions #####################################

gamma_list = Cauchy_BC.gamma_val_list()
a_list = Cauchy_BC.a_val_list()
Cauchy_facets = Cauchy_BC.facets()
Dirichlet_conditions = Dirichlet_BC.dirichlet_conditions()

solution = EquationAssembler(points,simplices,hull, k1_list, k2_list, rho_list,f_list,a_list,gamma_list,Dirichlet_conditions,Cauchy_facets)
sol = solution[1]

xx, yy = np.meshgrid(points[:,0], points[:,1])
fig = plt.figure()
plt.show()
plt.tricontourf(points[:,0], points[:,1], simplices, sol, 12, cmap="magma")
plt.tricontour(points[:,0], points[:,1], simplices, sol, 12, colors = "k", linewidths = 0.7)
plt.show()

xy = points
triangles = simplices
triang = mtri.Triangulation(xy[:,0], xy[:,1], triangles=triangles)
z = sol
    
fig, ax = plt.subplots(subplot_kw =dict(projection="3d"))
ax.plot_trisurf(triang, z, cmap = "magma")
ax.view_init(30, -50)
plt.show()

plot_streamline(solution[0],solution[1], solution[2], option = "cap")

