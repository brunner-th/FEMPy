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


#p_unitsquare = UnitSquareMesh(0.01)

mesh_path = r"C:\Users\brunn\Documents\GitHub\FEMPy\Mesh_files/DomainWithHoleFine.csv" #BoingAirfoil

p_unitsquare = CSVToMesh(mesh_path)

print("mesh generated/imported successfully")

points = p_unitsquare[0]
simplices = p_unitsquare[1]
hull = p_unitsquare[2]



############################## geometry handling ##############################



# unit square:
    
LeftSide =  BoundaryPointsInRectangle(-1,0.01,-2,2, points)
RightSide =  BoundaryPointsInRectangle(0.99,1.1,-2,2, points)
FullBoundaryCauchy = BoundaryFacetsOutsideRectangle(0.01, 0.99, 0.01, 0.99, points, hull)

#Heat kernel:
FullBoundaryDirichlet = BoundaryPointsOutsideRectangle(0.01, 0.99, 0.01, 0.99, points)

# flow:
SmallPatch = BoundaryPointsInRectangle(0.4,0.6,0.1,0.2, points)
Inlet = BoundaryFacetsInRectangle(-1,0.01,-2,2, points, hull)
Outlet = BoundaryFacetsInRectangle(0.99,1.1,-2,2, points, hull)

#capacitor:
    
Anode = BoundaryPointsInRectangle(0.3,0.35,0.3,0.7, points)
Kathode = BoundaryPointsInRectangle(0.65,0.7,0.3,0.7, points)
Dielectric = BoundaryPointsInRectangle(0.35,0.65,0.3,0.7, points)

############################## parameters #####################################



k1_list = np.zeros((len(points)))
k2_list = np.zeros((len(points)))
rho_list = np.zeros((len(points)))
f_list = np.zeros((len(points)))
a_list = np.zeros((len(points)))
gamma_list = np.zeros((len(points)))

k1_list.fill(1)
k2_list.fill(1)
#f_list.fill(1)

#f_list[300] = -10
#f_list = np.linspace(0,1,len(points))
gamma_list.fill(1)

#for point in Dielectric:
#    k1_list[point] = 1000
#    k2_list[point] = 1000
    

for facet in Outlet:
    p1 = facet[0]
    p2 = facet[1]
    gamma_list[p1]=-1
    gamma_list[p2]=-1
    
 
rho_list.fill(0)
#a_list.fill(0.5)

print(gamma_list)

################################ boundary conditions ##########################
#Dirichlet_conditions = [[FullBoundaryDirichlet, 0]]
#Dirichlet_conditions = [[SmallPatch, 10]]
#Dirichlet_conditions = [[Anode, 5],[Kathode, 0]]
#Dirichlet_conditions = [[RightSide, 0]]
#Dirichlet_conditions = [[LeftSide, 100]]
Dirichlet_conditions = []
#Dirichlet_conditions = [[FullBoundary, 0]]
#Cauchy_points = RightSideBoundary
Cauchy_points = []
Cauchy_points = np.concatenate((Inlet,Outlet))

################################ output to console ############################


#print("points:")
#print(points)
#print("hull:")
#print(hull)

#print("Dirichlet Boundary Points Index:")
#print(Dirichlet_conditions)
#print("Cauchy Boundary Points Index:")
print(Cauchy_points)


############################### functions #####################################

solution = EquationAssembler(points,simplices,hull, k1_list, k2_list, rho_list,f_list,a_list,gamma_list,Dirichlet_conditions,Cauchy_points)

sol = solution[1]


xx, yy = np.meshgrid(points[:,0], points[:,1])
fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
plt.show()
plt.tricontourf(points[:,0], points[:,1], simplices, sol, 12, cmap="magma")
plt.tricontour(points[:,0], points[:,1], simplices, sol, 12, colors = "k", linewidths = 0.7)
plt.show()
#ax2 = plt.figure().add_subplot(projection='3d')
#ax2.plot_trisurf(points[:,0], points[:,1], sol, linewidth=0.2,cmap = "magma", antialiased=False)
#ax2.view_init(40, 140)
#plt.show()

xy = points
triangles = simplices
triang = mtri.Triangulation(xy[:,0], xy[:,1], triangles=triangles)
z = sol
    
fig, ax = plt.subplots(subplot_kw =dict(projection="3d"))
ax.plot_trisurf(triang, z, cmap = "magma")
ax.view_init(30, -50)
plt.show()
    








plot_streamline(solution[0],solution[1], solution[2])



