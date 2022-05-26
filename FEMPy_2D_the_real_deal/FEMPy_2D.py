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
from equation_functions import *
import matplotlib.tri as mtri

############################## mesh generation / import #######################


#p_unitsquare = UnitSquareMesh(0.01)

mesh_path = r"C:\Users\brunn\Documents\GitHub\FEMPy\Mesh_files/DomainWithHole.csv"

p_unitsquare = CSVToMesh(mesh_path)

print("mesh generated/imported successfully")

points = p_unitsquare[0]
simplices = p_unitsquare[1]
hull = p_unitsquare[2]



############################## geometry handling ##############################



# unit square:
    
LeftSide =  BoundaryPointsInRectangle(-1,0.01,-2,2, points, OnlyBoundaryPoints=False)
RightSide =  BoundaryPointsInRectangle(0.99,1.1,-2,2, points, OnlyBoundaryPoints=False)
RightSideBoundary =  BoundaryPointsInRectangle(0.9,1.1,-2,2, points, index_list=hull, OnlyBoundaryPoints=True)
FullBoundary = BoundaryPointsOutsideRectangle(0.01, 0.99, 0.01, 0.99, points)


# domain with hole:
    
Inlet = BoundaryPointsInRectangle(-1,0.001,-2,2, points, index_list=hull, OnlyBoundaryPoints=True)
Outlet = BoundaryPointsInRectangle(0.99,1.1,-2,2, points, index_list=hull, OnlyBoundaryPoints=True)

############################## parameters #####################################
#print(hull)

Inlet = Inlet[:-1]

print(Inlet)
print(Outlet)

k1_list = np.zeros((len(points)))
k2_list = np.zeros((len(points)))
rho_list = np.zeros((len(points)))
f_list = np.zeros((len(points)))
a_list = np.zeros((len(points)))
gamma_list = np.zeros((len(points)))

k1_list.fill(1)
k2_list.fill(1)
#f_list.fill(10)
#f_list[78] = 10000
#f_list = np.linspace(0,100,len(points))
gamma_list.fill(1)
for point in Outlet:
    gamma_list[point]=-1
    
 
rho_list.fill(0)
#a_list.fill(0.5)


################################ boundary conditions ##########################

Dirichlet_conditions = []
#Dirichlet_conditions = [[LeftSide, 100]]

#Dirichlet_conditions = [[FullBoundary, 0]]
#Cauchy_points = RightSideBoundary
#Cauchy_points = []
Cauchy_points = np.concatenate((Inlet,Outlet))

################################ output to console ############################


#print("points:")
#print(points)
#print("hull:")
#print(hull)

#print("Dirichlet Boundary Points Index:")
#print(Dirichlet_conditions)
#print("Cauchy Boundary Points Index:")
#print(Cauchy_points)


############################### functions #####################################


def EquationAssembler(points, simplices, hull, k1_list, k2_list, rho_list, f_list, a_list, gamma_list, Dirichlet_points = [], Cauchy_points=[]):
    
    Master_Matrix = np.zeros((len(points),len(points)))
    Master_b = np.zeros((len(points)))
    
    for ind, element in enumerate(simplices): 
        
        k1_mean = (k1_list[element[0]]+k1_list[element[1]]+k1_list[element[2]])/3
        k2_mean = (k2_list[element[0]]+k2_list[element[1]]+k2_list[element[2]])/3
        rho_mean = (rho_list[element[0]]+rho_list[element[1]]+rho_list[element[2]])/3
        f_mean = (f_list[element[0]]+f_list[element[1]]+f_list[element[2]])/3
        
        M_mat = determine_M_mat(points[element[0]],points[element[1]],points[element[2]], k1_mean, k2_mean, rho_mean)
        B_vec = determine_B_vec(points[element[0]],points[element[1]],points[element[2]], f_mean)
        
        
        Master_Matrix[element[0],element[0]] += M_mat[0,0]
        Master_Matrix[element[1],element[1]] += M_mat[1,1]
        Master_Matrix[element[2],element[2]] += M_mat[2,2]
        Master_Matrix[element[0],element[1]] += M_mat[0,1]
        Master_Matrix[element[0],element[2]] += M_mat[0,2]
        Master_Matrix[element[1],element[0]] += M_mat[1,0]
        Master_Matrix[element[1],element[2]] += M_mat[1,2]
        Master_Matrix[element[2],element[0]] += M_mat[2,0]
        Master_Matrix[element[2],element[1]] += M_mat[2,1]
        Master_b[element[0]] += B_vec[0]
        Master_b[element[1]] += B_vec[1]
        Master_b[element[2]] += B_vec[2]
    
    
    def CauchyBC(Hull_Index):
        Equation_2 = determine_CauchyBC(points[hull[Hull_Index][0]],points[hull[Hull_Index][1]],a_list[Hull_Index], gamma_list[Hull_Index])
        return Equation_2
    

    for i in Cauchy_points:
        Equation2 = CauchyBC(i)
        Master_Matrix[hull[i][0],hull[i][0]] += Equation2[0][0,0]
        Master_Matrix[hull[i][0],hull[i][1]] += Equation2[0][0,1]
        Master_Matrix[hull[i][1],hull[i][0]] += Equation2[0][1,0]
        Master_Matrix[hull[i][1],hull[i][1]] += Equation2[0][1,1]        
        Master_b[hull[i][0]] += Equation2[1][0]
        Master_b[hull[i][1]] += Equation2[1][1]
    
    
    def DOF_Killer(Index, Value): 
        
        for num in range(len(points)):
            Master_b[num]=Master_b[num] - Master_Matrix[num][Index]*Value
            Master_Matrix[Index,num]=0.0
            Master_Matrix[num,Index]=0.0
        Master_Matrix[Index,Index]=1.0
        Master_b[Index]=Value
    
    
    def DirichletBC(Index_List, Value):
        for BC in Index_List:
            DOF_Killer(BC, Value)
    
    for ind, BC in enumerate(Dirichlet_points):
        DirichletBC(Dirichlet_conditions[ind][0], Dirichlet_conditions[ind][1])
    
    
    sol = np.linalg.solve(Master_Matrix,Master_b)
    
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


EquationAssembler(points,simplices,hull, k1_list, k2_list, rho_list,f_list,a_list,gamma_list,Dirichlet_conditions,Cauchy_points)

