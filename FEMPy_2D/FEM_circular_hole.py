# -*- coding: utf-8 -*-

#-----------------------------------------------------------------
#
# name: thomas brunner
# email: brunner.th@hotmail.com
# nr: 12018550
#
#-----------------------------------------------------------------

from math import *
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import numpy as np
from meshpy.tet import MeshInfo, build
from mpl_toolkits.mplot3d import Axes3D
import meshpy.triangle as triangle
from meshing import *
from geometry import *
from equation_functions import *

p = createTriangulation(500,6,True)
p2 = TriangleMesh([(10, 10), (-10, 10), (-10, -10), (10, -10)], 0.5)





def EquationAssembler(points, simplices, hull):
    
    Master_Matrix = np.zeros((len(points),len(points)))
    Master_b = np.zeros((len(points)))
    
    for element in simplices: #hoping indizes are ordered
        Equation = determine_M1_and_b1(points[element[0]],points[element[1]],points[element[2]], 2, 2, 3)
        Master_Matrix[element[0],element[0]] += Equation[0][0,0]
        Master_Matrix[element[1],element[1]] += Equation[0][1,1]
        Master_Matrix[element[2],element[2]] += Equation[0][2,2]
        Master_Matrix[element[0],element[1]] += Equation[0][0,1]
        Master_Matrix[element[0],element[2]] += Equation[0][0,2]
        Master_Matrix[element[1],element[0]] += Equation[0][1,0]
        Master_Matrix[element[1],element[2]] += Equation[0][1,2]
        Master_Matrix[element[2],element[0]] += Equation[0][2,0]
        Master_Matrix[element[2],element[1]] += Equation[0][2,1]
        
        Master_b[element[0]] += Equation[1][0]
        Master_b[element[1]] += Equation[1][1]
        Master_b[element[2]] += Equation[1][2]
    
   
    
    def CauchyBC(Hull_Index):
        Equation_2 = determine_CauchyBC(points[hull[Hull_Index][0]],points[hull[Hull_Index][1]],-0.5)
        return Equation_2
    
    Boundary_Outside = BoundaryPointsOutsideRectangle(-9,-9,9,9,points, index_list=hull, OnlyBoundaryPoints= True)
    
    
    Hull_Index_List = []#Boundary_Outside

    for i in Hull_Index_List:
        Equation2 = CauchyBC(i)
        
        Master_Matrix[hull[i][0],hull[i][0]] += Equation2[0][0,0]
        Master_Matrix[hull[i][0],hull[i][1]] += Equation2[0][0,1]
        Master_Matrix[hull[i][1],hull[i][0]] += Equation2[0][1,0]
        Master_Matrix[hull[i][1],hull[i][1]] += Equation2[0][1,1]        
        #b_vec stays 0
        
    
    
    def DOF_Killer(Index, Value): #quelle: sormann
        
        for np in range(len(points)):
            Master_b[np]=Master_b[np] - Master_Matrix[np][Index]*Value;
            Master_Matrix[Index,np]=0.0;
            Master_Matrix[np,Index]=0.0;
        Master_Matrix[Index,Index]=1.0;
        Master_b[Index]=Value;
    
    
    
    def DirichletBC(Index_List, Value):
        for BC in Index_List:
            DOF_Killer(BC, Value)
    
    
    InnerCircleIndizes = BoundaryPointsInCircle(4, 0, 0,points,hull)
    InRectangle =  BoundaryPointsInRectangle(-9.9,9.9,-9.9,9.9, points, OnlyBoundaryPoints=False)
    InsideRect = BoundaryPointsInRectangle(-7.5,7.5,-10,-8,points, OnlyBoundaryPoints= False)
    Boundary_Outside = BoundaryPointsOutsideRectangle(-9.9,-9.9,9.9,9.9,points, OnlyBoundaryPoints= True)
    
    
    
   
    
    
    DirichletBC(Boundary_Outside, 25)
    DirichletBC(InnerCircleIndizes, 100)
    #DirichletBC(InsideRect, 100)
    #DirichletBC([504], 25)
    
    
    sol = np.linalg.solve(Master_Matrix,Master_b)
    
    xx, yy = np.meshgrid(points[:,0], points[:,1])
    
    
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:,0], points[:,1], sol)
    plt.show()
    plt.tricontourf(points[:,0], points[:,1], simplices, sol, 12, cmap="magma")
    plt.tricontour(points[:,0], points[:,1], simplices, sol, 12, colors = "k", linewidths = 0.7)
    plt.show()
    ax2 = plt.figure().add_subplot(projection='3d')
    ax2.plot_trisurf(points[:,0], points[:,1], sol, linewidth=0.2,cmap = "magma", antialiased=False)
    ax2.view_init(40, 140)
    plt.show()
    

#EquationAssembler(p[0],p[1],p[2])
EquationAssembler(p2[0],p2[1],p2[2])
