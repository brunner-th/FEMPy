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
import matplotlib.tri as mtri

p_heatsink = HeatSinkMesh(0.1)




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
    
    InRectangle =  BoundaryPointsInRectangle(-9.9,9.9,-1,0.1, points, OnlyBoundaryPoints=False)
    OnEdge = BoundaryPointsInRectangle(-10,10,9.7,10, points,  OnlyBoundaryPoints=False)
    Boundary_Outside = BoundaryPointsInRectangle(-12,12,0.8,11,points, index_list=hull, OnlyBoundaryPoints= True)
    
    
    def CauchyBC(Hull_Index):
        Equation_2 = determine_CauchyBC(points[hull[Hull_Index][0]],points[hull[Hull_Index][1]],0.5, a5 = 0) #sollte nicht -0.5 sein??????
        return Equation_2
    
    Hull_Index_List = Boundary_Outside 

    for i in Hull_Index_List:
        Equation2 = CauchyBC(i)
        
        Master_Matrix[hull[i][0],hull[i][0]] += Equation2[0][0,0]
        Master_Matrix[hull[i][0],hull[i][1]] += Equation2[0][0,1]
        Master_Matrix[hull[i][1],hull[i][0]] += Equation2[0][1,0]
        Master_Matrix[hull[i][1],hull[i][1]] += Equation2[0][1,1]        
        Master_b[hull[i][0]] += Equation2[1][0]
        Master_b[hull[i][1]] += Equation2[1][1]
        
    
    
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
    
    
    #DirichletBC(InnerCircleIndizes, 100)
    #DirichletBC(Boundary_Outside, 25)
    DirichletBC(InRectangle, 100)
    #DirichletBC(OnEdge, 40)
    #DirichletBC([504], 25)
    
    
    sol = np.linalg.solve(Master_Matrix,Master_b)
    
    xx, yy = np.meshgrid(points[:,0], points[:,1])
    
   
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:,0], points[:,1], sol)
    ax.view_init(40,120)
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
    ax.view_init(20, 130)
    plt.show()
    
    return Master_Matrix

M = EquationAssembler(p_heatsink[0],p_heatsink[1],p_heatsink[2])
