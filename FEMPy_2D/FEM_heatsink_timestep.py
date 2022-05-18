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
import glob
from PIL import Image
import matplotlib.tri as mtri
import os



p_heatsink = HeatSinkMesh(0.5)
plt.rcParams["figure.dpi"] = 100



def StationaryEquationAssembler(points, simplices, hull, timesteps = 100, delta_t = 0.05):
    
    Master_Matrix = np.zeros((len(points),len(points)))
    Master_b = np.zeros((len(points)))
    B_mat = np.zeros_like(Master_Matrix)
    
    for element in simplices: #hoping indizes are ordered
        Equation = determine_M1_and_b1(points[element[0]],points[element[1]],points[element[2]], 2, 2, 3, a0=2)
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
    
        B_mat[element[0],element[0]] += Equation[2][0,0]
        B_mat[element[1],element[1]] += Equation[2][1,1]
        B_mat[element[2],element[2]] += Equation[2][2,2]
        B_mat[element[0],element[1]] += Equation[2][0,1]
        B_mat[element[0],element[2]] += Equation[2][0,2]
        B_mat[element[1],element[0]] += Equation[2][1,0]
        B_mat[element[1],element[2]] += Equation[2][1,2]
        B_mat[element[2],element[0]] += Equation[2][2,0]
        B_mat[element[2],element[1]] += Equation[2][2,1]
   
    
    InRectangle =  BoundaryPointsInRectangle(-9.9,9.9,-1,0.1, points, OnlyBoundaryPoints=False)
    OnEdge = BoundaryPointsInRectangle(-10,10,9.7,10, points,  OnlyBoundaryPoints=False)
    Boundary_Outside = BoundaryPointsInRectangle(-11,11,0.9,11,points, index_list=hull, OnlyBoundaryPoints= True)
    
    
    
    def CauchyBC(Hull_Index):
        Equation_2 = determine_CauchyBC(points[hull[Hull_Index][0]],points[hull[Hull_Index][1]],0.5, a5 = 0) #sollte nicht -0.5 sein??????
        return Equation_2
    
    Hull_Index_List = Boundary_Outside #!!!???

    for i in Hull_Index_List:
        Equation2 = CauchyBC(i)
        
        Master_Matrix[hull[i][0],hull[i][0]] += Equation2[0][0,0]
        Master_Matrix[hull[i][0],hull[i][1]] += Equation2[0][0,1]
        Master_Matrix[hull[i][1],hull[i][0]] += Equation2[0][1,0]
        Master_Matrix[hull[i][1],hull[i][1]] += Equation2[0][1,1]        
        #b_vec stays 0
        Master_b[hull[i][0]] += Equation2[1][0]
        Master_b[hull[i][1]] += Equation2[1][1]
        
    
    
    def DOF_Killer(Index, Value): #quelle: sormann
        
        for num in range(len(points)):
            Master_b[num]=Master_b[num] - Master_Matrix[num][Index]*Value
            Master_Matrix[Index,num]=0.0
            Master_Matrix[num,Index]=0.0
        Master_Matrix[Index,Index]=1.0
        Master_b[Index]=Value
    
    
    
    def DirichletBC(Index_List, Value):
        for BC in Index_List:
            DOF_Killer(BC, Value)
    
    
    #DirichletBC(InnerCircleIndizes, 100)
    #DirichletBC(Boundary_Outside, 25)
    #DirichletBC(InRectangle, 100)
    #DirichletBC(OnEdge, 40)
    #DirichletBC([504], 25)
    
    A_mat = Master_Matrix
    f = Master_b
    C_mat = np.add(A_mat, 2/delta_t*B_mat)
    D_mat = np.add(A_mat, -2/delta_t*B_mat)
    
    #dirichlet with timestep:
    
    for ind in InRectangle:
        value = 100
        for num in range(len(A_mat)):
            f[num]=f[num] - A_mat[num][ind]*value
            C_mat[ind][num]=0.0
            C_mat[num][ind]=0.0
            D_mat[ind][num]=0.0
            D_mat[num][ind]=0.0

        C_mat[ind][ind]=1.0
        D_mat[ind][ind]=-1.0
        f[ind]=0.0
        
    xy = points
    triangles = simplices
    triang = mtri.Triangulation(xy[:,0], xy[:,1], triangles=triangles)
    
    told = np.zeros_like(f)
    told.fill(20)
    animationcounter = 1000
    for t in range(timesteps):
        
        for ind in InRectangle:
            
            value = 100
            told[ind]=value
            
        RHS=(2*f-D_mat@told)
        sol = np.linalg.solve(C_mat,RHS)
        
        plt.tricontourf(points[:,0], points[:,1], simplices, told, 12, cmap="magma")
        plt.tricontour(points[:,0], points[:,1], simplices, told, 12, colors = "k", linewidths = 0.7)
        
        #fig, ax = plt.subplots(subplot_kw =dict(projection="3d"))
        #ax.view_init(20, 130)
        #ax.plot_trisurf(triang, told, cmap = "magma")
        
        
        plt.savefig(
            r"C:\Users\brunn\Desktop\FEM\Animation/Images/image_"
            + str(animationcounter))
        plt.show()
        told = sol
        animationcounter = animationcounter +1
    
    



StationaryEquationAssembler(p_heatsink[0],p_heatsink[1],p_heatsink[2], timesteps = 100, delta_t = 0.02)



# filepaths:
fp_in = (
    r"C:\Users\brunn\Desktop\FEM\Animation/Images/image_*.png"
)
fp_out = (
    r"C:\Users\brunn\Desktop\FEM\Animation/Output/animation.gif"
)
\

# this part is used to create the animation out of the pictures

img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
img.save(
    fp=fp_out, format="GIF", append_images=imgs, save_all=True, duration=100, loop=0
)
files_to_delete = glob.glob(fp_in)
for f in files_to_delete:
    os.remove(f)
