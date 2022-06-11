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
import matplotlib.tri as mtri
import scipy.interpolate as ip 
import boundary_condition_classes as BC



def EquationAssembler(points, simplices, hull, k1_list, k2_list, rho_list, 
                      f_list, a_list, gamma_list, Dirichlet_points = [], 
                      Cauchy_points=[]):
    
    Master_Matrix = np.zeros((len(points),len(points)))     # init
    Master_b = np.zeros((len(points)))
    
    for ind, element in enumerate(simplices):   #iterate over elements
        
        # calculate mean values of coeffs:
            
        k1_mean = (k1_list[element[0]]+k1_list[element[1]]+
                   k1_list[element[2]])/3
        k2_mean = (k2_list[element[0]]+k2_list[element[1]]+
                   k2_list[element[2]])/3
        rho_mean = (rho_list[element[0]]+rho_list[element[1]]+
                    rho_list[element[2]])/3
        f_mean = (f_list[element[0]]+f_list[element[1]]+
                  f_list[element[2]])/3
        
        
        # calculate M matrix and B vector using local functions:
            
        M_mat = determine_M_mat(points[element[0]],points[element[1]],
                                points[element[2]], k1_mean, k2_mean, rho_mean)
        B_vec = determine_B_vec(points[element[0]],points[element[1]],
                                points[element[2]], f_mean)
        
        
        # assemble master matrix/vec using the local matrix/vec:
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
    
    
    for facet in Cauchy_points:     # iterate over all cauchy facets
    
        # calculate mean value of coeffs:
        a_mean = (a_list[facet[0]]+a_list[facet[1]])/2
        gamma_mean = (gamma_list[facet[0]]+gamma_list[facet[1]])/2
        
        # calculate M matrix and B vector entries using local functions:
        Equation2 = determine_CauchyBC(points[facet[0]],points[facet[1]],
                                       a_mean, gamma_mean)
        
        # assemble master matrix/vec using the local matrix/vec:
        Master_Matrix[facet[0],facet[0]] += Equation2[0][0,0]
        Master_Matrix[facet[0],facet[1]] += Equation2[0][0,1]
        Master_Matrix[facet[1],facet[0]] += Equation2[0][1,0]
        Master_Matrix[facet[1],facet[1]] += Equation2[0][1,1]        
        Master_b[facet[0]] += Equation2[1][0]
        Master_b[facet[1]] += Equation2[1][1]
    
    
    def DOF_Killer(Index, Value): # reduce degrees of freedom for dirichlet BC
        
        for num in range(len(points)): # iterate over point index
            Master_b[num]=Master_b[num] - Master_Matrix[num][Index]*Value
            Master_Matrix[Index,num]=0.0
            Master_Matrix[num,Index]=0.0
        Master_Matrix[Index,Index]=1.0
        Master_b[Index]=Value
    
    
    def DirichletBC(Index_List, Value): # use DOF_Killer on all dirichlet points
        for BC in Index_List:
            DOF_Killer(BC, Value)
    
    
    for ind, BC in enumerate(Dirichlet_points): # iterate over all dirichlet BC
        DirichletBC(Dirichlet_points[ind][0], Dirichlet_points[ind][1])
    
    # solve equation using NumPy LU decomp:
    sol = np.linalg.solve(Master_Matrix,Master_b) 
    
    return points, sol, simplices


def plot_streamline(points, values, simplices, grid_dim=1000): 
    
    #create streamline/field line plot by calculating the gradient numerically
    
    x = np.linspace(np.min(points[:,0]), np.max(points[:,0]), grid_dim)
    y = np.linspace(np.min(points[:,1]), np.max(points[:,1]), grid_dim)

    grid_x, grid_y = np.meshgrid(x,y)
    grid_ip = ip.griddata(points, values, (grid_x, grid_y), method='linear')
    
    u_x = (grid_ip[0:-2,:]-grid_ip[1:-1,:])/grid_dim
    u_x = np.concatenate((u_x, u_x[-3:-1,:]),axis = 0) #np.zeros((2,grid_dim))
    
    u_y = (grid_ip[:,0:-2]-grid_ip[:,1:-1])/grid_dim
    u_y = np.concatenate((u_y, u_y[:,-3:-1]),axis = 1) #np.zeros((grid_dim,2))
    
    plt.tricontourf(points[:,0], points[:,1], simplices, values, 12, cmap="viridis")
    plt.streamplot(y, x, u_y, u_x, color = "k", linewidth = 1, density = 1)
    
    circle1 = plt.Circle((0.5, 0.5), 0.1, color='w',zorder=2)
    cap1 = plt.Rectangle((0.3, 0.3), 0.05, 0.4, color='firebrick',zorder=2)
    cap2 = plt.Rectangle((0.65, 0.3), 0.05,0.4, color='cornflowerblue',zorder=2)
    
    #plt.gca().add_patch(circle1)
    #plt.gca().add_patch(cap1)
    #plt.gca().add_patch(cap2)
