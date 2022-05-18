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


def BoundaryPointsDistance(point_p, point_q):
    return sqrt((point_p[0]-point_q[0])**2+(point_p[1]-point_q[1])**2)



def determine_M1_and_b1(p_i,p_j,p_k, a1, a2, h, a0=0):
    
    d_m = (p_j[0]-p_i[0])*(p_k[1]-p_i[1])-(p_k[0]-p_i[0])*(p_j[1]-p_i[1])
    e11 = (a1*(p_k[1]-p_i[1])**2+a2*(p_k[0]-p_i[0])**2)/d_m
    e12 = -(a1*(p_k[1]-p_i[1])*(p_j[1]-p_i[1])+a2*(p_k[0]-p_i[0])*(p_j[0]-p_i[0]))/d_m
    e22 = (a1*(p_j[1]-p_i[1])**2+a2*(p_j[0]-p_i[0])**2)/d_m
    M1 = np.zeros((3,3))
    B1 = np.zeros((3,3))
    M1[0,0]= e11+2*e12+e22
    M1[0,1]=-e11-e12
    M1[0,2]=-e12-e22
    M1[1,0]=-e11-e12
    M1[1,1]=e11
    M1[1,2]=e12
    M1[2,0]=-e12-e22
    M1[2,1]=e12
    M1[2,2]=e22
    B1[0,0]= a0*d_m/24*2
    B1[0,1]= a0*d_m/24
    B1[0,2]= a0*d_m/24
    B1[1,0]= a0*d_m/24
    B1[1,1]= a0*d_m/24*2
    B1[1,2]= a0*d_m/24
    B1[2,0]= a0*d_m/24
    B1[2,1]= a0*d_m/24
    B1[2,2]= a0*d_m/24*2
    B1 = np.multiply(M1, 0.5)
    b1 = np.multiply(np.full((3),1),h*d_m/6)
    
    return M1, b1, B1



def determine_CauchyBC(point_p,point_q,a4, a5 = 0): 
    
    d_pq = BoundaryPointsDistance(point_p, point_q)
    M_BC = np.zeros((2,2))
    M_BC[0,0]=2
    M_BC[0,1]=1
    M_BC[1,0]=1
    M_BC[1,1]=2
    M_BC = np.multiply(M_BC,a4*d_pq/6)
    a5_term = np.array([1,1])*a5*d_pq*1/2 
    b_BC = a5_term
    
    return M_BC, b_BC
















