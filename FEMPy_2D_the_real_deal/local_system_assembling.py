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


def determine_M_mat(p_i,p_j,p_k, k1, k2, rho):
    
    d_pq = (p_j[0]-p_i[0])*(p_k[1]-p_i[1])-(p_k[0]-p_i[0])*(p_j[1]-p_i[1])
    alpha = (k1*(p_k[1]-p_i[1])**2+k2*(p_k[0]-p_i[0])**2)/d_pq
    beta = -(k1*(p_k[1]-p_i[1])*(p_j[1]-p_i[1])+k2*(p_k[0]-p_i[0])*(p_j[0]-p_i[0]))
    beta = beta/d_pq
    gamma = (k1*(p_j[1]-p_i[1])**2+k2*(p_j[0]-p_i[0])**2)/d_pq
    M1 = np.zeros((3,3))
    M1[0,0] = alpha+2*beta+gamma
    M1[0,1] = -alpha-beta
    M1[0,2] = -beta-gamma
    M1[1,0] = -alpha-beta
    M1[1,1] = alpha
    M1[1,2] = beta
    M1[2,0] = -beta-gamma
    M1[2,1] = beta
    M1[2,2] = gamma
    M2 = np.zeros((3,3))
    M2.fill(1)
    M2[0,0] += 1
    M2[1,1] += 1
    M2[2,2] += 1
    M2 = np.multiply(M2*(-1),rho*d_pq/48)
    M_mat = np.add(M1,M2)
    
    return M_mat


def determine_B_vec(p_i,p_j,p_k, f):
    
    d_pq = (p_j[0]-p_i[0])*(p_k[1]-p_i[1])-(p_k[0]-p_i[0])*(p_j[1]-p_i[1])
    B_vec = np.zeros(3)
    B_vec.fill(f*d_pq/6)
    
    return B_vec


def determine_CauchyBC(point_p,point_q,a, gamma): 
    
    d_pq = BoundaryPointsDistance(point_p, point_q)
    M_BC = np.zeros((2,2))
    M_BC[0,0]=2
    M_BC[0,1]=1
    M_BC[1,0]=1
    M_BC[1,1]=2
    M_BC = np.multiply(M_BC,a*d_pq/12)
    gamma_term = np.array([1,1])*gamma*d_pq*1/2 
    b_BC = gamma_term
    
    return M_BC, b_BC
















