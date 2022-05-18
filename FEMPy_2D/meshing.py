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

plt.rcParams["figure.dpi"] = 300
#plt.rcParams["figure.figsize"] = (8.3*1.5, 11.7*1.5)

def createTriangulation(N_area, N_boundary,  CornerPoints=True):

    N_boundary_per_side = int(N_boundary/4)
    points_bottom = np.concatenate((np.random.rand(N_boundary_per_side, 1), np.zeros((N_boundary_per_side, 1))), axis = 1)
    points_left = np.concatenate((np.zeros((N_boundary_per_side, 1)),np.random.rand(N_boundary_per_side, 1)), axis = 1)
    points_top = np.concatenate((np.random.rand(N_boundary_per_side, 1), np.full((N_boundary_per_side, 1),1)), axis = 1)
    points_right = np.concatenate((np.full((N_boundary_per_side, 1),1),np.random.rand(N_boundary_per_side, 1)), axis = 1)
    random_points_area = np.random.rand(N_area,2)
    points = np.concatenate((random_points_area,points_left,points_top,points_right,points_bottom)) # + boundary points
    
    if CornerPoints==True:
        cornerPoints = np.array([[0,0],[1,1],[0,1],[1,0]])
        points = np.concatenate((points, cornerPoints), axis = 0)
            
    triangulation = Delaunay(points)
    hull = triangulation.convex_hull
    simplices = triangulation.simplices
    points = triangulation.points
    
    plt.triplot(points[:,0], points[:,1], triangulation.simplices)
    #print(triangulation.simplices)
    plt.show()
    
    return points, simplices, hull



##############################################################################





def round_trip_connect(start, end):
    return [(i, i + 1) for i in range(start, end)] + [(end, start)]


def TriangleMesh(point_list, max_vol = 1e-3):

    def round_trip_connect(start, end):
        result = []
        for i in range(start, end):
            result.append((i, i + 1))
        result.append((end, start))
        return result

    points = point_list
    facets = round_trip_connect(0, len(points) - 1)
    circ_start = len(points)
    points.extend(
        (3 * np.cos(angle), 3 * np.sin(angle))
        for angle in np.linspace(0, 2 * np.pi, 30, endpoint=False)
        )

    facets.extend(round_trip_connect(circ_start, len(points) - 1))
    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(facets)
    info.set_holes([(0, 0)])
    mesh = triangle.build(info, max_volume=max_vol, min_angle=25)
    mesh_points = np.array(mesh.points)
    mesh_tris = np.array(mesh.elements)

    #print(mesh_tris)
    plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, color ="k", linewidth = 1)
    plt.show()
    
    return mesh_points, mesh_tris, facets
    


def HeatSinkMesh(max_vol = 1):
    
    x_list = np.array([-5,-5,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,6,6])
    y_list = np.array([0,1,1,10,10,2,2,10,10,2,2,10,10,2,2,10,10,2,2,10,10,1,1,0])
    point_list = np.stack((x_list,y_list), axis = 1)

    def round_trip_connect(start, end):
        result = []
        for i in range(start, end):
            result.append((i, i + 1))
        result.append((end, start))
        return result
    
    points = point_list
    facets = round_trip_connect(0, len(points) - 1)
    
    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(facets)
    
    mesh = triangle.build(info, max_volume=max_vol, min_angle=25)
    mesh_points = np.array(mesh.points)
    mesh_tris = np.array(mesh.elements)
    
    #print(mesh_tris)
    #plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris)
    plt.show()
    
    return mesh_points, mesh_tris, facets


#######################################################################


def UnitSquareMesh(max_vol = 1):
    x_list = np.array([0,0,1,1])
    y_list = np.array([0,1,1,0])
    
    point_list = np.stack((x_list,y_list), axis = 1)
    
    
    def round_trip_connect(start, end):
        result = []
        for i in range(start, end):
            result.append((i, i + 1))
        result.append((end, start))
        return result

    points = point_list
    facets = round_trip_connect(0, len(points) - 1)
    
    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(facets)
    
    mesh = triangle.build(info, max_volume=max_vol, min_angle=30)

    mesh_points = np.array(mesh.points)
    mesh_tris = np.array(mesh.elements)
    
    plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, linewidth = 1)
    plt.show()
    
    return mesh_points, mesh_tris, facets


TriangleMesh([[0,10],[-10,10],[-10,-10],[10,-10], [10,0]],max_vol= 5)
UnitSquareMesh(max_vol = 0.01)
