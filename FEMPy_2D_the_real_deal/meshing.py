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
import csv

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






def TriangleMesh(point_list, max_vol = 1e-3, radius = 3):
    
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
        (radius * np.cos(angle), radius * np.sin(angle))
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

#def SnakeDomain():
#    edge_array_rise = np.linspace(0,0.2,30)
#    x_list = np.concatenate(np.zeros((30)),np.array(0.2,0.2,0.8,0.8,1),np.zeros((30)).fill(1), np.array(0.6,0.6,0.4,0.4))
#    y_list = np.concatenate(


def UnitSquareMesh(max_vol = 0.01, edge_num = 2):
    
    #edge_array_rise = np.linspace(0,1,edge_num)
    #edge_array_const = np.zeros_like(edge_array_rise)
    
    #x_list = np.concatenate((edge_array_const,edge_array_rise, (edge_array_const+1),np.flip(edge_array_rise)))
    #y_list = np.concatenate((edge_array_rise, (edge_array_const+1),np.flip(edge_array_rise), edge_array_const))
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
    
    #qqqqqqqqqqqqqqqqqqqqqqqqqqqqqq
    
    #circ_start = len(points)
    
    #data = np.genfromtxt('BoingPatentAirfoil.dat',
    #                 skip_header=3,
    #                 skip_footer=1,
    #                 dtype=None,
    #                 delimiter=',')
    #print(data)
    #x_aero = 0.5*data[:,0]+0.3
    
    #y_aero = 0.8*data[:,1]+0.5
    #points_extension = np.stack((x_aero,y_aero),axis = 1)
    #print(points_extension[64:,:])
    #points_extension[64:,:] = np.flip(points_extension[64:,:], axis=0)
    #points = np.concatenate((points, points_extension), axis = 0)
    #facets.extend(round_trip_connect(circ_start, len(points) - 1))
    #qqqqqqqqqqqqqqqqqqqqqqqqqqqqqq
    
    ###############
    #radius = 0.1
    #circ_start = len(points)
    
    
    #print(len(facets))
    
    #x_ext = 0.5+float(radius) * np.cos(np.linspace(0, 2 * np.pi, 30, endpoint=False))
    #y_ext = 0.5+float(radius) * np.sin(np.linspace(0, 2 * np.pi, 30, endpoint=False))
     
    #points_extension = np.stack((x_ext, y_ext), axis = 1)
    #print(np.shape(points_extension))
    #points = np.concatenate((points, points_extension), axis = 0)
    
    #facets.extend(round_trip_connect(circ_start, len(points) - 1))
    
    #print(len(facets))
    ######################
    
    info = triangle.MeshInfo()
    
    
    info.set_points(points)
    info.set_facets(facets)
    
    ########
    
    #info.set_holes([(0.5, 0.5)])
    
    #info.set_holes([(0.5, 0.52),(0.5,0.49)])
    
    ########
    
    mesh = triangle.build(info, max_volume=max_vol, min_angle=30)

    mesh_points = np.array(mesh.points)
    mesh_tris = np.array(mesh.elements)
    
    plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, linewidth = 1)
    plt.show()
    
    return mesh_points, mesh_tris, facets


def MeshToCSV(path, points, elements, facets):
    
    f = open(path, 'w')
    writer = csv.writer(f)
    print(len(facets))
    print(len(points))
    print(len(elements))
    for ind, element in enumerate(elements):
        
        
        if ind < len(facets):
            point = points[ind]
            row = [point[0], point[1], elements[ind][0],elements[ind][1],elements[ind][2], facets[ind][0],facets[ind][1]]
        elif ind < len(points):
            point = points[ind]
            row = [point[0], point[1], elements[ind][0],elements[ind][1],elements[ind][2], "x","x"]
        else:
            row = ["x", "x", elements[ind][0],elements[ind][1],elements[ind][2], "x","x"]
            
        writer.writerow(row)
    f.close()


def CSVToMesh(path):
    
    with open(path, 'r') as infile:
       
        #lines einlesen
        lines = infile.readlines()
        infile.close()
        
        #Werte Listen erstellen
        points = []
        elements = []
        facets = []
        
        
        for line in lines:
            
            words = line.split(',')
            
            words[-1] = str(words[-1]).rstrip()
            try:
                if words[0] == "x" and words[-1]=="x":
                    elements.append([int(str(words[2])),int((str(words[3]))),int((str(words[4])))])
                elif words[0] == "x":
                    elements.append([int(str(words[2])),int((str(words[3]))),int((str(words[4])))])
                    facets.append([float(str(words[5])),float((str(words[6])))])
                elif words[-1] == "x":
                    elements.append([int(str(words[2])),int((str(words[3]))),int((str(words[4])))])
                    points.append([float(str(words[0])),float((str(words[1])))])
                else:  
                    points.append([float(str(words[0])),float((str(words[1])))])
                    elements.append([int(str(words[2])),int((str(words[3]))),int((str(words[4])))])
                    facets.append([int(str(words[5])),int((str(words[6])))])
            except: pass #print("empty line")
                
                
    points = np.array(points)
    elements = np.array(elements)
    facets = np.array(facets)
    
    return points, elements, facets



#TriangleMesh([[0,10],[-10,10],[-10,-10],[10,-10], [10,0]],max_vol= 5)
#Test = TriangleMesh([[10,10],[-10,10],[-10,-10],[10,-10]],max_vol= 0.5)

Test = UnitSquareMesh(max_vol = 0.1,edge_num=50)

path = r"C:\Users\brunn\Documents\GitHub\FEMPy\Mesh_files/SnakeDomain.csv"

MeshToCSV(path, Test[0], Test[1], Test[2])
#Mesh = CSVToMesh(path)

#print(len(Mesh[2]))
