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


def BoundaryPointsInCircle(Radius, x, y, point_list, index_list, OnlyBoundaryPoints = True):
    
    IndexOfPointsInCircle = []
    for line in index_list:
        index = line[0]
        point = point_list[index]
        if (point[0]-x)**2+(point[1]-y)**2 <= Radius**2:
            IndexOfPointsInCircle.append(index)
            
    return IndexOfPointsInCircle


def BoundaryPointsInRectangle( x1,x2,y1,y2, point_list,index_list = [], OnlyBoundaryPoints = True):
    
    IndexOfPointsInRectangle = []
    if OnlyBoundaryPoints == True:
        for line in index_list:
            index = line[0]
            point = point_list[index]
            if x1 < point[0] and y1 < point[1] and   point[0]<x2 and point[1]<y2 :
                IndexOfPointsInRectangle.append(index)
    else:
        i = 0
        for point in point_list:
            if x1 < point[0] and y1 < point[1] and point[0]<x2 and point[1]<y2 :
                IndexOfPointsInRectangle.append(i)
            i = i + 1
                
    return IndexOfPointsInRectangle


def BoundaryPointsOutsideRectangle(x1,x2,y1,y2,point_list,index_list = [], OnlyBoundaryPoints = True):
    
    IndexOfPointsOutsideRectangle = []
    if OnlyBoundaryPoints == True:
        for line in index_list:
            index = line[0]
            point = point_list[index]
            if x2 < point[0] or y2 < point[1] or   point[0]<x1 or point[1]<y1 :
                IndexOfPointsOutsideRectangle.append(index)
    else:
        i = 0
        for point in point_list:
            if x2 < point[0] or y2 < point[1] or   point[0]<x1 or point[1]<y1 :
                IndexOfPointsOutsideRectangle.append(i)
            i = i + 1
                  
    return IndexOfPointsOutsideRectangle













