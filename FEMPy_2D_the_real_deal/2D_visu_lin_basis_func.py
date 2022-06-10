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
from global_assembler import *
import matplotlib.tri as mtri
import scipy.interpolate as ip 
import boundary_condition_classes as BC
plt.rcParams["figure.dpi"] = 300
plt.rcParams["figure.figsize"] = (10, 8)

mesh_path = r"C:\Users\brunn\Documents\GitHub\FEMPy\Mesh_files\UnitSquareNormal.csv"

p_unitsquare = CSVToMesh(mesh_path)

points = p_unitsquare[0]
simplices = p_unitsquare[1]
hull = p_unitsquare[2]

sol = np.zeros((len(points))) 
sol[4] = 1
xy = points
triangles = simplices
triang = mtri.Triangulation(xy[:,0], xy[:,1], triangles=triangles)
z = sol
    
fig, ax = plt.subplots(subplot_kw =dict(projection="3d"))
ax.plot_trisurf(triang, z, cmap = "Greys", linewidth = 1, edgecolors='gray')
ax.view_init(25, -40)

x_scale=1
y_scale=1
z_scale=0.3

scale=np.diag([x_scale, y_scale, z_scale, 1.0])
scale=scale*(1.0/scale.max())
scale[3,3]=1.0

def short_proj():
  return np.dot(Axes3D.get_proj(ax), scale)

ax.get_proj=short_proj

plt.show()