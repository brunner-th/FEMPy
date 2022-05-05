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
import scipy.integrate as integrate

plt.rcParams['figure.figsize'] = [5.0, 3.0]
plt.rcParams['figure.dpi'] = 200

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

plt.rc('text', usetex=True)
plt.rc('font', family='serif')



def inhom_function(x):
    f = 0
    #if  0.51 > x > 0.49:
    #    f = -50
    return f

def Integrand_1(var,j):
    val = inhom_function(var)
    return (var-x[j-1])/(x[j]-x[j-1])*val

def Integrand_2(var,j):
    val = inhom_function(var)
    return (x[j+1]-var)/(x[j+1]-x[j])*val


A = 0 # first boundary x-pos
B = 1 # second boundary x-pos
L = B-A
vertex_number = 10


leftBC = "dirichlet" #neumann or dirichlet
rightBC = "dirichlet"

a1 = 1 # value of left bc
a2 = 0 # value of right bc

k = 0
b = 2

Domain = [A, B]

x = np.linspace(Domain[0],Domain[1],vertex_number)
x_an = np.linspace(Domain[0],Domain[1], 1000)

A_mat = np.zeros((len(x),len(x)))
B_mat = np.zeros((len(x),len(x)))
C_mat = np.zeros((len(x),len(x)))
b_vec = np.zeros(len(x))



for j in range(len(x)):
    
    if j == 0 and leftBC == "neumann":
        print("left NBC")
        A_mat[j,j+1] = -1/(x[j+1]-x[j])
        B_mat[j,j+1] = (x[j+1]-x[j])/6
        C_mat[j,j+1] = 1/2
        
        A_mat[j,j] = 1/(x[j+1]-x[j])
        B_mat[j,j] = (x[j+1]-x[j])/3
        C_mat[j,j] = -1/2
    
    elif j == len(x)-1 and rightBC == "neumann":
        print("right NBC")
        A_mat[j,j-1] = -1/(x[j]-x[j-1])
        B_mat[j,j-1] = (x[j]-x[j-1])/6
        C_mat[j,j-1] = -1/2
        
        A_mat[j,j] = 1/(x[j]-x[j-1])
        B_mat[j,j] = (x[j]-x[j-1])/3
        C_mat[j,j] = 1/2
        
    elif j != len(x)-1 and j != 0:
        B_mat[j,j] = (-x[j-1])/3+(x[j+1])/3
        A_mat[j,j] = 1/(x[j]-x[j-1])+1/(x[j+1]-x[j]) 
        C_mat[j,j] = 0
        
        B_mat[j,j+1] = (x[j+1]-x[j])/6 
        A_mat[j,j+1] = -1/(x[j+1]-x[j]) 
        C_mat[j,j+1] = +1/2
        
        B_mat[j,j-1] = (x[j]-x[j-1])/6 
        A_mat[j,j-1] = -1/(x[j]-x[j-1])
        C_mat[j,j-1] = -1/2


mat_sum = -A_mat+k*B_mat+b*C_mat

for j in range(1,len(x)-1):
    Integral_1 = integrate.quad(lambda var: Integrand_1(var,j), x[j-1], x[j])[0]
    Integral_2 = integrate.quad(lambda var: Integrand_2(var,j), x[j], x[j+1])[0]
    b_vec[j] = Integral_1+Integral_2
    

if leftBC == "dirichlet":
    mat_sum[0,0] = 1
    b_vec[0] = a1
else:
    b_vec[0] += a1
    
if rightBC == "dirichlet":
    mat_sum[len(x)-1,len(x)-1] = 1
    b_vec[len(x)-1] = a2
else:
    b_vec[len(x)-1] -= a2
    
    
sol = np.linalg.solve(mat_sum, b_vec)

plt.grid()
plt.xlabel("distance $x$",fontsize = 12)
plt.ylabel("potential $\phi$",fontsize = 12)
plt.plot(x,sol, color = "firebrick", linestyle = "-", label = "Finite Element solution")
plt.legend()

# analytical solution in DBC and NBC file



def drawBasisFunction(ind1,ind2,ind3, col="k", solution = True):
    
    first_x_interval = np.linspace(x[ind1], x[ind2], 20)
    second_x_interval = np.linspace(x[ind2], x[ind3], 20)
    
   
    
    if solution == True:
        first_y_interval = np.linspace(0,sol[ind2],20)
        second_y_interval = np.linspace(sol[ind2],0,20)
    else:
        first_y_interval = np.linspace(0,1,20)
        second_y_interval = np.linspace(1,0,20)
    
    x_int = np.concatenate((first_x_interval, second_x_interval))
    y_int = np.concatenate((first_y_interval, second_y_interval))
    
    plt.plot(x_int, y_int, linestyle = "--", color = col)


drawBasisFunction(1, 0, 1)
drawBasisFunction(len(sol)-3, len(sol)-2, len(sol)-1)

for i in range(1,len(sol)-2):
    drawBasisFunction(i-1, i, i+1, solution=1)
    #print(" ")



