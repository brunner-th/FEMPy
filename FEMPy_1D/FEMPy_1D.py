# -*- coding: utf-8 -*-

#-----------------------------------------------------------------
#
# name: thomas brunner
# email: brunner.th@hotmail.com
# nr: 12018550
#
#-----------------------------------------------------------------

from math import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate as integrate
import functools as functools
import time

plt.rcParams['figure.figsize'] = [5.0, 3.0]
plt.rcParams['figure.dpi'] = 200

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

plt.rc('text', usetex=True)
plt.rc('font', family='serif')



#############################   parameters   ##################################


A = 0                           # first boundary x-pos
B = 250                        # second boundary x-pos

vertex_number = 20    # number of vertices


plot_hat_function = False       # plot with individual basis functions
constant_coeff = False      # use const coeff or coeff functions
plot_analytical_sol = False  # plot analytical solution
plot_legend = True


leftBC = "dirichlet"            # determine type of BC:
rightBC = "dirichlet"             # neumann or dirichlet

a1 = 5                          # value of left bc
a2 = 0                          # value of right bc

a = 0                        # parameter a if constant coeff is true
b = 0                   # parameter b if constant coeff is true


#############################   functions   ###################################


def inhom_function(x):          # inhomogeneous function (RHS)
    #f = -10**-9/(8.85*10**-6)
    
    f = 0
    if  100 > x > 0:
        f = +6*10**-9/(8.85*10**-6)
    elif 250 > x > 150:
        f = -3*10**-9/(8.85*10**-6)
    return f

def a_func(x):                  # function a(x) if constant coeff is false
    f = 0
    return f

def b_func(x):                  # function b(x) if constant coeff is false
    f = 0
    return f

def exact_solution_function(x):
    
    f = -0.0000564972*x**2 - 0.00587571*x + 5
    
    return f

################################   code   #####################################

L = B-A
Domain = [A, B]

x = np.linspace(Domain[0],Domain[1],vertex_number)
x_an = np.linspace(Domain[0],Domain[1], 1000)
n = len(x)

A_mat = np.zeros((len(x),len(x)))
B_mat = np.zeros((len(x),len(x)))
C_mat = np.zeros((len(x),len(x)))
b_vec = np.zeros(len(x))

def Integrand_1(var,j):
    val = inhom_function(var)
    return (var-x[j-1])/(x[j]-x[j-1])*val

def Integrand_2(var,j):
    val = inhom_function(var)
    return (x[j+1]-var)/(x[j+1]-x[j])*val


def phi_i(pos,i):
    f=0
    if i == 0:
        if x[0]<=pos<=x[1]:
            f = (x[1]-pos)/(x[1]-x[0])
    elif i == n-1:
        if x[n-2]<=pos<=x[n-1]:
            f = (pos-x[n-2])/(x[n-1]-x[n-2])
    else:
        if x[i-1]<=pos<=x[i]:
            f = (pos-x[i-1])/(x[i]-x[i-1])
        elif x[i]<=pos<=x[i+1]:
            f = (x[i+1]-pos)/(x[i+1]-x[i])
    return f

def d_phi_i(pos,i):
    f=0
    if i == 0:
        if x[0]<=pos<=x[1]:
            f = (-1)/(x[1]-x[0])
    elif i == n-1:
        if x[n-2]<=pos<=x[n-1]:
            f = (1)/(x[n-1]-x[n-2])
    else:
        if x[i-1]<=pos<=x[i]:
            f = (1)/(x[i]-x[i-1])
        elif x[i]<=pos<=x[i+1]:
            f = (-1)/(x[i+1]-x[i])
    return f

def plot_exact_sol(ex_func=exact_solution_function):
    
    y_an = np.zeros_like(x_an)
    for i, x_val in enumerate(x_an):
        y_an[i] = exact_solution_function(x_val)
    
    plt.plot(x_an,y_an, color = "skyblue", linestyle ="-.",label = "Analytical solution")
    
    return


def AssembleSystem(const_coeff = False):
    
    if const_coeff == False:
        for j in range(n):
            
            if j == 0 and leftBC == "neumann":
                
                A_mat[j,j+1] = -1/(x[j+1]-x[j])
                A_mat[j,j] = 1/(x[j+1]-x[j])
                
                B_mat[j,j] = integrate.quad(lambda var: b_func(var)*phi_i(var,j)**2, x[j], x[j+1])[0]
                B_mat[j,j+1] = integrate.quad(lambda var: b_func(var)*phi_i(var,j)*phi_i(var,j+1), x[j], x[j+1])[0]
                
                C_mat[j,j] = integrate.quad(lambda var: a_func(var)*d_phi_i(var,j)*phi_i(var,j), x[j], x[j+1])[0]
                C_mat[j,j+1] = -integrate.quad(lambda var: a_func(var)*d_phi_i(var,j)*phi_i(var,j+1), x[j], x[j+1])[0]
                # - sign is needed, not quite sure why yet
            
            elif j == n-1 and rightBC == "neumann":
                
                A_mat[j,j-1] = -1/(x[j]-x[j-1])
                A_mat[j,j] = 1/(x[j]-x[j-1])
                
                B_mat[j,j] = integrate.quad(lambda var: b_func(var)*phi_i(var,j)**2, x[j-1], x[j])[0]
                B_mat[j,j-1] = integrate.quad(lambda var: b_func(var)*phi_i(var,j)*phi_i(var,j-1), x[j-1], x[j])[0]
                
                C_mat[j,j] = integrate.quad(lambda var: a_func(var)*d_phi_i(var,j)*phi_i(var,j), x[j-1], x[j])[0]
                C_mat[j,j-1] = integrate.quad(lambda var: a_func(var)*d_phi_i(var,j-1)*phi_i(var,j), x[j-1], x[j])[0]
                 
                
            elif j != n-1 and j != 0:
                A_mat[j,j] = 1/(x[j]-x[j-1])+1/(x[j+1]-x[j])
                A_mat[j,j+1] = -1/(x[j+1]-x[j])
                A_mat[j,j-1] = -1/(x[j]-x[j-1])
                
                B_mat[j,j] = integrate.quad(lambda var: b_func(var)*phi_i(var,j)**2, x[j-1], x[j+1])[0]
                B_mat[j,j+1] = integrate.quad(lambda var: b_func(var)*phi_i(var,j)*phi_i(var,j+1), x[j], x[j+1])[0]
                B_mat[j,j-1] = integrate.quad(lambda var: b_func(var)*phi_i(var,j)*phi_i(var,j-1), x[j-1], x[j])[0]
                
                C_mat[j,j] = integrate.quad(lambda var: a_func(var)*d_phi_i(var,j)*phi_i(var,j), x[j-1], x[j+1])[0]
                
                
                C_mat[j,j+1] = integrate.quad(lambda var: a_func(var)*d_phi_i(var,j+1)*phi_i(var,j), x[j], x[j+1])[0]
                C_mat[j,j-1] = integrate.quad(lambda var: a_func(var)*d_phi_i(var,j-1)*phi_i(var,j), x[j-1], x[j])[0]
                
                Integral_1 = integrate.quad(lambda var: Integrand_1(var,j), x[j-1], x[j])[0]
                Integral_2 = integrate.quad(lambda var: Integrand_2(var,j), x[j], x[j+1])[0]
                b_vec[j] = Integral_1+Integral_2
                
        mat_sum = -A_mat+B_mat+C_mat

    else:
        
        for j in range(n): 
            
            if j == 0 and leftBC == "neumann":
                
                A_mat[j,j+1] = -1/(x[j+1]-x[j])
                B_mat[j,j+1] = (x[j+1]-x[j])/6
                C_mat[j,j+1] = 1/2 #error if - sign, will check calculations
                
                A_mat[j,j] = 1/(x[j+1]-x[j])
                B_mat[j,j] = (x[j+1]-x[j])/3
                C_mat[j,j] = -1/2
            
            elif j == n-1 and rightBC == "neumann":
                
                A_mat[j,j-1] = -1/(x[j]-x[j-1])
                B_mat[j,j-1] = (x[j]-x[j-1])/6
                C_mat[j,j-1] = -1/2
                
                A_mat[j,j] = 1/(x[j]-x[j-1])
                B_mat[j,j] = (x[j]-x[j-1])/3
                C_mat[j,j] = 1/2
                
            elif j != n-1 and j != 0:
                B_mat[j,j] = (-x[j-1])/3+(x[j+1])/3
                A_mat[j,j] = 1/(x[j]-x[j-1])+1/(x[j+1]-x[j]) 
                C_mat[j,j] = 0
                
                B_mat[j,j+1] = (x[j+1]-x[j])/6 
                A_mat[j,j+1] = -1/(x[j+1]-x[j]) 
                C_mat[j,j+1] = +1/2
                
                B_mat[j,j-1] = (x[j]-x[j-1])/6 
                A_mat[j,j-1] = -1/(x[j]-x[j-1])
                C_mat[j,j-1] = -1/2
                
                Integral_1 = integrate.quad(lambda var: Integrand_1(var,j), x[j-1], x[j])[0]
                Integral_2 = integrate.quad(lambda var: Integrand_2(var,j), x[j], x[j+1])[0]
                b_vec[j] = Integral_1+Integral_2
                
                
        mat_sum = -A_mat+b*B_mat+a*C_mat
   
    return mat_sum


if constant_coeff is True:
    mat_sum = AssembleSystem(const_coeff=True)
else:
    mat_sum = AssembleSystem()


if leftBC == "dirichlet":
    mat_sum[0,0] = 1
    b_vec[0] = a1
else:
    b_vec[0] += a1
    
if rightBC == "dirichlet":
    mat_sum[n-1,n-1] = 1
    b_vec[n-1] = a2
else:
    b_vec[n-1] -= a2
    

sol = np.linalg.solve(mat_sum, b_vec)


if plot_analytical_sol is True:
    plot_exact_sol()

#plt.grid(color='gainsboro')
plt.xticks(np.linspace(0,B,11))
plt.title("n = "+str(vertex_number))
plt.xlabel("$x$",fontsize = 12)
plt.ylabel("solution",fontsize = 12)
plt.plot(x,sol, color = "firebrick", linestyle = "-", label = "FE solution")


def drawBasisFunction(ind_list, col="k", solution = True):
    
    if len(ind_list) == 3:
        ind1 = ind_list[0]
        ind2 = ind_list[1]
        ind3 = ind_list[2]
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
        
    elif len(ind_list) == 2:
        ind1 = ind_list[0]
        ind2 = ind_list[1]
        x_interval = np.linspace(x[ind1], x[ind2], 20)
        
        if solution == True:
            y_interval = np.linspace(0,sol[ind2],20)
        else:
            y_interval = np.linspace(0,1,20)
  
        plt.plot(x_interval, y_interval, linestyle = "--", color = col)
        
if plot_hat_function == True:
    
    drawBasisFunction([1, 0])
    drawBasisFunction([len(sol)-2, len(sol)-1])
    for i in range(1,len(sol)-1):
        drawBasisFunction([i-1, i, i+1], solution=1)
    


if plot_legend is True:
    plt.legend()
