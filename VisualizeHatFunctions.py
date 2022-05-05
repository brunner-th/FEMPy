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

num_vert = 6

x = np.linspace(0,10,num_vert)
sol = np.zeros((num_vert))
sol.fill(1)

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
        plt.plot(x_int, y_int, linestyle = "--", color = col, linewidth = 1.3)
        
    elif len(ind_list) == 2:
        ind1 = ind_list[0]
        ind2 = ind_list[1]
        x_interval = np.linspace(x[ind1], x[ind2], 20)
        
        if solution == True:
            y_interval = np.linspace(0,sol[ind2],20)
        else:
            y_interval = np.linspace(0,1,20)
            
        
        plt.plot(x_interval, y_interval, linestyle = "--", color = col, linewidth = 1.3)
        
plt.title("Linear Basis Functions")
plt.xticks(np.linspace(0,10,11))
plt.yticks(np.linspace(0,1,6))
#plt.grid()
drawBasisFunction([1, 0])
drawBasisFunction([len(sol)-2, len(sol)-1])
for i in range(1,len(sol)-1):
    drawBasisFunction([i-1, i, i+1], solution=1)