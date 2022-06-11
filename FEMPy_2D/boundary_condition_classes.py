# -*- coding: utf-8 -*-

#-----------------------------------------------------------------
#
# name: thomas brunner
# email: brunner.th@hotmail.com
# nr: 12018550
#
#-----------------------------------------------------------------

from math import *
import numpy as np

class DirichletBoundaryCondition:
    
    def __init__(self, number_of_points):
        self.point_list = []
        self.val_list = np.zeros(number_of_points)
        self.condition_list = []  
         
    def add_dirichlet_points(self, point_list, dirichlet_val):
        for i, point in enumerate(point_list):
            self.point_list.append(point)
            self.val_list[point] = dirichlet_val
        self.condition_list.append([point_list, dirichlet_val])
            
    def points(self):
        return self.point_list
    
    def value(self):
        return self.val_list
    
    def dirichlet_conditions(self):
        return self.condition_list


class CauchyBoundaryCondition:
    
    def __init__(self, number_of_points):
        self.facet_list = []
        self.gamma_list = np.zeros((number_of_points))
        self.a_list = np.zeros((number_of_points))
        
    def add_cauchy_facets(self, facet_list, a_val, gamma_val):
        for i, facet in enumerate(facet_list):
            self.facet_list.append(facet)
            self.a_list[facet[0]]=a_val
            self.a_list[facet[1]]=a_val
            self.gamma_list[facet[0]]=gamma_val
            self.gamma_list[facet[1]]=gamma_val
    
    def facets(self):
        return self.facet_list
    
    def a_val_list(self):
        return self.a_list
    
    def gamma_val_list(self):
        return self.gamma_list

        
        
        
        
        
