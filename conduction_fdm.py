#!/usr/bin/env python
# coding: utf-8
import numpy as np

def initialize(x: int, y: int):
    '''This function initializes the domain for computation based on the length and bredth of the geometry
        INPUT : Length,Bredth
        OUTPUT: 2D array
        Aux INPUT : Mesh Element dimensions
        '''
    dx = input ("Enter the element size on X-axis")
    dy = input ("Enter the element size on Y-axis")
    domain = np.ones((x//dx,y//dy))
    domain = domain*25
    return domain

def boundary_condtions (domain: np.ndarray,t1:float=25, t2:float=25, t3:float=25, t4:float=25):
    '''Replaces the boundary values in the domain with the speccified values
        INPUT : Temperatures(t1,t2,t3,t4)
        |Takes Room Temperature if not specified|
        OUTPUT : Changes are made in-place at the array'''
    n1,n2 = domain.shape
    domain[:,0] = t1
    domain[0] = t2
    domain[:,n2-1] = t3
    domain[n1-1] = t4
    return domain

def fdm(domain: np.ndarray):
    n1,n2 = domain.shape
    for i in range(1,n1-1):
        for j in range(1,n2-1):
            domain[i][j] = domain[i-1][j] + domain [i][j+1] + domain[i+1][j] + domain[i][j-1]
