#!/usr/bin/env python
# coding: utf-8
import numpy as np

def initialize(x: int, y: int):
    
    '''This function initializes the domain for computation based on the length and breadth of the geometry
    
        INPUT : Length,Bredth
        OUTPUT: 2D array 
        
        Aux INPUT : Mesh Element dimensions dx ad dy
        '''
    while True:
        # the loop helps to restart the initialization procedure untill sucessful completion.
        try:
            # dx and dy expects a positive value which is less than the value x and y
            dx = float(input ("Enter the element size on X-axis "))
            dy = float(input ("Enter the element size on Y-axis "))
            if ((dx<=0 or dy<=0) or (dx>x or dy>y)):
                raise ValueError
            break
        except ValueError:
            print("Error: Recieved improper Value or datatype")
    print("Log: Domain created succesfully")
    return np.zeros((int(x//dx),int(y//dy)))

def boundary_condtions (domain: np.ndarray,t1:float=25, t2:float=25, t3:float=25, t4:float=25):
    
    '''Replaces the boundary values in the domain with the speccified values
        INPUT : Temperatures(t1,t2,t3,t4)
        note  : if parameters are not provided, the function assumes a deault 
                value of 25 (room temp in deg C)
        |Takes Room Temperature if not specified|
        OUTPUT : Changes are made in-place at the array'''
    
    n1,n2 = domain.shape
    domain[:,0] = t1
    domain[0] = t2
    domain[:,n2-1] = t3
    domain[n1-1] = t4
    return domain

def fdm(domain: np.ndarray,analysis:str,type:str="l_to_r"):
    ''' The function takes in the domain and modifies the element in place 
        to reflect the change as per the physics of the problem 
        
        INPUT : domain, analysis_Type(Conduction,Convection),
                march pattern as left to right, right to left etc
                
        OUTPUT : 5 point approximation of the current node'''
    
    
    n1,n2 = domain.shape
    # switch between equtions
    equations = {'conduction':lambda M,i,j: (M[i][j-1] + M[i-1][j] + M[i][j+1] +M[i+1][j])/4,
                 'convection':lambda domain,i,j:35}
    #governing model
    model = equations[analysis]
    
    #selection of march pattern
    method = {"l_to_r":(1,n1-1,1,1,n2-1,1),
              "r_to_l":(n1-2,0,-1,1,n2-1,1),
              "t_to_b":(1,n2-1,1,1,n1-1,1),
              "b_to_t":(n2-2,0,-1,1,n1-1,1)}
    p,q,g1,r,s,g2 = method[type]
    
    # dereferencing node positions
    for i in range(p,q,g1):
        for j in range(r,s,g2):
            print(i,j,domain[i][j])
            domain[i][j] = model(domain,i,j)
            
    return domain
