'''
Searches for root of optimal fill equation function, namely alpha0, by looking for change in sign
Requires knowing the function before hand 
!!!NOT GENERALLY APPLICABLE!!!
'''

import numpy as np
from sympy import *

def getalpha0(alpeqn):

    # y = symbols('y')
    # func = (y-sin(y)*cos(y))/(y**2) - alpeqn # Creates function
    # f = lambdify(y,func, "numpy")
    
    guess = 0.01 # intial guess of angle
    it = 0.0001
    # scale = np.arange(0,1.6,0.0001) # Iterates over this many
    # y = np.ones_like(scale)*guess + scale # scales y from 
    # ans = f(y) # inputs y
    y = guess
    scale = 1.6/it
    
    for i in range(0,1.6,scale):
        func = (y-sin(y)*cos(y))/(y**2) - alpeqn
        
        if func >= 0: # root is transition from negative to positive in this case
            alp0 = y - it # previous value before sign flip returned       
            break
        else:
            y += it
            continue
        print(alp0)
    return  alp0  

'''
Solves inverse lp equation, knowing lp, to allow determination of deformed angle
'''
def getalpha(A,lp,n,alpha0):
   
    y = symbols('y')
    func = sqrt((2*A*y**2)/(y-sin(y)*cos(y))) - lp
    f = lambdify(y,func, "numpy")
    
    guess = alpha0 # intial guess of angle previous
    scale = np.arange(0,1.6,0.0001) # Iterates over this possible range, possible redundant 
    y = np.ones_like(scale)*guess + scale # scales y from guess
    ans = f(y) # inputs y

    for i in range(0,len(scale)):
        if ans[i] <= 0: # root is transition from pos -> neg in this case
            alp = y[i-1]     
            break
        else:
            continue
    return alp