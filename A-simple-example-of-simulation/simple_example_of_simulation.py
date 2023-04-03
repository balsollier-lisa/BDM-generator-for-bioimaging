#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: lisabalsollier
"""

"""You will find in this file all the necessary programs to generate with 
the functions of process.py file, a simple example of BDM process."""


import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pickle
import csv
import math

import process
import draw


"""
Definition of auxiliary functions used in the main function process.proctotal
"""


"""this program defines the motion of the particles between jumps,
in this example, it is a Brownian motion"""
def move(t,x,y,r,c,d) :
    """
    Parameters
    ----------
    t : FLOAT
        duration of the motion that we want to generate.
    x : ARRAY
        vector that contains the x-coordinates of the points that form the initial condition.
    y : ARRAY
        vector that contains the y-coordinates of the points that form the initial condition.
    d : FLOAT
        pace discretization of the motion.

    Returns
    -------
    res : ARRAY
        table that contains all coordinates of the simulated points. 
        There are as many lines as there are discretizations of the duration time t. 
        Each line of the array contains all existing points at this time (in order x-coordinate, y-coordinate).
    ecarts : ARRAY
        vector that contains the time between two lines of resfinal.
    tps : ARRAY
        cumsum of ecarts 

    """
    P=len(x)
    if t>4*d :
        tps=np.concatenate((np.arange(0,t,d),np.array([t])))
        N=len(tps)
        ecarts=[d]*(len(tps)-2)+[t-(d*(len(tps)-2))]
    else :
        N=4
        tps=np.linspace(0,t,N)
        ecarts=[tps[i+1]-tps[i] for i in range(len(tps)-1)]
    res= np.zeros((N,4*P))
    for i in range (P) :
        normx=stats.norm.rvs(loc=0, scale=s1*np.sqrt(np.array(ecarts)))
        normy=stats.norm.rvs(loc=0, scale=s1*np.sqrt(np.array(ecarts)))
        xx=np.cumsum(np.concatenate((np.array([x[i]]),normx)))
        ide=np.where(np.abs(xx)>1)[0]
        while len(ide)!=0:
            k=ide[0]
            normx[k-1]=stats.norm.rvs(loc=0, scale=s1*np.sqrt(ecarts[k-1]))
            xx=np.cumsum(np.concatenate((np.array([x[i]]),normx)))
            ide=np.where(np.abs(xx)>1)[0]
        yy=np.cumsum(np.concatenate((np.array([y[i]]),normy)))
        idey=np.where(np.abs(yy)>1)[0]
        while len(idey)!=0:
            l=idey[0]
            normy[l-1]=stats.norm.rvs(loc=0, scale=s1*np.sqrt(ecarts[l-1]))
            yy=np.cumsum(np.concatenate((np.array([y[i]]),normy)))
            idey=np.where(np.abs(yy)>1)[0]
        res[:,4*i]=xx
        res[:,4*i+1]=yy
        cc=[1]*(N)
        rr=[0]*(N)
        res[:,4*i+2]=rr
        res[:,4*i+3]=cc
    return(res,ecarts,tps)


        
"""birth intensity"""
def beta(x) : 
    return(0.5)


"""death intensity"""
def delta(x): 
    if process.n(x)==1:
        return(0)
    else :
        return(0.5)
    
"""transition intensity"""
def tau(x) : 
    return(0)    
    
def alpha(x) :
    return(beta(x)+delta(x))  




"""birth kernel"""
def birthkernel(start): 
    X=stats.uniform.rvs(loc=-1,scale=2)
    Y=stats.uniform.rvs(loc=-1,scale=2)  
    return(X,Y,0,1)


"""death kernel"""
def deathkernel(start): 
    P=int(len(start)/4)
    i=stats.randint.rvs(0,P-1)
    tab=len(start)*[True]
    tab[4*i]=False
    tab[4*i+1]=False
    tab[4*i+2]=False
    tab[4*i+3]=False
    return(start[tab],0,1,i)

"""transition kernel"""
def transitionkernel(start): 
    return(None)



"""this program generates a random initial situation"""
def generatesituation(M):
    """
    Parameters
    ----------
    M : INTEGER
        number of particles of the configuration that will be returned.

    Returns
    -------
    X : ARRAY
        vector that contains the abscissae of the points that form the configuration that will be returned.
    Y : ARRAY
        vector that contains the ordinates of the points that form the configuration that will be returned.
    C : TYPE
        vector that contains the style of movement (1 for brownian, 2 for superdiffusive and 3 for confined) of the points that form the configuration that will be returned.
    """
    X=[]
    Y=[]
    C=[]
    R=[]
    for k in range(M):
        x=stats.uniform.rvs(loc=-1,scale=2)
        y=stats.uniform.rvs(loc=-1,scale=2)
        X+=[x]
        Y+=[y]
        C+=[1]
        R+=[0]
    track=[[i for i in range(1,M+1)]]
    cpttrack=M
    return(X,Y,R,C,track,cpttrack)



"""
Main: Simulation
"""

s1=0.02 # drift of the Brownian motion
T=10 # final time of the simulation
n_init= 10 # number of particles at initial time
d=0.1 # discretization step
Delta=1 # small time length for piecewise simulation



b=process.proctotal(T,n_init,d, Delta, generatesituation, move, beta, delta, tau, alpha, birthkernel, deathkernel, transitionkernel) 

(resfinal,TpsSauts,tabecarts,track,compteurs,tracknaissance, trackmort) =b   

"""
Post-processing: several plots concerning the Langerin channel
"""
 
#to draw the simulation as a movie
draw.drawsimpleex(b)


#extraction of the Langerin trajectories and other characteristics
intertps=0.15 # extraction each 0.15ms (intertps=d for full extraction)
(traj, tracktraj,trajmutantes,trajdec,tracktrajdec,couleurstrajdec, restronque,restottronque, tracktottronque)=draw.extract_traj(tabecarts, resfinal, track, 0, intertps)


#to draw Langerin trajectories
xcontour=[-1,-1,1,1,-1]
ycontour=[-1,1,1,-1,-1]
draw.traj(trajdec,couleurstrajdec, xcontour, ycontour)    

#to draw boxplots of the number of particles per frame
draw.boxplot_mean_particles(restronque,1)

#to draw histograms of the length of trajectories
draw.length_traj(traj, trajmutantes,1)

