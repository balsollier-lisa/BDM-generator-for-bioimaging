#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: lisabalsollier
"""

"""You will find in this file all the necessary programs to generate 
the process."""


import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pickle
import csv
import math


"""this program gives the regime of the particle x, 1 for brownian, 2 for superdiffusive and 3 for subdiffusive"""
def typee(z):
    """
    Parameters
    ----------
    z : FLOAT
        Fourth coordinate of a point of a configuration, the one that gives the way to move.

    Returns
    -------
    INTEGER
    1, 2 or 3.

    """
    if z==1 or z==3 :
        return(int(z)) 
    else :
        return(2)
    
    
"""number of particles in the configuration x"""
def n(x):
    """
    Parameters
    ----------
    x : ARRAY
        table with a process configuration in the order: abscissa, ordinate, type, speed.

    Returns
    -------
    INTEGER
    number of particules in the configuration x.
    """
    return(int(len(x)/4))


"""number of Rab11 particles in the configuration x"""
def nr(x):
    """
    Parameters
    ----------
    x : ARRAY
        table with a process configuration in the order: abscissa, ordinate, type, speed.

    Returns
    -------
    INTEGER
    number of Rab11 particules in the configuration x.
    """
    return(len(x[2::4][x[2::4]==1]))


"""number of Langerin particles in the configuration x"""
def nl(x):
    """
    Parameters
    ----------
    x : ARRAY
        table with a process configuration in the order: abscissa, ordinate, type, speed.

    Returns
    -------
    INTEGER
    number of Langerin particules in the configuration x.
    """
    return(len(x[2::4][x[2::4]==0]))



"""number of Langerin brownian particles in the configuration x"""
def nlb(x):
    choix=[]
    for ch in x[2::4]==0  :
        choix+=[ch]*4
    y=np.array(x)[np.array(choix)]
    return(len(y[3::4][y[3::4]==1]))


"""number of Langerin superdiffusive particles in the configuration x"""
def nlsp(x):
    choix=[]
    for ch in x[2::4]==0  :
        choix+=[ch]*4
    y=np.array(x)[np.array(choix)]
    #choix2=[typee(y[3::4][i])==2 for i in range(len(y[3::4]))]
    return(len(y[3::4][np.array([typee(y[3::4][i])==2 for i in range(len(y[3::4]))])]))


"""number of Langeirn subdiffusive particles in the configuration x"""
def nlsb(x):
    choix=[]
    for ch in x[2::4]==0  :
        choix+=[ch]*4
    y=np.array(x)[np.array(choix)]
    return(len(y[3::4][y[3::4]==3]))


"""number of Rab11 brownian particles in the configuration x"""
def nrb(x):
    choix=[]
    for ch in x[2::4]==1  :
        choix+=[ch]*4
    y=np.array(x)[np.array(choix)]
    return(len(y[3::4][y[3::4]==1]))


"""number of Rab11 superdiffusive particles in the configuration x"""
def nrsp(x):
    choix=[]
    for ch in x[2::4]==1  :
        choix+=[ch]*4
    y=np.array(x)[np.array(choix)]
    #choix2=[typee(y[3::4][i])==2 for i in range(len(y[3::4]))]
    return(len(y[3::4][np.array([typee(y[3::4][i])==2 for i in range(len(y[3::4]))])]))


"""number of Rab11 subdiffusive particles in the configuration x"""
def nrsb(x):
    choix=[]
    for ch in x[2::4]==1  :
        choix+=[ch]*4
    y=np.array(x)[np.array(choix)]
    return(len(y[3::4][y[3::4]==3]))




"""This programm below is used to compute the integral necessary to generate the jump time density"""
def integrale(res,ecarts, alpha) :
    s=0
    for i in range(len(ecarts)-1):
        s+= ((alpha(res[i,:])+alpha(res[i+1,:]))/2)*ecarts[i]
    s+=alpha(res[-1,:])*ecarts[-1]
    return(s)

"""This program gives us the value of the jump time density at point s"""
def densityf(p,res,ecarts,tps,s, alpha):
    i=len(tps[s>=tps])-1
    a=(1/(1-p))*alpha(res[i,:])*np.exp(-integrale(res[:i+1,:],ecarts[:i+1], alpha))
    if i<len(tps)-1 :
        b=(1/(1-p))*alpha(res[i+1,:])*np.exp(-integrale(res[:i+2,:],ecarts[:i+2], alpha))
        return((a+b)/2)
    else :
        return(a)
        
"""This program allows to simulate a value according to the jump time law."""
def drawlaw(p,res,ecarts,tps,tpsmax, alpha):
    tabf=[]
    for i in range(len(ecarts)-1):
        tabf+=[(1/(1-p))*alpha(res[i,:])*np.exp(-integrale(res[:i+1,:],ecarts[:i+1], alpha))]
    M=np.max(tabf)
    U=stats.uniform.rvs(0,tpsmax)
    V=stats.uniform.rvs(0,M)
    while V>=densityf(p,res,ecarts,tps,U, alpha) :
        U=stats.uniform.rvs(0,tpsmax)
        V=stats.uniform.rvs(0,M)
    return(U)


"""global program to generate the BDM process """
def proctotal(T,M,d, Delta, generatesituation, move, beta, delta, tau, alpha, birthkernel, deathkernel, transitionkernel):   
    """
    Parameters
    ----------
    T : FLOAT
        final time of the simulation.
    x : ARRAY
        vector that contains the abscissae of the points that form the initial condition.
    y : ARRAY
        vector that contains the ordinates of the points that form the initial condition.
    r : ARRAY
        vector that contains the type (0 for Langerin and 1 for Rab11 of the points that form the initial condition.
    c : TYPE
        vector that contains the style of movement (1 for brownian, 2 for superdiffusive and 3 for confined) of the points that form the initial condition.
    d : FLOAT
        pace discretization of the motion.
    move : FUNCTION, the default is move.
        Function that simulates the movement of particles between two jump times.
        The default is move.
    beta : FUNCTION, the default is beta.
        Function that gives the birth intensity.
    delta : FUNCTION, the default is delta.
        Function that gives the death intensity.
    tau : TYPE, the default is tau.
        Function that gives the transition intensity. 
    birthkernel : FUNCTION, the default is birthkernel.
        Function that gives the birth kernel.
    deathkernel : FUNCTION, the default is deathkernel.
        Function that gives the death kernel. 
    transitionkernel : FUNCTION, the default is transitionkernel.
        Function that gives the transition kernel. 

    Returns
    -------
    resfinal : ARRAY
        list of tables that contains all the coordinates of the points of the simulation between two jump times. Each lines of the array is all the points present.
    TpsSaut : ARRAY
        vector that contains the jump times.
    tabecarts : ARRAY
        list of list that contains the times between two lines of resfinal.
    """
    (x,y,r,c,track, cpttrack)=generatesituation(M)
    #print(c)
    t=0
    Tj=0
    TpsSauts=[0]
    resfinal=[]
    tabecarts=[]
    cptnlb=0
    cptnlsp=0
    cptnlsb=0
    cptnrb=0
    cptnrsp=0
    cptnrsb=0
    cptmlb=0
    cptmlsp=0
    cptmlsb=0
    cptmrb=0
    cptmrsp=0
    cptmrsb=0
    cptnaissances=0
    cptmorts=0
    trackmort=[]
    tracknaissance=[]
    while t<T :
        k=1
        #p=[1]
        res=[]
        #ecarts=[]
        while k!=0 :
            Deltaj=np.min([Delta,T-Tj-(k-1)*Delta])
            (resk,ecartsk,tpsk)=move(Deltaj,x,y,r,c,d)
            if k==1 :
                res=resk[:-1,:]
                ecarts=ecartsk
                tps=tpsk
            else:
                res=np.vstack((res,resk[:-1,:]))
                ecarts+=ecartsk
                tps=np.concatenate((tps,tps[-1]+tpsk[1:])) 
            p=np.exp(-integrale(resk,ecartsk, alpha))
            U=stats.uniform.rvs(loc=0,scale=1)
            if U<=p:
                if Deltaj==T-Tj-(k-1)*Delta:
                    resfinal+=[res]
                    tabecarts+=[ecarts]
                    t=T
                    k=0
                    #print('la')
                else:
                    #taujmax=np.min([(k+1)*taumax,T-Tj])
                    k+=1
                    x=resk[-1,0::4]
                    y=resk[-1,1::4]
                    r=resk[-1,2::4]
                    c=resk[-1,3::4]
                    #print('ici')
            else:
                #print('jesuisla')
                #print(Deltaj)
                tauj=(k-1)*Delta+drawlaw(p,resk,ecartsk,tpsk,Deltaj, alpha)#+(k-1)*Delta)
                Tj+=tauj
                #print(Tj)
                TpsSauts+=[Tj]
                i=len(tps[tauj>tps])
                resfinal+=[res[:i,:]]
                tabecarts+=[ecarts[:i-1]+[tauj-tps[tauj>tps][-1]]]
                U2=stats.uniform.rvs(0,1)
                if (U2<= (beta(res[i-1,:])/alpha(res[i-1,:]))) :
                    cptnaissances+=1
                    (X,Y,R,C)=birthkernel(res[i-1,:])
                    #print('naissance')
                    if R==0 :
                        if typee(C)==1 :
                            cptnlb+=1
                        elif typee(C)==3 :
                            cptnlsb+=1
                        else :
                            cptnlsp+=1
                    else :
                        if typee(C)==1 :
                            cptnrb+=1
                        elif typee(C)==3 :
                            cptnrsb+=1
                        else :
                            cptnrsp+=1
                    #print(cptnaissances,cptmorts,cptnlb,cptnlsp,cptnlsb,cptnrb,cptnrsp,cptnrsb,cptmlb,cptmlsp,cptmlsb,cptmrb,cptmrsp,cptmrsb)
                    newY=np.concatenate((res[i-1,:],np.array([X,Y,R,C])))
                    track+=[track[-1]+[cpttrack+1]]
                    cpttrack+=1
                    tracknaissance+=[cpttrack]
                elif (U2<=((beta(res[i-1,:])+tau(res[i-1,:]))/alpha(res[i-1,:]))):
                    #print(i-1)
                    newY=transitionkernel(res[i-1,:])
                    track+=[track[-1]]
                    #print('mutation')
                else:
                    #print('mort')
                    cptmorts+=1
                    newY,typemort,couleurmort,num=deathkernel(res[i-1,:])
                    #print(typemort,couleurmort)
                    if typemort==0 :
                        if couleurmort==1 :
                            cptmlb+=1
                        elif couleurmort==3 :
                            cptmlsb+=1
                        else :
                            cptmlsp+=1
                    else :
                        if couleurmort==1 :
                            cptmrb+=1
                        elif couleurmort==3 :
                            cptmrsb+=1
                        else :
                            cptmrsp+=1
                    trackmort+=[track[-1][num]]
                    track+=[track[-1][:num]+track[-1][num+1:]]
                    #print(cptnaissances,cptmorts,cptnlb,cptnlsp,cptnlsb,cptnrb,cptnrsp,cptnrsb,cptmlb,cptmlsp,cptmlsb,cptmrb,cptmrsp,cptmrsb)
                k=0  
                t+=tauj
                x=newY[0::4]
                y=newY[1::4]
                r=newY[2::4]
                c=newY[3::4]
            #print(k)
    compteurs=[cptnaissances,cptmorts,cptnlb,cptnlsp,cptnlsb,cptnrb,cptnrsp,cptnrsb,cptmlb,cptmlsp,cptmlsb,cptmrb,cptmrsp,cptmrsb]
    return(resfinal,TpsSauts,tabecarts,track,compteurs,tracknaissance, trackmort)    





