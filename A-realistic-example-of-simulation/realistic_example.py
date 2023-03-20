#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: lisabalsollier
"""

"""
You will find in this file all the necessary programs to generate with the 
functions of process.py file, a simulation of BDM which approaches the real data.
This file was used to generate the simulations presented in the article.
"""


import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pickle
import csv
import math

import process
import draw


#we import the data from the real dataset

#for Langerin particles
with open('estcoeflang.pickle', 'rb') as handle:
    brlang, splangdiff, splangdrift = pickle.load(handle)
    
with open('angledev.pickle', 'rb') as handle:
    angledev  = pickle.load(handle)

s3lang=0.78 #diffusion coefficient in the subdiffusive move
lalang=5.96 #  

# for rab11 particles
with open('estcoefrab.pickle', 'rb') as handle:
    brrab, sprabdiff, sprabdrift  = pickle.load(handle)
    
s3rab=1.19 #diffusion coefficient in the subdiffusive move    
larab=7.41 #coefficient in the Ornstein-Uhlenbeck process of subdiffusive move



prbirth=(1/2.98)*np.array([2.34,0.14,0.5]) #proportions for the birth of brownian, superdiffusive, subdiffusive Rab11 (in this order)
plbirth=(1/4.45)*np.array([3.59,0.47,0.39]) #proportions for the birth of brownian, superdiffusive, subdiffusive Langerin (in this order)
X0=np.array([97,141.2])



T=167.86 #final time of the simulation

d=0.01 #size of discretization pace between 2 images
Delta=0.2 # taumax of the alogrithm

theta=1.11 #parameter of variance in the birth kernel
q=0.069 #parameter that rules the uniform birth in the birth kernel



#we load the data of the cell tour

with open('maskc.pickle', 'rb') as handle:
    imgmod = pickle.load(handle)  

with open('xcontour.pickle', 'rb') as handle:
    xcontour = pickle.load(handle)  

with open('ycontour.pickle', 'rb') as handle:
    ycontour = pickle.load(handle)  


"""this program allows to say if a point is in the cell"""
def beincell(x):
    if np.floor(x[0])<=0 or np.floor(x[0])+1>=250:
        return(0)
    if np.floor(x[1])<=0 or np.floor(x[1])+1>=283:
        return(0)
    if np.floor(x[0])+0.5<x[0]:
        j=int(np.floor(x[0]))+1
    else :
        j=int(np.floor(x[0]))
    if np.floor(x[1])+0.5<x[1]:
        i=int(np.floor(x[1]))+1
    else :
        i=int(np.floor(x[1]))
    return(int(imgmod[i][j][0]))


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
    R : ARRAY
        vector that contains the type (0 for Langerin and 1 for Rab11 of the points that form the configuration that will be returned.
    C : TYPE
        vector that contains the style of movement (1 for brownian, 2 for superdiffusive and 3 for confined) of the points that form the configuration that will be returned.

    """
    X=[]
    Y=[]
    C=[]
    R=[]
    for k in range(M):
        x=stats.uniform.rvs(loc=0,scale=249)
        y=stats.uniform.rvs(loc=0,scale=282)
        while beincell([x,y])*True==False :
            x=stats.uniform.rvs(loc=0,scale=249)
            y=stats.uniform.rvs(loc=0,scale=282)
        X+=[x]
        Y+=[y]
        #r=0
        r=np.random.choice([0,1],p=[0.62,0.38])
        R+=[r]
        if r==0 :
            c=np.random.choice([1,2,3],p=plbirth)
        else :
            c=np.random.choice([1,2,3],p=prbirth)
        if c==2 :
            if r==1:
                C+=[stats.uniform.rvs(0,2*np.pi)]
            else :
                vect= -np.array([x,y])+X0
                C+=[math.atan2(vect[1],vect[0])+np.random.choice(angledev)]
        else:
            C+=[c]
        track=[[i for i in range(1,M+1)]]
        cpttrack=M
    return(X,Y,R,C,track,cpttrack)


"""this program generates a Brownian motion of diffusion coefficient s1lang 
or s1rab with as initial position the point (xi,yi,ri,ci)"""
def brownianmove(xi,yi,ri,ci,ecarts,N,s1lang,s1rab):
    if ri==0 :
        normx=stats.norm.rvs(loc=0, scale=s1lang*np.sqrt(np.array(ecarts)))
        normy=stats.norm.rvs(loc=0, scale=s1lang*np.sqrt(np.array(ecarts)))
        xx=np.cumsum(np.concatenate((np.array([xi]),normx)))
        yy=np.cumsum(np.concatenate((np.array([yi]),normy)))
        co=[[xx[m],yy[m]] for m in range(len(xx))]
        #print(co)
        ide=np.argwhere(np.array([beincell(co[p]) for p in range(len(co))])==0)
        #print('chgmt')
        #print(ide)
        while len(ide)!=0:
            k=ide[0][0]
            #print(i)
            normx[k-1]=stats.norm.rvs(loc=0, scale=s1lang*np.sqrt(ecarts[k-1]))
            normy[k-1]=stats.norm.rvs(loc=0, scale=s1lang*np.sqrt(ecarts[k-1]))
            #normx[k-1]=normx[k-1]-2*np.sign(xx[k])
            #print(normx)
            xx=np.cumsum(np.concatenate((np.array([xi]),normx)))
            yy=np.cumsum(np.concatenate((np.array([yi]),normy)))
            co=[[xx[m],yy[m]] for m in range(len(xx))]
            ide=np.argwhere(np.array([beincell(co[p]) for p in range(len(co))])==0)
            #print(ide)
    else :
        normx=stats.norm.rvs(loc=0, scale=s1rab*np.sqrt(np.array(ecarts)))
        normy=stats.norm.rvs(loc=0, scale=s1rab*np.sqrt(np.array(ecarts)))
        xx=np.cumsum(np.concatenate((np.array([xi]),normx)))
        yy=np.cumsum(np.concatenate((np.array([yi]),normy)))
        co=[[xx[m],yy[m]] for m in range(len(xx))]
        #print(co)
        ide=np.argwhere(np.array([beincell(co[p]) for p in range(len(co))])==0)
        #print('chgmt')
        #print(ide)
        while len(ide)!=0:
            k=ide[0][0]
            #print(i)
            normx[k-1]=stats.norm.rvs(loc=0, scale=s1rab*np.sqrt(ecarts[k-1]))
            normy[k-1]=stats.norm.rvs(loc=0, scale=s1rab*np.sqrt(ecarts[k-1]))
            #normx[k-1]=normx[k-1]-2*np.sign(xx[k])
            #print(normx)
            xx=np.cumsum(np.concatenate((np.array([xi]),normx)))
            yy=np.cumsum(np.concatenate((np.array([yi]),normy)))
            co=[[xx[m],yy[m]] for m in range(len(xx))]
            ide=np.argwhere(np.array([beincell(co[p]) for p in range(len(co))])==0)
            #print(ide)
    cc=[ci]*(N)
    rr=[ri]*(N)
    #print('fini')
    return(xx,yy,rr,cc)


"""this program generates a subdiffusive move with diffusion coefficient 
s3lang or s3rab and attraction coefficient lalang or larab 
with as initial position the point (xi,yi,ri,ci)"""
def subdiffusifmove(xi,yi,ri,ci,ecarts,tps,N,s3lang,s3rab, lalang, larab):
    if ri==0:
        normx=stats.norm.rvs(loc=0, scale=s3lang*np.sqrt(np.array(ecarts)))
        normy=stats.norm.rvs(loc=0, scale=s3lang*np.sqrt(np.array(ecarts)))
        xx=[xi]
        yy=[yi]
        n=0
        expla=[(np.exp(lalang*tps[i+1])+np.exp(lalang*tps[i]))/2 for i in range(len(tps)-1)]
        while n <= (len(ecarts)-1):
            nx=xi+np.exp(-lalang*tps[n+1])*(np.sum((expla*normx)[:n+1]))
            #print(nx)
            ny=yi+np.exp(-lalang*tps[n+1])*(np.sum((expla*normy)[:n+1]))
            #print(ny)
            if beincell([nx,ny])==1 :
                xx+=[nx]
                yy+=[ny]
                n=n+1
                #print(n)
            else:
                normx[n]=stats.norm.rvs(loc=0, scale=s3lang*np.sqrt(ecarts[n]))
                normy[n]=stats.norm.rvs(loc=0, scale=s3lang*np.sqrt(ecarts[n]))
            #else:
               # normx[n]=-normx[n]
               # normy[n]=-normy[n]
    else :
        normx=stats.norm.rvs(loc=0, scale=s3rab*np.sqrt(np.array(ecarts)))
        normy=stats.norm.rvs(loc=0, scale=s3rab*np.sqrt(np.array(ecarts)))
        xx=[xi]
        yy=[yi]
        n=0
        expla=[(np.exp(larab*tps[i+1])+np.exp(larab*tps[i]))/2 for i in range(len(tps)-1)]
        while n <= (len(ecarts)-1):
            nx=xi+np.exp(-larab*tps[n+1])*(np.sum((expla*normx)[:n+1]))
            #print(nx)
            ny=yi+np.exp(-larab*tps[n+1])*(np.sum((expla*normy)[:n+1]))
            #print(ny)
            if beincell([nx,ny])==1 :
                xx+=[nx]
                yy+=[ny]
                n=n+1
                #print(n)
            else:
                normx[n]=stats.norm.rvs(loc=0, scale=s3rab*np.sqrt(ecarts[n]))
                normy[n]=stats.norm.rvs(loc=0, scale=s3rab*np.sqrt(ecarts[n]))
            #else:
               # normx[n]=-normx[n]
               # normy[n]=-normy[n]
    cc=[ci]*(N)
    rr=[ri]*(N)
    return(xx,yy,rr,cc)


"""this program generates a subdiffusive move with diffusion coefficient 
s3lang or s3rab and drift coefficient drlang or drrab 
with as initial position the point (xi,yi,ri,ci)"""
def superdiffusifmove(xi,yi,ri,ci,ecarts,N,s2lang,s2rab,drlang, drrab):
    if ri==0:
        normx=stats.norm.rvs(loc=0, scale=s2lang*np.sqrt(np.array(ecarts)))
        normy=stats.norm.rvs(loc=0, scale=s2lang*np.sqrt(np.array(ecarts)))
        xx=[xi]
        yy=[yi]
        angle=ci
        cc=[ci]*(N)
        m=0
        #q=0
        while m <= (len(ecarts)-1):
            #if q==m+10:
             #   print('nul')
             #   return('nul')
            nx=xx[-1]+normx[m]+drlang*np.cos(angle)*ecarts[m]
            ny=yy[-1]+normy[m]+drlang*np.sin(angle)*ecarts[m]
            #print([nx,ny])
            if beincell([nx,ny])==1 :
                xx+=[nx]
                yy+=[ny]
                m=m+1
                #print(m)
            else :
                angle=stats.uniform.rvs(0,2*np.pi)
                normx[m]=stats.norm.rvs(loc=0, scale=s2lang*np.sqrt(ecarts[m]))
                normy[m]=stats.norm.rvs(loc=0, scale=s2lang*np.sqrt(ecarts[m]))
                cc[m:]=[angle]*len(cc[m:])
                #print('ici')
    else :
        normx=stats.norm.rvs(loc=0, scale=s2rab*np.sqrt(np.array(ecarts)))
        normy=stats.norm.rvs(loc=0, scale=s2rab*np.sqrt(np.array(ecarts)))
        xx=[xi]
        yy=[yi]
        angle=ci
        cc=[ci]*(N)
        m=0
        #q=0
        while m <= (len(ecarts)-1):
            #if q==m+10:
             #   print('nul')
             #   return('nul')
            nx=xx[-1]+normx[m]+drrab*np.cos(angle)*ecarts[m]
            ny=yy[-1]+normy[m]+drrab*np.sin(angle)*ecarts[m]
            #print([nx,ny])
            if beincell([nx,ny])==1 :
                xx+=[nx]
                yy+=[ny]
                m=m+1
                #print(m)
            else :
                angle=stats.uniform.rvs(0,2*np.pi)
                normx[m]=stats.norm.rvs(loc=0, scale=s2rab*np.sqrt(ecarts[m]))
                normy[m]=stats.norm.rvs(loc=0, scale=s2rab*np.sqrt(ecarts[m]))
                cc[m:]=[angle]*len(cc[m:])
                #print('ici')
    rr=[ri]*(N)
    return(xx,yy,rr,cc)


"""this program generates the motion of a particle up to time t according to its coordinate c"""
def move(t,x,y,r,c,d) :
    """
    Parameters
    ----------
    t : FLOAT
        time of the movement that we want to generate.
    x : ARRAY
        vector that contains the abscissae of the points that form the initial condition.
    y : ARRAY
        vector that contains the ordinates of the points that form the initial condition.
    r : ARRAY
        vector that contains the type (0 for Langerin and 1 for Rab11 of the points that form the initial condition.
    c : TYPE
        vector that contains the style of movement (1 for brownian, 2 for superdiffusive and 3 for confined) of the points that form the initial condition..
    d : FLOAT
        pace discretization of the motion.

    Returns
    -------
    res : ARRAY
        table that contains all the coordinates of the points of the simulation. There are as many lines as there are discretizations of time t. Each lines of the array contains all the points present (in the order abscissa, ordinate, type, regime).
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
    #print(tps)
    #print(ecarts)
    res= np.zeros((N,4*P))
    for i in range (P) :
        if c[i]==1 :
            #print('br1')
            xx,yy,rr,cc=brownianmove(x[i],y[i],r[i],c[i],ecarts,N, np.sqrt(np.random.choice(brlang)),np.sqrt(np.random.choice(brrab)))
            #xx,yy,rr,cc=brownianmove(x[i],y[i],r[i],c[i],ecarts,N, s1lang,s1rab)
            #print('br2')
        elif c[i]==3 :
            #print('sub1')
            #xx,yy,rr,cc=subdiffusifmove(x[i],y[i],r[i],c[i],ecarts,tps,N,np.sqrt(np.random.choice(sblangdiff)),np.sqrt(np.random.choice(sbrabdiff)),np.random.choice(sblanglb), np.random.choice(sblanglb))
            #xx,yy,rr,cc=subdiffusifmove(x[i],y[i],r[i],c[i],ecarts,tps,N,np.sqrt(np.random.choice(sblangdiff)),np.sqrt(np.random.choice(sbrabdiff)),lalang, larab)
            xx,yy,rr,cc=subdiffusifmove(x[i],y[i],r[i],c[i],ecarts,tps,N,s3lang/2,s3rab/2,lalang, larab)
            #print('sub2')
        else :
            #print('sup1')
            xx,yy,rr,cc=superdiffusifmove(x[i],y[i],r[i],c[i],ecarts,N,np.sqrt(np.random.choice(splangdiff))/2,np.sqrt(np.random.choice(sprabdiff))/2, np.random.choice(splangdrift), np.random.choice(sprabdrift))
            #xx,yy,rr,cc=superdiffusifmove(x[i],y[i],r[i],c[i],ecarts,N,np.sqrt(np.random.choice(splangdiff)),np.sqrt(np.random.choice(sprabdiff)), drlang, drrab)
            #xx,yy,rr,cc=superdiffusifmove(x[i],y[i],r[i],c[i],ecarts,N,s2lang,s2rab, drlang, drrab)
            #print('sup2')
        res[:,4*i]=xx
        res[:,4*i+1]=yy
        res[:,4*i+3]=cc
        res[:,4*i+2]=rr
        #print(i==P-1)
    return(res,ecarts,tps)



"""birth intensity"""
def beta(x) : 
    #return(0)
    return(4.45+2.98)#*n(x))


"""death intensity"""
def delta(x): 
    if process.n(x)==1:
        return(0)
    else :
        #return(0)
        #return(0.12*nl(x)+0.14*nr(x))
        return(0.15*process.nlb(x)+0.1*process.nlsp(x)+0.058*process.nlsb(x)+0.18*process.nrb(x)+0.22*process.nrsp(x)+0.07*process.nrsb(x))
    
    
"""mutation intensity"""
def tau(x) : 
    return(16/167.86)
    #return(0)
    
    
def alpha(x) :
    return(beta(x)+delta(x)+tau(x))  



"""birth kernel"""
def birthkernel(depart):
    """
    Parameters
    ----------
    depart : ARRAY
        Array that contains the coordinates of the point configuration present
        just before the new birth (in the order abscissa, ordinate, type, regime).

    Returns
    -------
    The coordinates of the new point (in the order abscissa, ordinate, type, regime).
    """
    r=stats.uniform.rvs(0,1)
    if r<(2.98/(2.98+4.45)):
        R=1
        X=stats.uniform.rvs(loc=0,scale=249)
        Y=stats.uniform.rvs(loc=0,scale=282)
        while beincell([X,Y])==0 :
            X=stats.uniform.rvs(loc=0,scale=249)
            Y=stats.uniform.rvs(loc=0,scale=282)
        CC=np.random.choice([1,2,3],p=prbirth)
        if CC==2 :
            C=stats.uniform.rvs(0,2*np.pi)
        else :
            C=CC
    else :
        R=0
        p=stats.uniform.rvs(0,1)
        if p<q and process.nr(depart)>0:
            P=process.nr(depart)
            i=stats.randint.rvs(1,P+1)
            m=0
            s=0
            while s<i and m<len(depart)/4:
                if depart[4*m+2]== 1 :
                    s+=1
                m+=1
            X=stats.norm.rvs(loc=depart[4*(m-1)],scale=theta)
            Y=stats.norm.rvs(loc=depart[4*(m-1)+1],scale=theta)
            while beincell([X,Y])==0:
                X=stats.norm.rvs(loc=depart[4*(m-1)],scale=theta)
                Y=stats.norm.rvs(loc=depart[4*(m-1)+1],scale=theta)
            C=depart[4*(m-1)+3]
        else :
            X=stats.uniform.rvs(loc=0,scale=249)
            Y=stats.uniform.rvs(loc=0,scale=282)
            while beincell([X,Y])==0 :
                X=stats.uniform.rvs(loc=0,scale=249)
                Y=stats.uniform.rvs(loc=0,scale=282)
            CC=np.random.choice([1,2,3],p=plbirth)
            if CC==2 :
                vect= -np.array([X,Y])+X0
                C=math.atan2(vect[1],vect[0])+np.random.choice(angledev)
                #C=stats.uniform.rvs(0,2*np.pi)
            else :
                C=CC
    #R=stats.randint.rvs(0,2)
    return(X,Y,R,C)


"""death kernel"""
def deathkernel(depart): 
    """
    Parameters
    ----------
    depart : ARRAY
        Array that contains the coordinates of the point configuration present
        just before the next death (in the order abscissa, ordinate, type, regime).

    Returns
    -------
    The same array as in the input but without the coordinates of the dead center.
    """
    r=stats.uniform.rvs(0,1)
    totder=1.15*(0.18*process.nrb(depart)+0.22*process.nrsp(depart)+0.07*process.nrsb(depart))
    totdel=1.15*(0.15*process.nlb(depart)+0.12*process.nlsp(depart)+0.058*process.nlsb(depart))
    if r<(totder/(totder+totdel)):
        R=1
        r2=stats.uniform.rvs(0,1)
        if r2< ((1.15*0.18*process.nrb(depart))/totder):
            P=process.nrb(depart)
            typechoisi=1
        elif r2< ((1.15*(0.18*process.nrb(depart)+0.22*process.nrsp(depart)))/totder):
            P=process.nrsp(depart)
            typechoisi=2
        else :
            P=process.nrsb(depart)
            typechoisi=3
            #print('ici')
        i=stats.randint.rvs(1,P+1)
        #print(P)
        #print(i)
        m=0
        s=0
        while s<i and m<len(depart)/4:
            if depart[4*m+2]== 1 and process.typee(depart[4*m+3])==typechoisi:
                s+=1
            m+=1
    else :
        R=0
        r2=stats.uniform.rvs(0,1)
        if r2< ((1.15*0.15*process.nlb(depart))/totdel):
            P=process.nlb(depart)
            typechoisi=1
        elif r2< ((1.15*(0.15*process.nlb(depart)+0.12*process.nlsp(depart)))/totdel):
            P=process.nlsp(depart)
            typechoisi=2
        else :
            P=process.nlsb(depart)
            typechoisi=3
            #print('icii')
        i=stats.randint.rvs(1,P+1)
        #print(P)
        #print(i)
        m=0
        s=0
        while s<i and m<len(depart)/4:
            if depart[4*m+2]== 0 and process.typee(depart[4*m+3])==typechoisi:
                s+=1
            m+=1
    C=depart[4*(m-1)+3]
    tab=len(depart)*[True]
    tab[4*(m-1)]=False
    tab[4*(m-1)+1]=False
    tab[4*(m-1)+2]=False
    tab[4*(m-1)+3]=False
    return(depart[tab],R,C,m-1)


"""transition kernel"""
def transitionkernel(depart): 
    """
    Parameters
    ----------
    depart : ARRAY
        Array that contains the coordinates of the point configuration present
        just before the next transition (in the order abscissa, ordinate, type, regime).

    Returns
    -------
    The same array as in the input but with the new coordinate "regime" of the point that has changed.
    """
    x=depart
    r=stats.uniform.rvs(0,1)
    if r<(1/16):
       P=process.nr(depart)
       i=stats.randint.rvs(1,P+1)
       m=0
       s=0
       while s<i and m<len(depart)/4:
           if depart[4*m+2]== 1 :
               s+=1
           m+=1 
    else:
        P=process.nl(depart)
        i=stats.randint.rvs(1,P+1)
        m=0
        s=0
        while s<i and m<len(depart)/4:
            if depart[4*m+2]== 0 :
                s+=1
            m+=1
    if depart[4*(m-1)+3]==1:
        nt=np.random.choice([2,3],p=[1/4,3/4])
    elif depart[4*(m-1)+3]==3:
        nt=1
    else:
        nt=np.random.choice([1,3],p=[1/2,1/2])
    if nt==2:
        nt=stats.uniform.rvs(0,2*np.pi)
    x[4*(m-1)+3]=nt
    return(x)



b=process.proctotal(T,60,d, Delta, generatesituation, move, beta, delta, tau, alpha, birthkernel, deathkernel, transitionkernel) 
(resfinal,TpsSauts,tabecarts,track,compteurs,tracknaissance, trackmort) =b 

(trajlang,tracktrajlang,trajmutantes,trajlangdec,tracktrajlangdec,couleurstrajlangdec, reslangtronque, restronque, tracktronque)=draw.trajlangtronque(tabecarts, resfinal, track)


#to draw Langerin trajectories
draw.dessintrajlang(trajlangdec,couleurstrajlangdec, xcontour, ycontour)    

#to draw the boxplot of the number of proteins per frame
draw.moyparframeboxplot(reslangtronque)

#to draw histograms of the length of trajectories
draw.longueurtraj(trajlang, trajmutantes)

#to draw the process as a video
#draw.draw(b,xcontour,ycontour)
