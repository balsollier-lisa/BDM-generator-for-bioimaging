#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 23:15:07 2023

@author: lisabalsollier
"""


import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pickle
import csv
import math

import process

rl=['v','.']
coul=['royalblue','orangered','limegreen']

intertps=0.14



with open('/Users/lisabalsollier/Documents/Thèse/données Lisa/maskc.pickle', 'rb') as handle:
    imgmod = pickle.load(handle)    


with open('/Users/lisabalsollier/Documents/Thèse/données Lisa/xcontour.pickle', 'rb') as handle:
    xcontour = pickle.load(handle)  

with open('/Users/lisabalsollier/Documents/Thèse/données Lisa/ycontour.pickle', 'rb') as handle:
    ycontour = pickle.load(handle)  


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



def drawsimpleex(a):
    (resfinal,TpsSauts,tabecarts,track,compteurs)=a
    fig,ax=plt.subplots()
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    plt.gca().set_aspect(aspect = 'equal')
    plt.show()
    for i in range(len(resfinal)):
        for j in range(len(resfinal[i])):
            ax.set_xlim([-1,1])
            ax.set_ylim([-1,1])
            plt.gca().set_aspect(aspect = 'equal')
            plt.scatter(resfinal[i][j,0::4],resfinal[i][j,1::4], marker='.')
            (l,c)=resfinal[i].shape
            nb=c//4
            ax.set_title('number of particles present : '+str(nb))
            plt.pause(tabecarts[i][j])
            ax.clear()
            


def draw(a):
    [resfinal, TpsSauts,tabecarts,track,compteurs]=a
    couleurs=[]
    for i in range(len(resfinal)):
        intercouleurs=[]
        for co in resfinal[i][0,3::4] :
            if co==1 or co==3 :
                intercouleurs+=[coul[int(co)-1]]
            else :
                intercouleurs+=['orangered']
        couleurs+=[intercouleurs]
    fig,ax=plt.subplots()
    ax.set_xlim([0,250])
    ax.set_ylim([0,283])
    plt.gca().set_aspect(aspect = 'equal')
    plt.plot(xcontour,ycontour,color='black')
    plt.axis('off')
    plt.show()
    for i in range(len(resfinal)):
        for j in range(len(resfinal[i])):
            ax.set_xlim([0,250])
            ax.set_ylim([0,283])
            plt.plot(xcontour,ycontour,color='black')
            plt.gca().set_aspect(aspect = 'equal')
            plt.axis('off')
            couleurs2=np.array(couleurs[i])
            plt.scatter(resfinal[i][j,0::4][resfinal[i][j,2::4]==0],resfinal[i][j,1::4][resfinal[i][j,2::4]==0], color=couleurs2[resfinal[i][j,2::4]==0], marker='.')
            plt.scatter(resfinal[i][j,0::4][resfinal[i][j,2::4]==1],resfinal[i][j,1::4][resfinal[i][j,2::4]==1], color=couleurs2[resfinal[i][j,2::4]==1], marker='v')
            cr=len(resfinal[i][j,0::4][resfinal[i][j,2::4]==1])
            cl=len(resfinal[i][j,0::4][resfinal[i][j,2::4]==0])
            (l,c)=resfinal[i].shape
            if i==0:
                ax.set_title(str(cr)+' Rab11 and '+str(cl) +' Langerin present')
            elif i>0 :
                (la,ca)=resfinal[i-1].shape
                if c==ca :
                    ax.set_title(str(cr)+' Rab11 and '+str(cl) +' Langerin present', color='red')
                else :
                    ax.set_title(str(cr)+' Rab11 and '+str(cl) +' Langerin present')
            plt.pause(tabecarts[i][j])
            ax.clear()




def trajlangtronqué(tabecarts, resfinal, track):
    s=0
    i=1
    lieux=[]
    for k in range(len(tabecarts)) :
        for l in range(len(tabecarts[k])) :
            if s+tabecarts[k][l]>i*intertps and s<=(i)*intertps :
                i+=1
                lieux+=[[k,l]]
            s+=tabecarts[k][l]
            
    
    reslangtronqué=[]
    tracklangtronqué=[]
    for tab in lieux :
        choix=[]
        for x in resfinal[tab[0]][tab[1],2::4]==0 :
            choix+=[x]*4
        reslangtronqué+=[resfinal[tab[0]][tab[1],:][np.array(choix)]]
        tracklangtronqué+=[np.array(track[tab[0]])[resfinal[tab[0]][tab[1],2::4]==0]]
    
    restronqué=[]
    tracktronqué=[]
    for tab in lieux :
        restronqué+=[resfinal[tab[0]][tab[1],:]]
        tracktronqué+=[np.array(track[tab[0]])]
    
    tracklangtronqué2=[]
    for i in range(len(tracklangtronqué)):
        if len(tracklangtronqué[i])!=0:
            tracklangtronqué2+=[tracklangtronqué[i]]
        
    nbtraj=np.max([np.max(tracklangtronqué2[i]) for i in range(len(tracklangtronqué2))])
    



    trajlang=[]
    tracktrajlang=[]
    for j in range(1,nbtraj+1):
        sstraj=[]
        for m in range(len(tracklangtronqué)):
            for o in range(len(tracklangtronqué[m])):
                if tracklangtronqué[m][o]==j :
                    sstraj+=[reslangtronqué[m][4*o], reslangtronqué[m][4*o+1],reslangtronqué[m][4*o+3]]
        if len(sstraj)!=0:
            trajlang+=[sstraj]
            tracktrajlang+=[j]
    trajmutantes=[]   
    trajlangdec=[]
    couleurstrajlangdec=[]
    for i in range(len(trajlang)):
        rep=0
        for j in range(1,len(trajlang[i])//3):
            if process.typee(trajlang[i][3*j+2])!=process.typee(trajlang[i][3*(j-1)+2]):
                trajlangdec+=[trajlang[i][3*rep:3*(j+1)+2]]
                trajmutantes+=[tracktrajlang[i]]
                couleurstrajlangdec+=[coul[process.typee(trajlang[i][3*(j-1)+2])-1]]
                rep=j
        trajlangdec+=[trajlang[i][3*(rep):]]
        couleurstrajlangdec+=[coul[process.typee(trajlang[i][-1])-1]]  
    return(trajlang, tracktrajlang,trajmutantes,trajlangdec,couleurstrajlangdec, reslangtronqué,restronqué, tracktronqué)




def dessintrajlang(trajlangdec,couleurstrajlangdec):
    fig,ax=plt.subplots()
    ax.set_xlim([0,250])
    ax.set_ylim([0,283])
    plt.plot(xcontour, ycontour, 'black')
    plt.axis('off')
    plt.gca().set_aspect(aspect = 'equal')
    plt.show()
    for i in range(len(trajlangdec)):
        plt.plot(trajlangdec[i][0::3],trajlangdec[i][1::3], color=couleurstrajlangdec[i])
    plt.show()
    


def moyparframe(reslangtronqué):
    resfinallang=reslangtronqué
    Cl11=[]
    Cl22=[]
    Cl33=[]
    for i in range(len(resfinallang)):
        cl1=0
        cl2=0
        cl3=0
        for j in range(len(resfinallang[i])//4):
            if resfinallang[i][4*j+3]==1:
                cl1+=1
            elif resfinallang[i][4*j+3]==3:   
                cl3+=1
            else:
                cl2+=1
        Cl11+=[cl1]
        Cl22+=[cl2]
        Cl33+=[cl3]

    fig, ax = plt.subplots()
    ax.plot([i for i in range(1,len(Cl11)+1)], Cl11, color='royalblue')
    ax.plot([i for i in range(1,len(Cl22)+1)], Cl22, color='orangered')
    ax.plot([i for i in range(1,len(Cl33)+1)], Cl33, color='limegreen')
    ax.plot([i for i in range(1,len(Cl11)+1)], np.array(Cl11)+np.array(Cl22)+np.array(Cl33), color='black')

    

def longueurtraj(trajlang, trajmutantes):
    Cr1b=[]
    Cr1sp=[]
    Cr1sb=[]
    for i in range(len(trajlang)):
        if i not in trajmutantes:
            if process.typee(trajlang[i][2])==1:
                if (len(trajlang[i])//3)<10 :
                    Cr1b+=[(len(trajlang[i])//3)+10] 
                else :
                    Cr1b+=[(len(trajlang[i])//3)] 
            if process.typee(trajlang[i][2])==2:
                if (len(trajlang[i])//3)<10 :
                    Cr1sp+=[(len(trajlang[i])//3)+10] 
                else :
                    Cr1sp+=[(len(trajlang[i])//3)] 
            if process.typee(trajlang[i][2])==3:
                if (len(trajlang[i])//3)<10 :
                    Cr1sb+=[(len(trajlang[i])//3)+10] 
                else :
                    Cr1sb+=[(len(trajlang[i])//3)] 

            
    fig, axs = plt.subplots(1,3)
    axs[0].hist(Cr1b, bins=30,range=(0,300), color='royalblue')
    axs[1].hist(Cr1sp, bins=30,range=(0,300),color='orangered')
    axs[2].hist(Cr1sb, bins=30,range=(0,300), color='limegreen')
    plt.show()
      
    
    
    
    
    
    
    