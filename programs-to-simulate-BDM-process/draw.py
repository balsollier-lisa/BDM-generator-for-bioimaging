#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: lisabalsollier
"""

"""You will find in this file all the necessary programs to plot 
the results of the simulations generated in the process.py file."""



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



#center of the cell
X0=np.array([97,141.2])




"""this program makes the simulation for the simple example in video form"""
def drawsimpleex(a):
    (resfinal,TpsSauts,tabecarts,track,compteurs,tracknaissance, trackmort) =a
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
            
"""this program makes the simulation for the other examples in video form"""
def draw(a, xcontour,ycontour):
    [resfinal, TpsSauts,tabecarts,compteurs]=a
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
                ax.set_title(str(cr)+' Rab11 et '+str(cl) +' Langerin pr??sents')
            elif i>0 :
                (la,ca)=resfinal[i-1].shape
                if c==ca :
                    ax.set_title(str(cr)+' Rab11 et '+str(cl) +' Langerin pr??sents', color='red')
                else :
                    ax.set_title(str(cr)+' Rab11 et '+str(cl) +' Langerin pr??sents')
            #plt.pause(tabecarts[i][j])
            #ax.clear()


"""this program allows to sort the particles by trajectories"""
def trajlangtronque(tabecarts, resfinal, track):
    s=0
    i=1
    lieux=[]
    for k in range(len(tabecarts)) :
        for l in range(len(tabecarts[k])) :
            #print(s+tabecarts[k][l] )
            if s+tabecarts[k][l]>i*intertps and s<=(i)*intertps :
                i+=1
                lieux+=[[k,l]]
            s+=tabecarts[k][l]
            
    
    reslangtronque=[]
    tracklangtronque=[]
    for tab in lieux :
        choix=[]
        for x in resfinal[tab[0]][tab[1],2::4]==0 :
            choix+=[x]*4
        reslangtronque+=[resfinal[tab[0]][tab[1],:][np.array(choix)]]
        tracklangtronque+=[np.array(track[tab[0]])[resfinal[tab[0]][tab[1],2::4]==0]]
    
    restronque=[]
    tracktronque=[]
    for tab in lieux :
        restronque+=[resfinal[tab[0]][tab[1],:]]
        tracktronque+=[np.array(track[tab[0]])]
    
    tracklangtronque2=[]
    for i in range(len(tracklangtronque)):
        if len(tracklangtronque[i])!=0:
            tracklangtronque2+=[tracklangtronque[i]]
        
    nbtraj=np.max([np.max(tracklangtronque2[i]) for i in range(len(tracklangtronque2))])
    #print('nbtraj='+str(nbtraj))
    
    trajlang=[]
    tracktrajlang=[]
    for j in range(1,nbtraj+1):
        sstraj=[]
        for m in range(len(tracklangtronque)):
            for o in range(len(tracklangtronque[m])):
                if tracklangtronque[m][o]==j :
                    sstraj+=[reslangtronque[m][4*o], reslangtronque[m][4*o+1],reslangtronque[m][4*o+3]]
        if len(sstraj)!=0:
            trajlang+=[sstraj]
            tracktrajlang+=[j]
    trajmutantes=[]   
    trajlangdec=[]
    couleurstrajlangdec=[]
    tracktrajlangdec=[]
    for i in range(len(trajlang)):
        rep=0
        for j in range(1,len(trajlang[i])//3):
            if process.typee(trajlang[i][3*j+2])!=process.typee(trajlang[i][3*(j-1)+2]):
                trajlangdec+=[trajlang[i][3*rep:3*(j+1)+2]]
                trajmutantes+=[tracktrajlang[i]]
                tracktrajlangdec+=[tracktrajlang[i]]
                #print('i='+str(i))
                #print('j='+str(j))
                #print(str(trajlang[i][3*(j-1)+2])+' et '+str(trajlang[i][3*(j)+2]))
                couleurstrajlangdec+=[coul[process.typee(trajlang[i][3*(j-1)+2])-1]]
                rep=j
        trajlangdec+=[trajlang[i][3*(rep):]]
        tracktrajlangdec+=[tracktrajlang[i]]
        couleurstrajlangdec+=[coul[process.typee(trajlang[i][-1])-1]]  
    return(trajlang, tracktrajlang,trajmutantes,trajlangdec,tracktrajlangdec,couleurstrajlangdec, reslangtronque,restronque, tracktronque)


"""this program allows to transform superdiffusive particles too small into Brownian particles"""
def transformproctotal(resfinal,TpsSauts,tabecarts,track,compteurs,tracknaissance, trackmort, nbframe,trajlang, tracktrajlang,trajmutantes,trajlangdec,tracktrajlangdec,couleurstrajlangdec, reslangtronque,restronque, tracktronque) :
    [cptnaissances,cptmorts,cptnlb,cptnlsp,cptnlsb,cptnrb,cptnrsp,cptnrsb,cptmlb,cptmlsp,cptmlsb,cptmrb,cptmrsp,cptmrsb]=compteurs
    resfinal2=resfinal
    nbchgmt=0
    for i in range(len(trajlangdec)):
        if process.typee(trajlangdec[i][2])==2:
            if len(trajlangdec[i])//3<nbframe :
                numtrack=tracktrajlangdec[i]
                for k in range(len(track)):
                    for l in range(len(track[k])):
                        if track[k][l]==numtrack:
                            for m in range(len(resfinal[k])):
                                if process.typee(resfinal[k][m,4*l+3])==2 :
                                    #print([resfinal[k][m,4*l],resfinal[k][m,4*l+1]])
                                    resfinal2[k][m,4*l+3]=1
                nbchgmt+=1
                if numtrack in tracknaissance:
                    cptnlsp-=1
                    cptnlb+=1
                if numtrack in trackmort :
                    cptmlsp-=1
                    cptmlb+=1
    #print('nbchgmt='+str(nbchgmt))
    cpt=[cptnaissances,cptmorts,cptnlb,cptnlsp,cptnlsb,cptnrb,cptnrsp,cptnrsb,cptmlb,cptmlsp,cptmlsb,cptmrb,cptmrsp,cptmrsb]
    return(resfinal2, TpsSauts, tabecarts, track, cpt, nbchgmt)


"""this program draws the trajectories of the Langerin or Rab11 particles of the simulation"""
def dessintrajlang(trajlangdec,couleurstrajlangdec, xcontour, ycontour):
    fig,ax=plt.subplots()
    #return(xinf,yinf,xmax,ymax)
    ax.set_xlim([0,250])
    ax.set_ylim([0,283])
    #plt.plot(np.array(xcontour), np.array(ycontour), linewidth=3,color='black')
    #plt.xticks([], [])
    #plt.yticks([], [])
    plt.axis('off')
    plt.gca().set_aspect(aspect = 'equal')
    #plt.imshow(tourc)
    #plt.plot(trajlangdec[1][0::3],trajlangdec[1][1::3], color=couleurstrajlangdec[1], label='Brownian')
    #plt.plot(trajlangdec[0][0::3],trajlangdec[0][1::3], color=couleurstrajlangdec[0], label='subdiffusive')
    #plt.plot(trajlangdec[75][0::3],trajlangdec[75][1::3], color=couleurstrajlangdec[75], label='superdiffusive')
    #plt.legend(fontsize=13)
    plt.show()
    #print(res)
    #tps=TpsSauts+[T]
    #ecart=np.array([tps[i+1]-tps[i] for i in range(len(tps)-1)])
    #print(len(ecart))
    #print(ecart)
    for i in range(len(trajlangdec)):
        plt.plot(trajlangdec[i][0::3],trajlangdec[i][1::3], color=couleurstrajlangdec[i])
    ycontourmod=[]
    for i in range(len(ycontour)):   
        if ycontour[i]<X0[1]:
            ycontourmod+=[ycontour[i]-2]
        else :
            ycontourmod+=[ycontour[i]+2] 
    xcontourmod=[]
    for i in range(len(xcontour)):   
        if xcontour[i]<X0[0]:
            xcontourmod+=[xcontour[i]-2]
        else :
            xcontourmod+=[xcontour[i]+2] 
    plt.plot(np.array(xcontourmod), np.array(ycontourmod), color='black')
    plt.show()        
 
    
"""this program draws the boxplot of the average number of particles in each frame of the simulation"""
def moyparframeboxplot(reslangtronque):
    resfinallang=reslangtronque
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
    
    bpl1=ax.boxplot(Cl11, positions=[1], whis=2,widths = [0.5])
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians',  'caps']:
            plt.setp(bpl1[element], color='royalblue',markeredgecolor='royalblue', linewidth=2)
    
    
    bpl2=ax.boxplot(Cl22,positions=[1.8], whis=2,widths = [0.5])
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(bpl2[element], color='orangered',markeredgecolor='orangered', linewidth=2)
    
    bpl3=ax.boxplot(Cl33,positions=[2.6], whis=2,widths = [0.5])
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(bpl3[element], color='limegreen',markeredgecolor='limegreen', linewidth=2)
    ax.xaxis.set_tick_params(bottom=False)
    #ax.xaxis.set_ticks(range(1,4))
    ax.xaxis.set_ticklabels([ ' ', ' ', ' '], rotation = 0, color = 'black', fontsize = 8, style = 'italic', verticalalignment = 'top') 
    plt.show()




"""this program plots the number of particles of each of the three types against the frame number"""
def moyparframe(reslangtronque):
    resfinallang=reslangtronque
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

    
"""this program plots the histograms of the length of the simulation's trjactories according to the three types"""
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

   
    