# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 12:22:20 2018

@author: julien
"""
#########################################################################
############################ Cibo code ##################################
############### Plot histogram of MAST result from : ####################
################# http://meme-suite.org/tools/mast ######################
#########################################################################

import xml.etree.ElementTree as ET
import numpy as np
from pylab import *
import matplotlib.patches as mpatches

tree = ET.parse('mast.xml')
root = tree.getroot()


h=0 #counter of hit
allmatch=[] #we put all macth in a list
for seq in root.findall('sequences'): #we search the sequences in root
    for a in seq.findall('sequence'): #we search the sequence in sequences
        #trouve = seq.find('hit').match
        #print (a.tag, a.attrib)
        name = a.get('name')#we get the name of the sequence
        #print(name)
        for b in a.findall('seg'): #we search the Seg
            for c in b.findall('hit'): #and in the seg the hits
                h+=1
                ID = c.get('idx') #we get the ID of the hit
                pos= c.get('pos') #and its posistion
                List=[] #we creat a list to store it
                List.append(name) #and fill it with the data
                List.append(h)
                List.append(ID)
                List.append(pos)
                #•print(name, h,ID,pos)
                #print(List)
                allmatch.append(List) #we creat a list of the characteristic of the match
#print(allmatch)  
allmotif=[] #we creat a list to put all motif
for seq in root.findall('motifs'): #we search the motifs in root
    for a in seq.findall('motif'): #in motifs we search the motif
        List2=[] #we creat a list to fill it with motif characteristic
        #print (a.tag, a.attrib)
        name = a.get('alt') #we get the name
        lenmot= a.get('length') #the lenth
        Shema= a.get('id') #and the id wich is te sghamma ex GGAG for SD
        List2.append(name)
        List2.append(Shema) #ad all in the list
        List2.append(lenmot)
        #print(name)
        allmotif.append(List2) #and put the list in the list of motifs
        
#print(allmotif)
matresult=np.zeros((len(allmotif),100)) #we creat an array for all motif in lenth 100 to count the number of happening, old version without histo, only used to have the average position 
for n in range(0, len(allmatch)): #parcour all match
    #print(allmatch[n])
    #print(allmatch[n][2])
    #print(n)
    for i in range(0,len(allmotif)): #parcour all motifs
        Numotif=int(allmatch[n][2]) #put the number of the motif of the hit in int
        #print(Numotif)
        if i==Numotif: #if motif of match = motif selected then
            locend=int(allmatch[n][3])+int(allmotif[i][2])-2 #calculate the location of thee end of the motif
            #♣print(locend)
            matresult[i,locend]=matresult[i,locend]+1 #ad one to the case coresponding
            
Listresult=[] #creat a list of list of the localisation of each motif on the sequecence, to creat histogram
for i in range(0,len(allmotif)): #parcour all motif
    listresmot=[] #creat a list for the localisation of this given motif
    for n in range(0, len(allmatch)): #parcour all hits
        Numotif=int(allmatch[n][2]) #put in int the number of the motif
        #print(Numotif)
        if i==Numotif: #if motif = selected one
            listresmot.append(int(allmatch[n][3])+int(allmotif[i][2])) #add the localisation of the end of the motif( loc start + lenth) to the list
    #print(listresmot)
    Listresult.append(listresmot) #add the formed list to the list of result
    
matmoy=np.zeros(len(allmotif))   #calculate the average localisation         
for i in range(0,len(allmotif)):
    n=0
    for t in range(0,100):
        if matresult[i,t]>1: #filter random noise
            n+=matresult[i,t]
            matmoy[i]+=matresult[i,t]*t
    matmoy[i]/=n
    matmoy[i]+=-100

#num_bins = 50
for plop in range(0,len(allmotif)): #for all motif print histo of same figure
    print(plop,allmotif[plop][1],"loc:",matmoy[plop])
    num_bins = int((max(Listresult[plop])-min(Listresult[plop]))/2.5)
    n, bins, patches = plt.hist(Listresult[plop], num_bins, alpha=0.5,label=allmotif[plop][1])
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, borderaxespad=0.)
plt.xlabel('Localisation on the sequences') 
plt.ylabel('Number of time we encounter a motif (bp)')
show()
# Et Voila !!!!
"""
 ▲
▲ ▲
"""