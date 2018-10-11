# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 10:12:25 2018

@author: julien
"""
#######################################################################
######################### Killy code ##################################
### Extract potential promoter sequences from intergenic sites ########
### To feed it to MEME: http://meme-suite.org/tools/meme ##############
#######################################################################

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import motifs
from Bio.SeqRecord import SeqRecord
import sys

file = sys.argv[1]#we read the name of the file given as argument
anotation=SeqIO.parse(file, "genbank") #we parse the file in a record named anotation
lenmax=100 #len maximunm of interested sequence 
lenin=24 #len min (joke) of interested sequence
for record in anotation:
    seq= record.seq
    N=0 #number of genes tested
    R=0 #number of gene strand +1 selected
    L=0 # number of gene strand -1 selected
    refus=0 #gene refused (to short)
    p=0 #  small accepted sequences
    Newrecords = [] #we prepare a list to put all our seqrecords
    for n in range(1,len(record.features)): #♥we parcour all the feature, not compting the first one, I don't remenber why

        if record.features[n].location.strand ==1 and record.features[n].type=="gene": #if gene sense 1 >> then
            loc=record.features[n].location.start #we take the begining of the gene
            #print(record.features[n].type)
            promlen= abs(loc-record.features[n-1].location.end) #and see the intergenic distance
            if promlen > lenin: #don't take it if too short
                if promlen >lenmax: #♦if biger than lenmax, put it at lenmax as that's what interest us
                    promlen= lenmax
                else: #we count the smal ones
                    p+=1
                seqpromhandle=seq[loc-(promlen+1):loc-1] #we put the sequence coresponding in an handle
                name="prom"
                name+=str(L+R)
                name+=str(record.features[n].qualifiers['locus_tag']) #♣we creat the name with the tag
                seqprom= SeqRecord(seqpromhandle,id=name,name=str(record.features[n].qualifiers['locus_tag']),description="")#we creat the seq record, with the sequence and the name, the rest is pointless as not writen in the end file, so no point....
                Newrecords.append(seqprom) #we add the seqprom to the list 
                R+=1 #we count the >>1 selected sequences
            else:
                refus+=1 #we count the refused ones
            
        if record.features[n].location.strand ==-1 and record.features[n].type=="gene": #if gene sense -1 << then
            loc=record.features[n].location.end #we take the end as sens -1
            promlen= abs(loc-record.features[n+1].location.start) #we calculate the intergenic distance
            #we check if intergenic distance not to small and (if not double reading or if double sense <200 nt) and if not after the end of genome
            if promlen > lenin and (record.features[n+1].location.strand ==-1 or promlen>200 ) and (record.features[n].location.end+100 < len(seq)):
                if promlen >lenmax: #☺if bigger than lenmax, we put at a lenmax
                    promlen= lenmax
                else:
                    p+=1 #☻count the smal ones
                seqpromhandle=seq[loc:loc+(promlen)] #we put the sequence in a handle
                seqpromhandle=seqpromhandle.reverse_complement() #we reverse complement it
                name="prom"
                name+=str(L+R)
                name+=str(record.features[n].location.strand) #we put -1 in the name as sens <<
                name+=str(record.features[n].qualifiers['locus_tag'])#we put the tag also
                seqprom= SeqRecord(seqpromhandle,id=name,name=str(record.features[n].qualifiers['locus_tag']),description="")#we creat the seq record, with the sequence and the name, the rest is pointless as not writen in the end file, so no point....
                Newrecords.append(seqprom) #we add it to teh list
                L+=1 #we count the << selected sequences
            else:
                refus+=1 #we count the refused ones
        if record.features[n].type=="gene":
            N+=1 #we count the genes 
    longname=0        
    for n in file: #check the lenth of the name of the file bifore .gb or .gbk
        if n != ".":
            longname +=1
        else: 
            break #I know dirty.....
    Nameout= file[0:longname]+"prom"+".fasta" #we creat the nam from the input minus the ".gb"+"prom" +".fasta"
    SeqIO.write(Newrecords ,Nameout, "fasta") #we write the file of sequenses
    
print("genes:",N,"1:",R,"-1:",L, "selected:",R+L, "refused: ",refus, "small:",p)#we print the charactéristic of the run, I love to do this.

"""
 ▲
▲ ▲
"""