#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:14:24 2019

@author: da39
"""

#### Making sql file for refSeqTable
'''
thisFile = 'uniprot_COMP533_2.csv'


#LINUX
path="/Users/da39/Desktop/CLASS_WORK/building_a_database/curated_csvs"
filePath = path + '/' + thisFile

openedFile = open(filePath)
allLines = openedFile.readlines()

for line in allLines:
    #print (line)
    uniprot = line.split(',')[0]
    refseq = line.split(',')[1]
    fullname = line.split(',')[2]
    proteinName = line.split(',')[3]
    geneName = line.split(',')[4]
    goTerms = line.split(',')[5]
    goprocesses = line.split(',')[6]
    mass = line.split(',')[7]
    length = line.split(',')[8]
    sequence = line.split(',')[9]
    
    #print (uniprot)
    
    #print (refseq)
    
    refseq_list = refseq[1:-2].split(";")
    #print (refseq_list)
    
    for i in refseq_list:
        if i != 'NUL' and i != 'NULL':
            i = i.split('[')
            print (uniprot + ", ", i[0])
 '''           
            



## How to make go term table
'''
thisFile = 'uniprot_COMP533_2.csv'

#LINUX
path="/Users/da39/Desktop/CLASS_WORK/building_a_database/curated_csvs"
filePath = path + '/' + thisFile

openedFile = open(filePath)
allLines = openedFile.readlines()

for line in allLines:
    #print (line)
    uniprot = line.split(',')[0]
    refseq = line.split(',')[1]
    fullname = line.split(',')[2]
    proteinName = line.split(',')[3]
    geneName = line.split(',')[4]
    goTerms = line.split(',')[5]
    goprocesses = line.split(',')[6]
    mass = line.split(',')[7]
    length = line.split(',')[8]
    sequence = line.split(',')[9]
    
    #print (uniprot)
    
    #print (refseq)
    
    goTermsl = goTerms[1:-1].split(";")
    #print (refseq_list)
    
    concat= ''
    gotermtracker = []
    for i in goTermsl:
        if "[" not in i and len(concat) < 1:
            concat = i 
        elif "[" not in i and len(concat) >= 1:
            concat = concat + ' ' + i
        elif i != 'NUL' and i != 'NULL':
            i = concat + ' ' + i
            golistsplit = i.split('[') #.split(']')
            gotermm = golistsplit[1][:-2]
            if gotermm not in gotermtracker:
                print (uniprot + ", ",  "'" + gotermm + "'" + ", ", 
                       "'" + golistsplit[0][1:-1] + "'")
                gotermtracker.append(gotermm)
            concat= ''
'''

## How to make go processes

thisFile = 'uniprot_COMP533_2.csv'

#LINUX
path="/Users/da39/Desktop/CLASS_WORK/building_a_database/curated_csvs"
filePath = path + '/' + thisFile

openedFile = open(filePath)
allLines = openedFile.readlines()

for line in allLines:
    #print (line)
    uniprot = line.split(',')[0]
    refseq = line.split(',')[1]
    fullname = line.split(',')[2]
    proteinName = line.split(',')[3]
    geneName = line.split(',')[4]
    goTerms = line.split(',')[5]
    goprocesses = line.split(',')[6]
    mass = line.split(',')[7]
    length = line.split(',')[8]
    sequence = line.split(',')[9]
    
    #print (uniprot)
    
    #print (refseq)
    
    goprocessesl = goprocesses[1:-1].split(";")
    #print (refseq_list)
    
    concat= ''
    gotermtracker = []
    for i in goprocessesl:
    
        if "[" not in i and len(concat) < 1:
            concat = i 
        elif "[" not in i and len(concat) >= 1:
            concat = concat + ' ' + i
        elif i != 'NUL' and i != 'NULL':
            i = concat + ' ' + i
            golistsplit = i.split('[') #.split(']')
            gotermm = golistsplit[1][:-2]
            if gotermm not in gotermtracker:
                print (uniprot + ", ",  "'" + gotermm + "'" + ", ", 
                       "'" + golistsplit[0][1:-1] + "'")
                gotermtracker.append(gotermm)
            concat= ''

    

'''
## How to make go processes

thisFile = 'uniprot_COMP533_2.csv'

#LINUX
path="/Users/da39/Desktop/CLASS_WORK/building_a_database/curated_csvs"
filePath = path + '/' + thisFile

openedFile = open(filePath)
allLines = openedFile.readlines()

for line in allLines:
    #print (line)
    uniprot = line.split(',')[0]
    refseq = line.split(',')[1]
    fullname = line.split(',')[2]
    proteinName = line.split(',')[3]
    geneName = line.split(',')[4]
    goTerms = line.split(',')[5]
    goprocesses = line.split(',')[6]
    mass = line.split(',')[7]
    length = line.split(',')[8]
    sequence = line.split(',')[9]

    print (
    uniprot + 
    ", "  
    + fullname + 
    ", "  
    + proteinName + 
    ", "  
    + geneName + 
    ", "  
    + mass + 
    ", "
    + length + 
    ", "  
    + sequence[])

''' 
   