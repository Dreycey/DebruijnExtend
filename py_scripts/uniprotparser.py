# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 17:40:54 2019

@author: Dreycey
"""

#### Making sql file for CRYSTALSTRUCTURESEQUENCE
'''
thisFile = 'pdb_data_seq.csv'

#LINUX
path="//mnt/c/Users/Dreycey/Desktop/building_a_database/protein-data-set"
filePath = path + '/' + thisFile


#WINDOWS
#path = "C:\\Users\\Dreycey\\Desktop\\building_a_database\\protein-data-set"
#filePath = path + '\\' + thisFile


#structureId	classification	experimentalTechnique	macromoleculeType	residueCount	resolution	structureMolecularWeight	crystallizationMethod	crystallizationTempK	densityMatthews	densityPercentSol	pdbxDetails	phValue	publicationYear


openedFile = open(filePath)
allLines = openedFile.readlines()

structureId_list = []
chainid_list = []
sequence_list = []
residueCount_list = []
macroMoleculeType_list = []


for line in allLines:
    
    structureId = line.split(',')[0]  
    chainid = line.split(',')[1]
    sequence = line.split(',')[2]
    residueCount = line.split(',')[3]
    macroMoleculeType = line.split(',')[4]

    
    #DEBUGGING
    #if uniprot == 'C5MK56': 
        #print (pdb_names, pdb_names, uniprot)
        #None
    
    # appending data to attributes
    if structureId == '':
        structureId_list.append('NULL')
    else:
        structureId_list.append(structureId)
    
    if chainid == '':
        chainid_list.append('NULL')
    else: 
        chainid_list.append(chainid)
    
    if sequence == '':
        sequence_list.append('NULL')
    else:
        sequence_list.append(sequence)
    
    if residueCount == '':
        residueCount_list.append('NULL')
    else:
        residueCount_list.append(residueCount)
    
    if macroMoleculeType == '':
        macroMoleculeType_list.append('NULL')
    else:
        macroMoleculeType_list.append(macroMoleculeType)


for entry_index in range( len(structureId_list) ):
    structureId_name = structureId_list[entry_index] 
    chainid_name = chainid_list[entry_index]
    sequence_name = sequence_list[entry_index]
    residueCount_name = residueCount_list[entry_index]
    macroMoleculeType_name = macroMoleculeType_list[entry_index]
    
    print (
        "'" 
        + structureId_name + 
        "'" + 
        "," + 
        "'" 
        + chainid_name + 
        "'" + 
        "," + 
        "'" 
        + sequence_name + 
        "'" + 
        "," + 
        residueCount_name + 
        "," + 
        "'" +
        macroMoleculeType_name[:-2] + 
        "'" )
'''

'''
#### Making sql file for CRYSTALSTRUCTURE

thisFile = 'pdb_data_no_dups_2.csv'

#LINUX
path="/Users/da39/Desktop/CLASS_WORK/building_a_database/protein-data-set"
filePath = path + '/' + thisFile


#WINDOWS
#path = "C:\\Users\\Dreycey\\Desktop\\building_a_database\\protein-data-set"
#filePath = path + '\\' + thisFile


#structureId	classification	experimentalTechnique	macromoleculeType	residueCount	resolution	structureMolecularWeight	crystallizationMethod	crystallizationTempK	densityMatthews	densityPercentSol	pdbxDetails	phValue	publicationYear


openedFile = open(filePath)
allLines = openedFile.readlines()

structureId_list = []
classification_list = []
experimentalTechnique_list = []
macromoleculeType_list = []
residueCount_list = []
resolution_list = []
structureMolecularWeight_list = []
crystallizationMethod_list = []
crystallizationTempK_list = []
densityMatthews_list = []
densityPercentSol_list = []
pdbxDetails_list = []
phValue_list = []
publicationYear_list = []


pdb_se = []

for line in allLines:
    
    structureId = line.split(',')[0]  
    classification_id = line.split(',')[1]
    experimentalTechnique = line.split(',')[2]
    macromoleculeType = line.split(',')[3]
    residueCount = line.split(',')[4]
    resolution = line.split(',')[5]
    structureMolecularWeight = line.split(',')[6]
    crystallizationMethod = line.split(',')[7]
    crystallizationTempK = line.split(',')[8]
    densityMatthews = line.split(',')[9]
    densityPercentSol = line.split(',')[10]
    pdbxDetails = line.split(',')[11]
    phValue = line.split(',')[12]
    publicationYear = line.split(',')[13][:-2]
    
    #DEBUGGING
    #if uniprot == 'C5MK56': 
        #print (pdb_names, pdb_names, uniprot)
        #None
        
    #print ('Hello_1')
    #print (pdb_se)
    
    if structureId not in pdb_se:
        #print ('hello')
        pdb_se.append(structureId)
    else: 
        #print ('AHA')
        continue
      
        
    
    # appending data to attributes
    if structureId == '':
        structureId_list.append('NULL')
    else:
        structureId_list.append(structureId)
    
    if classification_id == '':
        classification_list.append('NULL')
    else: 
        classification_list.append(classification_id)
    
    if experimentalTechnique == '':
        experimentalTechnique_list.append('NULL')
    else:
        experimentalTechnique_list.append(experimentalTechnique)
        
    if macromoleculeType == '':
        macromoleculeType_list.append('NULL')
    else:
        macromoleculeType_list.append(macromoleculeType)
    
    if residueCount == '':
        residueCount_list.append(-1)
    else:
        residueCount_list.append(residueCount)
    
    if resolution == '':
        resolution_list.append(-1)
    else:
        resolution_list.append(resolution)
        
    if structureMolecularWeight == '':
        structureMolecularWeight_list.append(-1)
    else:
        structureMolecularWeight_list.append(structureMolecularWeight)
    
    if crystallizationMethod == '':
        crystallizationMethod_list.append('NULL')
    else:
        crystallizationMethod_list.append(crystallizationMethod)
        
    if crystallizationTempK == '' or crystallizationTempK == 'NULL':
        crystallizationTempK_list.append(-1)
    else:
        crystallizationTempK_list.append(crystallizationTempK)
        
    if densityMatthews == '':
        densityMatthews_list.append(-1)
    else:
        densityMatthews_list.append(densityMatthews)
    
    if densityPercentSol == '':
        densityPercentSol_list.append(-1)
    else:
        densityPercentSol_list.append(densityPercentSol)
        
    if pdbxDetails == '':
        pdbxDetails_list.append('NULL')
    else:
        pdbxDetails_list.append(pdbxDetails)
    
    if phValue == '':
        phValue_list.append(-1)
    else:
        phValue_list.append(phValue)
    
    #print(len(publicationYear))
    if len(publicationYear) == 4:
        publicationYear_list.append(int(publicationYear))
    else: 
        publicationYear_list.append(-1)
        



for entry_index in range( len(structureId_list) ):
    
    structureId_name = structureId_list[entry_index] 
    classification_name = classification_list[entry_index]
    experimentalTechnique_name = experimentalTechnique_list[entry_index]
    macromoleculeType_name = macromoleculeType_list[entry_index]
    residueCount_name = residueCount_list[entry_index]
    resolution_name = resolution_list[entry_index]
    structureMolecularWeight_name = structureMolecularWeight_list[entry_index]
    crystallizationMethod_name = crystallizationMethod_list[entry_index]
    crystallizationTempK_name = crystallizationTempK_list[entry_index]
    densityMatthews_name = densityMatthews_list[entry_index]
    densityPercentSol_name = densityPercentSol_list[entry_index]
    pdbxDetails_name = pdbxDetails_list[entry_index]
    phValue_name = phValue_list[entry_index]
    publicationYear_name = publicationYear_list[entry_index]
    
    print ( 
        "'" 
        + structureId_name + 
        "'" + 
        "," + 
        "'" 
        + classification_name + 
        "'" + 
        "," + 
        "'" 
        + experimentalTechnique_name + 
        "'" + 
        "," + 
                "'" 
        + macromoleculeType_name + 
        "'" + 
        "," + 
        str(residueCount_name) + 
        "," + 
        str(resolution_name) + 
        "," + 
        str(structureMolecularWeight_name) + 
        "," + 
        "'" +
        crystallizationMethod_name + 
        "'" + 
        "," + 
        str(crystallizationTempK_name) + 
        "," + 
        str(densityMatthews_name) + 
        "," + 
        str(densityPercentSol_name) + 
        "," + 
        "'" +
        pdbxDetails_name + 
        "'" + 
        "," + 
        str(phValue_name) + 
        "," + 
        str(publicationYear_name) )

print(len(structureId_list))
'''


#### Making sql file for SECONDARYSTRUCTURE
'''
thisFile = '2018-06-06-ss.cleaned.csv'


#LINUX
path="//mnt/c/Users/Dreycey/Desktop/building_a_database/protein-secondary-structure"
filePath = path + '/' + thisFile


#WINDOWS
#path = "C:\\Users\\Dreycey\\Desktop\\building_a_database\\protein-secondary-structure"
#filePath = path + '\\' + thisFile


openedFile = open(filePath)
allLines = openedFile.readlines()

length_list = []
pdb_id_list = []
sequence_list = []
sst8_list = []
sst3_list = []
chaincode_list = []


for line in allLines:
    
    length = line.split(',')[5]  
    pdb_id = line.split(',')[0]
    sequence = line.split(',')[2]
    sst8 = line.split(',')[3]
    sst3 = line.split(',')[4]
    chaincode = line.split(',')[1]
    
    #DEBUGGING
    #if uniprot == 'C5MK56': 
        #print (pdb_names, pdb_names, uniprot)
        #None
    
    # appending data to attributes
    length_list.append(length)
    pdb_id_list.append(pdb_id)
    sequence_list.append(sequence)
    sst8_list.append(sst8)
    sst3_list.append(sst3)
    chaincode_list.append(chaincode)



for entry_index in range( len(pdb_id_list) ):
    length_name = length_list[entry_index] 
    pdb_id_name = pdb_id_list[entry_index]
    sequence_name = sequence_list[entry_index]
    sst8_name = sst8_list[entry_index]
    sst3_name = sst3_list[entry_index]
    chaincode_name = chaincode_list[entry_index]
    
    print ( 
         length_name +
        "," + 
        "'" 
        + pdb_id_name + 
        "'" + 
        "," + 
        "'" 
        + sequence_name + 
        "'" + 
        "," + 
                "'" 
        + sst8_name + 
        "'" + 
        "," + 
                "'" 
        + sst3_name + 
        "'" + 
        "," + 
                "'" 
        + chaincode_name + 
        "'" )

'''

#### Making sql file for UNIPROT
'''
thisFile = 'uniprot_COMP533_2.csv'


#LINUX
path="/Users/da39/Desktop/CLASS_WORK/building_a_database/curated_csvs"
filePath = path + '/' + thisFile


#WINDOWS
#path = "C:\\Users\\Dreycey\\Desktop\\building_a_database"
#filePath = path + '\\' + thisFile


openedFile = open(filePath)
allLines = openedFile.readlines()

uniprot_list = []
refseq_list = []
fullname_list = []
proteinName_list = []
geneName_list = []
goTerms_list = []
goprocesses_list = []
mass_list = []
length_list = []
sequence_list = []
location_list = []


for line in allLines:
      
    uniprot = line.split(',')[0].strip('"')
    refseq = line.split(',')[1].strip('"')
    fullname = line.split(',')[2].strip('"')
    proteinName = line.split(',')[3].strip('"')
    geneName = line.split(',')[4].strip('"')
    goTerms = line.split(',')[5].strip('"')
    goprocesses = line.split(',')[6].strip('"')
    mass = line.split(',')[7].strip('"')
    length = line.split(',')[8].strip('"')
    sequence = line.split(',')[9].strip('"')
    location = line.split(',')[10].strip('"')[:-2]
    
    #DEBUGGING
    if uniprot == 'C5MK56': 
        #print (pdb_names, pdb_names, uniprot)
        None
    
    # appending data to attributes
    if uniprot == '':
        uniprot_list.append('NULL')
    elif len(uniprot) > 50:
        uniprot_list.append(uniprot[:12])
        #print('GOTCHYA')
    else:
        uniprot_list.append(uniprot)

    if refseq == '':
        refseq_list.append('NULL')
    else:
        refseq_list.append(refseq)
 
    if fullname == '':
        fullname_list.append('NULL')
    else:
        fullname_list.append(fullname)

    if proteinName == '':
        proteinName_list.append('NULL')
    else:
        proteinName_list.append(proteinName)

    if geneName == '':
        geneName_list.append('NULL')
    else:
        geneName_list.append(geneName)

    if goTerms == '':
        goTerms_list.append('NULL')
    else:
        goTerms_list.append(goTerms)

    if goprocesses == '':
        goprocesses_list.append('NULL')
    else:
        goprocesses_list.append(goprocesses)
    
    if mass == '':
        mass_list.append(-1)
    else:
        mass_list.append(float(mass))

    if length == '':
        length_list.append(-1)
    else:
        length_list.append(float(length))
    
    if sequence == '' or len(sequence) < 6:
        sequence_list.append('NULL')
    else:
        sequence_list.append(sequence)

    if location == '': #or len(location) < 6:
        location_list.append('NULL')
    else:
        location_list.append(location)


for entry_index in range(len(uniprot_list)):
    uniprot_name = uniprot_list[entry_index] 
    refseq_name = refseq_list[entry_index]
    fullname_name = fullname_list[entry_index]
    proteinName_name = proteinName_list[entry_index]
    geneName_name = geneName_list[entry_index]
    goTerms_name = goTerms_list[entry_index]
    goprocesses_name = goprocesses_list[entry_index]
    mass_name = mass_list[entry_index]
    length_name = length_list[entry_index]
    sequence_name = sequence_list[entry_index]
    location_name = location_list[entry_index]

    
    print ("'" 
        + uniprot_name + 
        "'" + 
        "," + 
        "'" 
        + refseq_name + 
        "'" + 
        "," + 
        "'" 
        + fullname_name + 
        "'" + 
        "," + 
        "'" 
        + proteinName_name + 
        "'" + 
        "," + 
        "'" 
        + geneName_name + 
        "'" + 
        "," + 
        "'" 
        + goTerms_name.strip('"') + 
        "'" + 
        "," + 
        "'" 
        + goprocesses_name.strip('"') + 
        "'" + 
        "," 
        + str(mass_name) + 
        ","
        + str(length_name) + 
        "," + 
        "'" 
        + sequence_name + 
        "'" )

#uniprot = line.split(',')
#print (uniprot)
    #+ 
        #"," + 
        #"'" 
        #+ location_name + 
        #"'"
'''



#### Making sql file for hasproteins
    
'''
thisFile = 'mycsv_2.csv'


#LINUX
path="//mnt/c/Users/Dreycey/Desktop/building_a_database"
filePath = path + '/' + thisFile


#WINDOWS
#path = "C:\\Users\\Dreycey\\Desktop\\building_a_database"
#filePath = path + '\\' + thisFile


openedFile = open(filePath)
allLines = openedFile.readlines()



pdb_list = []
uniprot_list = []

for line in allLines:
    
    pdb_names = line.split(',')[1]
    pdb_names = line.strip('\n')
    uniprot = line.split(',')[0]
    
    #DEBUGGING
    if uniprot == 'C5MK56': 
        #print (pdb_names, pdb_names, uniprot)
        None
        
        
    #line.split('\t')
    pdb_list.append(pdb_names)
    uniprot_list.append(uniprot)
    
for pdb_index in range(0,len(pdb_list)):
    pdb_lis = pdb_list[pdb_index].split(',')
    pdb_li = []
    for i in pdb_lis:
        if ";" in i:
            spliti = i.split(';') # split pdb names
            for j in spliti:
                if len(i) > 2:        
                    pdb_li.append(j)
    #pdb_li = pdb_list[pdb_index].split(';')
    
    #IF DEBUGGING
    if uniprot_list[pdb_index] == 'C5MK56':
        #print("pdbliss: " + str(pdb_lis))
        #print ("pdb_list:" + pdb_list[pdb_index])
        #print ("pdb_lis: " +str(pdb_li))
        None
        
    #Used for printing all the values for each uniprot-pdb pair match
    for pdb in range(len(pdb_li) - 1):
        if len(pdb_li[pdb]) < 5 and len(pdb_li[pdb]) > 3:
            #print (uniprot_list[pdb_index], pdb)
            print ( "'"+ uniprot_list[pdb_index] +"'"+ "," +"'"+ pdb_li[pdb] + "'")
            
            #IF DEBUGGING
            if uniprot_list[pdb_index] == 'C5MK56': #'O19707':
                #print ("The pdb is: " + str(pdb_li[pdb]))
                None
'''


#print (pdb_list)


'''
with open(filePath) as tsv:
    for column in zip(*[line for line in csv.reader(tsv, dialect="excel-tab")]):
'''
