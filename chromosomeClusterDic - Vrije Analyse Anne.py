"""
Vrij analyse expressie, chromosoom
"""
import pandas as pd
import numpy as np
import urllib.request
import gzip
from collections import Counter
import matplotlib.pyplot as plt

def readAccession(filename_accessionnumbers, clusterresultdict):
    """
    Input: accessionnumbers.txt en dictionary van clusterresult
    Output: dataframe van accessionnumbers met hun bijhorende cluster nummer
    """
    #open het accessionnumbers.txt bestand
    infile_acc = open(filename_accessionnumbers)
    lines = infile_acc.read().splitlines()
    infile_acc.close()
    
    #data omzetten naar een dataframe    
    accessionlist = []
    for line in lines:
        linesplit = line.split('\t')
        accessionlist.append(linesplit)
    
    df_accessionnumbers = pd.DataFrame.from_records(accessionlist, columns = ['cloneID','EMBL(genbank)ID','unigeneID','ensemblegeneid'])
    
    #cluster nummer bij elke cloneID toevoegen
    df_accessionnumbers['cluster_nummer'] = df_accessionnumbers['cloneID'].map(clusterresultdict)
    #NaN in cluster_nummer weghalen
    df_accessionnumbers = df_accessionnumbers[df_accessionnumbers['cluster_nummer'].notna()].reset_index()
    
    return df_accessionnumbers

def resultLadenCluster(tekstbestand):
    """
    Input: resultaten clusteren.txt
    Output: dictionary met CloneID en clusternummer
    """
    
    #open het tekstbestand, lees het regel voor regel in en sluit het
    f = open(tekstbestand, "r")
    clusterresultlines = f.read().splitlines()
    f.close()
    
    #maak een lege dictionary aan
    clusterresultdict = {}
    #iedere regel wordt gesplitst op 3 spaties, zodat ze als key:value aan de dictionary toegevoegd kunnen worden
    for i in clusterresultlines:
        linesplit= i.split("   ")
        clusterresultdict[linesplit[0]]=linesplit[1]
    
    return clusterresultdict

def dicChromosomeExpress(url):
    '''
    Input: url van data van unigene ID
    Output: dictionary met key unigene ID en het chromosoom nummer als value;
            dictionary met key unigene ID en waar het tot expressie komt als value (string formaat)
    '''
    #laad de url in van het bestand (.data.gz bestand)
    with urllib.request.urlopen(url) as response:
        
        #lees het gzip bestand (word beschouwd als binaire bestand met bytes)
        #vb.: b(ID  Mm.297)
        with gzip.GzipFile(fileobj=response) as binary_file:
            
            #maak lege dictionaries aan voor output
            chromosomeDic = {}
            expressDic = {}
            
            #lees het binaire bestand en split de lijnen
            data = binary_file.read().splitlines()
            
            #voor elke lijn in het bestand
            for line in data:
                
                #zet bytes om naar string
                line = line.decode()
                #split lijn bij elke white space
                line = line.split()
                
                #als de lijn niet leeg is dan:
                if line != []:
                    #als het eerste woord in de lijn 'ID' is
                    #stel de variabele key gelijk aan de ID nummer (unigene ID)
                    if line[0] == 'ID':
                        key = str(line[1:][0])
                    #als het eerste woord in de lijn 'EXPRESS' is
                    elif line[0] == 'EXPRESS':
                        #voeg  de unigene ID als key en waar het tot expressie kan komen als value aan expressDic
                        expressDic.update({key: line[1:]})
                    #als het eerste woord in de lijn 'CHROMOSOME' is
                    elif line[0] == 'CHROMOSOME':
                        #voeg de unigene ID als key en bijhorende chromosoom nummer als value toe aan chromosomeDic
                        chromosomeDic.update({key: line[1]})
                        
    return chromosomeDic, expressDic

def creatingUnigeneDf(df_accessionnumbers):
    """
    Input: dataframe van accessionnumbers met hun clusternummers
    Output: dataframe van alleen unigeneID met bijhorende cluster nummer
    """
    #selecteer alleen cluster nummer en unigeneID kolommen uit df_accessionnumbers
    df_unigene = df_accessionnumbers[['cluster_nummer','unigeneID']].copy()
    
    #haal alle rijen uit df waar er geen unigeneID is
    df_unigene.replace('', np.nan, inplace=True)
    df_unigene.dropna(inplace=True)
    
    return df_unigene

def chromosomeCluster(chromosomeDic, df_unigene):
    """
    Input: dictionary van chromsome nummers; dataframe van unigeneID en bijhorende cluster nummer
    Output: dictionary met key: cluster nummer;
                           value: lijst met chromosoom nummer gesorteerd op het voorkomen in die cluster
            dictionary met key: cluster nummmer;
                           value: lijst met voorkomende chromosoom nummers, zonder de duplicates eruit te hebben gefiltert
    """
    
    chromosomeClusterDic = {}
    chromosomeClusterDicFull = {}
    
    #maak dictionary aan met key: clusternummer en value: df met alle rijen met dat clusternummer
    dicDfPerCluster = dict(tuple(df_unigene.groupby('cluster_nummer')))
    
    #lijst met alle keys van de dictionary wat de voorkomende cluster nummers zijn
    listClusterNr = list(dicDfPerCluster.keys())
    
    #voor elke cluster de voorkomende chromosoom nummers in een lijst gesorteerd op voorkomens
    for clusterNr in listClusterNr:
        df_per_cluster = dicDfPerCluster.get(clusterNr)             #krijg df van die cluster
        key = clusterNr                                             
        unigeneIDs = df_per_cluster['unigeneID'].tolist()           
        ##list comprehension voor alle unigeneIDs in die cluster hun chromosoom nummer te krijgen uit chromosomeDic
        valueList = [chromosomeDic.get(key) for key in unigeneIDs if key in chromosomeDic]
        value = [i[0] for i in Counter(valueList).most_common()]    #alle chromosoom nummers in die cluster gesorteerd op voorkomens
        chromosomeClusterDic.update({key: value})                   
        chromosomeClusterDicFull.update({key: valueList})
        
    return chromosomeClusterDic, chromosomeClusterDicFull
       
def chromosomePlot(chromosomeClusterDicFull):
    """
    Maak barplot van chromosoom nummers per cluster.
    Input: dictionary met key: cluster nummmer;
                          value: lijst met voorkomende chromosoom nummers, zonder de duplicates eruit te hebben gefiltert
    Output: barplots per cluster met y-as het aantal voorkomens en x-as de chromosoom nummers in die cluster
    """ 
    for chromosoomNrsPerCluster in list(chromosomeClusterDicFull.values()):
        count = {}
        for chromosoomNr in chromosoomNrsPerCluster: 
            count[chromosoomNr] = count.get(chromosoomNr, 0) + 1
        keys = count.keys()
        values = count.values()
        plt.bar(keys, values)
    

#laad definities
clusterresultdict = resultLadenCluster('Voorbeeld_clusterresult.txt')
df_accessionnumbers = readAccession('accessionnumbers.txt', clusterresultdict)
chromosomeDic, expressDic = dicChromosomeExpress('https://ftp.ncbi.nlm.nih.gov/repository/UniGene/Mus_musculus/Mm.data.gz')
df_unigene = creatingUnigeneDf(df_accessionnumbers)
chromosomeClusterDic, chromosomeClusterDicFull = chromosomeCluster(chromosomeDic, df_unigene)
chromosomePlot(chromosomeClusterDicFull)