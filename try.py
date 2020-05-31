"""
Informatie Ensemble gene ID van website: http://www.ensembl.org/index.html
Voornamelijk phenotype
"""
import requests, sys
import pandas as pd
import numpy as np
import urllib.request
import gzip
from collections import Counter

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

def creatingEnsembleDf(df_accessionnumbers):
    """
    Input: dataframe van accessionnumbers met hun clusternummers
    Output: dataframe van alleen unigeneID met bijhorende cluster nummer
    """
    #selecteer alleen cluster nummer en unigeneID kolommen uit df_accessionnumbers
    df_ensemble = df_accessionnumbers[['cluster_nummer','ensemblegeneid']].copy()
    
    #haal alle rijen uit df waar er geen ensemblegeneid is
    df_ensemble.replace('', np.nan, inplace=True)
    df_ensemble.dropna(inplace=True)
    
    return df_ensemble

def phenotypeEnsemble(df_ensemble):
    # phenotype + gene family
    phenotypeDic = {}
    for ensemblegeneid in df_ensemble['ensemblegeneid']:
        try:
            server = "https://rest.ensembl.org"
            ext = "/phenotype/gene/mus_musculus/" + ensemblegeneid + "?"
                    
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
            
            if not r.ok:
                r.raise_for_status()
                sys.exit()
            
            decoded = r.json() 
            key = ensemblegeneid
            value = [d['description'] for d in decoded]
            phenotypeDic.update({key: value})
        except requests.HTTPError:
            break
    return phenotypeDic

def phenotypePerCluster(df_ensemble, phenotypeDic):
    #maak lege dictionary aan voor de uiteindelijke keys en values eraan toe te voegen
    descriptionClusterDic = {}
    
    #maak dictionary aan met key: clusternummer en value: df met alle rijen met dat clusternummer
    dicDfPerCluster = dict(tuple(df_ensemble.groupby('cluster_nummer')))
    
    #lijst met alle keys van de dictionary wat de voorkomende cluster nummers zijn
    listClusterNr = list(dicDfPerCluster.keys())
    alleDes = []
    
    #voor elke cluster de voorkomende chromosoom nummers in een lijst gesorteerd op voorkomens
    for clusterNr in listClusterNr:
        key = clusterNr
        valueList = []
        df_per_cluster = dicDfPerCluster.get(clusterNr)             #krijg df van die cluster
        key = clusterNr                                             #stel key gelijk aan het cluster nummer
        ensemblegeneIDs = df_per_cluster['ensemblegeneid'].tolist()           #alle unigeneIDs van die cluster in een lijst
        ##list comprehension voor alle unigeneIDs in die cluster hun chromosoom nummer te krijgen uit chromosomeDic
        for ensemblegeneid in ensemblegeneIDs:
            if ensemblegeneid in phenotypeDic:
                alleDes += phenotypeDic.get(ensemblegeneid)
                valueList += phenotypeDic.get(ensemblegeneid)
                value = [i[0] for i in Counter(valueList).most_common()]    #alle chromosoom nummers in die cluster gesorteerd op voorkomens
                descriptionClusterDic.update({key: value})                   #voeg cluster nummer met bovenstaande gevormde lijst toe aan dictionary
    return descriptionClusterDic, alleDes

def checkIfDuplicates_3(listOfElems):
    ''' Check if given list contains any duplicates '''    
    for elem in listOfElems:
        if listOfElems.count(elem) > 1:
            return True
    return False

def phenotypeVoorkomen(df_ensemble, descriptionClusterDic, alleDes):
    voorkomende_des = list(dict.fromkeys(alleDes))
    voorkomen_des_cluster_dic = {}
    df_ensemble['cluster_nummer']=pd.to_numeric(df_ensemble['cluster_nummer'])
    maxClusterNr = df_ensemble['cluster_nummer'].max()
    for des in voorkomende_des:  #for loop met alle woorden in de beschrijvingen
        voorkomen_des_cluster_lijst = []    #lege lijst waar aantal voorkomen van woorden per cluster in komen
        for index_cluster in range(1,maxClusterNr+1):  #for loop met index van cluster
            telling = descriptionClusterDic.get(str(index_cluster), []).count(des)    #aantal voorkomen van woord in die cluster
            voorkomen_des_cluster_lijst.append(telling)  #voorkomen van woord in cluster wordt toegevoegd aan de lijst
        voorkomen_des_cluster_dic[des] = voorkomen_des_cluster_lijst  #uiteindelijke dictionary wordt gemaakt
        
    return voorkomen_des_cluster_dic, maxClusterNr

def phenotypeClusterCheck(voorkomen_des_cluster_dic, maxClusterNr):
    meerdereClusterVoorkomens = {}
    for item in voorkomen_des_cluster_dic.items():
        key = item[0]
        lijst = item[1]
        if lijst.count(0)<maxClusterNr-2:
            value = [i+1 for i,x in enumerate(lijst) if x != 0]
            meerdereClusterVoorkomens.update({key: value})
    
    return meerdereClusterVoorkomens
        

clusterresultdict = resultLadenCluster('clusterresult_K_means.txt')
df_accessionnumbers = readAccession('accessionnumbers.txt', clusterresultdict)
df_ensemble = creatingEnsembleDf(df_accessionnumbers)
phenotypeDic = phenotypeEnsemble(df_ensemble)
descriptionClusterDic, alleDes = phenotypePerCluster(df_ensemble, phenotypeDic)
result = checkIfDuplicates_3(alleDes)
if result:
    print('Yes, list contains duplicates')
else:
    print('No duplicates found in list')   
voorkomen_des_cluster_dic, maxClusterNr = phenotypeVoorkomen(df_ensemble, descriptionClusterDic, alleDes)
meerdereClusterVoorkomens = phenotypeClusterCheck(voorkomen_des_cluster_dic, maxClusterNr)
print(meerdereClusterVoorkomens)

