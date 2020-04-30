"""
1. Inladen van data uit tekstbestanden en die verwerken zodat het bruikbaar is. Data die wordt ingeladen:
   'voorbeeld_clusterdata.txt','voorbeeld_clusterresults.txt','GenDescription2.txt','GenDescription.txt','accessionnumbers.txt', 'CloneIdFamily.txt'
"""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns  # also improves the look of plots
sns.set()
plt.rcParams['figure.figsize'] = 10, 5  # default hor./vert. size of plots, in inches
plt.rcParams['lines.markeredgewidth'] = 1  # to fix issue with seaborn box plots; needed after import seaborn

def dataLadenCluster(tekstbestand):
    #input: data clusteren .txt
    #output: dataframe met index CloneID
    
    #open het tekstbestand, lees het regel voor regel in en sluit het
    f = open(tekstbestand, "r")
    clusterdatalines = f.read().splitlines()
    f.close()

    #maak een lege lijst aan
    clusterdatalist = []
    #iedere regel wordt gesplitst op 2 spaties, waardoor elke regel een lijst van values wordt
    for i in clusterdatalines:
        linesplit = i.split("  ")
        #de lijst van de regel wordt toegevoegd aan de lege lijst
        clusterdatalist.append(linesplit)
        
    #clusterdatalist is een lijst met lijsten, hiervan wordt een dataframe gemaakt met CloneID als index
    df_clusterdata = pd.DataFrame.from_records(clusterdatalist, columns = ['CloneID','1','2','3','4','5','6','7','8']).set_index('CloneID')
        
    return df_clusterdata
    
def resultLadenCluster(tekstbestand):
    #input: resultaten clusteren .txt
    #output: dictionary met CloneID en clusternummer
    
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

def loadDataGen2and1(locGen2, locGen):
    """
    Python code en txt files moeten in deze map gelocaliseerd zijn.
    Laad data in dat nodig is voor fase 3: de beschrijving van de genen als values in dictionary.
    """
    #lees txt bestanden
    infileGen2 = open(locGen2)
    lines2 = infileGen2.read().splitlines()
    infileGen2.close()
    
    infileGen = open(locGen)
    lines = infileGen.read().splitlines()
    infileGen.close()
    
    #extract keys
    locsID2=list(range(0,len(lines2),3))  #list of keys positions
    IDsGen2=[]
    for locID2 in locsID2:
        #add ID's to list
        ID2=lines2[locID2]
        IDsGen2.append(ID2)
    
    locsID=list(range(0,len(lines),3))  #list of keys positions
    IDsGen=[]
    for locID in locsID:
        #add ID's to list
        ID=lines[locID]
        IDsGen.append(ID)
    
    #extract values
    locsDescription2=list(range(1,len(lines2),3))    #list of values positions
    descriptionGen2=[]
    for locDescription2 in locsDescription2:
        #add descriptions to list
        description2=lines2[locDescription2]
        descriptionGen2.append(description2)
    
    locsDescription=list(range(1,len(lines),3))    #list of values positions
    descriptionGen=[]
    for locDescription in locsDescription:
        #add descriptions to list
        description=lines[locDescription]
        descriptionGen.append(description)
    
    #make dictionary
    tuplePair2 = zip(IDsGen2, descriptionGen2) #voeg de lijsten bij elkaar toe als tuples
    dictGenDes2 = dict(tuplePair2)  #maak dictionary van de ID's en hun beschrijving
    
    tuplePair = zip(IDsGen, descriptionGen) #voeg de lijsten bij elkaar toe als tuples
    dictGenDes = dict(tuplePair)  #maak dictionary van de ID's en hun beschrijving
    
    return dictGenDes2, dictGenDes

def readFileAccAndCloneId(filenameAccNumbers, filenameCloneId):
    """ Inladen van bestanden en het verwerken"""
    
    #Inladen van bestand
    infileAccNumbers = open(filenameAccNumbers)
    AccNumbers = infileAccNumbers.readlines()
    infileAccNumbers.close()
    
    # Dictionary maken met als key de cloneId en als values de emblID, unigeneID en ensemblegeneID.
    dicAccNumbersCloneId = {}
    for line in AccNumbers:
        lines = line.replace('\n','')
        lines = line.split('\t')
        key = lines[0]
        value = {
            "emblID": lines[1],
            "unigeneID": lines[2],
            "ensemblegeneID": lines[3]}
        dicAccNumbersCloneId[key] = value
    
    #Inladen van bestand
    infileCloneIdFamily = open(filenameCloneId)
    cloneIdFamily = infileCloneIdFamily.readlines()
    infileCloneIdFamily.close()
    
    # Dictionary maken met als key de cloneID en als value het familienummer
    dicCloneIdFamily = {}
    for cloneLine in cloneIdFamily:
        lines = cloneLine.split()
        key = lines[0]
        value = lines[1]
        dicCloneIdFamily[key] = value
    
    return dicAccNumbersCloneId, dicCloneIdFamily
        

#alle data laden mbv. def
dataLadenCluster("Voorbeeld_clusterdata.txt")
resultLadenCluster("Voorbeeld_clusterresult.txt")

loadDataGen2and1('GenDescription2.txt','GenDescription.txt')

readFileAccAndCloneId('accessionnumbers.txt', 'CloneIdFamily.txt')
