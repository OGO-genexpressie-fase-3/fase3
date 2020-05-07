import numpy as np
import pandas as pd

"""
Stap 2
"""

def stap2(tekstbestand, locGen2, locGen):
    """
    Python code en txt files moeten in deze map gelocaliseerd zijn.
    Laad data in dat nodig is voor de stap:
    
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
        clusterresultdict[linesplit[0]]=int(linesplit[1])
    
    #result is clusterresultdict
    
    """
    Laad verdere data in dat nodig is voor deze stap: 
    de beschrijving van de genen als values in dictionary.
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
    
    #result is dictGenDes2, dictGenDes
    
    """
    2.1: Maakt dictionary aan met key alle woorden in de descriptions en met value
    een lijst. Die lijst bestaat uit het aantal voorkomen van die woorden
    in de verschillende clusters, in de descriptions van de genen die 
    voorkomen in die clusters.
    """
    
    #voeg cluster nummers toe aan gen descripties voor gelijke CloneID's
    df_clusterresult = pd.DataFrame(list(clusterresultdict.items()), columns = ['CloneID','cluster_nummer'])
    df_genDes2 = pd.DataFrame(list(dictGenDes2.items()), columns = ['CloneID','beschrijving'])
    df_merge = pd.merge(df_genDes2, df_clusterresult, how='inner', on='CloneID')
    
    #creëert een lijst met lijsten bestaande uit de woorden in de gen beschrijvingen per cluster
    alle_woord_per_cluster = [] #lege lijst voor bovenstaande benoemde lijst
    
    alle_woorden = [] #lege lijst waar alle woorden inkomen
    
    for cluster_nr in range(1,df_merge['cluster_nummer'].max()+1):  #for loop met clusternummer
        df_per_cluster = df_merge[df_merge['cluster_nummer']==cluster_nr][['beschrijving']] #dataframe per cluster met beschrijving
        des_per_cluster = df_per_cluster['beschrijving'].tolist()   #zet kolom met beschrijving om in lijst
        
        woord_per_cluster = []  #lege lijst voor alle woorden per cluster
        
        for beschrijving in des_per_cluster:
            beschrijving = beschrijving.split(" ")
        
            #overbodige entries en tekens worden verwijderd
            if "" in beschrijving:
                beschrijving.remove("")
            for woord in beschrijving:
                beschrijving[beschrijving.index(woord)] = woord.strip("\(")
            for woord in beschrijving:
                beschrijving[beschrijving.index(woord)] = woord.strip(",")
            for woord in beschrijving:
                beschrijving[beschrijving.index(woord)] = woord.strip("\)")
            #beschrijving is de description van de gen gesplit in woorden in een lijst
            
            #lengt woord_per_cluster aan met nieuwe woorden van die cluster
            woord_per_cluster.extend(beschrijving)
            
            #voegt elk woord toe aan alle_woorden                
            alle_woorden.extend(beschrijving)
            
        #voegt lijst met woorden per cluster toe aan alle_woord_per_cluster
        alle_woord_per_cluster.append(woord_per_cluster)
    
    voorkomen_woorden_cluster_dict = {} #lege dictionary waar value de lijst met voorkomen woorden per cluster inkomen
    
    for w in alle_woorden:  #for loop met alle woorden in de beschrijvingen
        voorkomen_woorden_cluster_lijst = []    #lege lijst waar aantal voorkomen van woorden per cluster in komen
        for index_cluster in range(df_merge['cluster_nummer'].max()):  #for loop met index van cluster
            telling = alle_woord_per_cluster[index_cluster].count(w)    #aantal voorkomen van woord in die cluster
            voorkomen_woorden_cluster_lijst.append(telling)  #voorkomen van woord in cluster wordt toegevoegd aan de lijst
        voorkomen_woorden_cluster_dict[w] = voorkomen_woorden_cluster_lijst  #uiteindelijke dictionary wordt gemaakt
        
    return voorkomen_woorden_cluster_dict
    
    
    
dict = stap2('Voorbeeld_clusterresult.txt','GenDescription2.txt','GenDescription.txt')
print(dict)