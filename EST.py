import pandas as pd
import matplotlib.pyplot as plt

def loadData(GenDescription):
    """
    Python code en txt files moeten in deze map gelocaliseerd zijn.
    Laad data in dat nodig is voor fase 3 de beschrijving van de genen.
    Input: textbestand GenDescription
    Output: Dictionary
    """
    #lees txt bestanden
    infileGen = open(GenDescription)
    lines = infileGen.read().splitlines()
    infileGen.close()
    
    locsID=list(range(0,len(lines),3))  #list of keys positions
    IDsGen=[]
    for locID in locsID:
        #add ID's to list
        ID=lines[locID]
        IDsGen.append(ID)
    
    locsDescription=list(range(1,len(lines),3))    #list of values positions
    descriptionGen=[]
    for locDescription in locsDescription:
        #add descriptions to list
        description=lines[locDescription]
        descriptionGen.append(description)
    
    dictGenDes = {}
    for ID in IDsGen:
        dictGenDes[ID]=descriptionGen[IDsGen.index(ID)]
    
    return dictGenDes
loadData('GenDescription2.txt')

def resultladen(tekstbestand):
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
resultladen("Voorbeeld_clusterresult.txt")

def ESTs(df_clusterdata, dictGenDes):
    df_clusterdata['Description'] = df_clusterdata.CloneID.replace(dictGenDes)
    
    df_ESTs = df_clusterdata.loc[df_clusterdata['Description'] == 'ESTs']
    df_ESTs = df_ESTs[['CloneID', 'cluster_nummer', 'Description']]
    grouped_EST = df_ESTs.groupby(['cluster_nummer', 'Description']).size()
    grouped_EST = grouped_EST.to_frame() 
    grouped_EST = grouped_EST.rename(columns={0:'count'})
    grouped_EST = grouped_EST.reset_index()
    
    my_colors = list('rgbkymc')
    ax = grouped_EST.plot.bar(x='cluster_nummer', y='count', rot=0, color=my_colors)
    ax.set_xlabel('Cluster nummer')
    ax.set_ylabel('Aantal genen met een EST in het cluster')
    ax.get_legend().remove()
    plt.tight_layout()
