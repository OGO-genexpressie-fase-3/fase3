"""
Code versie 2.1
Fase 3 
Groep 1 OGO genexpressie
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import gzip
from collections import Counter
import operator
import requests, sys
import urllib.request
import seaborn as sns  # also improves the look of plots
sns.set()
plt.rcParams['figure.figsize'] = 10, 5  # default hor./vert. size of plots, in inches
plt.rcParams['lines.markeredgewidth'] = 1  # to fix issue with seaborn box plots; needed after import seaborn

def clusterDataLaden(tekstbestand):
    """
    Laad de data (clusterdata.txt) in in een dataframe
    Input: data clusteren .txt
    Output: dataframe met index CloneID
    """
    
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
    df_clusterdata = pd.DataFrame.from_records(clusterdatalist, columns = ['cloneID','1','2','4','7','14','21','45','90'])
        
    return df_clusterdata

def clusterResultLaden(tekstbestand):
    """
    Laad het tekstbestand met clusterresults in een dictionary
    Input: resultaten clusteren .txt
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
    
    #aantal clusters
    df_clusterresult = pd.DataFrame(clusterresultdict.items(), columns=['cloneID', 'cluster_nummer'])
    df_clusterresult_numeric = df_clusterresult.copy()
    df_clusterresult_numeric['cluster_nummer']=pd.to_numeric(df_clusterresult_numeric['cluster_nummer'])
    aantal_clusters = max(df_clusterresult_numeric['cluster_nummer'])+1
    
    return clusterresultdict, aantal_clusters, df_clusterresult

def genDescriptionLaden(GenDescription):
    """
    Laad de gen beschrijvingen in.
    Input: naam van gendescriptions txt bestand
    Output: Dictionary met beschrijving van de genen als values en de cloneID als key
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

def dataCloneIdFam(filenameCloneIdFam):
    """  
    Bestand openen en lezen en maakt een dictionairy met als keys de cloneID en als value het familienummer
    Input: cloneIDFamily.txt
    Output: dictionairy met cloneID en familienummer
    """
    
    #open het tekstbestand, lees het regel voor regel in en sluit het    
    infileCloneIdFamily = open(filenameCloneIdFam)
    cloneIdFamily = infileCloneIdFamily.readlines()
    infileCloneIdFamily.close()
    
    #dictionary maken met als key de cloneID en als value het familienummer
    dicCloneIdFamily = {}
    for cloneLine in cloneIdFamily[1:]:
        lines = cloneLine.split()
        dicCloneIdFamily[lines[0]] = int(lines[1])
    return dicCloneIdFamily

def readAccession(filename_accessionnumbers, clusterresultdict):
    """
    Laad de data van accesionnumbers in en zet daar bijhorende cluster nummer bij.
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
    #aparte df maken waarbij cluster nummers geen integer zijn
    df_accesionnumbers_nonnumeric = df_accessionnumbers.copy()
    #cluster nummers omzetten naar integers
    df_accessionnumbers['cluster_nummer']=pd.to_numeric(df_accessionnumbers['cluster_nummer'])
    
    return df_accessionnumbers, df_accesionnumbers_nonnumeric

def grafiekenClusters(df_clusterdata, clusterresultdict):
    """
    Maakt grafieken per cluster met hun relatieve expressiewaarden en de dagen.
    Input: dataframe van clusterdata, dictionary van clusterresult
    Output: grafieken met (dag, relatieve expressiewaarden van die dag) voor alle genen in elke cluster
    """
    #De bijbehorende cluster nummers zijn toegevoegd aan de CloneId's
    df_clusterdata.set_index('cloneID')
    df_clusterdata['cluster_nummer'] = df_clusterdata.cloneID.replace(clusterresultdict)
    
    #lijst aanmaken met nummers van de gebruikte dagen voor x-as
    dagen = ['1','2','4','7','14','21','45','90']
    
    #lege lijst waarin de gemiddelde expressiewaardes inkomen per dag
    gemPerDag = []
    
    #voeg de gemiddeldes per dag toe aan de lijst
    for dag in dagen:
        gemPerDag.append(df_clusterdata[dag].astype(float).mean())
    
    #voor elk cluster een apart dataframe
    df_cluster=[ df_clusterdata.loc[df_clusterdata['cluster_nummer']==val,:] for val in df_clusterdata['cluster_nummer'].unique() ]
    for df_per_cluster in df_cluster:
        df_per_cluster.loc[df_per_cluster.index[0],'cluster_nummer']
    
        # voor elk cluster een grafiek
        for i, row in df_per_cluster.iterrows():
            for i in range(8):
                dag = dagen[i]
                plt.scatter(dag, float(row[i+1]), c='purple')
                
        #plot dag gemiddeldes
        plt.scatter(dagen,gemPerDag, c='red')
        
        #grafiek en assen benoemen
        plt.suptitle('Grafiek van cluster ' + df_per_cluster.loc[df_per_cluster.index[0],'cluster_nummer'])
        plt.xlabel('Dagen')
        plt.ylabel('Relatieve expressiewaarde')
        plt.show()

def telWoorden(dict):
    """
    Telt voor elk woord hoevaak hij in totaal in alle genbeschrijvingen bij elkaar voorkomt
    Input: dictionary met CodeID's en descriptions (output van loadData)
    Output: dictionary met woord en hoe vaak het in de beschrijvingen voorkomt
    """
    
    #maak een lijst met de descriptions
    lijst = list(dict.values())
    #maak een lege dictionary, dit wordt de output
    woordendict = {}
    #voor elke beschrijving wordt een lijst van woorden gemaakt
    
    for beschrijving in lijst:
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
                
        #voor elk woord wordt 1 opgeteld bij het aantal keer dat het is voorgekomen in de woordendict   
        for woord in beschrijving:
            if woord in woordendict.keys():
                woordendict[woord] +=1
            else: #als het woord nog niet in woordendict staat, wordt het toegevoegd met value 1
                woordendict[woord] = 1
                
    return woordendict
            
def WoordenClusters(dictGenDes,clusterresultdict):
    """
    Kijkt per woord in de beschrijving in welke clusters het voorkomt
    Input: dictionary met CloneID: description dictGenDes (output van loadData), dictionary met CloneID: cluster, clusterresultdict (output van resultladen)
    Output: dictionary met key: woord en value: lijst met item voor elk cluster, bijv op index 1 een getal met hoevaak het woord voorkomt in cluster 1
    """
    
    #voeg cluster nummers toe aan gen descripties voor gelijke CloneID's
    df_clusterresult = pd.DataFrame(list(clusterresultdict.items()), columns = ['cloneID','cluster_nummer'])
    df_genDes = pd.DataFrame(list(dictGenDes.items()), columns = ['cloneID','beschrijving'])
    
    df_merge = pd.merge(df_genDes, df_clusterresult, how='inner', on='cloneID')
    #maak van de cluster nummers integers
    df_merge['cluster_nummer']=pd.to_numeric(df_merge['cluster_nummer'])
    
    
    #creëert een lijst met lijsten bestaande uit de woorden in de gen beschrijvingen per cluster
    alle_woord_per_cluster = [] #lege lijst voor bovenstaande benoemde lijst
    
    alle_woorden = [] #lege lijst waar alle woorden inkomen

    
    for cluster_nr in range(1,int(df_merge['cluster_nummer'].max())+1):  #for loop met clusternummer
        df_per_cluster = df_merge[df_merge['cluster_nummer']==cluster_nr]#dataframe per cluster met beschrijving
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

def woord_aantal_cluster(voorkomen_woorden_cluster_dict, aantal_clusters):
    """
    Geeft alle woorden die in 1 cluster voorkomen, in 2 of 3 clusters voorkomen, in 0 clusters voorkomen, in 4, of in geen clusters voorkomen
    input: output van WoordenClusters
    output: 4 lijsten voor verschillende hoeveelheid cluster, met de woorden die in zoveel clusters voorkomen
    """
    
    for woord in voorkomen_woorden_cluster_dict:
        #maak lege lijst aan, dit wordt de lijst met clusters waarin een woord voorkomt
        clusterlijst = []
        for cluster in voorkomen_woorden_cluster_dict[woord]:
            #voor elk cluster kan 1 of 0 staan, 1 betekent het woord zit erin
            if cluster >= 1:
                if voorkomen_woorden_cluster_dict[woord].index(cluster)+1 not in clusterlijst:
                    #voeg het nummer van de cluster toe aan de lijst met clusters
                    clusterlijst.append(voorkomen_woorden_cluster_dict[woord].index(cluster)+1)
        #als het woord in geen een cluster zit, zet dan een 0 neer
        if clusterlijst == []:
            clusterlijst= [0]
        #vervang de value van het woord met de lijst van clusters
        voorkomen_woorden_cluster_dict[woord]=clusterlijst

    woorden_1_cluster = {}
    for woord in voorkomen_woorden_cluster_dict:
        #als het aantal clusters waarin een woord voorkomt 1 is, voeg het woord dan toe aan de dict bij het juiste woord
        if len(voorkomen_woorden_cluster_dict[woord])==1:
            if voorkomen_woorden_cluster_dict[woord][0] not in woorden_1_cluster.keys():
                woorden_1_cluster[voorkomen_woorden_cluster_dict[woord][0]] = [woord]
            else:
                woorden_1_cluster[voorkomen_woorden_cluster_dict[woord][0]].append(woord)
    
    woorden_23_cluster = {}
    for woord in voorkomen_woorden_cluster_dict:
        #als het aantal clusters waarin een woord voorkomt 2 of 3 is, voeg het woord dan toe aan de dict bij het juiste woord
        if len(voorkomen_woorden_cluster_dict[woord])==2 or len(voorkomen_woorden_cluster_dict[woord])==3:
            for cluster in voorkomen_woorden_cluster_dict[woord]:
                if cluster not in woorden_23_cluster.keys():
                    woorden_23_cluster[cluster]=[woord]
                else:
                    woorden_23_cluster[cluster].append(woord)
            
    woorden_0_cluster = {}
    for woord in voorkomen_woorden_cluster_dict:
        #als het aantal clusters waarin een woord voorkomt 0 is, voeg het woord dan toe aan de dict bij het juiste woord
        if voorkomen_woorden_cluster_dict[woord]==[0]:
            if voorkomen_woorden_cluster_dict[woord][0] not in woorden_0_cluster.keys():
                woorden_0_cluster[voorkomen_woorden_cluster_dict[woord][0]] = [woord]
            else:
                woorden_0_cluster[voorkomen_woorden_cluster_dict[woord][0]].append(woord)

    woorden_all_cluster = []
    for woord in voorkomen_woorden_cluster_dict:
        #als het aantal clusters waarin een woord voorkomt 5 is, voeg het woord dan toe aan de dict bij het juiste woord
        if len(voorkomen_woorden_cluster_dict[woord])==aantal_clusters:
            woorden_all_cluster.append(woord)
        
    return woorden_1_cluster, woorden_23_cluster, woorden_0_cluster, woorden_all_cluster, voorkomen_woorden_cluster_dict

def lijkende_clusters(woorden_23_cluster,voorkomen_woorden_cluster_dict):
    """
    geeft in een dictionary weer hoe veel woorden overeenkomen tussen welke clusters
    input: woorden die in 2 of 3 clusters voorkomen, per cluster (woorden_23_cluster, output van woord_aantal_cluster)
    output: een dictionary met voor elk cluster een dictionary met de cluster waarmee ze gemeenschappelijke woorden hebben, en als value hoeveel woorden
    """
    
    gelijkenis = {}
    #maak een versie van voorkomen_woorden_cluster_dict waarin alleen de woorden staan die in 2 of 3 clusters voorkomen
    voorkomen_woorden_23_cluster={}
    for woord in voorkomen_woorden_cluster_dict:
        #als het woord de lengte van de lijst 2 of 3 is, voeg het toe aan voorkomen_woorden_23_cluster
        if len(voorkomen_woorden_cluster_dict[woord])==2 or len(voorkomen_woorden_cluster_dict[woord])==3:
            voorkomen_woorden_23_cluster[woord]=voorkomen_woorden_cluster_dict[woord]
            
            
    for woord1 in voorkomen_woorden_23_cluster:
        for cluster in voorkomen_woorden_23_cluster[woord1]:
            #per cluster dat in de lijst staat van een woord verwijder je het cluster
            voorkomen_woorden_23_cluster[woord1].remove(cluster)
            #van de clusters die overblijven doe je per cluster een koppel maken
            for cluster1 in voorkomen_woorden_23_cluster[woord1]:
                koppellijst = [cluster,cluster1]
                koppel = min(koppellijst),max(koppellijst)
                #bij dit koppel wordt 1 opgeteld als er een overeenkomstig woord is
                if koppel in gelijkenis.keys():
                    gelijkenis[koppel]+=1
                else:
                    gelijkenis[koppel]=1
    return gelijkenis

def grafiekenGenFamilies(df_clusterdata, dicCloneIdFamily):
    """
    Maakt grafieken aan per genfamilie met de relatieve expressiewaarden van de genen die voorkomen daarin.
    Input: df_clusterdata, dicCloneIdFamily
    Output: Grafieken per voorkomende genfamilie met dagen en relatieve expressiewaarden op de assen
    """
    #dicCloneIdFamily omzetten naar dataframe met CloneID en bijhorende familienummer
    df_clusterdata.set_index('cloneID')
    
    df_fam = pd.DataFrame(dicCloneIdFamily.items(), columns=['cloneID', 'Familienummer'])
    df_fam.drop(df_fam.index[0], inplace=True)
    
    #nieuwe dataframe bestaande uit samengevoegde df_fam en df_cluster met overeenkomende CloneID
    df_merge = pd.merge(df_fam, df_clusterdata, how='inner', on='cloneID')
    
    #lijst aanmaken met dataframes gegroepeerd volgens familienummer
    dagenFam = [df_merge.loc[df_merge['Familienummer']==nr,:] for nr in df_merge['Familienummer'].unique()]
    
    #lijst aanmaken met nummers van de gebruikte dagen voor x-as en gemiddelde
    dagen = ['1','2','4','7','14','21','45','90']
    
    #lege lijst waarin de gemiddelde expressiewaardes inkomen per dag
    gemPerDag = []
    
    #voeg de gemiddeldes per dag  toe aan de lijst
    for dag in dagen:
        gemPerDag.append(df_merge[dag].astype(float).mean())
        
    #dataframe aanmaken met de familienummer en bijhorende familie naam 
    namenFam = ['A Disintegrin-like and Metalloprotease with Thrombospondin Type 1 Motif Family','Adenylate Cyclase Family','Alpha Actinin Family','Annexin Family','Apolipoprotein C Family','Aquaporin Family','Calcium Channel, Voltage-dependent, Gamma Subunit Family','Calpain Family','Claudin Family','Collagen Family','Cyclic Nucleotide Phosphodiesterase Superfamily','E2f Transcription Factor Family','EH-domain Containing Family','Fatty Acid Coenzyme A Ligase Family','Gamma-Aminobutyric Acid Ionotropic Type A Receptor Family','Glutamate Receptor, Ionotropic Superfamily','Glycine Receptor, Alpha Family','Glypican Family','Iroquois Family','Leucine-rich Repeat LGI Family','Opsin Family','Phosphatidylinositol-4-phosphate 5-Kinase Family','Protocadherin Gamma Family', 'Secretory Carrier Membrane Protein Family','Synaptotagmin Family','Tripartite Motif Family']
    df_namen = pd.DataFrame(namenFam, columns=['Fam_naam'])
    
    #een grafiek voor elke voorkomende genfamilie bij experiment
    for fam in dagenFam:
        famNr = fam.loc[fam.index[0],'Familienummer']
        
        #(dag, relatieve expressiewaarden van die dag) plotten
        for i, row in fam.iterrows():
            for i in range(8):
                dag = dagen[i]
                plt.scatter(dag, float(row[i+1]), c='blue')
        
        #plot dag gemiddeldes
        plt.scatter(dagen,gemPerDag, c='red')
        
        #grafiek en assen benoemen
        plt.suptitle('Grafiek van genfamilie: ' + df_namen.loc[df_namen.index[int(famNr)-1], 'Fam_naam'])
        plt.xlabel('Dagen')
        plt.ylabel('Relatieve expressiewaarde')
        
        plt.show()


def barplotsGenFamilies(df_clusterdata, dicCloneIdFamily):
    """ Grafieken maken per gen familie en de clusters
    #input: gen families .txt, df_clusterdata
    #output: grafieken met (dag, relatieve expressiewaarden van die dag) voor elke gen familie die voorkomt bij het experiment
    """

    # De familienummers toegevoegd aan de bijbehorende cloneID in het dataframe   
    df_clusterdata.set_index('cloneID')
    df_clusterdata['Familienummer'] = df_clusterdata.cloneID.replace(dicCloneIdFamily)
    #dictNamen = {'0':'A Disintegrin-like and Metalloprotease with Thrombospondin Type 1 Motif Family' ,'1':'Adenylate Cyclase Family','2':'Alpha Actinin Family','3':'Annexin Family','4':'Apolipoprotein C Family','5':'Aquaporin Family','6':'Calcium Channel, Voltage-dependent, Gamma Subunit Family','7':'Calpain Family','8':'Claudin Family','9':'Collagen Family','10':'Cyclic Nucleotide Phosphodiesterase Superfamily','11':'E2f Transcription Factor Family','12':'EH-domain Containing Family','13':'Fatty Acid Coenzyme A Ligase Family','14':'Gamma-Aminobutyric Acid Ionotropic Type A Receptor Family','15':'Glutamate Receptor, Ionotropic Superfamily','16':'Glycine Receptor, Alpha Family','17':'Glypican Family','18':'Iroquois Family','19':'Leucine-rich Repeat LGI Family','20':'Opsin Family','21':'Phosphatidylinositol-4-phosphate 5-Kinase Family','22':'Protocadherin Gamma Family','23': 'Secretory Carrier Membrane Protein Family','24':'Synaptotagmin Family','25':'Tripartite Motif Family'}
    df_clusterdata = df_clusterdata[df_clusterdata.Familienummer < 27]
    
    #dataframe aanmaken met de familienummer en bijhorende familie naam 
    namenFam = ['A Disintegrin-like and Metalloprotease with Thrombospondin Type 1 Motif Family','Adenylate Cyclase Family','Alpha Actinin Family','Annexin Family','Apolipoprotein C Family','Aquaporin Family','Calcium Channel, Voltage-dependent, Gamma Subunit Family','Calpain Family','Claudin Family','Collagen Family','Cyclic Nucleotide Phosphodiesterase Superfamily','E2f Transcription Factor Family','EH-domain Containing Family','Fatty Acid Coenzyme A Ligase Family','Gamma-Aminobutyric Acid Ionotropic Type A Receptor Family','Glutamate Receptor, Ionotropic Superfamily','Glycine Receptor, Alpha Family','Glypican Family','Iroquois Family','Leucine-rich Repeat LGI Family','Opsin Family','Phosphatidylinositol-4-phosphate 5-Kinase Family','Protocadherin Gamma Family', 'Secretory Carrier Membrane Protein Family','Synaptotagmin Family','Tripartite Motif Family']
    df_namen = pd.DataFrame(namenFam, columns=['Fam_naam'])
    
    # Dataframe maken met als groepen: familienummer en cluster nummer
    grouped_cluster_familie = df_clusterdata.groupby(['Familienummer', 'cluster_nummer']).size()
    grouped_cluster_familie = grouped_cluster_familie.to_frame()
    grouped_cluster_familie = grouped_cluster_familie.rename(columns={0:'count'})
    
    # plots maken van de verschillende families met aantal per cluster
    fig, ax = plt.subplots(figsize=(10,20))

    for group in grouped_cluster_familie:
        grouped_cluster_familie.unstack(level=0).plot(kind='bar', subplots=True, ax=ax)
        #ax.set_suptitle('Grafiek van clusters per genfamilie: ' + df_namen.loc[df_namen.index[int(count)-1], 'Fam_naam'])
        ax.set_xlabel('Cluster nummer')
        ax.set_ylabel("Aantal CloneId's in een cluster")
        ax.set_xlim(0,6)
        ax.set_ylim(0,15)
        plt.tight_layout()

"""
Vrij analyse - EST (1)
"""
def ESTs(df_clusterresult, dictGenDes):
    """ Een bar plot maken met op de x-as de clusters en op de y-as het aantal genen
    met een EST in de beschijving van het gen. """
    # Een dataframe maken met waar de gen beschrijving EST bevast, de cloneIDs en de clusternummers per gen 
    df_clusterresult['Description'] = df_clusterresult.cloneID.replace(dictGenDes)
    df_ESTs = df_clusterresult.loc[df_clusterresult['Description'] == 'ESTs']
    df_ESTs = df_ESTs[['cloneID', 'cluster_nummer', 'Description']]
    
    # Tellen hoeveel genen met een EST beschijving er per cluster zijn
    grouped_EST = df_ESTs.groupby(['cluster_nummer', 'Description']).size()
    grouped_EST = grouped_EST.to_frame() 
    grouped_EST = grouped_EST.rename(columns={0:'count'})
    grouped_EST = grouped_EST.reset_index()
    
    # Een bar plot maken
    my_colors = list('rgbkymc')
    ax = grouped_EST.plot.bar(x='cluster_nummer', y='count', rot=0, color=my_colors)
    ax.set_title('Aantal genen met een EST in de beschrijving per cluster')
    ax.set_xlabel('Cluster nummer')
    ax.set_ylabel('Aantal genen')
    ax.get_legend().remove()
    plt.tight_layout()

"""
Vrij analyse - Tissuetype (2)
"""
def tissueType(df_accessionnumbers):
    """ function waar kan worden gezien hoeveel genen er uit een bepaald weefsel
    komen in de verschillende clusters
    input: df_accessionnumbers
    output: plots voor elk cluster met het aantal genen per soort weefsel
    """
    
    #selecteer alleen cluster nummer en unigeneID kolommen uit df_accessionnumbers
    df_cluster_EMBL = df_accessionnumbers[['cluster_nummer','EMBL(genbank)ID']].copy()
    
    #haal alle rijen uit df waar er geen unigeneID is
    df_cluster_EMBL.replace('', np.nan, inplace=True)
    df_cluster_EMBL.dropna(inplace=True)
    
    # Data inlezen van website met elk EMBL(genbank)ID
    website = "https://www.ebi.ac.uk/ena/browser/api/embl/"
    download = ".1?lineLimit=1000"
    tissue_type = {}
    for ID in df_cluster_EMBL['EMBL(genbank)ID']:
        url = website + ID + download
        rf = urllib.request.urlopen(url)
        data = rf.readlines()
        embldata = []
        for line in data:
            line = str(line)
            line = line[23:-4]
            line.split("\t")
            embldata.append(line)
            
            # Alle tissue types in een dictionairy zetten met CloneID
            if '/tissue_type=' in line:
                tissue_line = line[14:]
                tissue_type[ID]=tissue_line
    
    #df_cluster_EMBL['Tissue type'] = df_cluster_EMBL['EMBL(genbank)ID'].replace(tissue_type)
    df_cluster_EMBL['Tissue type'] = df_cluster_EMBL['EMBL(genbank)ID'].map(tissue_type)    
    df_cluster_EMBL = df_cluster_EMBL[df_cluster_EMBL['Tissue type'].notna()].reset_index()
    
    # Dataframe maken met als groepen, cluster nummer en tissue type
    grouped_cluster_tissue = df_cluster_EMBL.groupby(['cluster_nummer', 'Tissue type']).size()
    grouped_cluster_tissue = grouped_cluster_tissue.to_frame()
    grouped_cluster_tissue = grouped_cluster_tissue.rename(columns={0:'count'})
    
    # plots maken van de verschillende clusters met aantal tissue types
    fig, ax = plt.subplots(figsize=(10,30))

    for group in grouped_cluster_tissue:
        g = grouped_cluster_tissue.unstack(level=0).plot(kind='bar', subplots=True)
        g.xlabel('Cluster nummer')
        g.ylabel("Aantal CloneId's in een cluster per weefsel type")
        g.xlim(0,6)
        g.ylim(0,15)
        plt.tight_layout()
"""
Vrij Analyse - sequences (3)
"""

def embldata(ID):
    """
    laadt voor een EMBL(genbank)ID informatie uit de databank
    input: EMBL(genbank)ID
    output: list van regels van embldata (wat uit de database komt)
    """
    
    website = "https://www.ebi.ac.uk/ena/browser/api/embl/"
    path = ".1?lineLimit=1000"
    url = website + ID + path
    
    rf = urllib.request.urlopen(url)
    data = rf.readlines()
    
    embldata = []
    #maak van elke regel een string, gesplit per tab
    for line in data:
        line = str(line)
        line = line[2:-1]
        line.split("\t")
        embldata.append(line)
    
    return embldata

def sequence(embldata):
    """
    haalt uit de ruwe data van embldata de sequenceinformatie en de sequence zelf
    input: tekstbestand embldata (output van functie embldata)
    output: sequenceinfo: lijst met aantal baseparen, A, C, G, T, other als string
        sequencelist: lijst van de regels van de sequence
        sequence: een string van de sequence
    """
    
    #zoek waar in embldata de sequence zit
    for line in embldata:
        if "SQ" in line:
            sequenceindex = embldata.index(line)
    sequenceinfo = embldata[sequenceindex]
    sequencelist = embldata[sequenceindex+1:-1]
    sequence = ""
    #voeg elke regel van de sequence toe aan 'sequence' zonder het getal wat aangeeft op het hoeveelste basepaar je zit
    for line in sequencelist:
        sequence += line[:60]
    #pas de sequence zo aan dat er geen lange spaties inzitten
    while "  " in sequence:
        sequence = sequence.replace("  "," ")
        
        #pas de sequenceinfo aan zodat het een lijst wordt met het aantal baseparen, a, c,g,t en other, zonder de letters en andere omschrijvingen
    sequenceinfo = sequenceinfo.replace("\\n","")
    sequenceinfo = sequenceinfo.replace("SQ   Sequence ","")
    sequenceinfo = sequenceinfo.replace(" BP","")
    sequenceinfo = sequenceinfo.replace(" A","")
    sequenceinfo = sequenceinfo.replace(" C","")
    sequenceinfo = sequenceinfo.replace(" G","")
    sequenceinfo = sequenceinfo.replace(" T","")
    sequenceinfo = sequenceinfo.replace(" other","")
    sequenceinfo = sequenceinfo.replace(" ","")   
    sequenceinfo = sequenceinfo.split(';')
    sequenceinfo = sequenceinfo[:-1]

    return sequenceinfo, sequencelist, sequence

def df_accession_cluster_sequenceinfo(df_accessionnumbers, clusterresult):
    """
    voegt de sequence informatie en sequencetoe aan df_accessionnumbers
    input accessionnumbers dataframe (output accessioninladen)
    output accessionnumbers dataframe aangevuld met clusternummer, sequence en baseparen, a,c,g,t,other
    """

    #geef df andere naam
    df_merge = df_accessionnumbers.copy()
    
    bp = {}
    a = {}
    c = {}
    g = {}
    t = {}
    other = {}
    sequences = {}
    
    # Voor elk ID in de dataframe, probeer de embldata te verkrijgen, tenzij de error HTTPError komt (en het gen dus niet in de database staat)
    for ID in df_merge['EMBL(genbank)ID']:
        try:
            embldatas=embldata(ID)
            #roep de functie sequence aan op de data om sequenceinfo etc. te verkrijgen
            sequenceinfo, sequencelist, sequenced = sequence(embldatas)
            #voeg de info toe aan de respectievelijke dictionaries, zodat bijv. komt: a = {A642732 : 157}
            bp[ID]=sequenceinfo[0]
            a[ID]=sequenceinfo[1]
            c[ID]=sequenceinfo[2]
            g[ID]=sequenceinfo[3]
            t[ID]=sequenceinfo[4]
            other[ID]=sequenceinfo[5]
            sequences[ID]=sequenced
        except urllib.error.HTTPError:
            break

    #maak de kolommen baseparen (bp), a, c, g, t, other in een dataframe
    df_bp = pd.DataFrame(bp.items(), columns=['EMBL(genbank)ID','bp'])
    df_a = pd.DataFrame(a.items(), columns=['EMBL(genbank)ID','a'])
    df_c = pd.DataFrame(c.items(), columns=['EMBL(genbank)ID','c'])
    df_g = pd.DataFrame(g.items(), columns=['EMBL(genbank)ID','g'])
    df_t = pd.DataFrame(t.items(), columns=['EMBL(genbank)ID','t'])
    df_other = pd.DataFrame(other.items(), columns=['EMBL(genbank)ID','other'])
    
    #maak een dataframe van de sequences en voeg deze toe aan de dataframe
    df_sequence = pd.DataFrame(sequences.items(), columns = ['EMBL(genbank)ID','sequence'])
    
    #voeg alle gemaakte dataframes toe aan df_merge
    df_merge = pd.merge(df_merge, df_bp, how='inner', on='EMBL(genbank)ID')
    df_merge = pd.merge(df_merge, df_a, how='inner', on='EMBL(genbank)ID')
    df_merge = pd.merge(df_merge, df_c, how='inner', on='EMBL(genbank)ID')
    df_merge = pd.merge(df_merge, df_g, how='inner', on='EMBL(genbank)ID')
    df_merge = pd.merge(df_merge, df_t, how='inner', on='EMBL(genbank)ID')
    df_merge = pd.merge(df_merge, df_other, how='inner', on='EMBL(genbank)ID')
    df_merge = pd.merge(df_merge, df_sequence, how='inner', on='EMBL(genbank)ID')
    
    #zet kolommen om in numeric values
    df_merge['bp'] = pd.to_numeric(df_merge['bp'])
    df_merge['a'] = pd.to_numeric(df_merge['a'])
    df_merge['c'] = pd.to_numeric(df_merge['c'])
    df_merge['g'] = pd.to_numeric(df_merge['g'])
    df_merge['t'] = pd.to_numeric(df_merge['t'])
    df_merge['other'] = pd.to_numeric(df_merge['other'])
    
    return df_merge

def sequence_compare(df_merge,cluster):
    """
    Kijkt voor elk cluster wat per base de 'range' is, in hoe vaak het in de sequence voorkomt.
    Input: df_merge (de output uit de functie df_accession_cluster_sequenceinfo)
    Output:voor elke cluster een dictionary van de range van a,c,g,t etc.
    """
    
    #sequenceinfo = ['bp','a','c','g','t','other']
    range_cluster = {}
    #maak een dataframe aan van df_merge met maar 1 cluster
    cluster = df_merge[df_merge['cluster_nummer']==cluster]
    #bereken voor elke nucleotide de 'range' oftewel het verschil in maximale waarde en minimale waarde en voeg deze toe aan de dictionary
    range_cluster['bp']= (cluster['bp'].max()-cluster['bp'].min())
    range_cluster['a']= (cluster['a'].max()-cluster['a'].min())
    range_cluster['c']= (cluster['c'].max()-cluster['c'].min())
    range_cluster['g']= (cluster['g'].max()-cluster['g'].min())
    range_cluster['t']= (cluster['t'].max()-cluster['t'].min())

    return range_cluster


def plot_cluster(df_merge, base, clusternr):
    """
    Maakt een density plot
    Input: df_merge (output van functie df_accession_cluster_sequenceinfo), en de base die je wil plotten (bv 'a', 'bp')
    Output: density plot van de base van een cluster
    """
    
    #maak een dataframe van elk apart cluster
    df = df_merge[df_merge['cluster_nummer']==clusternr]
    
    #plot voor elk cluster de gewenste base
    plot = df[base].plot(kind='density')
    plot.set_xlabel('Aantal '+ base + ' in sequentie')
    plot.set_title('Densityplot voor aantal '+ ' in sequentie per cluster')

    return plot
    
def plot_base_per_cluster(df_merge, base, aantalclusters):
    """
    Plot voor alle clusters de density van de base.
    Input: df_merge, de gewenste base in een string, het aantal clusters als numeric value
    Output: densityplot van een base met elke cluster
    """
    
    #roep de functie plot_cluster af op elk cluster om de plot te krijgen
    for cluster in range(aantalclusters):
        plot_cluster(df_merge,base,cluster)


"""
Vrije Analyse - expressie (4)
"""
def creatingUnigeneDf(df_accessionnumbers):
    """
    Maakt dataframe aan met alleen unigeneIDs.
    Input: dataframe van accessionnumbers met hun clusternummers
    Output: dataframe van alleen unigeneID met bijhorende cluster nummer
    """
    #selecteer alleen cluster nummer en unigeneID kolommen uit df_accessionnumbers
    df_unigene = df_accessionnumbers[['cluster_nummer','unigeneID']].copy()
    
    #haal alle rijen uit df waar er geen unigeneID is
    df_unigene.replace('', np.nan, inplace=True)
    df_unigene.dropna(inplace=True)
    
    return df_unigene

def dicChromosomeExpress(url):
    '''
    Inladen van unigeneID data via internet in dictionaries.
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

def expressieClusters(expressDic,clusterresultdict, df_unigene):
    """
    Kijkt voor de verschillende plekken waar een gen tot expressie komt in welke clusters dat hoevaak voorkomt.
    input: expressDic,clusterresultdict, df_accessionnumbers
    output: dictionary
    """
    #voeg cluster nummers toe aan gen descripties voor gelijke CloneID's

    df_express = pd.DataFrame(list(expressDic.items()), columns = ['unigeneID','expressie']).set_index('unigeneID')
    df_merge = pd.merge(df_unigene, df_express, how='inner', on='unigeneID')

  
    
    #creëert een lijst met lijsten bestaande uit de woorden in de gen beschrijvingen per cluster
    alle_woord_per_cluster = [] #lege lijst voor bovenstaande benoemde lijst
    
    alle_woorden = [] #lege lijst waar alle woorden inkomen

    
    for cluster_nr in range(0,int(df_merge['cluster_nummer'].max())+1):  #for loop met clusternummer
        df_per_cluster = df_merge[df_merge['cluster_nummer']==cluster_nr]#dataframe per cluster met beschrijving
        exp_per_cluster = df_per_cluster['expressie']   #zet kolom met beschrijving om in lijst

        woord_per_cluster = []  #lege lijst voor alle woorden per cluster
        
        for woord in exp_per_cluster:
            
            woord_per_cluster.extend(woord)
            
            #voegt elk woord toe aan alle_woorden                
            alle_woorden.extend(woord)
            
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

"""
Vrije Analyse - chromosoom nummer (5)
"""

def chromosomeCluster(chromosomeDic, df_unigene):
    """
    De chromosoom nummers die voorkomen per cluster in een dictionary.
    Input: dictionary van chromosoom nummers; dataframe van unigeneID en bijhorende cluster nummer
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
    for key, lijst in chromosomeClusterDicFull.items():
        cntDic = Counter(lijst).most_common()
        pd.DataFrame({key: dict(cntDic)}).plot(kind='bar')
        plt.title('Chromosoom nummer barplot')
        plt.xlabel('Chromosoom nummer')
        plt.ylabel('Aantal voorkomens')
        plt.legend()
        plt.show()

"""
Vrije analyse - ziekte beschrijving (6)
"""

def creatingEnsembleDf(df_accessionnumbers):
    """
    Maakt een aparte dataframe aan met alleen de voorgekozen ID.
    Input: dataframe van accessionnumbers met hun clusternummers
    Output: dataframe van alleen ensemble gene ID met bijhorende cluster nummer
    """
    #selecteer alleen cluster nummer en unigeneID kolommen uit df_accessionnumbers
    df_ensemble = df_accessionnumbers[['cluster_nummer','ensemblegeneid']].copy()
    
    #haal alle rijen uit df waar er geen ensemblegeneid is
    df_ensemble.replace('', np.nan, inplace=True)
    df_ensemble.dropna(inplace=True)
    
    return df_ensemble

def phenotypeEnsemble(df_ensemble):
    """
    Haalt ziektebeschrijvingen met de ensemble gene ID van het internet.
    Input: dataframe van alleen ensemble gene ID met bijhorende cluster nummer
    Output: dictionary met key: ensemblegeneid; value: beschrijving van ziektes
            lijst met alle voorkomende ziektebeschrijvingen
    """
    phenotypeDic = {}
    alleDes = []
    #de ziektebeschrijvingen per ensemblegeneid van het internet halen
    #wanneer ensemblegeneid niet gevonden kan worden komt er geen error en probeerd het voor een nieuwe id
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
            alleDes += value
        except requests.HTTPError:
            break
    return phenotypeDic, alleDes

def phenotypePerCluster(df_ensemble, phenotypeDic):
    """
    Maakt dictionary aan dat de ziektebeschrijvingen sorteert per cluster.
    Input: dataframe van alleen ensemble gene ID met bijhorende cluster nummer
           dictionary met key: ensemblegeneid; value: beschrijving van ziektes
    Output: dictionary met key: cluster nummer; value: lijst met alle voorkomende ziektebeschrijvingen
    """
    descriptionClusterDic = {}
    
    #maak dictionary aan met key: clusternummer en value: df met alle rijen met dat clusternummer
    dicDfPerCluster = dict(tuple(df_ensemble.groupby('cluster_nummer')))
    
    #lijst met alle keys van de dictionary wat de voorkomende cluster nummers zijn
    listClusterNr = list(dicDfPerCluster.keys())
    
    #voor elke cluster de voorkomende ziektebeschrijvingen in een lijst gesorteerd op voorkomens
    for clusterNr in listClusterNr:
        key = clusterNr
        valueList = []
        df_per_cluster = dicDfPerCluster.get(clusterNr)             #krijg df van die cluster
        key = clusterNr
        ensemblegeneIDs = df_per_cluster['ensemblegeneid'].tolist()           #alle ensemblegeneids van die cluster in een lijst
        ##list comprehension voor alle ensemblegeneid in die cluster hun ziektebeschrijving te krijgen uit phenotypeDic
        for ensemblegeneid in ensemblegeneIDs:
            if ensemblegeneid in phenotypeDic:
                valueList += phenotypeDic.get(ensemblegeneid)
                value = [i[0] for i in Counter(valueList).most_common()]    #alle ziektebeschrijvingen in die cluster gesorteerd op voorkomens
                descriptionClusterDic.update({key: value})                   #voeg cluster nummer met bovenstaande gevormde lijst toe aan dictionary
    
    return descriptionClusterDic



def phenotypeVoorkomen(df_ensemble, descriptionClusterDic, alleDes, aantal_clusters):
    """
    Maakt dictionary dat de voorkomens van ziektebeschrijvingen beschrijft per cluster.
    Input: dataframe van alleen ensemble gene ID met bijhorende cluster nummer
           dictionary met key: cluster nummer; value: lijst met alle voorkomende ziektebeschrijvingen
           lijst met alle voorkomende ziektebeschrijvingen
    Output: dictionary met key: ziektebeschrijving; value: lijst met aantal voorkomens per cluster
    
    """
    #haal duplicates eruit
    voorkomende_des = list(dict.fromkeys(alleDes))
    voorkomen_des_cluster_dic = {}
    
    for des in voorkomende_des:  #for loop met alle woorden in de beschrijvingen
        voorkomen_des_cluster_lijst = []        
        #for index_cluster in range(aantal_clusters):  #for loop met index van cluster
            #telling = descriptionClusterDic.get(str(index_cluster)).count(des)    #aantal voorkomen van woord in die cluster
            #print(type(descriptionClusterDic.get(str(index_cluster))))
            #voorkomen_des_cluster_lijst.append(telling)
        voorkomen_des_cluster_dic[des] = voorkomen_des_cluster_lijst
        
    return voorkomen_des_cluster_dic

def phenotypeClusterCheck(voorkomen_des_cluster_dic, aantal_clusters):
    """
    Dictionaries van ziektebeschrijving die in minstens 3 cluster of maar in 1 cluster voorkomen.
    Input:
        dictionary met key: cluster nummer; value: lijst met alle voorkomende ziektebeschrijvingen
    Output:
        dictionary met key: ziektebeschrijving; value: lijst met alle voorkomende clusters (min. 3 cluster)
        dictionary met key: ziektebeschrijving; value: voorkomende cluster (1 cluster)
    """
    meerdereClusterVoorkomens = {}
    eenClusterVoorkomen = {}
    #for loop met key, value paar van de dictionary
    for item in voorkomen_des_cluster_dic.items():
        key = item[0]
        lijst = item[1]
        #als het in min. 3 clusters voorkomt
        if lijst.count(0)<aantal_clusters-2:
            #zet de voorkomens om in de clusternummer het voorkomt
            value = [i+1 for i,x in enumerate(lijst) if x != 0]
            meerdereClusterVoorkomens.update({key: value})
        #als het in 1 cluster voorkomt
        if lijst.count(0) == aantal_clusters:
             #zet de voorkomens om in de clusternummer het voorkomt
            value = [i+1 for i,x in enumerate(lijst) if x != 0]
            eenClusterVoorkomen.update({key: value})
            
    return meerdereClusterVoorkomens, eenClusterVoorkomen
        
def meestVoorkomendeZiekteBeschrijving(voorkomen_des_cluster_dic):
    """
    Welke ziektebeschrijvingen komen het meeste voor.
    Input: voorkomen_des_cluster_dic met key: beschrijving ziekte; value: voorkomens in lijst per cluster
    Output: meestVoorkomendeDes met key: aantal voorkomens (hoogste en hoogste-1); value: beschrijvingen die zoveel keer voorkomen (lijst)
    """
    #maak dictionary met key: beschrijving; value: aantal voorkomens in alle cluster
    voorkomendeDes = {}
    for des, voorkomens in voorkomen_des_cluster_dic.items():
        key = des
        value = sum(voorkomens)
        voorkomendeDes.update({key: value})
     #sorteer de dictionary op dalende value waarde
    voorkomendeDes = dict(sorted(voorkomendeDes.items(), key=operator.itemgetter(1),reverse=True))
    
    #eerste value uit dictionary halen, het hoogste aantal voorkomens
    eerste_value = next(iter(voorkomendeDes.values()))
    
    #maak een nieuwe dictionary aan waarin de hoogste aantal voorkomens in komen en de hoogste aantal voorkomens -1
    meestVoorkomendeDes = {}
    hoogsteVoorkomen = []
    opEenNaHoogsteVoorkomen = []
    
    for beschrijving, aantal in voorkomendeDes.items():
        if aantal == eerste_value:
            hoogsteVoorkomen.append(beschrijving)
        elif aantal == eerste_value-1:
            opEenNaHoogsteVoorkomen.append(beschrijving)
    meestVoorkomendeDes.update({eerste_value: hoogsteVoorkomen})
    meestVoorkomendeDes.update({eerste_value-1: opEenNaHoogsteVoorkomen})
    
    return meestVoorkomendeDes

"""
alle data laden mbv. def
"""
df_clusterdata = clusterDataLaden('Voorbeeld_clusterdata.txt')
clusterresultdict, aantal_clusters, df_clusterresult = clusterResultLaden('Voorbeeld_clusterresult.txt')
dictGenDes = genDescriptionLaden('GenDescription2.txt')
dicCloneIdFamily = dataCloneIdFam('CloneIdFamily.txt')
df_accessionnumbers, df_accesionnumbers_nonnumeric = readAccession('accessionnumbers.txt', clusterresultdict)

#grafiekenClusters(df_clusterdata, clusterresultdict)

woordendict = telWoorden(dictGenDes)
voorkomen_woorden_cluster_dict = WoordenClusters(dictGenDes,clusterresultdict)
woorden_1_cluster, woorden_23_cluster, woorden_0_cluster, woorden_all_cluster, voorkomen_woorden_cluster_dict = woord_aantal_cluster(voorkomen_woorden_cluster_dict, aantal_clusters)
gelijkenis = lijkende_clusters(woorden_23_cluster,voorkomen_woorden_cluster_dict)

#grafiekenGenFamilies(df_clusterdata, dicCloneIdFamily)
#barplotsGenFamilies(df_clusterdata,dicCloneIdFamily)

#Vrij Analyse
ESTs(df_clusterresult, dictGenDes)

#tissueType(df_accessionnumbers)
df_cluster = df_accession_cluster_sequenceinfo(df_accessionnumbers, clusterresultdict)
plot_base_per_cluster(df_cluster, 'a', aantal_clusters)
for cluster in range(aantal_clusters):
    sequence_compare(df_cluster, cluster)

df_unigene = creatingUnigeneDf(df_accessionnumbers)
chromosomeDic, expressDic = dicChromosomeExpress('https://ftp.ncbi.nlm.nih.gov/repository/UniGene/Mus_musculus/Mm.data.gz')
voorkomen_woorden_cluster_dict = expressieClusters(expressDic,clusterresultdict, df_unigene)

chromosomeClusterDic, chromosomeClusterDicFull = chromosomeCluster(chromosomeDic, df_unigene)
#chromosomePlot(chromosomeClusterDicFull)

df_ensemble = creatingEnsembleDf(df_accesionnumbers_nonnumeric)
phenotypeDic, alleDes = phenotypeEnsemble(df_ensemble)
descriptionClusterDic = phenotypePerCluster(df_ensemble, phenotypeDic)
voorkomen_des_cluster_dic = phenotypeVoorkomen(df_ensemble, descriptionClusterDic, alleDes, aantal_clusters)
meerdereClusterVoorkomens = phenotypeClusterCheck(voorkomen_des_cluster_dic, aantal_clusters)
meestVoorkomendeDes = meestVoorkomendeZiekteBeschrijving(voorkomen_des_cluster_dic)
"""
Mogelijke verbeterpunten:
    - inladen van definities op efficientere manier
    - Nienke haar code zoals tissueType, EST, barplots laten werken met dezelfde inlaad definities
    - Plots van tissueType en barplots van gen families verbeteren
    - Misschien per cluster de gevonden informatie laten zien (print)
"""