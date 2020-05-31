"""
Code versie 1.2
Fase 3 
Groep 1 OGO genexpressie
"""
import pandas as pd

import matplotlib.pyplot as plt

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
    df_clusterdata = pd.DataFrame.from_records(clusterdatalist, columns = ['CloneID','1','2','4','7','14','21','45','90'])
        
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

def grafiekenClusters(df_clusterdata, clusterresultdict):
    #input: dataframe van clusterdata, dictionary van clusterresult
    #output: grafieken met (dag, relatieve expressiewaarden van die dag) voor alle genen in elke cluster
    
    #De bijbehorende cluster nummers zijn toegevoegd aan de CloneId's
    df_clusterdata['cluster_nummer'] = df_clusterdata.CloneID.replace(clusterresultdict)
    
    #lijst aanmaken met nummers van de gebruikte dagen voor x-as
    dagen = ['1','2','4','7','14','21','45','90']
    
    #lege lijst waarin de gemiddelde expressiewaardes inkomen per dag
    gemPerDag = []
    
    #voeg de gemiddeldes per dag  toe aan de lijst
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
    #input: dictionary met CodeID's en descriptions
    #output: dictionary met woord en hoe vaak het in de beschrijvingen voorkomt
    
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


            
def woordenClusters(dictGenDes,clusterresultdict):
    
    df_clusterdata.set_index('CloneID')
    
    #voeg cluster nummers toe aan gen descripties voor gelijke CloneID's
    df_clusterresult = pd.DataFrame(list(clusterresultdict.items()), columns = ['CloneID','cluster_nummer'])
    df_genDes = pd.DataFrame(list(dictGenDes.items()), columns = ['CloneID','beschrijving'])
    df_merge = pd.merge(df_genDes, df_clusterresult, how='inner', on='CloneID')
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
    print(alle_woord_per_cluster)
    for w in alle_woorden:  #for loop met alle woorden in de beschrijvingen
        voorkomen_woorden_cluster_lijst = []    #lege lijst waar aantal voorkomen van woorden per cluster in komen
        for index_cluster in range(df_merge['cluster_nummer'].max()):  #for loop met index van cluster
            telling = alle_woord_per_cluster[index_cluster].count(w)    #aantal voorkomen van woord in die cluster
            voorkomen_woorden_cluster_lijst.append(telling)  #voorkomen van woord in cluster wordt toegevoegd aan de lijst
        voorkomen_woorden_cluster_dict[w] = voorkomen_woorden_cluster_lijst  #uiteindelijke dictionary wordt gemaakt
            
        
    return voorkomen_woorden_cluster_dict
    
def woord_aantal_cluster(voorkomen_woorden_cluster_dict):
    
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

    woorden_1_cluster = {1:[],2:[],3:[],4:[],5:[]}
    for woord in voorkomen_woorden_cluster_dict:
        #als het aantal clusters waarin een woord voorkomt 1 is, voeg het woord dan toe aan de dict bij het juiste woord
        if len(voorkomen_woorden_cluster_dict[woord])==1:
            woorden_1_cluster[voorkomen_woorden_cluster_dict[woord][0]].append(woord)
    
    woorden_23_cluster = {1:[],2:[],3:[],4:[],5:[]}
    for woord in voorkomen_woorden_cluster_dict:
        #als het aantal clusters waarin een woord voorkomt 2 of 3 is, voeg het woord dan toe aan de dict bij het juiste woord
        if len(voorkomen_woorden_cluster_dict[woord])==2 or len(voorkomen_woorden_cluster_dict[woord])==3:
            for cluster in voorkomen_woorden_cluster_dict[woord]:
                woorden_23_cluster[cluster].append(woord)
            
    woorden_0_cluster = {1:[],2:[],3:[],4:[],5:[]}
    for woord in voorkomen_woorden_cluster_dict:
        #als het aantal clusters waarin een woord voorkomt 0 is, voeg het woord dan toe aan de dict bij het juiste woord
        if voorkomen_woorden_cluster_dict[woord]==[0]:
            woorden_0_cluster[voorkomen_woorden_cluster_dict[woord][0]].append(woord)

    woorden_all_cluster = []
    for woord in voorkomen_woorden_cluster_dict:
        #als het aantal clusters waarin een woord voorkomt 5 is, voeg het woord dan toe aan de dict bij het juiste woord
        if len(voorkomen_woorden_cluster_dict[woord])==5:
            woorden_all_cluster.append(woord)
        
    return woorden_1_cluster, woorden_23_cluster, woorden_0_cluster, woorden_all_cluster, voorkomen_woorden_cluster_dict

#deze functie is niet af
def lijkende_clusters(woorden_23_cluster,voorkomen_woorden_cluster_list,voorkomen_woorden_cluster_dict):
    #maak een lege dictionary aan, de key wordt een tuple van de twee clusters die hetzelfde woord bevatten, de value wordt hoevaak dat voorkomt
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
    #dicCloneIdFamily omzetten naar dataframe met CloneID en bijhorende familienummer
    df_clusterdata.set_index('CloneID')
    
    df_fam = pd.DataFrame(dicCloneIdFamily.items(), columns=['CloneID', 'Familienummer'])
    df_fam.drop(df_fam.index[0], inplace=True)
    
    #nieuwe dataframe bestaande uit samengevoegde df_fam en df_cluster met overeenkomende CloneID
    df_merge = pd.merge(df_fam, df_clusterdata, how='inner', on='CloneID')
    
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
        
        #geef plot weer
        plt.show()


def barPlotGenFamilies(df_clusterdata,clusterresultdict,dicCloneIdFamily):
    
    #De bijbehorende cluster nummers zijn toegevoegd aan de CloneId's
    df_clusterdata['cluster_nummer'] = df_clusterdata.CloneID.replace(clusterresultdict)

    # De familienummers toegevoegd aan de dataframe    
    df_clusterdata['Familienummer'] = df_clusterdata.CloneID.replace(dicCloneIdFamily)
    
    df_clusterdata['Familienummer'] = df_clusterdata['Familienummer'].astype(int)
    
    df_clusterdata = df_clusterdata[df_clusterdata.Familienummer < 27]
    
    
    # Dataframe maken met als groepen, familienummer en cluster nummer
    grouped_cluster_familie = df_clusterdata.groupby(['Familienummer', 'cluster_nummer']).size()
    grouped_cluster_familie = grouped_cluster_familie.to_frame()
    grouped_cluster_familie = grouped_cluster_familie.rename(columns={0:'count'})
    #print(grouped_cluster_familie)
    
    # plots maken van de verschillende families met aantal per cluster
    fig, ax = plt.subplots(figsize=(10,20))

    for group in grouped_cluster_familie:
        grouped_cluster_familie.unstack(level=0).plot(kind='bar', subplots=True, ax=ax)
        #ax.set_suptitle('Barplot van aantal per clusters van genfamilie: ' + df_namen.loc[df_namen.index[int(count)-1], 'Fam_naam'])
        ax.set_xlabel('Cluster nummer')
        ax.set_ylabel("Aantal CloneId's in een cluster")
        ax.set_xlim(0,6)
        ax.set_ylim(0,15)
        plt.tight_layout()

"""
alle data laden mbv. def
"""

df_clusterdata = dataLadenCluster("Voorbeeld_clusterdata.txt")
clusterresultdict = resultLadenCluster("Voorbeeld_clusterresult.txt")

dictGenDes2, dictGenDes = loadDataGen2and1('GenDescription2.txt','GenDescription.txt')

dicAccNumbersCloneId, dicCloneIdFamily = readFileAccAndCloneId('accessionnumbers.txt', 'CloneIdFamily.txt')

#grafiekenClusters(df_clusterdata,clusterresultdict)

voorkomen_woorden_cluster_dict = woordenClusters(dictGenDes2,clusterresultdict)

woorden_1_cluster, woorden_23_cluster, woorden_0_cluster, woorden_all_cluster, voorkomen_woorden_cluster_list = woord_aantal_cluster(voorkomen_woorden_cluster_dict)
lijkende_clusters(woorden_23_cluster,voorkomen_woorden_cluster_list,voorkomen_woorden_cluster_dict)
print(voorkomen_woorden_cluster_list)
#grafiekenGenFamilies(df_clusterdata, dicCloneIdFamily)

#barPlotGenFamilies(df_clusterdata, clusterresultdict, dicCloneIdFamily)
