import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import urllib.request
import seaborn as sns  # also improves the look of plots
sns.set()
plt.rcParams['figure.figsize'] = 10, 5  # default hor./vert. size of plots, in inches
plt.rcParams['lines.markeredgewidth'] = 1  # to fix issue with seaborn box plots; needed after import seaborn

def dataladen(tekstbestand):
    """
    laadt de data (clusterdata.txt) in in een dataframe
    input: data clusteren .txt
    output: dataframe met index CloneID
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
    df_clusterdata = pd.DataFrame.from_records(clusterdatalist, columns = ['CloneID','1','2','3','4','5','6','7','8']).set_index('CloneID')
        
    return df_clusterdata

def resultladen(tekstbestand):
    """
    Laadt het tekstbestand met clusterresults in een dictionary
    input: resultaten clusteren .txt
    output: dictionary met CloneID en clusternummer
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

def loadData(GenDescription):
    """
    Python code en txt files moeten in deze map gelocaliseerd zijn.
    Laadt data in dat nodig is voor fase 3 de beschrijving van de genen.
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

def accessioninladen(tekstbestand):
    """
    laadt de data uit accessionnumbers.txt in in een dataframe
    input: tekstbestand accessionnumbers.txt
    output: dataframe met alle accessionnumbers en CloneID als index
    """
    
    #open het bestand en lees het regel voor regel in
    f = open(tekstbestand,"r")
    accessionlines = f.read().splitlines()
    f.close()
    
    accessionlines.pop(0)
    accessionlist = []
    #voeg per regel de verschillende waarden toe aan accessionlist, gesplitst op een tab
    for i in accessionlines:
        linesplit = i.split("\t")
        accessionlist.append(linesplit)
    
        #maak een dataframe van de lijst, met de juiste kolomnamen en index
    df_accessionnumbers = pd.DataFrame.from_records(accessionlist, columns = ['CloneID','EMBL(genbank)ID','unigeneID','ensemblegeneid']).set_index('CloneID')
    
    return df_accessionnumbers

def Telwoorden(dict):
    """
    Telt voor elk woord hoevaak hij in totaal in alle genbeschrijvingen bij elkaar voorkomt
    input: dictionary met CodeID's en descriptions (output van loadData)
    output: dictionary met woord en hoe vaak het in de beschrijvingen voorkomt
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
    input: dictionary met CloneID: description dictGenDes (output van loadData), dictionary met CloneID: cluster, clusterresultdict (output van resultladen)
    output: dictionary met key: woord en value: lijst met item voor elk cluster, bijv op index 1 een getal met hoevaak het woord voorkomt in cluster 1
    """
    
    #voeg cluster nummers toe aan gen descripties voor gelijke CloneID's
    df_clusterresult = pd.DataFrame(list(clusterresultdict.items()), columns = ['CloneID','cluster_nummer'])
    df_genDes = pd.DataFrame(list(dictGenDes.items()), columns = ['CloneID','beschrijving'])
    print(df_genDes)
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
            print(woord_per_cluster)
            
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

def woord_aantal_cluster(voorkomen_woorden_cluster_dict):
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

    woorden_1_cluster = {1:[],2:[],3:[],4:[],5:[]}
    for woord in voorkomen_woorden_cluster_dict:
        #als het aantal clusters waarin een woord voorkomt 1 is, voeg het woord dan toe aan de dict bij het juiste cluster
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
        
    return woorden_1_cluster, woorden_23_cluster, woorden_0_cluster, woorden_all_cluster

def lijkende_clusters(woorden_23_cluster,voorkomen_woorden_cluster_list,voorkomen_woorden_cluster_dict):
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
    #print(embldata)
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

def df_accession_cluster_sequenceinfo(df_accessionnumbers):
    """
    voegt de sequence informatie, sequence en clusternummer toe aan df_accessionnumbers
    input accessionnumbers dataframe (output accessioninladen)
    output accessionnumbers dataframe aangevuld met clusternummer, sequence en baseparen, a,c,g,t,other
    """
    
    #maak een dataframe met accessionnumbers en clusternummer voor elk gen
    df_clusterresult = pd.DataFrame(list(clusterresult.items()), columns = ['CloneID','cluster_nummer'])
    df_merge = pd.merge(df_accessionnumbers, df_clusterresult, how='inner', on='CloneID')
    #maak van de cluster nummers integers
    df_merge['cluster_nummer']=pd.to_numeric(df_merge['cluster_nummer'])
    
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

def sequence_compare(df_merge):
    """
    kijkt voor elk cluster wat per base de 'range' is in hoevaak het in de sequence voorkomt
    input: df_merge (de output uit de functie df_accession_cluster_sequenceinfo)
    output:voor elke cluster een dictionary van de range van a,c,g,t etc.
    """
    
    #sequenceinfo = ['bp','a','c','g','t','other']
    range_1 = {}
    #maak een dataframe aan van df_merge met maar 1 cluster
    cluster_1 = df_merge[df_merge['cluster_nummer']==1]
    #bereken voor elke nucleotide de 'range' oftewel het verschil in maximale waarde en minimale waarde en voeg deze toe aan de dictionary
    range_1['bp']= (cluster_1['bp'].max()-cluster_1['bp'].min())
    range_1['a']= (cluster_1['a'].max()-cluster_1['a'].min())
    range_1['c']= (cluster_1['c'].max()-cluster_1['c'].min())
    range_1['g']= (cluster_1['g'].max()-cluster_1['g'].min())
    range_1['t']= (cluster_1['t'].max()-cluster_1['t'].min())
    
    #doe voor elke cluster hetzelfde als voor cluster 1
    range_2 = {}
    
    cluster_2 = df_merge[df_merge['cluster_nummer']==2]
    range_2['bp']= (cluster_2['bp'].max()-cluster_2['bp'].min())
    range_2['a']= (cluster_2['a'].max()-cluster_2['a'].min())
    range_2['c']= (cluster_2['c'].max()-cluster_2['c'].min())
    range_2['g']= (cluster_2['g'].max()-cluster_2['g'].min())
    range_2['t']= (cluster_2['t'].max()-cluster_2['t'].min())
    
    range_3 = {}
    
    cluster_3 = df_merge[df_merge['cluster_nummer']==3]
    range_3['bp']= (cluster_3['bp'].max()-cluster_3['bp'].min())
    range_3['a']= (cluster_3['a'].max()-cluster_3['a'].min())
    range_3['c']= (cluster_3['c'].max()-cluster_3['c'].min())
    range_3['g']= (cluster_3['g'].max()-cluster_3['g'].min())
    range_3['t']= (cluster_3['t'].max()-cluster_3['t'].min())
    
    range_4 = {}
    
    cluster_4 = df_merge[df_merge['cluster_nummer']==4]
    range_4['bp']= (cluster_4['bp'].max()-cluster_4['bp'].min())
    range_4['a']= (cluster_4['a'].max()-cluster_4['a'].min())
    range_4['c']= (cluster_4['c'].max()-cluster_4['c'].min())
    range_4['g']= (cluster_4['g'].max()-cluster_4['g'].min())
    range_4['t']= (cluster_4['t'].max()-cluster_4['t'].min())
    
    range_5 = {}
    
    cluster_5 = df_merge[df_merge['cluster_nummer']==5]
    range_5['bp']= (cluster_5['bp'].max()-cluster_5['bp'].min())
    range_5['a']= (cluster_5['a'].max()-cluster_5['a'].min())
    range_5['c']= (cluster_5['c'].max()-cluster_5['c'].min())
    range_5['g']= (cluster_5['g'].max()-cluster_5['g'].min())
    range_5['t']= (cluster_5['t'].max()-cluster_5['t'].min())

    return range_1,range_2,range_3,range_4,range_5

def plot(df_merge, base):
    """
    Plot de density van de hoeveelheid van een base in een sequence, per cluster
    input: df_merge (output van functie df_accession_cluster_sequenceinfo), en de base die je wil plotten (bv 'a', 'bp')
    output: density plot van de base per cluster
    """
    
    #maak een dataframe van elk apart cluster
    df_1 = df_merge[df_merge['cluster_nummer']==1]
    df_2 = df_merge[df_merge['cluster_nummer']==2]
    df_3 = df_merge[df_merge['cluster_nummer']==3]
    df_4 = df_merge[df_merge['cluster_nummer']==4]
    df_5 = df_merge[df_merge['cluster_nummer']==5]
    
    #plot voor elk cluster de gewenste base
    df_1[base].plot(kind='density')
    df_2[base].plot(kind='density')
    df_3[base].plot(kind='density')
    df_4[base].plot(kind='density')
    df_5[base].plot(kind='density')
    
def expressieClusters(expressDic,clusterresultdict, df_accessionnumbers):
    """
    kijkt voor de verschillende plekken waar een gen tot expressie komt in welke clusters dat hoevaak voorkomt
    input: expressDic,clusterresultdict, df_accessionnumbers
    output: dictionary
    """
    #voeg cluster nummers toe aan gen descripties voor gelijke CloneID's

    df_express = pd.DataFrame(list(expressDic.items()), columns = ['unigeneID','expressie']).set_index('unigeneID')
    df_merge = pd.merge(df_accessionnumbers, df_express, how='inner', on='unigeneID')
    df_merge['cluster_nummer']= pd.to_numeric(df_merge['cluster_nummer'])

  
    
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
    
    
df_accessionnumbers = accessioninladen('accessionnumbers.txt')
df_merge = df_accession_cluster_sequenceinfo(df_accessionnumbers)
plot(df_merge, 'a')


