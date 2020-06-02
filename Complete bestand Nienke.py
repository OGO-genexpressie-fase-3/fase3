import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics as st
import urllib.request

def dataGenDescription(GenDescription):
    """
    Python code en txt files moeten in deze map gelocaliseerd zijn.
    Laad data in dat nodig is voor fase 3 de beschrijving van de genen.
    Input: textbestand GenDescription
    Output: Dictionary
    """
    #open het tekstbestand, lees het regel voor regel in en sluit het    
    infileGen = open(GenDescription)
    lines = infileGen.read().splitlines()
    infileGen.close()
    
    # Een lijst maken met alle CloneIds als key position 
    locsID=list(range(0,len(lines),3))  
    IDsGen=[]
    for locID in locsID:
        ID=lines[locID]
        IDsGen.append(ID)

    # Een lijst maken met alle descriptions als value position 
    locsDescription=list(range(1,len(lines),3))    
    descriptionGen=[]
    for locDescription in locsDescription:
        #add descriptions to list
        description=lines[locDescription]
        descriptionGen.append(description)
    
    # dictionairy maken met de keys en values
    dictGenDes = {}
    for ID in IDsGen:
        dictGenDes[ID]=descriptionGen[IDsGen.index(ID)]
    
    return dictGenDes

dataGenDescription('GenDescription2.txt')


def dataClusternummer(tekstbestand):
    """
    Bestand openen en lezen en een dictionairy maken met als key het cloneId 
    nummer en als value het bijbehorende clusternummer
    input: resultaten clusteren .txt
    output: dictionary met CloneID en clusternummer
    """
    
    #open het tekstbestand, lees het regel voor regel in en sluit het
    f = open(tekstbestand, "r")
    clusterresultlines = f.read().splitlines()
    f.close()
    
    # Een dictionairy maken met als keys de cloneIDs en als value het bijbehorende clusternummer
    clusterresultdict = {}
    for i in clusterresultlines:
        linesplit= i.split("   ")
        clusterresultdict[linesplit[0]]=linesplit[1]
    
    return clusterresultdict

dataClusternummer("Voorbeeld_clusterresult.txt")


def dataClusterdata(file_clusterData, clusterresultdict):
    """ Een functie die een dataframe maakt met de clusterdata en clusternummer per cloneID
    input: data clusteren .txt
    output: 
    """
    
    #open het tekstbestand, lees het regel voor regel in en sluit het
    f = open(file_clusterData, "r")
    clusterdatalines = f.read().splitlines()
    f.close()

    # Een lijst maken met de clusterdata en deze omzetten naar een dataframe emt het CloneID als index
    clusterdatalist = []
    for i in clusterdatalines:
        linesplit = i.split("  ")
        clusterdatalist.append(linesplit)
    df_clusterdata = pd.DataFrame.from_records(clusterdatalist, columns = ['CloneID','1','2','3','4','5','6','7','8'])
    
    #De bijbehorende cluster nummers zijn toegevoegd aan de CloneId's
    df_clusterdata['cluster_nummer'] = df_clusterdata.CloneID.replace(clusterresultdict)
    
    return df_clusterdata

dataClusterdata("Voorbeeld_clusterdata.txt")


def plotClusters(df_clusterdata):
    """ Er wordt een scatterplot gemaakt voor elk cluster met het clustergemiddelde erbij
    input: df_clusterdata
    output: voor elk cluster een plot
    """
    # dagen voor op de x-as
    dagen = ['1','2','4','7','14','21','45','90']

    #voor elk cluster een apart dataframe maken
    df_cluster=[ df_clusterdata.loc[df_clusterdata['cluster_nummer']==val,:] for val in df_clusterdata['cluster_nummer'].unique() ]
    for df_per_cluster in df_cluster:
        df_per_cluster.loc[df_per_cluster.index[0],'cluster_nummer']
    
        #creëer lege lijsten voor het clustergemiddelde te berekenen
        x=[]
        y=[]
        
        # voor elk cluster een grafiek
        for i, row in df_per_cluster.iterrows():
            for i in range(8):
                dag = dagen[i]
                plt.scatter(dag, float(row[i+1]), c='orange')
                x.append(i)
                y.append(float(row[i+1]))
        
        #plot clustergemiddelde en de assen benoemen
        plt.scatter(float(st.mean(x)),float(st.mean(y)), c='red')
        plt.suptitle('A scatterplot for cluster ' + df_per_cluster.loc[df_per_cluster.index[0],'cluster_nummer'])
        plt.xlabel('Dag')
        plt.ylabel('Relatieve expressiewaarde')
        plt.show()
    
def dataCloneIdFam(filenameCloneIdFam):
    """  Bestand openen en lezen en maakt een dictionairy met als keys de cloneID en als value het familienummer
    input: cloneIDFamily.tx
    output: dictionairy met cloneID en familienummer
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

dataCloneIdFam('CloneIdFamily.txt')


def grafiekenGenFamilies(dicCloneIdFamily, df_clusterdata):
    """ Grafieken maken per gen familie en de clusters
    #input: gen families .txt, df_clusterdata
    #output: grafieken met (dag, relatieve expressiewaarden van die dag) voor elke gen familie die voorkomt bij het experiment
    """

    # De familienummers toegevoegd aan de bijbehorende cloneID in het dataframe   
    df_clusterdata['Familienummer'] = df_clusterdata.CloneID.replace(dicCloneIdFamily)
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

def dataAccessionnumbers(filename_EMBL_ID):
    """ Bestand openen en lezen en maakt een dataframe met de data
    input: accessionnumbers.txt
    output: dataframe van de data
    """
    
    #open het tekstbestand, lees het regel voor regel in en sluit het    
    infile_EMBL = open(filename_EMBL_ID)
    lines = infile_EMBL.readlines()
    infile_EMBL.close()
    
    # data omzetten naar een dataframe    
    accessionlist = []
    for line in lines:
        linesplit = line.split("\t")
        accessionlist.append(linesplit)
    
    # columns toevoegen aan het dataframe
    df_accessionnumbers = pd.DataFrame.from_records(accessionlist, columns = ['CloneID','EMBL(genbank)ID','unigeneID','ensemblegeneid'])
    df_accessionnumbers = df_accessionnumbers.drop(df_accessionnumbers.index[0])
    return df_accessionnumbers

dataAccessionnumbers('accessionnumbers.txt')

    
def tissueType(df_accessionnumbers, df_clusterdata):
    """ function waar kan worden gezien hoeveel genen er uit een bepaald weefsel
    komen in de verschillende clusters
    input: df_accessionnumbers, df_clusterdata
    output: plots voor elk cluster met het aantal genen per soort weefsel
    """
    
    # dataframe df_cluster en df_accessionnumbers samenvoegen
    df_cluster_EMBL = pd.merge(df_clusterdata, df_accessionnumbers[["CloneID", "EMBL(genbank)ID"]], on="CloneID", how="left")
    df_cluster_EMBL['EMBL(genbank)ID'].replace('', np.nan, inplace=True)
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
        grouped_cluster_tissue.unstack(level=0).plot(kind='bar', subplots=True)
        ax.set_xlabel('Cluster nummer')
        ax.set_ylabel("Aantal CloneId's in een cluster per weefsel type")
        ax.set_xlim(0,6)
        ax.set_ylim(0,15)
        plt.tight_layout()
    

def ESTs(df_clusterdata, dictGenDes):
    """ Een bar plot maken met op de x-as de clusters en op de y-as het aantal genen
    met een EST in de beschijving van het gen. """
    # Een dataframe maken met waar de gen beschrijving EST bevast, de cloneIDs en de clusternummers per gen 
    df_clusterdata['Description'] = df_clusterdata.CloneID.replace(dictGenDes)
    df_ESTs = df_clusterdata.loc[df_clusterdata['Description'] == 'ESTs']
    df_ESTs = df_ESTs[['CloneID', 'cluster_nummer', 'Description']]
    
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

