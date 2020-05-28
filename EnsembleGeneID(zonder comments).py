"""
Informatie Ensemble gene ID van website: http://www.ensembl.org/index.html
Definities maken een dictionary aan met verkregen info
"""
import requests, sys
import urllib.request
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def phenotypeEnsemble(filename_EMBL_ID, filenameClusterData, file_clusterResult):

    # open het tekstbestand, lees het regel voor regel in en sluit het
    f = open(filenameClusterData, "r")
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
    df_clusterdata = df_clusterdata.drop(['1','2','4','7','14','21','45','90'], axis=1)
    
    #open het tekstbestand, lees het regel voor regel in en sluit het
    f = open(file_clusterResult, "r")
    clusterresultlines = f.read().splitlines()
    f.close()
    
    #maak een lege dictionary aan
    clusterresultdict = {}
    #iedere regel wordt gesplitst op 3 spaties, zodat ze als key:value aan de dictionary toegevoegd kunnen worden
    for i in clusterresultlines:
        linesplit = i.split('   ')
        clusterresultdict[linesplit[0]]=linesplit[1]

    #De bijbehorende cluster nummers zijn toegevoegd aan de CloneId's
    df_clusterdata['cluster_nummer'] = df_clusterdata.CloneID.replace(clusterresultdict)
    
    # open file
    infile_accessionnumbers = open(filename_EMBL_ID)
    lines = infile_accessionnumbers.readlines()
    infile_accessionnumbers.close()
    
    # data omzetten naar een dataframe    
    accessionlist = []
    for line in lines:
        linesplit = line.split("\t")
        accessionlist.append(linesplit)
    
    df_accessionnumbers = pd.DataFrame.from_records(accessionlist, columns = ['CloneID','EMBL(genbank)ID','unigeneID','ensemblegeneid'])
    df_accessionnumbers = df_accessionnumbers.drop(df_accessionnumbers.index[0])
    
    df_cluster_accessionnumbers = pd.merge(df_clusterdata, df_accessionnumbers[["CloneID", "EMBL(genbank)ID", 'ensemblegeneid']], on="CloneID", how="left")
    
    df_cluster_accessionnumbers['ensemblegeneid'].replace('', np.nan, inplace=True)
    df_cluster_accessionnumbers.dropna(inplace=True)
    
    
    #https://www.ebi.ac.uk/ena/browser/api/embl/AI552275.1?lineLimit=1000 openen en in een lijst zetten
    website = "https://rest.ensembl.org/phenotype/gene/mus_musculus/"
    download = "?include_associated=1"
    
    description = {}
    for ID in df_cluster_accessionnumbers['ensemblegeneid']:
        url = website + ID + download
        r = requests.get(url, headers={ "Content-Type" : "application/json"})
    
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.json()

        embldata = []
        for line in data:
            line = str(line)
            line = line[23:-4]
            line.split("\t")
            embldata.append(line)
            
            # Alle tissue types in een dictionairy zetten met CloneID
            if '/tissue_type=' in line:
                tissue_line = line[14:]
                description[ID]=tissue_line
                print(tissue_line)

"""
    #df_cluster_EMBL['Tissue type'] = df_cluster_EMBL['EMBL(genbank)ID'].replace(tissue_type)
    df_cluster_accessionnumbers['Tissue type'] = df_cluster_accessionnumbers['EMBL(genbank)ID'].map(tissue_type)    
    df_cluster_accessionnumbers = df_cluster_accessionnumbers[df_cluster_accessionnumbers['Tissue type'].notna()].reset_index()
    
    

    # Dataframe maken met als groepen, cluster nummer en tissue type
    grouped_cluster_tissue = df_cluster_EMBL.groupby(['cluster_nummer', 'Tissue type']).size()
    grouped_cluster_tissue = grouped_cluster_tissue.to_frame()
    grouped_cluster_tissue = grouped_cluster_tissue.rename(columns={0:'count'})
    
    # plots maken van de verschillende clusters met aantal tissue types
    fig, ax = plt.subplots(figsize=(10,30))

    for group in grouped_cluster_tissue:
        grouped_cluster_tissue.unstack(level=0).plot(kind='bar', subplots=True, ax=ax)
        ax.set_xlabel('Cluster nummer')
        ax.set_ylabel("Aantal CloneId's in een cluster per weefsel type")
        ax.set_xlim(0,6)
        ax.set_ylim(0,15)
        plt.tight_layout()
    """
phenotypeEnsemble('accessionnumbers.txt', 'Voorbeeld_clusterdata.txt', 'Voorbeeld_clusterresult.txt')
