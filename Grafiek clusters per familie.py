import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def grafiekenGenFamilies(filenameClusterData, filenameCloneIdFam, file_clusterResult):
    #input: data clusteren .txt en gen families .txt
    #output: grafieken met (dag, relatieve expressiewaarden van die dag) voor elke gen familie die voorkomt bij het experiment
    
    #open het tekstbestand, lees het regel voor regel in en sluit het
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
    
    #inladen van tweede bestand 
    infileCloneIdFamily = open(filenameCloneIdFam)
    cloneIdFamily = infileCloneIdFamily.readlines()
    infileCloneIdFamily.close()
    
    #dictionary maken met als key de cloneID en als value het familienummer
    dicCloneIdFamily = {}
    for cloneLine in cloneIdFamily[1:]:
        lines = cloneLine.split()
        dicCloneIdFamily[lines[0]] = int(lines[1])
        
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

    # De familienummers toegevoegd aan de     
    df_clusterdata['Familienummer'] = df_clusterdata.CloneID.replace(dicCloneIdFamily)
    #dictNamen = {'0':'A Disintegrin-like and Metalloprotease with Thrombospondin Type 1 Motif Family' ,'1':'Adenylate Cyclase Family','2':'Alpha Actinin Family','3':'Annexin Family','4':'Apolipoprotein C Family','5':'Aquaporin Family','6':'Calcium Channel, Voltage-dependent, Gamma Subunit Family','7':'Calpain Family','8':'Claudin Family','9':'Collagen Family','10':'Cyclic Nucleotide Phosphodiesterase Superfamily','11':'E2f Transcription Factor Family','12':'EH-domain Containing Family','13':'Fatty Acid Coenzyme A Ligase Family','14':'Gamma-Aminobutyric Acid Ionotropic Type A Receptor Family','15':'Glutamate Receptor, Ionotropic Superfamily','16':'Glycine Receptor, Alpha Family','17':'Glypican Family','18':'Iroquois Family','19':'Leucine-rich Repeat LGI Family','20':'Opsin Family','21':'Phosphatidylinositol-4-phosphate 5-Kinase Family','22':'Protocadherin Gamma Family','23': 'Secretory Carrier Membrane Protein Family','24':'Synaptotagmin Family','25':'Tripartite Motif Family'}
    df_clusterdata = df_clusterdata[df_clusterdata.Familienummer < 27]

    
    #dataframe aanmaken met de familienummer en bijhorende familie naam 
    namenFam = ['A Disintegrin-like and Metalloprotease with Thrombospondin Type 1 Motif Family','Adenylate Cyclase Family','Alpha Actinin Family','Annexin Family','Apolipoprotein C Family','Aquaporin Family','Calcium Channel, Voltage-dependent, Gamma Subunit Family','Calpain Family','Claudin Family','Collagen Family','Cyclic Nucleotide Phosphodiesterase Superfamily','E2f Transcription Factor Family','EH-domain Containing Family','Fatty Acid Coenzyme A Ligase Family','Gamma-Aminobutyric Acid Ionotropic Type A Receptor Family','Glutamate Receptor, Ionotropic Superfamily','Glycine Receptor, Alpha Family','Glypican Family','Iroquois Family','Leucine-rich Repeat LGI Family','Opsin Family','Phosphatidylinositol-4-phosphate 5-Kinase Family','Protocadherin Gamma Family', 'Secretory Carrier Membrane Protein Family','Synaptotagmin Family','Tripartite Motif Family']
    df_namen = pd.DataFrame(namenFam, columns=['Fam_naam'])
    
    # Dataframe maken met als groepen, familienummer en cluster nummer
    grouped_cluster_familie = df_clusterdata.groupby(['Familienummer', 'cluster_nummer']).size()
    grouped_cluster_familie = grouped_cluster_familie.to_frame()
    grouped_cluster_familie = grouped_cluster_familie.rename(columns={0:'count'})
    print(grouped_cluster_familie)
    
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
        
#data runnen
grafiekenGenFamilies('Voorbeeld_clusterdata.txt','CloneIdFamily.txt', 'Voorbeeld_clusterresult.txt')