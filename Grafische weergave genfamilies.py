import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

def grafiekenGenFamilies(filenameClusterData, filenameCloneIdFam):
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
    df_clusterdata = pd.DataFrame.from_records(clusterdatalist, columns = ['CloneID','1','2','4','7','14','21','45','90']).set_index('CloneID')
    
    #inladen van tweede bestand 
    infileCloneIdFamily = open(filenameCloneIdFam)
    cloneIdFamily = infileCloneIdFamily.readlines()
    infileCloneIdFamily.close()
    
    #dictionary maken met als key de cloneID en als value het familienummer
    dicCloneIdFamily = {}
    for cloneLine in cloneIdFamily:
        lines = cloneLine.split()
        key = lines[0]
        value = lines[1]
        dicCloneIdFamily[key] = value
    
    #dicCloneIdFamily omzetten naar dataframe met CloneID en bijhorende familienummer
    df_fam = pd.DataFrame(dicCloneIdFamily.items(), columns=['CloneID', 'Familienummer'])
    df_fam.drop(df_fam.index[0], inplace=True)
    
    #nieuwe dataframe bestaande uit samengevoegde df_fam en df_cluster met overeenkomende CloneID
    df_merge = pd.merge(df_fam, df_clusterdata, how='inner', on='CloneID')
    
    #lijst aanmaken met dataframes gegroepeerd volgens familienummer
    dagenFam = [df_merge.loc[df_merge['Familienummer']==nr,:] for nr in df_merge['Familienummer'].unique()]
    
    #lijst aanmaken met nummers van de gebruikte dagen voor x-as
    dagen = ['1','2','4','7','14','21','45','90']
    
    #dataframe aanmaken met de familienummer en bijhorende familie naam 
    namenFam = ['A Disintegrin-like and Metalloprotease with Thrombospondin Type 1 Motif Family','Adenylate Cyclase Family','Alpha Actinin Family','Annexin Family','Apolipoprotein C Family','Aquaporin Family','Calcium Channel, Voltage-dependent, Gamma Subunit Family','Calpain Family','Claudin Family','Collagen Family','Cyclic Nucleotide Phosphodiesterase Superfamily','E2f Transcription Factor Family','EH-domain Containing Family','Fatty Acid Coenzyme A Ligase Family','Gamma-Aminobutyric Acid Ionotropic Type A Receptor Family','Glutamate Receptor, Ionotropic Superfamily','Glycine Receptor, Alpha Family','Glypican Family','Iroquois Family','Leucine-rich Repeat LGI Family','Opsin Family','Phosphatidylinositol-4-phosphate 5-Kinase Family','Protocadherin Gamma Family', 'Secretory Carrier Membrane Protein Family','Synaptotagmin Family','Tripartite Motif Family']
    df_namen = pd.DataFrame(namenFam, columns=['Fam_naam'])
    
    #een grafiek voor elke voorkomende genfamilie bij experiment
    for fam in dagenFam:
        famNr = fam.loc[fam.index[0],'Familienummer']
        
        #(dag, relatieve expressiewaarden van die dag) plotten
        for x, row in fam.iterrows():
            for x in range(8):
                dag = dagen[x]
                plt.scatter(dag, float(row[x+1]))
        
        #grafiek en assen benoemen
        plt.suptitle('Grafiek van genfamilie: ' + df_namen.loc[df_namen.index[int(famNr)-1], 'Fam_naam'])
        plt.xlabel('Dagen')
        plt.ylabel('Relatieve expressiewaarde')
        plt.show()
        
#data runnen
grafiekenGenFamilies('Voorbeeld_clusterdata.txt','CloneIdFamily.txt')