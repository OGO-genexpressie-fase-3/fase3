import urllib.request
import pandas as pd

def tissueType(filename_EMBL_ID):
    # open file
    infile_EMBL = open(filename_EMBL_ID)
    lines = infile_EMBL.readlines()
    infile_EMBL.close()
    
    # data omzetten naar een dataframe    
    accessionlist = []
    for line in lines:
        linesplit = line.split("\t")
        accessionlist.append(linesplit)
    
    df_accessionnumbers = pd.DataFrame.from_records(accessionlist, columns = ['CloneID','EMBL(genbank)ID','unigeneID','ensemblegeneid'])
    
    #variable_EMBL_ID = df_accessionnumbers[['EMBL(genbank)ID']]
    df_variable_EMBL_ID = df_accessionnumbers['EMBL(genbank)ID'].copy()
    df_variable_EMBL_ID.replace('', np.nan, inplace=True)
    df_variable_EMBL_ID.dropna( inplace=True)
    
    #https://www.ebi.ac.uk/ena/browser/api/embl/AI552275.1?lineLimit=1000
    website = "https://www.ebi.ac.uk/ena/browser/api/embl/"
    download = ".1?lineLimit=1000"
    for ID in df_variable_EMBL_ID[1:]:
        url = website + ID + download
        rf = urllib.request.urlopen(url)
        data = rf.readlines()
        embldata = []
        for line in data:
            line = str(line)
            line = line[23:-3]
            line.split("\t")
            embldata.append(line)
            
            tissue_type = []
            if '/tissue_type=' in line:
                line = line[14:]
                tissue_type.append(line)
    print(tissue_type)

tissueType('accessionnumbers.txt')
