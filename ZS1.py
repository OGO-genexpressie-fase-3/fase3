
"""
Inladen van data: gendescriptions als dictionary met keys de ID's en de
bij horende beschrijving als values
"""
def loadData(locGen2, locGen):
    """
    Python code en txt files moeten in deze map gelocaliseerd zijn.
    Laad data in dat nodig is voor fase 3 de beschrijving van de genen.
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

loadData('GenDescription2.txt','GenDescription.txt')
