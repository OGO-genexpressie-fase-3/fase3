
def readFile(filenameAccNumbers, filenameCloneId):
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
    print(dicAccNumbersCloneId)
    
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
    print(dicCloneIdFamily)


readFile('accessionnumbers.txt', 'CloneIdFamily.txt')
