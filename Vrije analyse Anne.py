"""
Vrij analyse expressie
"""
import urllib.request
import gzip


# Read the first 64 bytes of the file inside the .gz archive located at `url`
url ='https://ftp.ncbi.nlm.nih.gov/repository/UniGene/Mus_musculus/Mm.data.gz'

with urllib.request.urlopen(url) as response:
    with gzip.GzipFile(fileobj=response) as binary_file:
        
        chromosomeDic = {}
        expressDic = {}
        
        data = binary_file.read().splitlines()
        for line in data:
            line = line.decode()
            line = line.split()
            
            if line != []:
                
                if line[0] == 'ID':
                    key = str(line[1:])
                    
                elif line[0] == 'EXPRESS':
                    expressDic.update({key: line[1:]})
    
                elif line[0] == 'CHROMOSOME':
                    chromosomeDic.update({key: tuple(line[1:])})
        
        
        
        print(chromosomeDic)
        print(expressDic)