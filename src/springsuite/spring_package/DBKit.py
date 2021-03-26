from os.path import isfile
import gzip

class DBKit:
    def __init__(self, indexFile, databaseFile):
        if not isfile(indexFile):
            raise Exception("Index file not found: %s." % indexFile)    
        if not isfile(databaseFile):
            raise Exception("Database file not found: %s." % databaseFile)    
        self.databaseFile = databaseFile
        self.index = dict()
        with open(indexFile) as file:
            for line in file:
                cols = line.split()
                try:
                    identifier = cols[0]
                    start = int(cols[1])
                    size = int(cols[2])
                    self.index[identifier] = [start, size]
                except Exception:
                    raise Exception("Invalid DBKit Index file format: %s." % line)

    def createFile(self, identifier, outputName, zipped=None):
        if zipped:
            identifier = "%s.%s" % (identifier, zipped)
        if identifier in self.index:
            entry = self.index[identifier]
            start = entry[0]
            size = entry[1]
            with open(self.databaseFile) as file:
                file.seek(start)
                content = file.read(size)
                outputFile = open(outputName, "w")
                outputFile.write(content)
                outputFile.close()
                if zipped:
                    outputFile = gzip.open(outputName, 'rb')
                    content = outputFile.read()
                    outputFile.close()
                    outputFile = open(outputName, "w")
                    outputFile.write(content)
                    outputFile.close()
            return True
        else:
            return False
