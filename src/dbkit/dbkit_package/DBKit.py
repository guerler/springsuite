from os.path import getsize, isfile


class DBKit:
    def __init__(self, indexFile, databaseFile):
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

    def createFile(self, identifier, outputName):
        if identifier in self.index:
            entry = self.index[identifier]
            start = entry[0]
            size = entry[1]
            with open(self.databaseFile, "rb") as file:
                file.seek(start)
                content = file.read(size)
                outputFile = open(outputName, "wb")
                outputFile.write(content)
                outputFile.close()
            return True
        else:
            return False

    def getIndex(self):
        return self.index


def writeEntry(identifier, fileName, outputIndex, outputDatabase):
    if isfile(outputDatabase):
        currentSize = getsize(outputDatabase)
    else:
        currentSize = 0
    if isfile(fileName):
        entrySize = getsize(fileName)
    else:
        entrySize = 0
    if entrySize > 0:
        outputIndexFile = open(outputIndex, "a+")
        outputIndexFile.write("%s\t%s\t%s\n" % (identifier, currentSize, entrySize))
        tempFile = open(fileName, "rb")
        databaseFile = open(outputDatabase, "ab+")
        databaseFile.write(tempFile.read())
        databaseFile.close()
        tempFile.close()
        outputIndexFile.close()
        return True
    else:
        return False
