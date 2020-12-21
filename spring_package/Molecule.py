class Molecule:
    def __init__(self, fileName=None):
        self.calpha = dict()
        self.biomol = dict()
        self.atoms = list()
        if fileName is not None:
            self.fromFile(fileName)

    def fromFile(self, fileName):
        biomolNumber = 0
        biomolChains = list()
        self.biomol[0] = None
        with open(fileName) as file:
            for line in file:
                key = line[0:6].strip()
                if key == "ATOM":
                    atom = line[12:16]
                    atomNumber = line[6:11]
                    chainName = line[21:22]
                    if chainName not in self.calpha:
                        self.calpha[chainName] = dict()
                    x = self.toFloat(line[30:38])
                    y = self.toFloat(line[38:46])
                    z = self.toFloat(line[46:54])
                    occupancy = self.toFloat(line[54:60], optional=True)
                    temperature = self.toFloat(line[54:60], optional=True)
                    residue = line[17:20]
                    residueNumber = self.toInt(line[22:26])
                    atomNumber = self.toInt(line[6:11])
                    atomName = line[12:16]
                    atomDict = dict(x=x, y=y, z=z,
                                    residue=residue,
                                    occupancy=occupancy,
                                    temperature=temperature,
                                    atomNumber=atomNumber,
                                    atomName=atomName,
                                    residueNumber=residueNumber,
                                    chainName=chainName)
                    if atom.strip() == "CA":
                        self.calpha[chainName][residueNumber] = atomDict
                    self.atoms.append(atomDict)
                biokey = "REMARK 350 BIOMOLECULE:"
                if line[0:len(biokey)] == biokey:
                    biomolNumber = self.toInt(line[len(biokey):])
                    if biomolNumber == 0:
                        raise Exception("Invalid biomolecule identifier [0].")
                    biokey = "REMARK 350 APPLY THE FOLLOWING TO CHAINS:"
                    nextLine = next(file)
                    while nextLine[:len(biokey)] != biokey:
                        nextLine = next(file)
                    biomolChains = nextLine[len(biokey):].split(",")
                    biomolChains = list(map(lambda x: x.strip(), biomolChains))
                    biokey = "REMARK 350                    AND CHAINS:"
                    nextLine = next(file)
                    while nextLine[:len(biokey)] == biokey:
                        moreChains = nextLine[len(biokey):].split(",")
                        moreChains = list(map(lambda x: x.strip(), moreChains))
                        biomolChains = biomolChains + moreChains
                        nextLine = next(file)
                    biokey = "REMARK 350   BIOMT"
                    if nextLine[:len(biokey)] == biokey:
                        biomolMatId1, biomolMat1 = self.getFloats(nextLine)
                        nextLine = next(file)
                        biomolMatId2, biomolMat2 = self.getFloats(nextLine)
                        nextLine = next(file)
                        biomolMatId3, biomolMat3 = self.getFloats(nextLine)
                        if biomolMatId1 != biomolMatId2 or biomolMatId1 != biomolMatId3:
                            raise Exception("Invalid rotation matrix format [%s]." % biomolMatId1)
                        matrix = [biomolMat1, biomolMat2, biomolMat3]
                        biomolChains = [c for c in biomolChains if c]
                        if biomolNumber not in self.biomol:
                            self.biomol[biomolNumber] = list()
                        self.biomol[biomolNumber].append(dict(chains=biomolChains, matrix=matrix))
        removeChains = []
        for chainName in self.calpha:
            if len(self.calpha[chainName]) == 0:
                removeChains.append(chainName)
        for chainName in removeChains:
            del self.calpha[chainName]
        if not self.calpha:
            raise Exception("Molecule has no atoms.")

    def getFloats(self, nextLine):
        matId = self.toInt(nextLine[20:23])
        matLine = nextLine[23:].split()
        matLine = list(map(lambda x: self.toFloat(x), matLine))
        return matId, matLine

    def createUnit(self, biomolNumber=0):
        if biomolNumber == 0:
            return self
        else:
            molecule = Molecule()
            chainCount = dict()
            for matrixDict in self.biomol[biomolNumber]:
                for chain in matrixDict["chains"]:
                    if chain in self.calpha:
                        chainCopy = dict()
                        for residue in self.calpha[chain]:
                            chainCopy[residue] = self.calpha[chain][residue].copy()
                        for atomNumber in chainCopy:
                            atom = chainCopy[atomNumber]
                            rotmat = matrixDict["matrix"]
                            self.applyMatrix(atom, rotmat)
                        if chain in chainCount:
                            chainCount = chainCount[chain]
                            chainName = "%s_%d" % (chain, chainCount)
                            chainCount[chain] = chainCount + 1
                        else:
                            chainName = chain
                            chainCount[chain] = 0
                        molecule.calpha[chainName] = chainCopy
            return molecule

    def getSequence(self, chainName):
        seq = ""
        if chainName not in self.calpha:
            raise Exception("Chain identifier not found [%s]" % chainName)
        chainDict = self.calpha[chainName]
        for residueNumber in sorted(chainDict):
            seq = seq + self.toSingleAmino(chainDict[residueNumber]["residue"])
        return seq

    def applyMatrix(self, atom, rotmat):
        newx = atom["x"] * rotmat[0][0] + atom["y"] * rotmat[0][1] + atom["z"] * rotmat[0][2] + rotmat[0][3]
        newy = atom["x"] * rotmat[1][0] + atom["y"] * rotmat[1][1] + atom["z"] * rotmat[1][2] + rotmat[1][3]
        newz = atom["x"] * rotmat[2][0] + atom["y"] * rotmat[2][1] + atom["z"] * rotmat[2][2] + rotmat[2][3]
        atom["x"] = newx
        atom["y"] = newy
        atom["z"] = newz
        return atom

    def toFloat(self, x, optional=False):
        try:
            return float(x)
        except Exception:
            if not optional:
                raise Exception("Invalid float conversion [%s]." % x)
        return 0.0

    def toInt(self, x):
        try:
            return int(x)
        except Exception:
            raise Exception("Invalid integer conversion [%s]." % x)

    def toSingleAmino(self, seq):
        code = dict(GLY="G", ALA="A", VAL="V", LEU="L", ILE="I", MET="M", PHE="F", PRO="P", TYR="Y", TRP="W",
                    LYS="K", SER="S", CYS="C", ASN="N", GLN="Q", HIS="H", THR="T", GLU="E", ASP="D", ARG="R")
        return code[seq] if seq in code else "X"

    def saveChain(self, chainName, outputName):
        f = open(outputName, "w")
        for residueNumber in sorted(self.calpha[chainName].keys()):  
            ca = self.calpha[chainName][residueNumber]
            if ca["residue"] is not None:
                f.write(self.atomString(ca))
        f.close()

    def save(self, outputName, append=False, chainName=None, payload=None):
        fileFlag = "+a" if append else "w"
        f = open(outputName, fileFlag)
        if payload:
            f.write("%s\n" % payload)
        for atom in self.atoms:
            atom["chainName"] = chainName if chainName else atom["chainName"]
            f.write(self.atomString(atom))
        f.write("TER\n")
        f.close()

    def atomString(self, atom):
        return "ATOM  %5d %s %s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (atom["atomNumber"], atom["atomName"], atom["residue"], atom["chainName"], atom["residueNumber"],
                                                                          atom["x"], atom["y"], atom["z"], atom["occupancy"], atom["temperature"])
