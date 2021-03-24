import math
from os.path import dirname, realpath

NTYPE = 21
NDIST = 20
NSCALE = 2.0


class Energy:
    def __init__(self):
        self.dfire = list()
        dirPath = dirname(realpath(__file__))
        with open("%s/Energy.data" % dirPath) as file:
            for line in file:
                self.dfire.append(float(line))

    def get(self, residuesA, residuesB):
        result = 0
        for atomA in residuesA:
            indexA = self.toResCode(atomA["alignedResidue"])
            for atomB in residuesB:
                indexB = self.toResCode(atomB["alignedResidue"])
                dist2 = ((atomA["x"] - atomB["x"]) ** 2 +
                         (atomA["y"] - atomB["y"]) ** 2 +
                         (atomA["z"] - atomB["z"]) ** 2)
                dist = int((math.sqrt(dist2) * NSCALE))
                if dist < NDIST:
                    index = indexA * NTYPE * NDIST + indexB * NDIST + dist
                    result = result + self.dfire[index]
        return result

    def getClashes(self, moleculeA, moleculeB, minDist=5.0):
        minDist = minDist ** 2
        clashes = 0
        chainA = list(moleculeA.calpha.keys())[0]
        chainB = list(moleculeB.calpha.keys())[0]
        calphaA = moleculeA.calpha[chainA]
        calphaB = moleculeB.calpha[chainB]
        lenA = len(calphaA.keys())
        lenB = len(calphaB.keys())
        if lenA > lenB:
            temp = calphaB
            calphaB = calphaA
            calphaA = temp
            lenA = len(calphaA.keys())
            lenB = len(calphaB.keys())
        for i in calphaA:
            atomA = calphaA[i]
            for j in calphaB:
                atomB = calphaB[j]
                dist2 = ((atomA["x"] - atomB["x"]) ** 2 +
                         (atomA["y"] - atomB["y"]) ** 2 +
                         (atomA["z"] - atomB["z"]) ** 2)
                if dist2 < minDist:
                    clashes = clashes + 1
                    break
        return clashes / float(lenA)

    def toResCode(self, seq):
        code = dict(A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9, M=10,
                    N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19)
        return code[seq] if seq in code else 20
