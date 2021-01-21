#! /usr/bin/env python3
import argparse
from os import mkdir, system
from os.path import isfile, isdir

from spring_package.DBKit import DBKit
from spring_package.Molecule import Molecule
from spring_package.Utilities import getId, getChain, getName


def getPDB(line, pdbDatabase):
    pdb = getName(line)
    pdbChain = getChain(line)
    pdbFile = "temp/temp.pdb"
    pdbDatabaseId = "%s.pdb" % pdb
    pdbDatabase.createFile(pdbDatabaseId, pdbFile)
    return pdbFile, pdbChain


def getSequences(fileName):
    sequences = dict()
    with open(fileName) as file:
        for line in file:
            if line.startswith(">"):
                name = getId(line.split()[0][1:])
                nextLine = next(file)
                sequences[name] = nextLine
    return sequences


def findMatch(identifier, templates, databaseFile, pdbDatabase, evalue=0.0):
    if identifier in templates:
        return identifier
    resultSub = identifier[:2]
    fastaFile = "temp/%s/%s.fasta" % (resultSub, identifier)
    resultFile = "%s.result" % fastaFile
    if not isfile(resultFile):
        resultDir = "temp/%s" % resultSub
        if not isdir(resultDir):
            mkdir(resultDir)
        pdbFile, pdbChain = getPDB(identifier, pdbDatabase)
        mol = Molecule(pdbFile)
        seq = mol.getSequence(pdbChain)
        with open(fastaFile, "w") as fasta:
            fasta.write(">%s\n" % identifier)
            fasta.write("%s" % seq)
        system("psiblast -query %s -db %s -out %s" % (fastaFile, databaseFile, resultFile))
    maxMatch = None
    try:
        with open(resultFile) as file:
            for i in range(38):
                line = next(file)
            columns = line.split()
            maxMatch = getId(columns[0])
            maxScore = float(columns[2])
            if maxScore > evalue:
                return None
    except Exception:
        return None
    return maxMatch


def main(args):
    if not isdir("temp"):
        mkdir("temp")
    logFile = open(args.log, "w")
    templateSequenceFile = "temp/templates.fasta"
    pdbDatabase = DBKit(args.index, args.database)
    templates = set()
    with open(args.list) as file:
        for rawId in file:
            templateId = getId(rawId)
            templates.add(templateId)
    if not isfile(templateSequenceFile):
        templateSequences = open(templateSequenceFile, "w")
        with open(args.list) as file:
            for rawId in file:
                templateId = getId(rawId)
                pdbFile, pdbChain = getPDB(templateId, pdbDatabase)
                try:
                    templateMol = Molecule(pdbFile)
                    templateSeq = templateMol.getSequence(pdbChain)
                    templateSequences.write(">%s\n" % templateId)
                    templateSequences.write("%s\n" % templateSeq)
                except Exception:
                    logFile.write("Warning: File not found [%s].\n" % templateId)
        templateSequences.close()
        system("makeblastdb -in %s -dbtype prot" % templateSequenceFile)
    else:
        logFile.write("Using existing sequences for templates [%s].\n" % templateSequenceFile)
    logFile.write("Found %s template entries from `%s`.\n" % (len(templates), args.list))
    logFile.flush()

    crossReference = list()
    with open(args.cross) as file:
        for line in file:
            cols = line.split()
            if len(cols) != 2:
                raise Exception("Invalid line in crossreference [%s]." % line)
            crossReference.append(dict(core=cols[0], partner=cols[1]))
    logFile.write("Loaded crossreference with %d entries.\n" % len(crossReference))
    logFile.flush()

    for refEntry in crossReference:
        coreId = refEntry["core"]
        logFile.write("Processing %s.\n" % coreId)
        coreMatch = findMatch(coreId, templates, templateSequenceFile, pdbDatabase, evalue=args.evalue)
        partnerId = refEntry["partner"]
        logFile.write("Processing %s.\n" % partnerId)
        partnerMatch = findMatch(partnerId, templates, templateSequenceFile, pdbDatabase, evalue=args.evalue)
        if partnerMatch is None or coreMatch is None:
            logFile.write("Warning: Failed alignment [%s, %s].\n" % (coreId, partnerId))
        else:
            logFile.write("Found matching entries [%s, %s].\n" % (coreMatch, partnerMatch))
            refEntry["coreMatch"] = coreMatch
            refEntry["partnerMatch"] = partnerMatch
        logFile.flush()

    entryCount = 0
    with open(args.output, 'w') as output_file:
        for refEntry in crossReference:
            coreId = refEntry["core"]
            partnerId = refEntry["partner"]
            if "coreMatch" in refEntry and "partnerMatch" in refEntry:
                entry = "%s\t%s\t%s\t%s\n" % (refEntry["coreMatch"], refEntry["partnerMatch"], refEntry["core"], refEntry["partner"])
                output_file.write(entry)
                entryCount = entryCount + 1
            else:
                logFile.write("Warning: Skipping failed missing partner match [%s, %s].\n" % (coreId, partnerId))
    logFile.write("Found %s cross reference entries.\n" % entryCount)
    logFile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Maps binding partners to template library')
    parser.add_argument('-l', '--list', help='List of template entries `PDB_CHAIN`', required=True)
    parser.add_argument('-i', '--index', help='PDB Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='PDB Database files (ffdata)', required=True)
    parser.add_argument('-c', '--cross', help='Cross reference (unmapped)', required=True)
    parser.add_argument('-o', '--output', help='Cross reference', required=True)
    parser.add_argument('-g', '--log', help='Log File', required=True)
    parser.add_argument('-e', '--evalue', help='e-Value threshold', type=float, default=0.0)
    args = parser.parse_args()
    main(args)
