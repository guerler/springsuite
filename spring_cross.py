#! /usr/bin/env python3
import argparse
from os import mkdir
from os.path import isdir

from spring_package.DBKit import DBKit
from spring_package.Molecule import Molecule
from spring_package.Utilities import getName


def main(args):
    logFile = open(args.log, "w")
    if not isdir("temp"):
        mkdir("temp")
    pdbCount = 0
    partnerList = set()
    entries = set()
    with open(args.index) as file:
        for line in file:
            entries.add(getName(line))
    logFile.write("Found %s template entries.\n" % len(entries))
    pdbDatabase = DBKit(args.index, args.database)
    for pdb in sorted(entries):
        print("Processing %s" % pdb)
        pdbFile = "temp/temp.pdb"
        pdbDatabaseId = "%s.pdb" % pdb
        pdbDatabase.createFile(pdbDatabaseId, pdbFile)
        try:
            mol = Molecule(pdbFile)
        except Exception as e:
            logFile.write("Warning: Entry '%s' not found. %s.\n" % (pdbDatabaseId, str(e)))
            continue
        pdbCount = pdbCount + 1
        for pdbChain in mol.calpha.keys():
            logFile.write("Processing %s, chain %s.\n" % (pdb, pdbChain))
            logFile.write("Found %d biomolecule(s).\n" % len(mol.biomol.keys()))
            for biomolNumber in mol.biomol:
                if biomolNumber == 0:
                    logFile.write("Processing biomolecule.\n")
                    bioMolecule = mol
                else:
                    logFile.write("Processing biomolecule %d.\n" % biomolNumber)
                    bioMolecule = mol.createUnit(biomolNumber)
                nChains = len(bioMolecule.calpha.keys())
                if nChains > 1 and pdbChain in bioMolecule.calpha:
                    for bioChain in bioMolecule.calpha:
                        if bioChain == pdbChain:
                            continue
                        corePdbChain = "%s_%s" % (pdb.upper(), pdbChain[:1])
                        partnerPdbChain = "%s_%s" % (pdb.upper(), bioChain[:1])
                        partnerList.add("%s\t%s" % (corePdbChain, partnerPdbChain))
                else:
                    logFile.write("Skipping: Chain not found or single chain [%s].\n" % pdbChain)
            logFile.flush()
    with open(args.output, 'w') as output_file:
        for entry in sorted(partnerList):
            output_file.write("%s\n" % entry)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='List filtering.')
    parser.add_argument('-i', '--index', help='PDB Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='PDB Database files (ffdata)', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    parser.add_argument('-g', '--log', help='Log File', required=True)
    args = parser.parse_args()
    main(args)
