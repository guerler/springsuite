#! /usr/bin/env python3
import argparse
from os import mkdir, remove
from os.path import isdir, isfile

from spring_package.DBKit import DBKit
from spring_package.Modeller import createModel


class ModelArguments:
    def __init__(self, args):
        self.log = args.log
        self.index = args.index
        self.database = args.database
        self.cross = args.cross
        self.wenergy = args.wenergy
        self.minscore = args.minscore
        self.maxtries = args.maxtries
        self.maxclashes = args.maxclashes
        self.showtemplate = args.showtemplate

    def set(self, a_hhr, b_hhr, output):
        self.a_hhr = a_hhr
        self.b_hhr = b_hhr
        self.output = output


def main(args):
    modelArgs = ModelArguments(args)
    outPath = args.outputpath.rstrip("/")
    if not isdir(outPath):
        mkdir(outPath)
    if not isdir("temp"):
        mkdir("temp")
    dbkit = DBKit(args.hhr_index, args.hhr_database)
    with open(args.pairs, "r") as file:
        for line in file:
            param = line.split()
            aIdentifier = param[0]
            bIdentifier = param[1]
            aFile = "temp/%s" % aIdentifier
            bFile = "temp/%s" % bIdentifier
            if not dbkit.createFile(aIdentifier, aFile):
                print("Failed to retrieve entry %s." % aIdentifier)
                continue
            if not dbkit.createFile(bIdentifier, bFile):
                print("Failed to retrieve entry %s." % bIdentifier)
                continue
            output = "%s/%s.%s.pdb" % (outPath, aIdentifier, bIdentifier)
            modelArgs.set(a_hhr=aFile, b_hhr=bFile, output=output)
            createModel(modelArgs)
            if isfile(aFile):
                remove(aFile)
            if isfile(bFile):
                remove(bFile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create 3D models from HH-search results.')
    parser.add_argument('-p', '--pairs', help='Interaction table e.g. from min-Z evaluation (2-columns)', required=True)
    parser.add_argument('-ih', '--hhr_index', help='HHR Index database file (ffindex)', required=True)
    parser.add_argument('-dh', '--hhr_database', help='HHR Database file (ffdata)', required=True)
    parser.add_argument('-i', '--index', help='PDB Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='PDB Database file (ffdata)', required=True)
    parser.add_argument('-c', '--cross', help='PDB Cross Reference', required=True)
    parser.add_argument('-g', '--log', help='Log file', required=True)
    parser.add_argument('-o', '--outputpath', help='Path to output directory', required=True)
    parser.add_argument('-we', '--wenergy', help='Weight Energy term', type=float, default=-0.01, required=False)
    parser.add_argument('-ms', '--minscore', help='Minimum min-Z score threshold', type=float, default=10.0, required=False)
    parser.add_argument('-mt', '--maxtries', help='Maximum number of templates', type=int, default=20, required=False)
    parser.add_argument('-mc', '--maxclashes', help='Maximum fraction of clashes', type=float, default=0.1, required=False)
    parser.add_argument('-sr', '--showtemplate', help='Add reference template to model structure', required=False, default="true")
    args = parser.parse_args()
    main(args)
