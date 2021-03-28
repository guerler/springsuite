#!/usr/bin/env python3
import argparse
from os.path import getsize
from shutil import copyfile

from dbkit_package.DBKit import DBKit, writeEntry


def main(args):
    logFile = open(args.log, "w")
    outputIndex = args.outputindex
    outputDatabase = args.outputdatabase
    if getsize(args.firstindex) > getsize(args.secondindex):
        firstIndex = args.firstindex
        firstData = args.firstdata
        secondIndex = args.secondindex
        secondData = args.seconddata
    else:
        firstIndex = args.secondindex
        firstData = args.seconddata
        secondIndex = args.firstindex
        secondData = args.firstdata
    copyfile(firstIndex, outputIndex)
    copyfile(firstData, outputDatabase)
    firstEntries = set()
    with open(firstIndex, "r") as f:
        for line in f:
            name = line.split()[0]
            firstEntries.add(name)
    logFile.write("Detected %s entries.\n" % len(firstEntries))
    secondEntries = list()
    with open(secondIndex, "r") as f:
        for line in f:
            name = line.split()[0]
            secondEntries.append(name)
    fileName = "temp.dat"
    count = 0
    dbkit = DBKit(secondIndex, secondData)
    for secondKey in secondEntries:
        if secondKey not in firstEntries:
            dbkit.createFile(secondKey, fileName)
            writeEntry(secondKey, fileName, outputIndex, outputDatabase)
            count = count + 1
        else:
            logFile.write("Skipping existing entry %s.\n" % secondKey)
    logFile.write("Added %s entries.\n" % count)
    logFile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DBKit - Merge database pair.')
    parser.add_argument('-i', '--firstindex', help='First Index file', required=True)
    parser.add_argument('-d', '--firstdata', help='First Data file', required=True)
    parser.add_argument('-si', '--secondindex', help='Second Index file', required=True)
    parser.add_argument('-sd', '--seconddata', help='Second Data file', required=True)
    parser.add_argument('-oi', '--outputindex', help='Output Index file', required=True)
    parser.add_argument('-od', '--outputdatabase', help='Output Data file', required=True)
    parser.add_argument('-log', '--log', help='Log file', required=True)
    args = parser.parse_args()
    main(args)
