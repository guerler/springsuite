#!/usr/bin/env python3
import argparse

from dbkit_package.DBKit import DBKit, writeEntry


def main(args):
    logFile = open(args.log, "w")
    outputIndex = args.outputindex
    outputDatabase = args.outputdatabase
    entries = list()
    with open(args.list, "r") as f:
        for line in f:
            name = line.split()[0]
            entries.append(name)
    logFile.write("Detected %s entries.\n" % len(entries))
    fileName = "temp.dat"
    count = 0
    dbkit = DBKit(args.index, args.database)
    for entry in sorted(entries):
        success = dbkit.createFile(entry, fileName)
        if success:
            writeEntry(entry, fileName, outputIndex, outputDatabase)
            count = count + 1
        else:
            logFile.write("Entry %s not found.\n" % entry)
    logFile.write("Extracted %s entries.\n" % count)
    logFile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DBKit - Merge database pair.')
    parser.add_argument('-l', '--list', help='List of entries to be extracted', required=True)
    parser.add_argument('-i', '--index', help='Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='Database Data file (ffdata)', required=True)
    parser.add_argument('-oi', '--outputindex', help='Output Index file', required=True)
    parser.add_argument('-od', '--outputdatabase', help='Output Data file', required=True)
    parser.add_argument('-g', '--log', help='Log file', required=True)
    args = parser.parse_args()
    main(args)
