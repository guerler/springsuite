#!/usr/bin/env python3
import argparse
from os import remove
from os.path import getsize, isfile

from dbkit_package.DBKit import writeEntry

try:
    import wget
except Exception:
    wget = None


def getIdentifiers(args):
    entries = set()
    with open(args.list) as file:
        for line in file:
            entry = line.split()[0]
            idLength = int(args.idlength)
            if idLength > 0:
                entry = entry[:idLength]
            if args.idcase == "lower":
                entry = entry.lower()
            elif args.idcase == "upper":
                entry = entry.upper()
            if args.idextension:
                entry = "%s%s" % (entry, args.idextension)
            if args.idprefix:
                entry = "%s%s" % (args.idprefix, entry)
            entries.add(entry)
    return sorted(entries)


def main(args):
    entries = getIdentifiers(args)
    logFile = open(args.log, "w")
    logFile.write("Found %s entries.\n" % len(entries))
    outputIndex = args.index
    outputDatabase = args.database
    if isfile(outputDatabase):
        remove(outputDatabase)
    for entryId in entries:
        logFile.write("Loading %s.\n" % entryId)
        if args.url:
            fileName = wget.download("%s%s" % (args.url, entryId))
        else:
            pathName = args.path.rstrip("/")
            fileName = "%s/%s" % (pathName, entryId)
        if isfile(fileName):
            entrySize = getsize(fileName)
            if entrySize == 0:
                logFile.write("Entry `%s` not found.\n" % entryId)
            else:
                writeEntry(entryId, fileName, outputIndex, outputDatabase)
        else:
            logFile.write("Content not found: %s.\n" % fileName)
        logFile.flush()
    logFile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DBKit - Download and Merge files into a single file.')
    parser.add_argument('-l', '--list', help='List of entries', required=True)
    parser.add_argument('-u', '--url', help='Source Url', required=False)
    parser.add_argument('-p', '--path', help='Path to files', required=False)
    parser.add_argument('-il', '--idlength', help='Format Identifier Length (integer)', required=False, default="0")
    parser.add_argument('-ic', '--idcase', help='Format Identifier Case (lower, upper)', required=False, default=None)
    parser.add_argument('-ie', '--idextension', help='Format Identifier Suffix', required=False, default=None)
    parser.add_argument('-ip', '--idprefix', help='Format Identifier Prefix', required=False, default=None)
    parser.add_argument('-o', '--index', help='Output Database Index', required=True)
    parser.add_argument('-d', '--database', help='Output Database', required=True)
    parser.add_argument('-g', '--log', help="Log file", required=True)
    args = parser.parse_args()
    main(args)
