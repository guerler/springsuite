#!/usr/bin/env python
import argparse
import math
import re
from os.path import isfile

import pandas as pd

METHODS = ["Biochemical Activity",
           "Co-fractionation",
           "Co-localization",
           "Far Western",
           "FRET",
           "PCA",
           "Co-crystal Structure",
           "Co-purification",
           "Two-hybrid",
           "Affinity Capture-MS"]


def getIds(rawIds):
    return rawIds.split("|")


def getCenterId(rawId):
    elements = rawId.split("|")
    if len(elements) > 1:
        return elements[1]
    return rawId


def getOrganism(rawId):
    elements = rawId.split("_")
    return elements[-1]


def getKey(a, b):
    if a > b:
        name = "%s_%s" % (a, b)
    else:
        name = "%s_%s" % (b, a)
    return name


def getPercentage(rate, denominator):
    if denominator > 0:
        return 100.0 * rate / denominator
    return 0.0


def getFilter(filterName):
    print("Loading target organism(s)...")
    filterSets = dict()
    with open(filterName) as filterFile:
        for line in filterFile:
            columns = line.split()
            for colIndex in [0, 1]:
                if colIndex >= len(columns):
                    break
                colEntry = columns[colIndex]
                id = getCenterId(colEntry)
                organism = getOrganism(colEntry)
                if organism not in filterSets:
                    filterSets[organism] = set()
                filterSets[organism].add(id)
    print("Organism(s) in set: %s." % filterSets.keys())
    return filterSets


def getReference(fileName, filterA=None, filterB=None, minScore=None, aCol=0,
                 bCol=1, scoreCol=-1, separator=None,
                 skipFirstLine=False, filterValues=None):
    if not filterValues:
        filterValues = list()
    index = dict()
    count = 0
    with open(fileName) as fp:
        line = fp.readline()
        if skipFirstLine:
            line = fp.readline()
        while line:
            ls = line.split(separator)
            skipEntry = False
            if separator is not None:
                aList = getIds(ls[aCol])
                bList = getIds(ls[bCol])
            else:
                aId = getCenterId(ls[aCol])
                bId = getCenterId(ls[bCol])
                aList = [aId]
                bList = [bId]
            if not skipEntry:
                validEntry = False
                for a in aList:
                    for b in bList:
                        skip = False
                        if a == "-" or b == "-":
                            skip = True
                        if filterA is not None and filterB is not None:
                            skip = not ((a in filterA and b in filterB) or (a in filterB and b in filterA))
                        for f in filterValues:
                            if len(ls) > f[0]:
                                columnEntry = ls[f[0]].lower()
                                searchEntry = f[1].lower()
                                if columnEntry.find(searchEntry) == -1:
                                    skip = True
                        if not skip:
                            name = getKey(a, b)
                            if name not in index:
                                validEntry = True
                                if scoreCol >= 0 and len(ls) > scoreCol:
                                    score = float(ls[scoreCol])
                                    skip = False
                                    if minScore is not None:
                                        if minScore > score:
                                            return index, count
                                    if not skip:
                                        index[name] = score
                                else:
                                    index[name] = 1.0
                if validEntry:
                    count = count + 1
            line = fp.readline()
    return index, count


def getMCC(prediction, positive, positiveCount, negative):
    sortedPrediction = sorted(prediction.items(), key=lambda x: x[1],
                              reverse=True)
    positiveTotal = positiveCount
    negativeTotal = len(negative)
    x = list([0])
    y = list([0])
    xMax = 0
    topCount = 0
    topMCC = 0.0
    topFP = 0.0
    topTP = 0.0
    topScore = 0.0
    tp = 0
    fp = 0
    count = 0
    for (name, score) in sortedPrediction:
        found = False
        if name in positive:
            found = True
            tp = tp + 1
        if name in negative:
            found = True
            fp = fp + 1
        fn = positiveTotal - tp
        tn = negativeTotal - fp
        denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
        yValue = getPercentage(tp, tp + fn)
        xValue = getPercentage(fp, fp + tn)
        if denom > 0.0:
            mcc = (tp * tn - fp * fn) / math.sqrt(denom)
            if mcc >= topMCC:
                topMCC = mcc
                topScore = score
                topCount = count
                topFP = xValue
                topTP = yValue
        if found:
            y.append(yValue)
            x.append(xValue)
            xMax = max(xValue, xMax)
        count = count + 1
    if len(sortedPrediction) > 0:
        print("Top ranking prediction %s." % str(sortedPrediction[0]))
        print("Total count of prediction set: %s (tp=%1.2f, fp=%1.2f)." % (topCount, topTP, topFP))
        print("Total count of positive set: %s." % len(positive))
        print("Total count of negative set: %s." % len(negative))
        print("Matthews-Correlation-Coefficient: %s at Score >= %s." % (round(topMCC, 2), topScore))
    return topMCC


def getNegativeSet(args, filterA, filterB, negativeRequired, jSize=5):
    # determine negative set
    print("Identifying non-interacting pairs...")
    negative = set()
    if args.negative and isfile(args.negative):
        # load from explicit file
        with open(args.negative) as file:
            for line in file:
                cols = line.split()
                nameA = cols[0]
                nameB = cols[1]
                key = getKey(nameA, nameB)
                if key not in negative:
                    negative.add(key)
    else:
        if args.region_a and args.region_b:
            locations = dict()
            regionA = args.region_a.lower()
            regionB = args.region_b.lower()
            locations[regionA] = list()
            locations[regionB] = list()
            regions = [regionA, regionB]
            print("Filtering regions %s" % str(regions))
            with open(args.locations) as locFile:
                for line in locFile:
                    searchKey = "SUBCELLULAR LOCATION"
                    searchPos = line.find(searchKey)
                    if searchPos != -1:
                        uniId = line.split()[0]
                        if uniId not in filterA and uniId not in filterB:
                            continue
                        locStart = searchPos + len(searchKey) + 1
                        locId = line[locStart:]
                        locId = re.sub(r"\s*{.*}\s*", "", locId)
                        locId = locId.replace(".", ",")
                        locId = locId.strip().lower()
                        filter_pos = locId.find("note=")
                        if filter_pos >= 0:
                            locId = locId[:filter_pos]
                        filter_pos = locId.find(";")
                        if filter_pos >= 0:
                            locId = locId[:filter_pos]
                        if locId:
                            locId = list(map(lambda x: x.strip(), locId.split(",")))
                            if (regionA in locId and regionB not in locId):
                                locations[regionA].append(uniId)
                            elif (regionA not in locId and regionB in locId):
                                locations[regionB].append(uniId)
            filterAList = sorted(locations[regionA])
            filterBList = sorted(locations[regionB])
        else:
            filterAList = sorted(filterA)
            filterBList = sorted(filterB)
        for i, j in randomPairs(len(filterAList), len(filterBList), jSize):
            nameA = filterAList[i]
            nameB = filterBList[j]
            key = getKey(nameA, nameB)
            if key not in negative:
                negative.add(key)
                negativeRequired = negativeRequired - 1
            if negativeRequired == 0:
                break
    return negative


def randomPairs(iLen, jLen, jSize):
    i = 0
    jStart = 0
    while i < iLen:
        jMax = min(jStart + jSize, jLen)
        for j in range(jStart, jMax):
            yield i, j
        i = i + 1
        if i == iLen and jMax < jLen:
            i = 0
            jStart = jStart + jSize + 1


def main(args):
    # load source files
    filterSets = getFilter(args.input)
    filterKeys = list(filterSets.keys())
    filterA = filterSets[filterKeys[0]]
    if len(filterKeys) > 1:
        filterB = filterSets[filterKeys[1]]
    else:
        filterB = filterA

    # identify biogrid filter options
    performance = dict()
    for methodReference in METHODS:

        # process biogrid database
        print("Loading positive set from BioGRID file (%s)..." % methodReference)
        filterValues = [[11, methodReference]]
        positive, positiveCount = getReference(args.biogrid, aCol=23, bCol=26,
                                               separator="\t", filterA=filterA,
                                               filterB=filterB, skipFirstLine=True,
                                               filterValues=filterValues)

        # estimate negative set
        negative = getNegativeSet(args, filterA, filterB, positiveCount)

        # evaluate other methods
        yValues = list()
        for method in METHODS:
            if methodReference != method:
                print("Method: %s" % method)
                filterValues = [[11, method]]
                prediction, _ = getReference(args.biogrid, aCol=23, bCol=26,
                                             separator="\t", filterA=filterA,
                                             filterB=filterB, skipFirstLine=True,
                                             filterValues=filterValues)
                mcc = getMCC(prediction, positive, positiveCount, negative)
                yValues.append(mcc)
            else:
                yValues.append(0.0)

        # add results to performance dication
        performance[methodReference] = yValues

        # get and append prediction results
        print("Loading prediction file...")
        prediction, _ = getReference(args.input, scoreCol=2, minScore=0.0)
        mcc = getMCC(prediction, positive, positiveCount, negative)
        performance[methodReference].append(mcc)

    # build yTicks
    yTicks = METHODS[:]
    yTicks.append("SPRING")

    # create plot
    print("Producing plot data...")
    print("Total count in prediction file: %d." % len(prediction))
    print("Total count in positive file: %d." % len(positive))
    df = pd.DataFrame(performance, index=yTicks)
    ax = df.plot.barh()
    ax.set_title(args.experiment)
    ax.set_xlabel("Matthews-Correlation Coefficient (MCC)")
    plt = ax.get_figure()
    plt.tight_layout()
    plt.savefig(args.output, format="png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create ROC plot.')
    parser.add_argument('-i', '--input', help='Input prediction file (2-columns).', required=True)
    parser.add_argument('-b', '--biogrid', help='BioGRID interaction database file', required=True)
    parser.add_argument('-l', '--locations', help='UniProt export table with subcellular locations', required=False)
    parser.add_argument('-ra', '--region_a', help='First subcellular location', required=False)
    parser.add_argument('-rb', '--region_b', help='Second subcellular location', required=False)
    parser.add_argument('-n', '--negative', help='Negative set (2-columns)', required=False)
    parser.add_argument('-t', '--throughput', help='Throughput (low/high)', required=False)
    parser.add_argument('-e', '--experiment', help='Experiment Title', required=False, default="Results")
    parser.add_argument('-o', '--output', help='Output (png)', required=True)
    args = parser.parse_args()
    main(args)
