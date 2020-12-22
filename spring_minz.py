#! /usr/bin/env python3
import argparse

from spring_package.Utilities import getCrossReference, getTemplates


def main(args):
    minScore = args.minscore
    logFile = open(args.log, 'w')
    targets = list()
    targetPath = args.targetpath.rstrip("/")
    hhrResults = dict()
    with open(args.targetlist) as file:
        for line in file:
            name = line.strip()
            targets.append(name)
    print("Loaded %s target names from `%s`." % (len(targets), args.targetlist))
    for targetName in targets:
        targetFile = "%s/%s" % (targetPath, targetName)
        hhrResults[targetName] = getTemplates(targetFile, minScore)
    if args.inputlist:
        inputs = list()
        inputPath = args.inputpath.rstrip("/")
        with open(args.inputlist) as file:
            for line in file:
                name = line.strip()
                inputs.append(name)
        print("Loaded %s input names from `%s`." % (len(inputs), args.inputlist))
        for inputName in inputs:
            if inputName not in hhrResults:
                inputFile = "%s/%s" % (inputPath, inputName)
                hhrResults[inputName] = getTemplates(inputFile, minScore)
    else:
        inputs = targets
    print("Loaded hhr results for %s entries." % len(hhrResults.keys()))
    crossReference = getCrossReference(args.cross)
    print("Loaded cross reference from `%s`." % args.cross)
    interactions = dict()
    for targetName in targets:
        matchScores(hhrResults=hhrResults,
                    targetName=targetName,
                    inputs=inputs,
                    crossReference=crossReference,
                    minScore=args.minscore,
                    logFile=logFile,
                    interactions=interactions)
    interactions = sorted(interactions.values(), key=lambda item: item["minZ"], reverse=True)
    with open(args.output, 'w') as output_file:
        for entry in interactions:
            output_file.write("%s\t%s\t%s\t%s\n" % (entry["targetName"],
                              entry["inputName"], entry["minZ"],
                              entry["minInfo"]))
    logFile.close()


def matchScores(hhrResults, targetName, inputs, crossReference, minScore, logFile, interactions):
    if targetName not in hhrResults:
        print("Target not found `%s`" % targetName)
    else:
        targetTop, targetHits = hhrResults[targetName]
        print("Evaluating %s." % targetName)
        logFile.write("Evaluating %s.\n" % targetName)
        logFile.flush()
        for inputName in inputs:
            if inputName in hhrResults:
                inputTop, inputHits = hhrResults[inputName]
                minZ = 0
                minInfo = ""
                for t in targetHits:
                    if t in crossReference:
                        partners = crossReference[t]["partners"]
                        for p in partners:
                            if p in inputHits:
                                score = min(targetHits[t], inputHits[p])
                                if score > minZ:
                                    minZ = score
                                    minInfo = "%s\t%s\t%s\t%s" % (targetTop, inputTop, t, p)
                if minZ > minScore:
                    if targetName > inputName:
                        interactionKey = "%s_%s" % (targetName, inputName)
                    else:
                        interactionKey = "%s_%s" % (inputName, targetName)
                    if interactionKey in interactions:
                        if interactions[interactionKey]["minZ"] >= minZ:
                            continue
                    interactions[interactionKey] = dict(targetName=targetName,
                                                        inputName=inputName,
                                                        minZ=minZ, minInfo=minInfo)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script identifies interactions by detecting matching HH-search results.')
    parser.add_argument('-tl', '--targetlist', help='Text file containing identifiers.', required=True)
    parser.add_argument('-tp', '--targetpath', help='Directory containing `hhr` files', required=True)
    parser.add_argument('-il', '--inputlist', help='Text file containing identifiers.', required=False)
    parser.add_argument('-ip', '--inputpath', help='Directory containing `hhr` files', required=False)
    parser.add_argument('-c', '--cross', help='PDB Cross Reference', required=True)
    parser.add_argument('-o', '--output', help='Output file containing min-Z scores', required=True)
    parser.add_argument('-l', '--log', help='Log file', required=True)
    parser.add_argument('-m', '--minscore', help='min-Z score threshold', type=int, default=25)
    args = parser.parse_args()
    main(args)
