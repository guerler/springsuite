import subprocess
from os import mkdir
from os.path import basename, isdir

from spring_package.Alignment import Alignment
from spring_package.DBKit import DBKit
from spring_package.Energy import Energy
from spring_package.Molecule import Molecule
from spring_package.Utilities import getChain, getCrossReference, getName, getTemplates


def createPDB(identifier, pdbDatabase, outputName, zipped=None):
    pdb = getName(identifier)
    pdbDatabaseId = "pdb%s.ent" % pdb
    return pdbDatabase.createFile(pdbDatabaseId, outputName, zipped=zipped)


def createMonomer(resultFile, identifier, pdbDatabase, outputName, zipped=None):
    print("Building model with: %s." % identifier)
    if not createPDB(identifier, pdbDatabase, outputName, zipped=zipped):
        print("Template not found in database [%s]" % identifier)
        return False
    template = Molecule(outputName)
    pdbChain = getChain(identifier)
    if pdbChain not in template.calpha:
        print("Chain not found in template [%s]" % pdbChain)
        return False
    chain = template.calpha[pdbChain]
    alignment = Alignment(resultFile)
    alignment.createModel(chain)
    template.saveChain(pdbChain, outputName)
    try:
        subprocess.run(["pulchra", outputName], check=True)
    except subprocess.CalledProcessError as e:
        print(str(e))
        return False
    return True


def TMalign(fileA, fileB, tmName="temp/tmalign"):
    try:
        tmResult = open("%s.out" % tmName, "w")
        subprocess.run(["TMalign", fileA, fileB, "-m", "%s.mat" % tmName], check=True, stdout=tmResult)
        tmResult.close()
    except subprocess.CalledProcessError as e:
        raise Exception(str(e))
    rotmat = list()
    with open("%s.mat" % tmName) as file:
        line = next(file)
        line = next(file)
        for i in range(3):
            line = next(file)
            rotmatLine = line.split()
            rotmatLine = list(map(lambda x: float(x), rotmatLine))
            rotmatLine = [rotmatLine[2], rotmatLine[3], rotmatLine[4], rotmatLine[1]]
            rotmat.append(rotmatLine)
    with open("%s.out" % tmName) as file:
        for i in range(18):
            line = next(file)
        try:
            tmscore = float(line[9:17])
            line = next(file)
            tmscore = max(tmscore, float(line[9:17]))
        except Exception:
            raise Exception("TMalign::Failed to retrieve TMscore.")
    molecule = Molecule(fileA)
    for atom in molecule.atoms:
        molecule.applyMatrix(atom, rotmat)
    return tmscore, molecule


def TMalignAlignment(bioMolecule, templateChain, tmName="temp/tmalign"):
    aligned = list()
    with open("%s.out" % tmName) as file:
        for i in range(23):
            line = next(file)
        try:
            modelAlign = line
            line = next(file)
            alignment = line
            line = next(file)
            templateAlign = line
        except Exception:
            raise Exception("TMalign::Failed to retrieve TMalign results.")
    templateResidues = sorted(bioMolecule.calpha[templateChain].values(), key=lambda item: item["residueNumber"])
    templateIndex = 0
    for i in range(len(alignment)):
        t = templateAlign[i]
        if alignment[i] in [":", "."]:
            templateResidue = templateResidues[templateIndex]
            templateResidue["alignedResidue"] = modelAlign[i]
            aligned.append(templateResidue)
        if t != "-":
            templateIndex = templateIndex + 1
    return aligned


def getFrameworks(aTemplates, bTemplates, crossReference, minScore, maxTries):
    templateHits = list()
    for aTemplate in aTemplates:
        if aTemplate in crossReference:
            partners = crossReference[aTemplate]["partners"]
            templates = crossReference[aTemplate]["templates"]
            for pIndex in range(len(partners)):
                pTemplate = partners[pIndex]
                templatePair = templates[pIndex]
                if pTemplate in bTemplates:
                    minZ = min(aTemplates[aTemplate], bTemplates[pTemplate])
                    templateHits.append(dict(templatePair=templatePair, score=minZ))
    templateList = sorted(templateHits, key=lambda item: item["score"], reverse=True)
    print("Found %d templates." % len(templateList))
    for templateHit in templateList:
        if templateHit["score"] < minScore or maxTries == 0:
            break
        maxTries = maxTries - 1
        yield templateHit["templatePair"], templateHit["score"]


def createModel(args):
    print("SPRING - Complex Model Creation")
    aName = basename(args.a_hhr)
    bName = basename(args.b_hhr)
    print("Sequence A: %s" % aName)
    print("Sequence B: %s" % bName)
    aTops, aTemplates = getTemplates(args.a_hhr)
    bTops, bTemplates = getTemplates(args.b_hhr)
    if not isdir("temp"):
        mkdir("temp")
    outputName = args.output
    pdbDatabase = DBKit(args.index, args.database)
    crossReference = getCrossReference(args.cross)
    interfaceEnergy = Energy()
    maxScore = -9999
    maxInfo = None
    minScore = float(args.minscore)
    maxTries = int(args.maxtries)
    for aTop in aTops:
        if not createMonomer(args.a_hhr, aTop, pdbDatabase, "temp/monomerA.pdb", zipped=args.zipped):
            print("Warning: Failed to determine monomer model for %s." % args.a_hhr)
            continue
        for bTop in bTops:
            if not createMonomer(args.b_hhr, bTop, pdbDatabase, "temp/monomerB.pdb", zipped=args.zipped):
                print("Warning: Failed to determine monomer model for %s." % args.b_hhr)
                continue
            for [aTemplate, bTemplate], zscore in getFrameworks(aTemplates, bTemplates, crossReference, minScore=minScore, maxTries=maxTries):
                print("Evaluating Complex Template: %s." % aTemplate)
                templateFile = "temp/template.pdb"
                createPDB(aTemplate, pdbDatabase, templateFile, zipped=args.zipped)
                templateMolecule = Molecule(templateFile)
                aTemplateChain = getChain(aTemplate)
                bTemplateChain = getChain(bTemplate)
                if aTemplateChain == bTemplateChain:
                    bTemplateChain = "%s_0" % bTemplateChain
                print("Evaluating chain %s and %s..." % (aTemplate, bTemplate))
                biomolFound = False
                for biomolNumber in templateMolecule.biomol:
                    bioMolecule = templateMolecule.createUnit(biomolNumber)
                    if (len(bioMolecule.calpha.keys()) > 1
                    and aTemplateChain in bioMolecule.calpha
                    and bTemplateChain in bioMolecule.calpha):
                        print("Evaluating biomolecule %i..." % biomolNumber)
                        bioMolecule.saveChain(aTemplateChain, "temp/template_0.pdb")
                        bioMolecule.saveChain(bTemplateChain, "temp/template_1.pdb")
                        try:
                            coreScore, coreMolecule = TMalign("temp/monomerA.rebuilt.pdb", "temp/template_0.pdb")
                            coreAligned = TMalignAlignment(bioMolecule, aTemplateChain)
                            partnerScore, partnerMolecule = TMalign("temp/monomerB.rebuilt.pdb", "temp/template_1.pdb")
                            partnerAligned = TMalignAlignment(bioMolecule, bTemplateChain)
                        except Exception as e:
                            print("Warning: Failed TMalign [%s]." % bTemplateChain)
                            print(str(e))
                            continue
                        biomolFound = True
                        print("  zscore:\t%5.2f" % zscore)
                        tmscore = min(coreScore, partnerScore)
                        print("  tmscore:\t%5.2f" % tmscore)
                        energy = -interfaceEnergy.get(coreAligned, partnerAligned)
                        print("  energy:\t%5.2f" % energy)
                        clashes = interfaceEnergy.getClashes(coreMolecule, partnerMolecule)
                        print("  clashes:\t%5.2f" % clashes)
                        springscore = tmscore + energy * args.wenergy
                        print("  springscore:\t%5.2f" % springscore)
                        if springscore > maxScore and clashes < args.maxclashes:
                            maxScore = springscore
                            maxInfo = dict(aTemplate=aTemplate, bTemplate=bTemplate, springscore=springscore, tmscore=tmscore, energy=energy, clashes=clashes, zscore=zscore)
                            coreMolecule.save(outputName, chainName="0")
                            partnerMolecule.save(outputName, chainName="1", append=True)
                            if args.showtemplate == "true":
                                bioMolecule.save(outputName, append=True)
                    if biomolFound:
                        break
    if maxInfo is not None:
        print("Final Model:")
        for key in maxInfo:
            print("  %s:\t%s" % (key, maxInfo[key]))
        print("Completed.")
    else:
        print("Warning: Failed to determine model.")
    return maxInfo
