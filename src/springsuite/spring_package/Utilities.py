from os.path import isfile


def validateIdentifier(identifier):
    if len(identifier) < 6 or identifier[4:5] != "_":
        raise Exception("Invalid list entry (`PDB_CHAIN`): %s." % identifier)


def getId(identifier):
    identifier = identifier.strip()
    validateIdentifier(identifier)
    return identifier[:4].upper() + identifier[4:6]


def getChain(identifier):
    validateIdentifier(identifier)
    pdbChain = identifier[5:6]
    return pdbChain


def getName(identifier):
    if identifier.startswith("pdb"):
        pdb = identifier[3:7]
    else:
        pdb = identifier[:4]
    return pdb.lower()


def getCrossReference(crossReferenceFile, allPartners=False):
    crossReference = dict()
    crossCount = 0
    with open(crossReferenceFile) as file:
        for line in file:
            columns = line.split()
            if len(columns) < 2:
                raise Exception("Invalid Cross Reference Entry %s." % line)
            core = columns[0]
            partner = columns[1]
            if len(columns) < 4:
                templates = [core, partner]
            else:
                templates = [columns[2], columns[3]]
            if core not in crossReference:
                crossReference[core] = dict(partners=list(), templates=list())
            if allPartners or partner not in crossReference[core]["partners"]:
                crossReference[core]["partners"].append(partner)
                crossReference[core]["templates"].append(templates)
                crossCount = crossCount + 1
    print("Identified %s reference interactions." % crossCount)
    return crossReference


def getTemplates(hhrFile, minScore=10, nTemplates=5):
    result = dict()
    topTemplates = list()
    if isfile(hhrFile):
        with open(hhrFile) as file:
            for index, line in enumerate(file):
                if index > 8:
                    if not line.strip():
                        break
                    templateId = line[4:10]
                    templateScore = float(line[57:63])
                    if templateScore > minScore:
                        if len(topTemplates) < nTemplates:
                            topTemplates.append(templateId)
                        result[templateId] = templateScore
    return topTemplates, result
