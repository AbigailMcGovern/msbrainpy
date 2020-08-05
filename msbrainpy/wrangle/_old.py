import re
import pandas as pd
from ._structure_helpers import writeAllStructures, getStructDF, tree
# ------------------------------------------------- Deprecated? --------------------------------------------------------

def addAcronymCol_name(df, namesCol, tree, header='acronym', loc=0):
    """
    this was specifically designed for the ero et al., 2018 data set, in which ','s were removed
    from the names of regions. This means that the names cannot be used to directly query the
    ontology structure tree to get the acronyms. BTW a function exists to get IDs but depends on the
    presence of acronyms. This module typically uses acronyms unless IDs are required (e.g., decends from?).
    Therefore, it is better to have both. Additionally, I have officially used a comical number of characters
    in this docstring without technically defining the function, its arguments, or its output properly.
    Also, the wording is less than eloquent.
    Please enjoy:
    """
    acronyms = []
    structs_all = writeAllStructures(tree)
    for name in df[namesCol]:
        try:
            structName = getSearchWord(name, structs_all, sep=r', ')
        except:
            acronym = None
        try:
            structDict = tree.get_structures_by_name([structName])[0]
        except KeyError:
            acronym = None
        acronym = structDict['acronym']
        acronyms.append(acronym)

    df.insert(loc, header, acronyms)
    return df

def addAcronymCol_id(df, idCol, tree, header='structure_acronym', loc=0):
    acronyms = []
    for ID in df[idCol]:
        try:
            structDict = tree.get_structures_by_id([ID])[0]
        except KeyError:
            acronym = None
        acronym = structDict['acronym']
        acronyms.append(acronym)

    df.insert(loc, header, acronyms)
    return df

def addIDCol_acronym(df, acronymCol, tree, header='structure_id', loc=0):
    IDs = []
    for acronym in df[acronymCol]:
        try:
            structDict = tree.get_structures_by_acronym([acronym])[0]
        except KeyError:
            ID = None
        ID = structDict['id']
        IDs.append(ID)

    df.insert(loc, header, IDs)
    return df

def addIDCol_name(df, nameCol, tree, header='structure_id', loc=0):
    IDs = []
    for name in df[nameCol]:
        try:
            structName = getSearchWord(name, structs_all, sep=r', ')
        except:
            ID = None
        try:
            structDict = tree.get_structures_by_name([structName])[0]
            ID = structDict['id']
        except KeyError:
            ID = None

        IDs.append(ID)

    df.insert(loc, header, IDs)
    return df

def addNameCol_acronym(df, acronymCol, tree, header='structure_name', loc=0):
    names = []
    for acronym in df[acronymCol]:
        try:
            structDict = tree.get_structures_by_acronym([acronym])[0]
            name = structDict['id']
        except KeyError:
            name = None

        names.append(name)

    df.insert(loc, header, names)
    return df

def addNameCol_id(df, idCol, tree, header='structure_name', loc=0):
    names = []
    for ID in df[idCol]:
        try:
            structDict = tree.get_structures_by_id([ID])[0]
            name = structDict['acronym']
        except KeyError:
            name = None

        names.append(name)

    df.insert(loc, header, names)
    return df

def getScatterPoints(df_x, df_y, reHeader_x, reHeader_y, meanHeader_x, meanHeader_y, xHeader, yHeader,
                     saveName=None, outHeader='parent_acronym', regex='[A-Z]*[a-z]*'):
    '''
    FUNCTION: y is first obtained (i.e., measurment points) and equivalents from x are added to the dataframe
    '''
    y = getAverage_re(df_y, reHeader_y, meanHeader_y, regex=regex, outHeader=outHeader)
    x = getAverage_re(df_x, reHeader_x, meanHeader_x, regex=regex, outHeader=outHeader)

    xlist = []
    ylist = []
    parents = []
    for j in range(len(y[outHeader])):
        parent = y[outHeader][j]
        for i in range(len(x[outHeader])):
            if x[outHeader][i] == parent:
                xlist.append(x[meanHeader_x][i])
                ylist.append(y[meanHeader_y][j])
                parents.append(parent)

    out = pd.DataFrame()
    out[outHeader] = parents
    out[xHeader] = xlist
    out[yHeader] = ylist

    if saveName is not None:
        out.to_csv(saveName, sep='\t')
    return out


def getAverage_re(df, reHeader, meanHeader, regex='[A-Z]*[a-z]*', outHeader='parent_acronym'):
    parent = []
    pattern = re.compile(regex)
    for string in df[reHeader]:
        matchObj = pattern.match(string)
        parent.append(matchObj[0])
    df[outHeader] = parent
    means = df.groupby(outHeader)[meanHeader].mean()
    out = pd.DataFrame()
    parents0 = []
    for index in means.index:
        parents0.append(index)

    out[outHeader] = parents0
    out[meanHeader] = means.values
    return out

def writeStructDFsandMeans(df, structs, prefix, groupby='section_data_set_id'):
    started = False
    for struct in structs:
        structDF = getStructDF(df, struct, tree, IDheader='structure_id')
        saveName = prefix + '_' + struct + '.txt'
        structDF.to_csv(saveName, sep='\t')
        meanStruct = structDF.groupby([groupby]).mean()
        saveName = prefix + '_' + struct + '_mean.txt'
        meanStruct.to_csv(saveName, sep='\t')
        meanStruct = pd.read_csv(saveName, sep='\t')
        names = []
        for _ in range(len(meanStruct[groupby])):
            names.append(struct)
        meanStruct['region'] = names
        if not started:
            out = meanStruct
            started = True
        if started:
            out = out.append(meanStruct)
    saveName = prefix + '_structure_means.txt'
    out.to_csv(saveName, sep='\t')
    return out