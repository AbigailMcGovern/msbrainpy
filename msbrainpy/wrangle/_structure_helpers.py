import re
import numpy as np
import pandas as pd
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree


oapi = OntologiesApi()
structure_graph = oapi.get_structures_with_sets([1])
structure_graph = StructureTree.clean_structures(structure_graph)
tree = StructureTree(structure_graph)
del oapi
name_map = tree.get_name_map()

def writeAllStructures(tree):
    nameMap = tree.get_name_map()
    structures = [nameMap[i] for i in tree.descendant_ids([997])[0]]
    return structures

def writeAllCorticalStructures(tree):
    """B fulcher allensdk repo"""
    name_map = tree.get_name_map()
    cortexStructures = [name_map[i] for i in tree.descendant_ids([315])[0]]
    return cortexStructures

def writeSubStructures(structName, tree, filename):
    structDict = tree.get_structures_by_name([structName])
    strid = structDict[0]['id']
    name_map = tree.get_name_map()
    structures = [name_map[i] for i in tree.descendant_ids([strid])[0]]
    ids = [i for i in tree.descendant_ids([strid])[0]]
    df = pd.DataFrame()
    df['structure_name'] = structures
    df['structure_id'] = ids
    df.to_csv(filename, sep='\t')
    return df

def multiSubStids(queryStructs, tree):
    query_stids = []
    for struct in queryStructs:
        strids = subStructures_id(struct, tree)
        for strid in strids:
            query_stids.append(strid)
    return query_stids

def subStructures_id(structName, tree):
    structDict = tree.get_structures_by_name([structName])
    strid = structDict[0]['id']
    ids = [i for i in tree.descendant_ids([strid])[0]]
    return ids

def getMsTree():
    oapi = OntologiesApi()
    # The 1 refers to the adult mouse brain atlas
    structure_graph = oapi.get_structures_with_sets([1])
    # clean_structures() removes unused fields
    structure_graph = StructureTree.clean_structures(structure_graph)
    # a class with methods for aceessing and using ontologies data
    tree = StructureTree(structure_graph)
    return tree

def getSubStrids(strids, strid, verbose=False):
    out = []
    for strida in strids:
        if verbose == True:
            is_desc = '' if tree.structure_descends_from(strida, strid) else ' not'
            print('{0} is{1} in {2}'.format(name_map[strida], is_desc, name_map[strid]))
        if tree.structure_descends_from(strida, strid) == True:
            out.append(strida)
    return out

def getSearchWord(string, toSearch, sep=r', '):
    """
    The Ero et al., 2018, Front. Neuroinfo. data has no commas in structure names and this
    returns an acceptable search-term when fed a list of Allen Atlas structure names
    """
    arr = []
    for op in toSearch:
        searchList = re.split(sep, op)
        successCount = 0
        for item in searchList:
            if item in string:
                successCount += 1
        arr.append(successCount)
    arr = np.array(arr)
    i = np.argmax(arr)
    return toSearch[i]

def getStructDF(df, structName, tree, IDheader='id', exclude=None):  # use with np.s_[:]
    """
    FUNCTION: obtain decendents of a particular structure from a DataFrame
    ARGUMENTS:
        df = data (pandas.DataFrame)
        structName = name of the parent structure for which data should be obtained (str)
        tree = instance of StructureTree (allensdk.core.structure_tree.StructureTree())
        IDheader = header of the id column (str)
        exclude = ids to exclude. id col is treated as list and spliced using np.s_[...] (np.s_[:])
    DEPENDENCIES: (1) getSubStrids(..); (2) dropBool(...)
        1) obtain list of structures that decend from a given parent structure from a list of structures (function)
        2) removes specified bool from DataFrame when supplied list of bool of len(df) (function)
    RETURNS: DataFrame
    """
    structDict = tree.get_structures_by_name([structName])
    strid = structDict[0]['id']

    strids = df[IDheader]
    if exclude is not None:
        strids = strids[exclude]
    strid_list = getSubStrids(strids, strid, verbose=False)
    bools = []
    for ID in strids:
        if ID in strid_list:
            bools.append(True)
        else:
            bools.append(False)
    out = df
    out['bools'] = bools
    out = out[out['bools'] == True]
    out = out.drop(['bools'], axis=1)
    return out

def addLayer(df, acronymHeader):
    out = df.copy()
    regex = r'[A-Z]*.*[a-z]*\d.*'
    pattern = re.compile(regex)
    d = r'\d.*'
    digit = re.compile(d)
    layers = []
    for string in df[acronymHeader]:
        try:
            match = pattern.match(string)[0]
            layer = digit.findall(match)[0]
            layers.append(layer)
        except TypeError:
            print('A layer could not be identified in {}'.format(string))
            layers.append('NA')
    out['cortical_layer'] = layers
    return out