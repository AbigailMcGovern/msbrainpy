import numpy as np
import pandas as pd
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
from ._structure_helpers import getSubStrids, OntologiesApi, tree, getStructDF
from ._helpers import dropBool
from ._re_based import add_reColumn, drop_reRows


# ------------- Wrangling DFs on the basis of atlas structures ----------------

def writeStructDFs(df, structName_list, filename, IDheader='id', exclude=None):
    oapi = OntologiesApi()
    # The 1 refers to the adult mouse brain atlas
    structure_graph = oapi.get_structures_with_sets([1])
    # clean_structures() removes unused fields
    structure_graph = StructureTree.clean_structures(structure_graph)
    # a class with methods for aceessing and using ontologies data
    tree = StructureTree(structure_graph)
    for structName in structName_list:
        towrite = getStructDF(df, structName, tree, IDheader=IDheader, exclude=exclude)
        writename = structName + filename
        towrite.to_csv(writename, sep='\t')



# ------------------------------------ Adding Structure tree info to DF ------------------------------------------------

def addFromStructureTree(df, useCol, toAdd, toUse, tree):
    """
     write a docstring. I mean it.
    """
    if toAdd == toUse:
        message = 'The toAdd ({}) argument cannot be the same as the toUse argument({})'.format(toAdd, toUse)
        print('''Please choose the type of column to add (acronym, id, name, ...) \n 
            and the name (useCol) and type of column containing (acronym, id, or name) \n 
            the information with which to query the tree''')
        raise ValueError(message)
    header = 'structure_' + toAdd
    if toUse == 'id':
        out = addCol_id(df, useCol, toAdd, toUse, tree, header, loc=0)
    if toUse == 'name':
        out = addCol_name(df, useCol, toAdd, toUse, tree, header, loc=0)
    if toUse == 'acronym':
        out = addCol_acronym(df, useCol, toAdd, toUse, tree, header, loc=0)
    return out

def addCol_id(df, idCol, toAdd, toUse, tree, header, loc=0):
    additions = []
    for ID in df[idCol]:
        try:
            structDict = tree.get_structures_by_id([ID])[0]
            addition = structDict[toAdd]
        except KeyError:
            addition = None
        additions.append(addition)
    df.insert(loc, header, additions)
    return df

def addCol_acronym(df, acronymCol, toAdd, toUse, tree, header, loc=0):
    additions = []
    for acronym in df[acronymCol]:
        try:
            structDict = tree.get_structures_by_acronym([acronym])[0]
            addition = structDict[toAdd]
        except KeyError:
            addition = None
        additions.append(addition)
    df.insert(loc, header, additions)
    return df

def addCol_name(df, nameCol, toAdd, toUse, tree, header, loc=0):
    additions = []
    for name in df[nameCol]:
        try:
            structName = getSearchWord(name, structs_all, sep=r', ')
        except:
            addition = None
        try:
            structDict = tree.get_structures_by_name([structName])[0]
            addition = structDict[toAdd]
        except KeyError:
            addition = None

        additions.append(addition)

    df.insert(loc, header, additions)
    return df

def dropParents(df, tree, IDheader='structure_id'):
    """
     FUNCTION: remove any partent structures from a dataframe using the structure_id column.
         is dependent on the dropBool() function.
     ARGS:
         df : data to clean (pandas.DataFrame)
         tree : instance of StructureTree (allensdk.core.structure_tree.StructureTree())
         IDheader : header of the id column (str)
     DEPENDENCIES: (1) dropBool(...)
         1) removes specified bool from DataFrame when supplied list of bool of len(df) (function)
     RETURNS: copy df without the parent structures (pandas.DataFrame)
     """
    parents = []
    uniqueElem = np.unique([strid for strid in df[IDheader]])
    print('there are {} unique IDs in the data frame'.format(len(uniqueElem)))

    for strid_a in uniqueElem:
        values = []
        matches = []
        for strid_b in uniqueElem:
            if strid_a != strid_b:
                parent = tree.structure_descends_from(strid_b, strid_a)
                values.append(parent)
                if parent:
                    matches.append(strid_b)
        if len(matches) != 0:
            print('the children of {} are as follows'.format(strid_a))
            print(matches)
        a = True in values
        if a:
            parents.append(strid_a)
    parentBool = []
    for strid in df[IDheader]:
        b = strid in parents
        parentBool.append(b)
    out = dropBool(df, parentBool, todrop=True)
    return out


def saveStructFiles(df, prefix, structName='Isocortex', regex=r'[A-Z]*[a-z]*\d.*',
                    reNames=['sanslayers', 'layers'], IDheader='structure_id',
                    tree=tree, acronymHeader='structure_acronym', exclude=None):
    name = prefix + '_' + structName
    name0 = name + '.csv'
    struct = getStructDF(df, structName, tree, IDheader=IDheader)
    struct.to_csv(name0)
    if regex != None:
        name1 = name + '_' + reNames[0]
        name2 = name + '_' + reNames[1]
        struct_sre = drop_reRows(struct, acronymHeader, regex=r'[A-Z]*[a-z]*\d.*', remove=True)
        struct_sre = dropParents(struct_sre, tree, IDheader=IDheader)
        struct_sre.to_csv(name1)
        struct_wre = drop_reRows(struct, acronymHeader, regex=r'[A-Z]*[a-z]*\d.*', remove=False)
        struct_wre = dropParents(struct_wre, tree, IDheader=IDheader)
        struct_wre.to_csv(name2)
        return struct_sre, struct_wre
    else:
        return struct

# ------------------------------------ Adding acronym-associated info to DF --------------------------------------------

def addICtxGroups(df, ac_header='structure_acronym',
                  addParent_colNames=('ctx_subregion', 'ctx_region'), writeOut=None):
    groups = []
    for ac in df[ac_header]:
        ac = Acronym(ac)
        groups.append(ac.group)
    df['group'] = groups
    if addParent_colNames is not None:
        add_reColumn(df, ac_header, '[A-Z]*[a-z]*', addParent_colNames[0], loc=0)
        add_reColumn(df, ac_header, '[A-Z]*', addParent_colNames[1], loc=0)
    if writeOut != None:
        df.to_csv(writeOut, sep='\t')
    return df

class Acronym:
    def __init__(self, string):
        """
        For shits and giggles - slightly evil class
        ARNUMENTS:
            string: acronym of choice
        """
        self.string = string
        self.group = self.getGroup()
        # self.otherInfo?

    def getGroup(self):
        ac = self.string
        somatomotor = ['MOp', 'SS']
        medial = ['PTLp', 'VISam', 'VISpm', 'RSP']
        temporal = ['AUD', 'PERI', 'TE', 'ECT']
        visual = ['VISal', 'VISl', 'VISp', 'VISpl']
        anterolateral = ['VISC', 'GU', 'AI']
        prefrontal = ['MOs', 'FRP', 'ACA', 'ORB', 'PL', 'ILA']
        groupDefs = {'somatomotor': somatomotor, 'medial': medial, 'temporal': temporal,
                     'anterolateral': anterolateral, 'prefrontal': prefrontal, 'visual': visual}
        group = 'undefined'
        for key in groupDefs.keys():
            inGroup = False
            for value in groupDefs[key]:
                if ac.find(value) != -1:
                    inGroup = True
            if inGroup:
                group = key
        return group
