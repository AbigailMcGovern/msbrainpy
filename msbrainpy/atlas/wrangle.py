import os
import json
import nrrd
import re
import numpy as np
import pandas as pd
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
from skimage.measure import regionprops

# --------------------------------------------------- Globals ----------------------------------------------------------

oapi = OntologiesApi()
structure_graph = oapi.get_structures_with_sets([1])
structure_graph = StructureTree.clean_structures(structure_graph)
tree = StructureTree(structure_graph)
del oapi
name_map = tree.get_name_map()

# ------------------------------------------- Obtain CSV from json  ----------------------------------------------------

def parseExpressionToDF(fileName, queryList, saveName=None):
    """
    FUNCTION: pull a specified pd.DataFrame object out of a jsonn file obtained from an AllenSDK:
        RMA api-based AGEA structure unionize query.
    ARGUMENTS:
        fileName = json file path (str).
            Forces you to save the direct database output to json.
            Please keep all of this (regardless of how much you choose to obtain)!
        queryList = contains (1) the keys from which you wish to obtain data for each row.
            (2) if the data you wish to acess is nested provide a list of indexes and or keys by which the
            data can be accessed (ascending depth of key/index).
    DEPENDENCIES: Packages/modules: json, pandas as pd
    RETURNS: pd.DataFrame
    """
    with open(fileName, 'r') as jsonFile:
        expression = json.load(jsonFile)
    df = pd.DataFrame()
    rows = []
    for i in range(len(expression)):
        row = []
        for query in queryList:
            if type(query) == str:
                row.append(expression[i][query])
            if type(query) == list:
                val = expression[i][query[0]]
                for ii in range(len(query) - 1):
                    val = val[query[ii + 1]]
                row.append(val)
        rows.append(row)
    for j in range(len(queryList)):
        query = queryList[j]
        if type(query) == list:
            started = False
            for k in range(len(query)):
                if not started:
                    colName = str(query[0])
                    started = True
                elif started:
                    colName = colName + '_' + str(query[k])
        if type(query) == str:
            colName = query
        values = []
        for i in range(len(rows)):
            row = rows[i]
            value = row[j]
            values.append(value)
        df.insert(j, colName, values)
    if saveName != None:
        df.to_csv(saveName)
    return df

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
    structs_all = writeAllStructures(tree)
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


def saveStructFiles(df, prefix, structName='Isocortex', regex='[A-Z]*[a-z]*\d.*',
                    reNames=['sanslayers', 'layers'], IDheader='structure_id',
                    tree=tree, acronymHeader='structure_acronym', exclude=None):
    """
    save dataframes containing only information about specified structure,
    which should be specifed according to the structure ontology name.

    Parameters
    ----------
    df: pd.DataFrame
        Rows correspond to structures or brain sample x structure
    prefix: str
        Prefix with which to name 
    """
    name = prefix + '_' + structName
    name0 = name + '.csv'
    struct = getStructDF(df, structName, tree, IDheader=IDheader)
    struct.to_csv(name0)
    if regex != None:
        name1 = name + '_' + reNames[0]
        name2 = name + '_' + reNames[1]
        struct_sre = drop_reRows(struct, acronymHeader, regex='[A-Z]*[a-z]*\d.*', remove=True)
        struct_sre = dropParents(struct_sre, tree, IDheader=IDheader)
        struct_sre.to_csv(name1)
        struct_wre = drop_reRows(struct, acronymHeader, regex='[A-Z]*[a-z]*\d.*', remove=False)
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
        For shits and giggles or planned expansion and use, either way
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

# ------------------------------------------- Regex-based DF wrangling  ------------------------------------------------

def add_reColumn(df, strColHeader, regex, reColHeader, loc=0):
    out = df
    pattern = re.compile(regex)
    reColumn = []
    for string in out[strColHeader]:
        matchObj = pattern.match(string)
        reColumn.append(matchObj[0])
    out.insert(loc, reColHeader, reColumn)
    return out

def add_reColumn_exc(df, strColHeader, regex, reColHeader, loc=0, exc=['4', '5', '6']):
    out = df
    pattern = re.compile(regex)
    reColumn = []
    for string in out[strColHeader]:
        matchObj = pattern.match(string)
        exclude = False
        for cond in exc:
            if cond in string:
                exclude = True
        if not exclude:
            reColumn.append(matchObj[0])
        if exclude:
            reColumn.append('NA')
    out.insert(loc, reColHeader, reColumn)
    return out

def drop_reRows(df, strColHeader, regex='[A-Z]*[a-z]*\d.*', remove=True):
    """
    When applied to acronyms, this function will remove all acronyms that
    belong to a cortical layer. Invert to keep only layers (technically with # in acronym)
    """
    matches = []
    pattern = re.compile(regex)
    for string in df[strColHeader]:
        match = pattern.findall(string)
        if len(match) > 0:
            matches.append(True)
        if len(match) == 0:
            matches.append(False)
    out = dropBool(df, matches, todrop=remove)
    return out

def addLayer(df, acronymHeader):
    out = df.copy()
    regex = '[A-Z]*.*[a-z]*\d.*'
    pattern = re.compile(regex)
    d = '\d.*'
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

# --------------------------------------------- Statistical wrangling --------------------------------------------------

def meansSEMsSort(df, savename=None, groupby='structure_acronym'):
    """
    FUNCTION: Obtain a data frame containing means and standard errors for groups (e.g., anatomical structures)
    ARGUMENTS:
        df: pd.DataFrame or dictionary with save names as keys and dfs as values
        savename: if a df is supplied, name with which to save output
        groupby: column header of categorical variable for which to compute means
    DEPENDENCIES: Packages/modules/etc: pandas as pd
    Native: (1) amalgamateDFs(df_dict)
    1) joins two data frames with corresponding rows
    RETURNS: pd.DataFrame
    """
    if type(df) == dict:
        df_dict = df
        keys = list(df_dict.keys())
        assert len(savenames) == len(keys)
        out = {}
        for i in range(len(keys)):
            df = df_dict[keys[i]]
            name = key
            means = df.groupby([groupby]).mean()
            means.to_csv(name)
            means = pd.read_csv(name)
            sems = df.groupby([groupby]).sem()
            sems.to_csv(name)
            sems = pd.read_csv(name)
            sems = sems.drop([groupby], axis=1)
            result = amalgamateDFs({'mean': means, 'sem': sems})
            result.sort_values(groupby)
            result.to_csv(name)
            out[keys[i]] = result
    else:
        means = df.groupby([groupby]).mean()
        means.to_csv(savename)
        means = pd.read_csv(savename)
        sems = df.groupby([groupby]).sem()
        sems.to_csv(savename)
        sems = pd.read_csv(savename)
        sems = sems.drop([groupby], axis=1)
        out = amalgamateDFs({'mean': means, 'sem': sems})
        out.sort_values(groupby)
        out.to_csv(savename)
        print('You can supply a dictionary of save_name: df pairs to expedite bulk processing')
    return out

# ----------------------------------------- Row matching & column addition ---------------------------------------------

def matchAndAmalgamate(df0, df1, prefixes, groups=('structure_acronym', 'structure_acronym')):
    df0, df1 = matchRows(df0, df1, groups=groups)
    df = amalgamateDFs({prefixes[0]: df0, prefixes[1]: df1})
    return df

def matchRows(df0_, df1_, groups=('structure_acronym', 'structure_acronym')):
    """
    FUNCTION: Checks if two data frames contain identical rows on the basis of categories/IDs (i.e., groups).
        If not, unnecessary rows are recursively removed from longer data frame.
        Checks that rows are in correct sequence.
            NB: for DFs of means use 'structure_acronym' else create col of unique concatenated strings
            (strid_structSetID) - but this is a less likely scenario
    ARGUMENTS:
        df0_ = pd.DataFrame
        df1_ = pd.DataFrame
        groups = columns containing the categories or IDs that must match
    DEPENDENCIES: Packages/modules/etc: numpy as np, pandas as pd
    RETURNS: pd.DataFrame, pd.DataFrame
    """
    df0 = df0_.copy().sort_values(groups[0])
    df1 = df1_.copy().sort_values(groups[1])
    df_dict = {True: {'df': df0, 'index': 0}, False: {'df': df1, 'index': 1}}
    if set(df1[groups[1]]) == set(df0[groups[0]]):
        if len(df0) == len(df1):
            errors = 0
            vals_0 = [i for i in df0[groups[0]]]
            vals_1 = [i for i in df1[groups[0]]]
            for i in range(len(vals_0)):
                try:
                    assert vals_0[i] == vals_1[i]
                except AssertionError:
                    errors += 1
                    print('{} does not match {}'.format(df0[groups[0]][i], df1[groups[1]][i]))
            message = 'there were {} mismatches'.format(errors)
            if errors > 0:
                raise AssertionError(message)
            else:
                return df_dict[True]['df'], df_dict[False]['df']
    else:
        zeroLonger = len(set(df0[groups[0]])) >= len(set(df1[groups[1]]))
        zeroShorter = not zeroLonger
        bools = [ac in set(df_dict[zeroShorter]['df'][groups[df_dict[zeroShorter]['index']]])
                 for ac in df_dict[zeroLonger]['df'][groups[df_dict[zeroLonger]['index']]]]
        df_dict[zeroLonger]['df'] = dropBool(df_dict[zeroLonger]['df'], bools, todrop=False)
        number = 0
        rows = []
        for i in range(len(bools)):
            b = bools[i]
            if not b:
                number += 1
                rows.append(i)
        print('{} rows were removed from df {}'.format(number, df_dict[zeroLonger]['index']))
        print('The rows removed were: {}'.format(rows))
        df0, df1 = df_dict[True]['df'], df_dict[False]['df']
        return matchRows(df0, df1, groups)

def amalgamateDFs(df_dict):
    for key in df_dict.keys():
        col_names = list(df_dict[key].columns)
        new_names = []
        for col in col_names:
            new_names.append(col + '_' + key)
        dictionary = dict(zip(col_names, new_names))
        df_dict[key] = df_dict[key].rename(columns=dictionary)
    result = pd.concat([df_dict[key] for key in df_dict.keys()], axis=1, sort=False)
    return result


# ----------------------------------- Get Absolute File and Directory Paths -------------------------------------------
def get_image_addresses(path, age_id=None):
    
    """
    FUNCTION:Return a dictionary with all of the absolute file paths for images belonging to a gene directory in
        the format produced by msbrainpy.atlas.query.download_gene_images().
        {'entrez_id-<gene id>_<Acronym>' : {'plane_of_section-<1 or 2>' : {'age_id-<age id>_id-<id>' :
            [<file path 0>, <file path 1>, ...]}
    :param path: Absolute file path of the directory containing the gene directories (str)
    :param age_id: age id (15 = P56) (int or None)
    :return: dictionary of image paths (dict)

     e.g., image_address_dict = get_image_addresses('/Users/amcg0011/Data/InSituData', age_id=None)
    """ 
    gene_directories = get_gene_paths(path)
    image_address_dictionary = {}
    for gene in gene_directories:
        key = re.compile(r'entrez_id\d*_\w*').findall(gene)[0]
        section_directories = get_gene_section_paths(gene, age_id=age_id)
        image_address_dictionary[key] = section_directories
    return image_address_dictionary


def get_gene_paths(path):
    path_list = get_directory_paths(path, r'entrez_id\d*_\w*')
    return path_list


def get_gene_section_paths(gene_section_path, age_id=None):
    plane_of_section_directories = get_directory_paths(gene_section_path, r'plane_of_section-[1-2]')
    sections_dictionary = {}
    if age_id is not None:
        string = 'age_id-{}_id-\d*'.format(age_id)
    else:
        string = 'age_id-\d*_id-\d*'.format(age_id)
    pattern = re.compile(r'' + string)
    for plane_of_section in plane_of_section_directories:
        key = re.compile(r'plane_of_section-[1-2]').findall(plane_of_section)[0]
        paths = get_section_data_set_paths(plane_of_section, age_id=None)
        sections_dictionary[key] = {}
        for path in paths:
            key_1 = pattern.findall(path)[0]
            image_paths = get_directory_paths(path, r'image_id-\d*\.jpeg')
            sections_dictionary[key][key_1] = image_paths
    return sections_dictionary


def get_section_data_set_paths(gene_plane_section_path, age_id=None):
    if age_id is not None:
        string = 'age_id-{}*_id-\d*'.format(age_id)
        path_list = get_directory_paths(gene_plane_section_path, r'' + string)
    else:
        path_list = get_directory_paths(gene_plane_section_path, r'age_id-\d*_id-\d*')
    return path_list


def get_directory_paths(path, pattern):
    paths = []
    pattern = re.compile(pattern)
    for file in os.listdir(path):
        if pattern.fullmatch(file) is not None:
            full_path = os.path.join(path, file)
            paths.append(full_path)
    return sorted(paths)


# ----------------------------------------------- Basic functions ------------------------------------------------------

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

def removeRows(df, rm_list, rm_header):
    toDrop = []
    for i in range(len(df[rm_header])):
        val = df[rm_header][i]
        if val in rm_list:
            toDrop.append(i)
    df = df.drop(toDrop)
    return df

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

def rowMatches(df, toMatch, matchCol, todrop=False):
    matchBool = []
    for elem in df[matchCol]:
        b = elem in toMatch
        matchBool.append(b)
    out = dropBool(df, matchBool, todrop=todrop)
    return out

def dropBool(df, bools, todrop=False):
    out = df.copy()
    out['bools'] = bools
    out = out[out['bools'] != todrop]
    out = out.drop(['bools'], axis=1)
    return out


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
        for i in range(len(meanStruct[groupby])):
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


