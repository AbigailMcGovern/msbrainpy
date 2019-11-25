import os
import json
import nrrd
import numpy as np
import pandas as pd
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
from skimage.measure import regionprops

### for gene expression data (anatomical gene expression atlas Lein et al., 2014, structure unionize)
def parseExpressionToDF(fileName, queryList, saveName = None):
    with open(fileName, 'r') as jsonFile:
        expression = json.load(jsonFile)

    df = pd.DataFrame()   
    rows = []
    for i in range(len(expression['msg'])):
        row = []
        for query in queryList:
            if type(query) == str:
                row.append(expression['msg'][i][query])
            if type(query) == list:
                row.append(expression['msg'][i][query[0]][query[1]])
        rows.append(row)
        
    for j in range(len(queryList)):
        query = queryList[j]
        if type(query) ==list:
            colName = query[0]+'_'+query[1]
        if type(query) == str:
            colName = query
        values = []
        for i in range(len(rows)):
            row = rows[i]
            value = row[j]
            values.append(value)
        df.insert(j, colName, values)

        
    if saveName != None:
        df.to_csv(saveName, sep='\t')

    return df


### wrangling dataframes with structures ################################################################
def writeStructDFs(df, structName_list, filename, IDheader = 'id', exclude = None):
    oapi = OntologiesApi()
    # The 1 refers to the adult mouse brain atlas
    structure_graph = oapi.get_structures_with_sets([1])
    # clean_structures() removes unused fields
    structure_graph = StructureTree.clean_structures(structure_graph)
    # a class with methods for aceessing and using ontologies data
    tree = StructureTree(structure_graph)
    for structName in structName_list:
        towrite = getStructDF(df, structName, tree, IDheader = IDheader, exclude = exclude)
        writename = structName+filename
        towrite.to_csv(writename, sep='\t')

def writeStructDFsandMeans(df, structs, prefix, groupby = 'section_data_set_id'):
    started = False
    for struct in structs:
        structDF = getStructDF(df, struct, tree, IDheader = 'structure_id')
        saveName = prefix+'_'+struct+'.txt'
        structDF.to_csv(saveName, sep = '\t')
        meanStruct = structDF.groupby([groupby]).mean()
        saveName = prefix+'_'+struct+'_mean.txt'
        meanStruct.to_csv(saveName, sep = '\t')
        meanStruct = pd.read_csv(saveName, sep = '\t')
        names = []
        for i in range(len(meanStruct[groupby])):
            names.append(struct)
        meanStruct['region'] = names
        if started == False:
            out = meanStruct
            started = True
        if started == True:
            out = out.append(meanStruct)
    saveName = prefix+'_structure_means.txt'
    out.to_csv(saveName, sep = '\t')
    return out

def getStructDF(df, structName, tree, IDheader = 'id', exclude = None): # use with np.s_[:]
    structDict = tree.get_structures_by_name([structName])
    strid = structDict[0]['id']
    
    strids = df[IDheader]
    if exclude != None:
        strids = strids[exclude]
    
    strid_list = getSubStrids(strids, strid, verbose = False)
    
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

def getSubStrids(strids, strid, verbose = False):
    out = []
    for strida in strids:
        if verbose == True:
            is_desc = '' if tree.structure_descends_from(strida, strid) else ' not'
            print( '{0} is{1} in {2}'.format(name_map[strida], is_desc, name_map[strid]))
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

## Dealing with allen atlas labled data via acronym

def addICtxGroups(df, ac_header = 'structure_acronym', 
                  addParent_colNames = ['ctx_subregion', 'ctx_region'], writeOut = None):
    groups = []
    for ac in df[ac_header]:
        ac = Acronym(ac)
        groups.append(ac.group)
    df['group'] = groups
    
    if addParent_colNames != None: 
        add_reColumn(df, ac_header, '[A-Z]*[a-z]*', addParent_colNames[0], loc=0) 
        add_reColumn(df, ac_header, '[A-Z]*', addParent_colNames[1], loc=0)
        
    if writeOut != None:
        df.to_csv(writeOut, sep = '\t')
        
    return df
        
        
class Acronym:
    def __init__(self, string):
        self.string = string
        self.group = self.getGroup()
        #self.parents
        #self.children
        #self.structureSets
         
    def getGroup(self):
        ac = self.string
        
        somatomotor = ['MOp', 'SS']
        medial = ['PTLp', 'VISam', 'VISpm', 'RSP']
        temporal = ['AUD', 'PERI', 'TE', 'ECT']
        visual = ['VISal', 'VISl', 'VISp', 'VISpl']
        anterolateral = ['VISC', 'GU', 'AI']
        prefrontal = ['MOs', 'FRP', 'ACA', 'ORB', 'PL', 'ILA']
        
        groupDefs = {'somatomotor' : somatomotor, 'medial' : medial, 'temporal' : temporal, 
                     'anterolateral' : anterolateral, 'prefrontal' : prefrontal, 'visual' : visual}
        
        group = 'undefined'
        for key in groupDefs.keys():
            inGroup = False
            for value in groupDefs[key]:
                if ac.find(value) != -1:
                    inGroup = True
            if inGroup == True:
                group = key
            
        return group


## Regular expressions for dataframe wrangling

def add_reColumn(df, strColHeader, regex, reColHeader, loc=0):
    out = df
    pattern = re.compile(regex)
    reColumn = []
    for string in out[strColHeader]:
        matchObj = pattern.match(string)
        reColumn.append(matchObj[0])
    out.insert(loc, reColHeader, reColumn)
    return out

def add_reColumn_exc(df, strColHeader, regex, reColHeader, loc=0, exc = ['4', '5', '6']):
    out = df
    pattern = re.compile(regex)
    reColumn = []
    for string in out[strColHeader]:
        matchObj = pattern.match(string)
        exclude = False
        for cond in exc:
            if cond in string:
                exclude = True
        if exclude == False:
            reColumn.append(matchObj[0])
        if exclude == True:
            reColumn.append('NA')
    out.insert(loc, reColHeader, reColumn)
    return out



### base functions ###################################################################################
def writeAllCorticalStructures(tree):
    '''B fulcher allensdk repo'''
    nameMap = tree.get_name_map()
    cortexStructures = [name_map[i] for i in tree.descendant_ids([315])[0]]
    return cortexStructures

def writeSubStructures(structName, tree, filename):
    structDict = tree.get_structures_by_name([structName])
    strid = structDict[0]['id']
    nameMap = tree.get_name_map()
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



#### probably not especially useful #################################################################
def getScatterPoints(df_x, df_y, reHeader_x, reHeader_y, meanHeader_x, meanHeader_y, xHeader, yHeader, 
                     saveName = None, outHeader = 'parent_acronym', regex = '[A-Z]*[a-z]*'):
    '''
    FUNCTION: y is first obtained (i.e., measurment points) and equivalents from x are added to the dataframe
    '''
    y = getAverage_re(df_y, reHeader_y, meanHeader_y, regex = regex, outHeader = outHeader) 
    x = getAverage_re(df_x, reHeader_x, meanHeader_x, regex = regex, outHeader = outHeader)
    
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

    if saveName != None:
        out.to_csv(saveName, sep='\t')
    return out

def getAverage_re(df, reHeader, meanHeader, regex = '[A-Z]*[a-z]*', outHeader = 'parent_acronym'):
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