import os
import numpy as np
import pandas as pd
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
from collections import defaultdict
#from msbrainpy import #### PLEASE ADD A WRAPPER FOR RESAMPLING & MULTIPLE TRANSFORMATIONS -> AMBA

### Density analysis (i.e., average voxel density img @ 25 um with annotation file)
#def getAreaDensity - should this be here or in amba
#def mapDensity():
#def 

def getProp(props, IDs, key):    
    out = []
    for ID in IDs:
        found = False
        for prop in props:
            if prop.label == ID: 
                out.append(prop[key])
                found = True
        if found == False:
            out.append('NA')
    return out

#### Old code 
def writeStridCells(points, annotation, prefix, directory = None):
    tree = getStructTree()
    stridDict = getStridCells(points, annotation, tree, prefix, directory)
    
    df = pd.DataFrame()
    strids = []
    strgos = []
    strnms = []
    stracs = []
    strccs = []
    for key in stridDict.keys():
        strids.append(stridDict[key]['structure_id'])
        strgos.append(stridDict[key]['structure_graph_order'])
        strnms.append(stridDict[key]['structure_name'])
        stracs.append(stridDict[key]['structure_acronym'])
        strccs.append(stridDict[key]['cellCount'])
    df['structure_id'] = strids
    df['structure_graph_order'] = strgos
    df['structure_name'] = strnms
    df['structure_acronym'] = stracs
    df['cellCount'] = strccs
    
    dfName = prefix+'_regionalCellCounts.txt'
    if directory != None:
        dfName = os.path.join(directory, dfName)
    df.to_csv(dfName, sep='\t')
    print('Regional cell counts were saved to {}'.format(dfName))
    return df

def getStridCells(points, annotation, tree, prefix, directory):
    levels = np.unique(annotation)
    stridDict = defaultdict()
    name_map = tree.get_name_map()
    for strid in levels:
        try:
            IDdict = tree.get_structures_by_name([name_map[strid]]) # by name in order to raise KeyError if NA
            stridDict[strid] = defaultdict()
            stridDict[strid]['structure_name'] = IDdict[0]['name']
            stridDict[strid]['structure_graph_order'] = IDdict[0]['graph_order']
            stridDict[strid]['structure_id'] = strid
            stridDict[strid]['structure_acronym'] = IDdict[0]['acronym']
            stridDict[strid]['cellCount'] = 0
        except KeyError:
            print('No structure corresponds to the ID {}'.format(strid))
    
    for point in points:
        oor = 0
        nir = 0
        unassigned = []
        try:
            strid = annotation[point[0], point[1], point[2]]
        except IndexError:
            oor += 1
            unassigned.append(point)
        try:
            stridDict[strid]['cellCount'] += 1
        except KeyError:
            nir += 1
            unassigned.append(point)
    print('{} cell indexes were out of range\n{} indexes were not assigned to a named structure'.format(oor, nir))
    unassigned = np.concatenate(unassigned)
    unassignedName = prefix+'_unassigned.txt'
    if directory != None:
        unassignedName = os.path.join(directory, unassignedName)
    np.savetxt(unassignedName, unassigned)
    print('All unassigned cells were saved to {}'.format(unassignedName))
    return stridDict


