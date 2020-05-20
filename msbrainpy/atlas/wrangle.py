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

oapi = ontologies_api()
structure_graph = oapi.get_structures_with_sets([1])
structure_graph = structure_tree.clean_structures(structure_graph)
tree = structure_tree(structure_graph)
del oapi

# ------------------------------------------- Obtain CSV from json  ----------------------------------------------------

def parse_expression_to_df(file_name, query_list, save_name=None):
    """
    FUNCTION: pull a specified pd.data_frame object out of a jsonn file obtained from an allen_sdk:
        RMA api-based AGEA structure unionize query.
    ARGUMENTS:
        file_name = json file path (str).
            Forces you to save the direct database output to json.
            Please keep all of this (regardless of how much you choose to obtain)!
        query_list = contains (1) the keys from which you wish to obtain data for each row.
            (2) if the data you wish to acess is nested provide a list of indexes and or keys by which the
            data can be accessed (ascending depth of key/index).
    DEPENDENCIES: Packages/modules: json, pandas as pd
    RETURNS: pd.data_frame
    """
    with open(file_name, 'r') as json_file:
        expression = json.load(json_file)
    df = pd.data_frame()
    rows = []
    for i in range(len(expression)):
        row = []
        for query in query_list:
            if type(query) == str:
                row.append(expression[i][query])
            if type(query) == list:
                val = expression[i][query[0]]
                for ii in range(len(query) - 1):
                    val = val[query[ii + 1]]
                row.append(val)
        rows.append(row)
    for j in range(len(query_list)):
        query = query_list[j]
        if type(query) == list:
            started = False
            for k in range(len(query)):
                if not started:
                    col_name = str(query[0])
                    started = True
                elif started:
                    col_name = col_name + '_' + str(query[k])
        if type(query) == str:
            col_name = query
        values = []
        for i in range(len(rows)):
            row = rows[i]
            value = row[j]
            values.append(value)
        df.insert(j, col_name, values)
    if save_name != None:
        df.to_csv(save_name)
    return df

# ----------------------------- Wrangling DFs on the basis of atlas structures -----------------------------------------

def write_struct_d_fs(df, struct_name_list, filename, IDheader='id', exclude=None):
    oapi = ontologies_api()
    # The 1 refers to the adult mouse brain atlas
    structure_graph = oapi.get_structures_with_sets([1])
    # clean_structures() removes unused fields
    structure_graph = structure_tree.clean_structures(structure_graph)
    # a class with methods for aceessing and using ontologies data
    tree = structure_tree(structure_graph)
    for struct_name in struct_name_list:
        towrite = get_struct_df(df, struct_name, tree, IDheader=IDheader, exclude=exclude)
        writename = struct_name + filename
        towrite.to_csv(writename, sep='\t')

def get_struct_df(df, struct_name, tree, IDheader='id', exclude=None):  # use with np.s_[:]
    """
    FUNCTION: obtain decendents of a particular structure from a data_frame
    ARGUMENTS:
        df = data (pandas.data_frame)
        struct_name = name of the parent structure for which data should be obtained (str)
        tree = instance of structure_tree (allensdk.core.structure_tree.structure_tree())
        IDheader = header of the id column (str)
        exclude = ids to exclude. id col is treated as list and spliced using np.s_[...] (np.s_[:])
    DEPENDENCIES: (1) get_sub_strids(..); (2) drop_bool(...)
        1) obtain list of structures that decend from a given parent structure from a list of structures (function)
        2) removes specified bool from data_frame when supplied list of bool of len(df) (function)
    RETURNS: data_frame
    """
    struct_dict = tree.get_structures_by_name([struct_name])
    strid = struct_dict[0]['id']

    strids = df[IDheader]
    if exclude is not None:
        strids = strids[exclude]
    strid_list = get_sub_strids(strids, strid, verbose=False)
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

def add_from_structure_tree(df, use_col, to_add, to_use, tree):
    """
     write a docstring. I mean it.
    """
    if to_add == to_use:
        message = 'The to_add ({}) argument cannot be the same as the to_use argument({})'.format(to_add, to_use)
        print('''Please choose the type of column to add (acronym, id, name, ...) \n 
            and the name (use_col) and type of column containing (acronym, id, or name) \n 
            the information with which to query the tree''')
        raise value_error(message)
    header = 'structure_' + to_add
    if to_use == 'id':
        out = add_col_id(df, use_col, to_add, to_use, tree, header, loc=0)
    if to_use == 'name':
        out = add_col_name(df, use_col, to_add, to_use, tree, header, loc=0)
    if to_use == 'acronym':
        out = add_col_acronym(df, use_col, to_add, to_use, tree, header, loc=0)
    return out

def add_col_id(df, id_col, to_add, to_use, tree, header, loc=0):
    additions = []
    for ID in df[id_col]:
        try:
            struct_dict = tree.get_structures_by_id([ID])[0]
            addition = struct_dict[to_add]
        except key_error:
            addition = None
        additions.append(addition)
    df.insert(loc, header, additions)
    return df

def add_col_acronym(df, acronym_col, to_add, to_use, tree, header, loc=0):
    additions = []
    for acronym in df[acronym_col]:
        try:
            struct_dict = tree.get_structures_by_acronym([acronym])[0]
            addition = struct_dict[to_add]
        except key_error:
            addition = None
        additions.append(addition)
    df.insert(loc, header, additions)
    return df

def add_col_name(df, name_col, to_add, to_use, tree, header, loc=0):
    additions = []
    for name in df[name_col]:
        try:
            struct_name = get_search_word(name, structs_all, sep=r', ')
        except:
            addition = None
        try:
            struct_dict = tree.get_structures_by_name([struct_name])[0]
            addition = struct_dict[to_add]
        except key_error:
            addition = None

        additions.append(addition)

    df.insert(loc, header, additions)
    return df

def drop_parents(df, tree, IDheader='structure_id'):
    """
     FUNCTION: remove any partent structures from a dataframe using the structure_id column.
         is dependent on the drop_bool() function.
     ARGS:
         df : data to clean (pandas.data_frame)
         tree : instance of structure_tree (allensdk.core.structure_tree.structure_tree())
         IDheader : header of the id column (str)
     DEPENDENCIES: (1) drop_bool(...)
         1) removes specified bool from data_frame when supplied list of bool of len(df) (function)
     RETURNS: copy df without the parent structures (pandas.data_frame)
     """
    parents = []
    unique_elem = np.unique([strid for strid in df[IDheader]])
    print('there are {} unique IDs in the data frame'.format(len(unique_elem)))

    for strid_a in unique_elem:
        values = []
        matches = []
        for strid_b in unique_elem:
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
    parent_bool = []
    for strid in df[IDheader]:
        b = strid in parents
        parent_bool.append(b)
    out = drop_bool(df, parent_bool, todrop=True)


def save_struct_files(df, prefix, struct_name='Isocortex', regex='[A-Z]*[a-z]*\d.*',
                    re_names=['sanslayers', 'layers'], IDheader='structure_id',
                    tree=tree, acronym_header='structure_acronym', exclude=None):
    name = prefix + '_' + struct_name
    name0 = name + '.csv'
    struct = get_struct_df(df, struct_name, tree, IDheader=IDheader)
    struct.to_csv(name0)
    if regex != None:
        name1 = name + '_' + re_names[0]
        name2 = name + '_' + re_names[1]
        struct_sre = drop_re_rows(struct, acronym_header, regex='[A-Z]*[a-z]*\d.*', remove=True)
        struct_sre = drop_parents(struct_sre, tree, IDheader=IDheader)
        struct_sre.to_csv(name1)
        struct_wre = drop_re_rows(struct, acronym_header, regex='[A-Z]*[a-z]*\d.*', remove=False)
        struct_wre = drop_parents(struct_wre, tree, IDheader=IDheader)
        struct_wre.to_csv(name2)
        return struct_sre, struct_wre
    else:
        return struct

# ------------------------------------ Adding acronym-associated info to DF --------------------------------------------

def add_i_ctx_groups(df, ac_header='structure_acronym',
                  add_parent_col_names=('ctx_subregion', 'ctx_region'), write_out=None):
    groups = []
    for ac in df[ac_header]:
        ac = Acronym(ac)
        groups.append(ac.group)
    df['group'] = groups
    if add_parent_col_names is not None:
        add_re_column(df, ac_header, '[A-Z]*[a-z]*', add_parent_col_names[0], loc=0)
        add_re_column(df, ac_header, '[A-Z]*', add_parent_col_names[1], loc=0)
    if write_out != None:
        df.to_csv(write_out, sep='\t')
    return df

class Acronym:
    def __init__(self, string):
        """
        For shits and giggles or planned expansion and use, either way
        ARNUMENTS:
            string: acronym of choice
        """
        self.string = string
        self.group = self.get_group()
        # self.other_info?

    def get_group(self):
        ac = self.string
        somatomotor = ['MOp', 'SS']
        medial = ['PTLp', 'VISam', 'VISpm', 'RSP']
        temporal = ['AUD', 'PERI', 'TE', 'ECT']
        visual = ['VISal', 'VISl', 'VISp', 'VISpl']
        anterolateral = ['VISC', 'GU', 'AI']
        prefrontal = ['MOs', 'FRP', 'ACA', 'ORB', 'PL', 'ILA']
        group_defs = {'somatomotor': somatomotor, 'medial': medial, 'temporal': temporal,
                     'anterolateral': anterolateral, 'prefrontal': prefrontal, 'visual': visual}
        group = 'undefined'
        for key in group_defs.keys():
            in_group = False
            for value in group_defs[key]:
                if ac.find(value) != -1:
                    in_group = True
            if in_group:
                group = key
        return group

# ------------------------------------------- Regex-based DF wrangling  ------------------------------------------------

def add_re_column(df, str_col_header, regex, re_col_header, loc=0):
    out = df
    pattern = re.compile(regex)
    re_column = []
    for string in out[str_col_header]:
        match_obj = pattern.match(string)
        re_column.append(match_obj[0])
    out.insert(loc, re_col_header, re_column)
    return out

def add_re_column_exc(df, str_col_header, regex, re_col_header, loc=0, exc=['4', '5', '6']):
    out = df
    pattern = re.compile(regex)
    re_column = []
    for string in out[str_col_header]:
        match_obj = pattern.match(string)
        exclude = False
        for cond in exc:
            if cond in string:
                exclude = True
        if not exclude:
            re_column.append(match_obj[0])
        if exclude:
            re_column.append('NA')
    out.insert(loc, re_col_header, re_column)
    return out

def drop_re_rows(df, str_col_header, regex='[A-Z]*[a-z]*\d.*', remove=True):
    """
    When applied to acronyms, this function will remove all acronyms that
    belong to a cortical layer. Invert to keep only layers (technically with # in acronym)
    """
    matches = []
    pattern = re.compile(regex)
    for string in df[str_col_header]:
        match = pattern.findall(string)
        if len(match) > 0:
            matches.append(True)
        if len(match) == 0:
            matches.append(False)

    out = drop_bool(df, matches, todrop=remove)
    return out

def add_layer(df, acronym_header):
    out = df.copy()
    regex = '[A-Z]*.*[a-z]*\d.*'
    pattern = re.compile(regex)
    d = '\d.*'
    digit = re.compile(d)

    layers = []
    for string in df[acronym_header]:
        try:
            match = pattern.match(string)[0]
            layer = digit.findall(match)[0]
            layers.append(layer)
        except type_error:
            print('A layer could not be identified in {}'.format(string))
            layers.append('NA')
    out['cortical_layer'] = layers
    return out

# --------------------------------------------- Statistical wrangling --------------------------------------------------

def means_se_ms_sort(df, savename=None, groupby='structure_acronym'):
    """
    FUNCTION: Obtain a data frame containing means and standard errors for groups (e.g., anatomical structures)
    ARGUMENTS:
        df: pd.data_frame or dictionary with save names as keys and dfs as values
        savename: if a df is supplied, name with which to save output
        groupby: column header of categorical variable for which to compute means
    DEPENDENCIES: Packages/modules/etc: pandas as pd
    Native: (1) amalgamate_d_fs(df_dict)
    1) joins two data frames with corresponding rows
    RETURNS: pd.data_frame
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
            result = amalgamate_d_fs({'mean': means, 'sem': sems})
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
        out = amalgamate_d_fs({'mean': means, 'sem': sems})
        out.sort_values(groupby)
        out.to_csv(savename)
        print('You can supply a dictionary of save_name: df pairs to expedite bulk processing')
    return out

# ----------------------------------------- Row matching & column addition ---------------------------------------------

def match_and_amalgamate(df0, df1, prefixes, groups=('structure_acronym', 'structure_acronym')):
    df0, df1 = match_rows(df0, df1, groups=groups)
    df = amalgamate_d_fs({prefixes[0]: df0, prefixes[1]: df1})
    return df

def match_rows(df0_, df1_, groups=('structure_acronym', 'structure_acronym')):
    """
    FUNCTION: Checks if two data frames contain identical rows on the basis of categories/IDs (i.e., groups).
        If not, unnecessary rows are recursively removed from longer data frame.
        Checks that rows are in correct sequence.
            NB: for DFs of means use 'structure_acronym' else create col of unique concatenated strings
            (strid_struct_set_id) - but this is a less likely scenario
    ARGUMENTS:
        df0_ = pd.data_frame
        df1_ = pd.data_frame
        groups = columns containing the categories or IDs that must match
    DEPENDENCIES: Packages/modules/etc: numpy as np, pandas as pd
    RETURNS: pd.data_frame, pd.data_frame
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
                except assertion_error:
                    errors += 1
                    print('{} does not match {}'.format(df0[groups[0]][i], df1[groups[1]][i]))
            message = 'there were {} mismatches'.format(errors)
            if errors > 0:
                raise assertion_error(message)
            else:
                return df_dict[True]['df'], df_dict[False]['df']
    else:
        zero_longer = len(set(df0[groups[0]])) >= len(set(df1[groups[1]]))
        zero_shorter = not zero_longer
        bools = [ac in set(df_dict[zero_shorter]['df'][groups[df_dict[zero_shorter]['index']]])
                 for ac in df_dict[zero_longer]['df'][groups[df_dict[zero_longer]['index']]]]
        df_dict[zero_longer]['df'] = drop_bool(df_dict[zero_longer]['df'], bools, todrop=False)
        number = 0
        rows = []
        for i in range(len(bools)):
            b = bools[i]
            if not b:
                number += 1
                rows.append(i)
        print('{} rows were removed from df {}'.format(number, df_dict[zero_longer]['index']))
        print('The rows removed were: {}'.format(rows))
        df0, df1 = df_dict[True]['df'], df_dict[False]['df']
        return match_rows(df0, df1, groups)

def amalgamate_d_fs(df_dict):
    for key in df_dict.keys():
        col_names = list(df_dict[key].columns)
        new_names = []
        for col in col_names:
            new_names.append(col + '_' + key)
        dictionary = dict(zip(col_names, new_names))
        df_dict[key] = df_dict[key].rename(columns=dictionary)
    result = pd.concat([df_dict[key] for key in df_dict.keys()], axis=1, sort=False)
    return result

# ----------------------------------------------- Basic functions ------------------------------------------------------

def write_all_structures(tree):
    name_map = tree.get_name_map()
    structures = [name_map[i] for i in tree.descendant_ids([997])[0]]
    return structures

def write_all_cortical_structures(tree):
    """B fulcher allensdk repo"""
    name_map = tree.get_name_map()
    cortex_structures = [name_map[i] for i in tree.descendant_ids([315])[0]]
    return cortex_structures

def write_sub_structures(struct_name, tree, filename):
    struct_dict = tree.get_structures_by_name([struct_name])
    strid = struct_dict[0]['id']
    name_map = tree.get_name_map()
    structures = [name_map[i] for i in tree.descendant_ids([strid])[0]]
    ids = [i for i in tree.descendant_ids([strid])[0]]
    df = pd.data_frame()
    df['structure_name'] = structures
    df['structure_id'] = ids
    df.to_csv(filename, sep='\t')
    return df

def multi_sub_stids(query_structs, tree):
    query_stids = []
    for struct in query_structs:
        strids = sub_structures_id(struct, tree)
        for strid in strids:
            query_stids.append(strid)
    return query_stids

def sub_structures_id(struct_name, tree):
    struct_dict = tree.get_structures_by_name([struct_name])
    strid = struct_dict[0]['id']
    ids = [i for i in tree.descendant_ids([strid])[0]]
    return ids

def get_ms_tree():
    oapi = ontologies_api()
    # The 1 refers to the adult mouse brain atlas
    structure_graph = oapi.get_structures_with_sets([1])
    # clean_structures() removes unused fields
    structure_graph = structure_tree.clean_structures(structure_graph)
    # a class with methods for aceessing and using ontologies data
    tree = structure_tree(structure_graph)
    return tree

def get_sub_strids(strids, strid, verbose=False):
    out = []
    for strida in strids:
        if verbose == True:
            is_desc = '' if tree.structure_descends_from(strida, strid) else ' not'
            print('{0} is{1} in {2}'.format(name_map[strida], is_desc, name_map[strid]))
        if tree.structure_descends_from(strida, strid) == True:
            out.append(strida)
    return out

def remove_rows(df, rm_list, rm_header):
    to_drop = []
    for i in range(len(df[rm_header])):
        val = df[rm_header][i]
        if val in rm_list:
            to_drop.append(i)
    df = df.drop(to_drop)
    return df

def get_search_word(string, to_search, sep=r', '):
    """
    The Ero et al., 2018, Front. Neuroinfo. data has no commas in structure names and this
    returns an acceptable search-term when fed a list of Allen Atlas structure names
    """
    arr = []
    for op in to_search:
        search_list = re.split(sep, op)
        success_count = 0
        for item in search_list:
            if item in string:
                success_count += 1
        arr.append(success_count)
    arr = np.array(arr)
    i = np.argmax(arr)
    return to_search[i]

def row_matches(df, to_match, match_col, todrop=False):
    match_bool = []
    for elem in df[match_col]:
        b = elem in to_match
        match_bool.append(b)
    out = drop_bool(df, match_bool, todrop=todrop)
    return out

def drop_bool(df, bools, todrop=False):
    out = df.copy()
    out['bools'] = bools
    out = out[out['bools'] != todrop]
    out = out.drop(['bools'], axis=1)
    return out


# ------------------------------------------------- Deprecated? --------------------------------------------------------

def add_acronym_col_name(df, names_col, tree, header='acronym', loc=0):
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
    structs_all = write_all_structures(tree)
    for name in df[names_col]:
        try:
            struct_name = get_search_word(name, structs_all, sep=r', ')
        except:
            acronym = None
        try:
            struct_dict = tree.get_structures_by_name([struct_name])[0]
        except key_error:
            acronym = None
        acronym = struct_dict['acronym']
        acronyms.append(acronym)

    df.insert(loc, header, acronyms)
    return df

def add_acronym_col_id(df, id_col, tree, header='structure_acronym', loc=0):
    acronyms = []
    for ID in df[id_col]:
        try:
            struct_dict = tree.get_structures_by_id([ID])[0]
        except key_error:
            acronym = None
        acronym = struct_dict['acronym']
        acronyms.append(acronym)

    df.insert(loc, header, acronyms)
    return df

def add_id_col_acronym(df, acronym_col, tree, header='structure_id', loc=0):
    IDs = []
    for acronym in df[acronym_col]:
        try:
            struct_dict = tree.get_structures_by_acronym([acronym])[0]
        except key_error:
            ID = None
        ID = struct_dict['id']
        IDs.append(ID)

    df.insert(loc, header, IDs)
    return df

def add_id_col_name(df, name_col, tree, header='structure_id', loc=0):
    IDs = []
    for name in df[name_col]:
        try:
            struct_name = get_search_word(name, structs_all, sep=r', ')
        except:
            ID = None
        try:
            struct_dict = tree.get_structures_by_name([struct_name])[0]
            ID = struct_dict['id']
        except key_error:
            ID = None

        IDs.append(ID)

    df.insert(loc, header, IDs)
    return df

def add_name_col_acronym(df, acronym_col, tree, header='structure_name', loc=0):
    names = []
    for acronym in df[acronym_col]:
        try:
            struct_dict = tree.get_structures_by_acronym([acronym])[0]
            name = struct_dict['id']
        except key_error:
            name = None

        names.append(name)

    df.insert(loc, header, names)
    return df

def add_name_col_id(df, id_col, tree, header='structure_name', loc=0):
    names = []
    for ID in df[id_col]:
        try:
            struct_dict = tree.get_structures_by_id([ID])[0]
            name = struct_dict['acronym']
        except key_error:
            name = None

        names.append(name)

    df.insert(loc, header, names)
    return df

def get_scatter_points(df_x, df_y, re_header_x, re_header_y, mean_header_x, mean_header_y, x_header, y_header,
                     save_name=None, out_header='parent_acronym', regex='[A-Z]*[a-z]*'):
    '''
    FUNCTION: y is first obtained (i.e., measurment points) and equivalents from x are added to the dataframe
    '''
    y = get_average_re(df_y, re_header_y, mean_header_y, regex=regex, out_header=out_header)
    x = get_average_re(df_x, re_header_x, mean_header_x, regex=regex, out_header=out_header)

    xlist = []
    ylist = []
    parents = []
    for j in range(len(y[out_header])):
        parent = y[out_header][j]
        for i in range(len(x[out_header])):
            if x[out_header][i] == parent:
                xlist.append(x[mean_header_x][i])
                ylist.append(y[mean_header_y][j])
                parents.append(parent)

    out = pd.data_frame()
    out[out_header] = parents
    out[x_header] = xlist
    out[y_header] = ylist

    if save_name is not None:
        out.to_csv(save_name, sep='\t')
    return out


def get_average_re(df, re_header, mean_header, regex='[A-Z]*[a-z]*', out_header='parent_acronym'):
    parent = []
    pattern = re.compile(regex)
    for string in df[re_header]:
        match_obj = pattern.match(string)
        parent.append(match_obj[0])
    df[out_header] = parent
    means = df.groupby(out_header)[mean_header].mean()
    out = pd.data_frame()
    parents0 = []
    for index in means.index:
        parents0.append(index)

    out[out_header] = parents0
    out[mean_header] = means.values
    return out

def write_struct_d_fsand_means(df, structs, prefix, groupby='section_data_set_id'):
    started = False
    for struct in structs:
        struct_df = get_struct_df(df, struct, tree, IDheader='structure_id')
        save_name = prefix + '_' + struct + '.txt'
        struct_df.to_csv(save_name, sep='\t')
        mean_struct = struct_df.groupby([groupby]).mean()
        save_name = prefix + '_' + struct + '_mean.txt'
        mean_struct.to_csv(save_name, sep='\t')
        mean_struct = pd.read_csv(save_name, sep='\t')
        names = []
        for i in range(len(mean_struct[groupby])):
            names.append(struct)
        mean_struct['region'] = names
        if not started:
            out = mean_struct
            started = True
        if started:
            out = out.append(mean_struct)
    save_name = prefix + '_structure_means.txt'
    out.to_csv(save_name, sep='\t')
    return out


# ----------------------------------- Get Absolute File and Directory Paths -------------------------------------------
def get_image_addresses(path, age_id=None):
    
    """
    FUNCTION:Return a dictionary with all of the absolute file paths for images belonging to a gene directory in
        the format produced by msbrainpy.atlas.query.download_gene_images().
        Approx. {'entrez_id-<gene id>_<Acronym>' : {'plane_of_section-<1 or 2>' : {'age_id-<age id>_id-<id>' :
            [<file path 0>, <file path 1>, ...]}
    :param path: Absolute file path of the directory containing the gene directories (str)
    :param age_id: age id (15 = P56) (int or None)
    :return: dictionary of image paths (dict)

     e.g., image_address_dict = get_image_addresses('/Users/amcg0011/Data/in_situ_data', age_id=None)
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
