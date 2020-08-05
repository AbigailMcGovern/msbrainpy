import pandas as pd
from ._helpers import dropBool


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
        out = {}
        for i, key in enumerate(keys):
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

