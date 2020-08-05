import json
import pandas as pd

def parseExpressionToDF(fileName, queryList, saveName=None):
    """
    Extract a specified pd.DataFrame object out of a 
        jsonn file obtained from an AllenSDK: RMA api-based AGEA structure 
        unionize query.

    Parameters
    ----------
    fileName: str
        path to the json file
    queryList: 
        contains (1) the keys for which you wish to obtain data
            (2) if nested, a list of indexes and or keys by which the
            data can be accessed (ascending depth of key/index).

    Returns
    -------
    df: pd.DataFrame
    
    Notes
    -----
    DEPENDENCIES: Packages/modules: json, pandas as pd
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

def removeRows(df, rm_list, rm_header):
    toDrop = []
    for i in range(len(df[rm_header])):
        val = df[rm_header][i]
        if val in rm_list:
            toDrop.append(i)
    df = df.drop(toDrop)
    return df

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