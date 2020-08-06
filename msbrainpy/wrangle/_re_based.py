import re
import pandas as pd
from ._helpers import dropBool

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

def drop_reRows(df, strColHeader, regex=r'[A-Z]*[a-z]*\d.*', remove=True):
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

