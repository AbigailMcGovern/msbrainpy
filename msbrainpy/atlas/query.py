import json
import numpy as np
import pandas as pd
from allensdk.api.queries.rma_api import RmaApi


# ---------------------------------------- Structure unionize query ----------------------------------------------------

def rmaStructUniQuery(structureList, geneIDs, graph_id=1, product_id=1,
                      options='[only$eq''genes.entrez_id,data_sets.id'']', blockSize=2000,
                      include='structure,section_data_set(genes)', startRow=0,
                      verbose=True, writeOut=False, maxLen=100):
    """
    FUNCTION: structure unionize model RMA query of AGEA
    ARGUMENTS:
        structureList = list of structure ids for atlas structutre (list)
        geneIDs = entrez ids (int, list, or None)
        graph_id = defalts to adult mouse: 1 (int)
        product_id = defaults to 1 (int)
        options = options for RmaApi().model_query() via rmaStructUniQuery() (string)
        include = include for RmaApi().model_query() via rmaStructUniQuery() (string)
        verbose = all I can think of is a text-based rpg
        writeOut = name for file if data is to be saved as json (str or None)
        maxLen = max number of strucutres to use in any api.model_query(). If too big requests throws error
    DEPENDENCIES: (1) generateQueryRows(...) / (2) RmaApi().model_query();
        1) retrives sequence of rows from serial RMA query (generator)
        2) allensdk front end for data access (class method)
    RETURNS = rows from all structure unionizes identified for the query terms (list)
    """
    model = 'StructureUnionize'
    allRows = []
    tooLong = len(structureList) > maxLen
    if not tooLong:
        criteria = writeStructUniCriteria(structureList, geneIDs, graph_id, product_id)
        rowGenObj = generateQueryRows(model, criteria, options=options,
                                      include=include, startRow=startRow, blockSize=blockSize, verbose=verbose)
        for rows in rowGenObj:
            allRows += rows
            if writeOut:
                with open(writeOut, 'w') as outfile:
                    json.dump(allRows, outfile)
                if verbose:
                    print('data was saved at {}'.format(writeOut))
    if tooLong:
        criteria = []
        num = int((np.ceil((len(structureList) / maxLen))))
        print('The query may be too big.')
        print('Structure unionizes will be obtained in {} groups of up to {} structures'.format(num, maxLen))
        if verbose:
            print('Structure group size can be adjusted by editing the maxLen argument of rmaStructUniQuery(...) ')
            print('maxLen default = 100, this is appropriate for genes with large numbers of experimental replicates')
            print('This value has no empirical basis, it was the first one that worked')
        counter = 0
        for i in range(num):
            if counter < (maxLen * (num - 1)):
                structs = structureList[counter:counter + maxLen]
            else:
                structs = structureList[counter:]
            crit = writeStructUniCriteria(structs, geneIDs, graph_id, product_id)
            criteria.append(crit)
            counter += len(structs)
        num_crit = 0
        for crit in criteria:
            print('Criteria set {}'.format(num_crit))
            num_crit += 1
            rowGenObj = generateQueryRows(model, crit, options=options,
                                          include=include, startRow=startRow, blockSize=blockSize, verbose=verbose)
            for rows in rowGenObj:
                allRows += rows
                if writeOut:
                    with open(writeOut, 'w') as outfile:
                        json.dump(allRows, outfile)
                    if verbose:
                        print('data was saved at {}'.format(writeOut))
    return allRows


def writeStructUniCriteria(structureList, geneIDs, graph_id=1, product_id=1):
    """
    FUNCTION: write a structure unionise query criteria string
    ARGUMENTS:
        structureList: structure IDs of desired structures (list: [int, ...])
        geneIDs: entrez ID of the gene/s in question (int, list, or None)
        graph_id: 1 = adult mouse (P56)
        product_id: 1 = also mouse?
    DEPENDENCES: NA
    RETURNS: str
    """
    struct = 'structure[id$in{}]'.format(''.join(str(structureList))[1:-1])
    rest0 = ' [graph_id$eq{}],section_data_set[failed$eqfalse]'.format(graph_id)
    if type(geneIDs) == int:
        genes = ' (genes[entrez_id$eq{}],products[id$eq{}])'.format(geneIDs, product_id)
    if type(geneIDs) == list:
        genes = ' (genes[entrez_id$in{}],products[id$eq{}])'.format(geneIDs, product_id)
    if type(geneIDs) is None:
        genes = '(products[id$eq{}])'.format(product_id)
    else:
        raise ValueError('Please enter an argument for geneIDs of type int, list, or None')
    string = struct + rest0 + genes
    return string


# --------------------------------------------- General RMA query ------------------------------------------------------

def generateQueryRows(model, criteria, options='[only$eq''genes.entrez_id,data_sets.id'']',
                      include='structure,section_data_set(genes)', startRow=0, blockSize=2000, verbose=True):
    """
    GENERATOR: obtain object for iteratively obtaining rows from an AllenSDK Rma query.
    ARGUMENTS:
        model: the type of query (str)
        criteria: criteria for the query (str). Where a function to automatically obtain this string does not exist,
            please write one.
        options: options string for query
        include: include string for query
        startRow: at which row of the query result would you like to start retriving data? (int)
        blockSize: how many rows would you like to retrive with each iteration? (int)
        verbose: additional print out, as per usual (bool)
    DEPENDENCIES: Packages/modules/etc: RmaApi (allensdk.api.queries.rma_api)
        (1) RmaApi().model_query() ;
        1) allensdk front end for data access (class method)
    YEILDS: rows of data ([{}, {}...]; json.dump() compatable)
    """
    # [only$eq''genes.entrez_id,data_sets.id'']
    if verbose:
        print('Obtaining generator object for {} query'.format(model))
        print('Criteria: {}'.format(criteria))
        print('Options: {}'.format(options))
        print('Include: {}'.format(include))
    done = False
    api = RmaApi()
    while not done:
        if verbose:
            print('Querying rows {} - {}'.format(startRow, startRow + blockSize))
        queryResult = api.model_query(model=model, criteria=criteria, options=options,
                                      include=include, start_row=startRow, num_rows=blockSize)
        if type(queryResult) == str:
            raise ValueError('Query returned an error message. RMA query likely contains a syntactic error')
        retrieved = len(queryResult)
        if verbose:
            print('{} rows were found'.format(retrieved))
        startRow += retrieved
        print('new start row: {}'.format(startRow))
        done = retrieved == 0
        yield queryResult
