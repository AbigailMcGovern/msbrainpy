import json
import os
import requests
import numpy as np
import pandas as pd
from skimage import io
from zipfile import ZipFile
from allensdk.api.queries.rma_api import RmaApi
from ._Allen_insitu_images import download_from_section_data_set

# -----------------------------------------------------------------------------
# Structure Unionize Data 
# -----------------------

def rmaStructUniQuery(structureList, 
                      geneIDs, 
                      graph_id=1, 
                      product_id=1,
                      options='[only$eq''genes.entrez_id,data_sets.id'']', 
                      blockSize=2000, 
                      include='structure,section_data_set(genes)', 
                      startRow=0,
                      verbose=True, 
                      writeOut=None, 
                      maxLen=100):
    """
    structure unionize model RMA query to obtain 

    Parameters
    ----------
    structureList: list of int
        list of structure ids for atlas structutre 
    geneIDs: int, list, or None
        entrez ids (int, list, or None)
    graph_id: int
        defalts to adult mouse brain: 1 
    product_id: int
        defaults to mouse: 1 
    options: str
        options strngfor RmaApi().model_query() via rmaStructUniQuery() 
    include:
        include string for RmaApi().model_query() via rmaStructUniQuery() 
    verbose: bool
        all I can think of is a text-based rpg
    writeOut: str or None
        name for file if data is to be saved as json
    maxLen: int
        max number of strucutres to use in any api.model_query(). 
        If too big requests throws error

    Returns
    -------
    all_rows: list
        rows from all structure unionizes identified for the query terms 

    Notes
    -----
    DEPENDENCIES: (1) generateQueryRows(...) / (2) RmaApi().model_query();
        1) retrives sequence of rows from serial RMA query (generator)
        2) allensdk front end for data access (class method)
    """
    model = 'StructureUnionize'
    allRows = []
    tooLong = len(structureList) > maxLen
    if not tooLong:
        criteria = writeStructUniCriteria(structureList, 
                                          geneIDs, 
                                          graph_id, 
                                          product_id)
        rowGenObj = generateQueryRows(model, 
                                      criteria, 
                                      options=options,
                                      include=include, 
                                      startRow=startRow, 
                                      blockSize=blockSize,
                                      verbose=verbose)
        for rows in rowGenObj:
            allRows += rows
            if writeOut is not None:
                with open(writeOut, 'w') as outfile:
                    json.dump(allRows, outfile)
                if verbose:
                    print('data was saved at {}'.format(writeOut))
    if tooLong:
        criteria = []
        num = int((np.ceil((len(structureList) / maxLen))))
        print('The query may be too big.')
        m0 = 'Structure unionizes will be obtained in '
        m1 = f'{num} groups of up to {maxLen} structures'
        print(m0 + m1)
        if verbose:
            m0 = 'Structure group size can be adjusted by editing the '
            m1 = 'maxLen argument of rmaStructUniQuery(...) '
            print(m0 + m1)
            m0 = 'maxLen default = 100, this is appropriate for genes '
            m1 = 'with large numbers of experimental replicates'
            print(m0 + m1)
            m0 = 'This value has no empirical basis, '
            m1 = 'it was the first one that worked'
            print(m0 + m1)
        counter = 0
        for _ in range(num):
            if counter < (maxLen * (num - 1)):
                structs = structureList[counter:counter + maxLen]
            else:
                structs = structureList[counter:]
            crit = writeStructUniCriteria(structs, geneIDs, 
                                          graph_id, product_id)
            criteria.append(crit)
            counter += len(structs)
        num_crit = 0
        for crit in criteria:
            print('Criteria set {}'.format(num_crit))
            num_crit += 1
            rowGenObj = generateQueryRows(model, 
                                          crit, 
                                          options=options,
                                          include=include, 
                                          startRow=startRow, 
                                          blockSize=blockSize, 
                                          verbose=verbose)
            for rows in rowGenObj:
                allRows += rows
                if writeOut:
                    with open(writeOut, 'w') as outfile:
                        json.dump(allRows, outfile)
                    if verbose:
                        print('data was saved at {}'.format(writeOut))
    return allRows


def writeStructUniCriteria(structureList, 
                           geneIDs, 
                           graph_id=1, 
                           product_id=1):
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
    strids = ''.join(str(structureList))[1:-1]
    struct = f'structure[id$in{strids}]'
    rest0 = f' [graph_id$eq{graph_id}],section_data_set[failed$eqfalse]'
    if type(geneIDs) == int:
        genes = f' (genes[entrez_id$eq{geneIDs}],products[id$eq{product_id}])'
    if type(geneIDs) == list:
        genes = f' (genes[entrez_id$in{geneIDs}],products[id$eq{product_id}])'
    if type(geneIDs) is None:
        genes = f'(products[id$eq{product_id}])'
    else:
        m = 'Please enter an argument for geneIDs of type int, list, or None'
        raise ValueError(m)
    string = struct + rest0 + genes
    return string


# -----------------------------------------------------------------------------
# Generic AllenSDK query 
# ----------------------

def generateQueryRows(model, 
                      criteria, 
                      options='[only$eq''genes.entrez_id,data_sets.id'']',
                      include='structure,section_data_set(genes)', 
                      startRow=0, 
                      blockSize=2000, 
                      verbose=True):
    """
    obtain object for iteratively obtaining rows from an AllenSDK Rma query.

    Parameters
    ----------
    model: str
        the type of query
    criteria: str
        RMA criteria string for the query. 
    options: str
        options string for query
    include: str
        include string for query
    startRow: int
        at which row of the query result would you 
        like to start retriving data? 
    blockSize: int
        how many rows would you like to retrive with each iteration? 
        verbose: additional print out, as per usual (bool)

    Yeilds
    ------
    queryResult: list of dict
        rows of data ([{}, {}...]; json.dump() compatable

    Notes
    -----
    DEPENDENCIES: Packages/modules/etc: RmaApi (allensdk.api.queries.rma_api)
        (1) RmaApi().model_query() ;
        1) allensdk front end for data access (class method)
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
            m0 = 'Query returned an error message. '
            m1 = 'RMA query likely contains a syntactic error'
            m = m0 + m1
            raise ValueError(m)
        retrieved = len(queryResult)
        if verbose:
            print('{} rows were found'.format(retrieved))
        startRow += retrieved
        print('new start row: {}'.format(startRow))
        done = retrieved == 0
        yield queryResult


# -----------------------------------------------------------------------------
# Gridded Data 
# ------------

def get_gridded_data(entrez_id, 
                     directory_path=None, 
                     age_id=None, 
                     product_id=1, 
                     verbose=False):
    download_from_section_data_set(get_section_grid, 
                                   entrez_id, 
                                   directory_path=directory_path, 
                                   age_id=age_id, 
                                   product_id=product_id, 
                                   verbose=verbose)


def get_section_grid(section_id, out_dir):
    # this should obtain energy, intensity, and density 
    # but only obtains energy. Not sure why.
    url0 = 'http://api.brain-map.org/grid_data/download/'
    url1 = f'{section_id}?include=energy,density,intensity'
    url = url0 + url1
    response = requests.get(url)
    filename = out_dir + '.zip'
    with open(filename, 'wb') as file:
        file.write(response.content)
    with ZipFile(filename, 'r') as zip_obj:
        zip_obj.extractall(out_dir)


# Some examples ------------------------------------------
# Gap43 entrez id is 14432 -------------------------------
# gap43_section_data_set = section_data_set_query(14432)
# gap43_section0_ims_json = section_image_query(gap43_section_data_set['msg'][0]['id'])
# img15 = download_section_images(gap43_section0_ims_json, 'gap43_79669511', return_img_n=15)
# download_gene_images(14432, age_id=15, product_id=1, verbose=True)
# get_gridded_data(14432, age_id=None, product_id=1, verbose=False) # need to make obtain intensity/density

# Bdnf entrez id is 12064 --------------------------------
# download_gene_images(12064, age_id=None, product_id=1, verbose=True)

