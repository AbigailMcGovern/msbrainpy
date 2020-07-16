import json
import os
import requests
import numpy as np
import pandas as pd
from skimage import io
from zipfile import ZipFile
from allensdk.api.queries.rma_api import RmaApi


# ---------------------------------------- Structure unionize query ----------------------------------------------------

def rmaStructUniQuery(structureList, geneIDs, graph_id=1, product_id=1,
                      options='[only$eq''genes.entrez_id,data_sets.id'']', blockSize=2000,
                      include='structure,section_data_set(genes)', startRow=0,
                      verbose=True, writeOut=None, maxLen=100):
    """
    FUNCTION: structure unionize model RMA query of AGEA
    ARGUMENTS:
        structureList = list of structure ids for atlas structutre (list)
        geneIDs = entrez ids (int, list, or None)
        graph_id = defalts to adult mouse brain: 1 (int)
        product_id = defaults to mouse: 1 (int)
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
            if writeOut is not None:
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



# ----------------------------------------------- Get Gene Images ------------------------------------------------------
def download_gene_images(entrez_id, directory_path=None, age_id=15, product_id=1, verbose=False):
    download_from_section_data_set(get_section_images, entrez_id, directory_path=directory_path, age_id=age_id, 
                                   product_id=product_id, verbose=verbose)


def get_section_images(data_set_id, out_dir, verbose=False):
    section_images = section_image_query(data_set_id, verbose=verbose)
    download_section_images(section_images, out_dir, check_img_exists=True, return_img_n=None)


# function for retrieving lists of images for a given data set
def section_image_query(data_set_id, verbose=True):
    url_0 = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::SectionImage,rma::criteria,"
    url_1 = '[data_set_id$eq{}]'.format(data_set_id)
    url = url_0 + url_1
    response = requests.get(url)
    image_json = response.json()
    if verbose:
        print(f'{image_json["num_rows"]} images matching the criteria were found')
    return image_json


def download_section_images(section_image_json, out_dir_name, check_img_exists=False, return_img_n=None, **kwargs):
    image_ids = [row['id'] for row in section_image_json['msg']]
    if not os.path.exists(out_dir_name):
        os.mkdir(out_dir_name)
    for i in range(len(image_ids)):
        img = image_ids[i]
        filename = f'image_id-{str(img)}.jpeg'
        filepath = os.path.join(out_dir_name, filename)
        if not check_img_exists or not os.path.exists(filepath):
            if i == return_img_n:
                img_n = download_image(img, filepath, return_image=True)
            else:
                download_image(img, filepath)
    if return_img_n is not None:
        return img_n


def download_image(image_id, filename, return_image=False):
    url = f'http://api.brain-map.org/api/v2/image_download/{image_id}'
    response = requests.get(url)
    with open(filename, 'wb') as file:
        file.write(response.content)
    if return_image:
        img = io.imread(filename)
        return img


# ------------------------------------------------ Get Gridded Data ----------------------------------------------------

def get_gridded_data(entrez_id, directory_path=None, age_id=None, product_id=1, verbose=False):
    download_from_section_data_set(get_section_grid, entrez_id, directory_path=directory_path, age_id=age_id, 
                                   product_id=product_id, verbose=verbose)


def get_section_grid(section_id, out_dir):
    # this should obtain energy, intensity, and density but only obtains energy. Not sure why.
    url = f'http://api.brain-map.org/grid_data/download/{section_id}?include=energy,density,intensity'
    response = requests.get(url)
    filename = out_dir + '.zip'
    with open(filename, 'wb') as file:
        file.write(response.content)
    with ZipFile(filename, 'r') as zip_obj:
        zip_obj.extractall(out_dir


# ----------------------------------------------- Find Data Sets -------------------------------------------------------
def download_from_section_data_set(func, entrez_id, directory_path=None, age_id=15, product_id=1, verbose=False):
    section_data_set = section_data_set_query(entrez_id, plane_of_section=None, product_id=product_id,
                                              verbose=verbose)
    dirname = 'entrez_id_' + str(entrez_id) + '_' + section_data_set['msg'][0]['genes'][0]['acronym']
    if directory_path is not None:
        dirname = os.path.join(directory_path, dirname)
    os.makedirs(dirname, exist_ok=True)
    plane_of_section_ids = [row['plane_of_section_id'] for row in section_data_set['msg']]
    for plane_of_section in list(set(plane_of_section_ids)):
        dirname_plane_of_section_id = 'plane_of_section-' + str(plane_of_section)
        dirname_plane_of_section_id = os.path.join(dirname, dirname_plane_of_section_id)
        os.makedirs(dirname_plane_of_section_id, exist_ok=True)
        for row in section_data_set['msg']:
            if row['plane_of_section_id'] == plane_of_section:
                if age_id is None or row['reference_space']['age_id'] == age_id:
                    dirname_section_id = 'age_id-' + str(row['reference_space']['age_id']) + '_id-' + str(row['id'])
                    os.makedirs(dirname_section_id, exist_ok=True)
                    out_dir = os.path.join(dirname_plane_of_section_id, dirname_section_id)
                    func(row['id'], out_dir)


# function for retrieving section data sets for a gene or genes
def section_data_set_query(entrez_ids, plane_of_section=None, product_id=1, write_out=None, verbose=True):
    criteria = section_data_set_criteria(entrez_ids, plane_of_section, product_id=product_id)
    url_0 = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::SectionDataSet,rma::criteria,"
    url_1 = criteria
    url_3 = ',rma::include,genes,reference_space'
    url = url_0 + url_1 + url_3
    response = requests.get(url)
    sections = response.json()
    if verbose:
        print('{} samples matching the criteria were found'.format(sections['num_rows']))
    if write_out is not None:
        with open(write_out, 'w') as outfile:
            json.dump(sections, outfile)
        if verbose:
            print('data was saved at {}'.format(write_out))
    return sections


# function for composing criteria string for experiment search
def section_data_set_criteria(entrez_ids, plane_of_section_id=None, product_id=1):
    """
    FUNCTION: This returns a criteria string with which to query the allen rma api
    ARGS:
        :param entrez_ids: gene id or ids (int, list, None)
        :param plane_of_section_id: 1 -> coronal, 2 -> saggital (None, int)
        :param product_id: Mouse -> 1, default=1 (int)
    RETURNS:
        :return: rma criteria (str)

    Example of successful criteria:
    rma::criteria,[failed$eq'false'],products[id$eq1],plane_of_section[id$eq1],genes[entrez_id$eq14432]
        TIP: Play with the RMA Query builder to find out about query options and variable results
    """
    criteria_0 = '[failed$eq\'false\'],products[id$eq{}],'.format(product_id)
    if type(plane_of_section_id) == int:
        criteria_1 = 'plane_of_section[id$eq{}],'.format(plane_of_section_id)
    elif plane_of_section_id is None:
        criteria_1 = ''
    raise_an_error = True
    if type(entrez_ids) == int:
        criteria_2 = 'genes[entrez_id$eq{}]'.format(entrez_ids)
        raise_an_error = False
    elif type(entrez_ids) == list:
        criteria_2 = 'genes[entrez_id$in{}]'.format(entrez_ids)
        raise_an_error = False
    elif type(entrez_ids) is None:
        criteria_2 = ''
        raise_an_error = False
    elif raise_an_error:
        raise TypeError('Please enter an argument for entrez_ids of type int, list, or None')
    criteria = criteria_0 + criteria_1 + criteria_2
    return criteria

# Some examples ------------------------------------------
# Gap43 entrez id is 14432 -------------------------------
# gap43_section_data_set = section_data_set_query(14432)
# gap43_section0_ims_json = section_image_query(gap43_section_data_set['msg'][0]['id'])
# img15 = download_section_images(gap43_section0_ims_json, 'gap43_79669511', return_img_n=15)
# download_gene_images(14432, age_id=15, product_id=1, verbose=True)
# get_gridded_data(14432, age_id=None, product_id=1, verbose=False) # need to make obtain intensity/density

# Bdnf entrez id is 12064 --------------------------------
# download_gene_images(12064, age_id=None, product_id=1, verbose=True)

