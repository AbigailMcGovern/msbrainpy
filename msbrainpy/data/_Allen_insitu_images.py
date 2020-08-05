import json
import os
import requests
from skimage import io

# Allen Institute in situ data (jpeg)
def download_gene_images(entrez_id, 
                         directory_path=None, 
                         age_id=15, 
                         product_id=1, 
                         verbose=False):
    download_from_section_data_set(get_section_images, 
                                   entrez_id, 
                                   directory_path=directory_path, 
                                   age_id=age_id, 
                                   product_id=product_id, 
                                   verbose=verbose)


def get_section_images(data_set_id, 
                       out_dir, 
                       verbose=False):
    section_images = section_image_query(data_set_id, 
                                         verbose=verbose)
    download_section_images(section_images, 
                            out_dir, 
                            check_img_exists=True, 
                            return_img_n=None)


# function for retrieving lists of images for a given data set
def section_image_query(data_set_id, verbose=True):
    url_0 = "http://api.brain-map.org/api/v2/data/query.json?"
    url_1 = "criteria=model::SectionImage,rma::criteria,"
    url_2 = f'[data_set_id$eq{data_set_id}]'
    url = url_0 + url_1 + url_2
    response = requests.get(url)
    image_json = response.json()
    if verbose:
        m = f'{image_json["num_rows"]} images matching the criteria were found'
        print(m)
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


def download_from_section_data_set(func, 
                                   entrez_id, 
                                   directory_path=None, 
                                   age_id=15, 
                                   product_id=1, 
                                   verbose=False):
    section_data_set = section_data_set_query(entrez_id, 
                                              plane_of_section=None, 
                                              product_id=product_id,
                                              verbose=verbose)
    gene_name = section_data_set['msg'][0]['genes'][0]['acronym']
    d = f'entrez_id_{entrez_id}_'
    dirname = d + gene_name
    if directory_path is not None:
        dirname = os.path.join(directory_path, dirname)
    os.makedirs(dirname, exist_ok=True)
    plane_of_section_ids = [row['plane_of_section_id'] 
                            for row in section_data_set['msg']]
    for plane_of_section in list(set(plane_of_section_ids)):
        d = f'plane_of_section-{plane_of_section}'
        dirname_plane_of_section_id = os.path.join(dirname, d)
        os.makedirs(dirname_plane_of_section_id, exist_ok=True)
        for row in section_data_set['msg']:
            if row['plane_of_section_id'] == plane_of_section:
                age = row['reference_space']['age_id']
                if age_id is None or age == age_id:
                    d = f'age_id-{age}'
                    id_ = row['id']
                    dirname_section_id = d + f'_id-{id_}'
                    os.makedirs(dirname_section_id, exist_ok=True)
                    out_dir = os.path.join(dirname_plane_of_section_id, dirname_section_id)
                    func(row['id'], out_dir)


# function for retrieving section data sets for a gene or genes
def section_data_set_query(entrez_ids, 
                           plane_of_section=None, 
                           product_id=1, 
                           write_out=None, 
                           verbose=True):
    criteria = section_data_set_criteria(entrez_ids, 
                                         plane_of_section, 
                                         product_id=product_id)
    url_ = "http://api.brain-map.org/api/v2/data/query.json?"
    url_0 = "criteria=model::SectionDataSet,rma::criteria,"
    url_1 = criteria
    url_3 = ',rma::include,genes,reference_space'
    url = url_ + url_0 + url_1 + url_3
    response = requests.get(url)
    sections = response.json()
    if verbose:
        num_rows = sections['num_rows']
        m = f'{num_rows} samples matching the criteria were found'
        print(m)
    if write_out is not None:
        with open(write_out, 'w') as outfile:
            json.dump(sections, outfile)
        if verbose:
            print('data was saved at {}'.format(write_out))
    return sections


# function for composing criteria string for experiment search
def section_data_set_criteria(entrez_ids, 
                              plane_of_section_id=None, 
                              product_id=1):
    """
    get a criteria string with which to query the allen rma api

    Parameters
    ----------
    entrez_ids: int, list, or None
        gene id or ids 
    plane_of_section_id: int, None
        1 -> coronal, 2 -> saggital 
    product_id: int
        Mouse -> 1, default=1
    
    Returns
    -------
        :return: rma criteria (str)

    Notes
    -----
    Example of successful criteria:
    rma::criteria,[failed$eq'false'],products[id$eq1],plane_of_section[id$eq1]
        ,genes[entrez_id$eq14432]
    TIP: Play with the RMA Query builder to find out about query options 
    and variable results
    """
    criteria_0 = f'[failed$eq\'false\'],products[id$eq{product_id}],'
    if type(plane_of_section_id) == int:
        criteria_1 = 'plane_of_section[id$eq{plane_of_section_id}],'
    elif plane_of_section_id is None:
        criteria_1 = ''
    raise_an_error = True
    if type(entrez_ids) == int:
        criteria_2 = 'genes[entrez_id$eq{entrez_ids}]'
        raise_an_error = False
    elif type(entrez_ids) == list:
        criteria_2 = 'genes[entrez_id$in{entrez_ids}]'
        raise_an_error = False
    elif type(entrez_ids) is None:
        criteria_2 = ''
        raise_an_error = False
    elif raise_an_error:
        m = 'Please enter an argument for entrez_ids of int, list, or None'
        raise TypeError(m)
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
