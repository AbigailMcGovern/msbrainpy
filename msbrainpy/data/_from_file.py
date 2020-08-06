import os
import re

def get_image_addresses(path, age_id=None):
    
    """
    FUNCTION:Return a dictionary with all of the absolute file paths for images belonging to a gene directory in
        the format produced by msbrainpy.atlas.query.download_gene_images().
        {'entrez_id-<gene id>_<Acronym>' : {'plane_of_section-<1 or 2>' : {'age_id-<age id>_id-<id>' :
            [<file path 0>, <file path 1>, ...]}
    :param path: Absolute file path of the directory containing the gene directories (str)
    :param age_id: age id (15 = P56) (int or None)
    :return: dictionary of image paths (dict)

     e.g., image_address_dict = get_image_addresses('/Users/amcg0011/Data/InSituData', age_id=None)
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
        string = r'age_id-{}_id-\d*'.format(age_id)
    else:
        string = r'age_id-\d*_id-\d*'
    pattern = re.compile(string)
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
        string = r'age_id-{}*_id-\d*'.format(age_id)
        path_list = get_directory_paths(gene_plane_section_path, string)
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