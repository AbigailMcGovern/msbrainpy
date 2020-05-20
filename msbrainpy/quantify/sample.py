import os
import pandas as pd
from msbrainpy.chain import Chain, make_chain_template_dict
from msbrainpy.io import generate_from_diretory
from msbrainpy.quantify.processing import nuclear_detection_list, img_file_count, extract_img_file_counts, \
    extract_chunk_corrected


# ------------------------------------------ Purpose-specific functions ------------------------------------------------

def sample_nuclear_detection(directory_path, prefix, name, selem=6, min_sigma=0.25,
                           max_sigma=5, threshold=0.04):
    func_list = nuclear_detection_list(selem=selem, min_sigma=min_sigma, max_sigma=max_sigma,
                                    threshold=threshold, directory_path=directory_path, name=name, write_out=True)
    df = process_directory(directory_path, func_list, prefix)
    return df


# to be continued ...
# specifically, may want to extend this process to include segmentation and record other measurements
#   (e.g., size and maybe average intensity) once I am confident in staining and imaging (i.e., that these would be
#   correspond to some biologically meaningful measure - e.g., nuclear size/promoter expression rather than
#   tissue/depth).


# ------------------------------------------- Chain related functions --------------------------------------------------

def process_directory(directory_path, function_dict_list, prefix, also_return_chain=True, record_name='count',
                     in_method=generate_from_diretory, in_kwargs='generate_from_directory', record_method=img_file_count,
                     record_kwargs=None, write_info_recordkw='file_', write_info_exists=True,
                     extract_method=extract_img_file_counts,
                     extract_kwargs='extract_img_file_counts'):
    """
    FUNCTION: A quantification method should be applied to a sample of volumes from a full data set (or multiple
        similar data sets) prior to use on full data sets. Here, a quantification method refers to a series of functions
        with specific parameters. By implementing the chain method with a specific
    ARGUMENTS: see Chain class for more.
        directory_path: path of the directory that should be processed (str)
        prefix: identifier for the images in the directory. Will be used to find correct files and in saving
            the output (str).
        record_name: a string which can be used as the header for output records in the output data frame in the
            case that extract_img_file_counts is used as the output function.
        also_return_chain: Should a chain object be a second returned object.
    DEPENDENCIES: Pkg/mod/etc: pandas as pd, numpy as np
        Native: (1) generate_from_directory(), (2) img_file_count(), (3) extract_img_file_counts
                via msbrainpy.chain.Chain()
        1) creates generator object which yields images from a directory of files. Only reads files with the
            designated prefix (arg) and the suffix .tif.
        2) returns a list containing the length of an array and some write_info ([record, file_name]).
        3) takes a list of lists of the form [record, file_name] and writes then returns a data frame with
            record and file columns.
    RETURNS: usually a data frame (pd.data_frame) and chain object (Chain), optionally just a data frame.
    """
    if 'generate_from_directory' in in_kwargs:
        in_kwargs = {'prefix': prefix}
    if 'extract_img_file_counts' in extract_kwargs:
        extract_kwargs = {'prefix': prefix, 'directory_path': directory_path, 'record_name': record_name}
    chain = make_process_dir_chain(function_dict_list, in_method=in_method, in_kwargs=in_kwargs, record_method=record_method,
                                record_kwargs=record_kwargs, write_info_recordkw=write_info_recordkw,
                                write_info_exists=write_info_exists, extract_method=extract_method,
                                extract_kwargs=extract_kwargs)
    result = chain.execute(data=directory_path)
    if also_return_chain:
        return result, chain
    else:
        return result


def make_process_dir_chain(function_dict_list, in_method=generate_from_diretory, in_kwargs=None,
                        record_method=img_file_count, record_kwargs=None, write_info_recordkw='file_', write_info_exists=True,
                        extract_method=extract_img_file_counts, extract_kwargs=None):
    template = process_directory_dict(in_method=in_method, in_kwargs=in_kwargs, record_method=record_method,
                                    record_kwargs=record_kwargs, write_info_recordkw=write_info_recordkw,
                                    write_info_exists=write_info_exists, extract_method=extract_method,
                                    extract_kwargs=extract_kwargs)
    template['function_dict_list'] = function_dict_list
    if in_kwargs is None and in_method is img_file_count:
        print('prior to executing chain, please ensure that in_kwargs is given an appropriate value')
    if extract_kwargs is None and extract_method is extract_img_file_counts:
        print('prior to executing chain, please ensure that extract_kwargs is given an appropriate value')
    chain = Chain(**template)
    return chain


def process_directory_dict(in_method=generate_from_diretory, in_kwargs=None, record_method=img_file_count,
                         record_kwargs=None, write_info_recordkw='file_', write_info_exists=True,
                         extract_method=None, extract_kwargs=None):
    template = make_chain_template_dict()
    template['in_method'] = in_method
    template['in_kwargs'] = in_kwargs
    template['record_method'] = record_method
    template['record_kwargs'] = record_kwargs
    template['write_info_exists'] = write_info_exists
    template['write_info_recordkw'] = write_info_recordkw
    template['extract_method'] = extract_method
    template['extract_kwargs'] = extract_kwargs
    return template


# ------------------------------------------ Old unnecessary functions -------------------------------------------------

def sample_cell_counts(samples_path, selem, min_sigma=0.25, max_sigma=5, threshold=0.1,
                     log_name='cell_counts__lo_g.txt', out_prefix=None):
    samples = os.listdir(samples_path)
    cell_counts = []
    files = []
    for _file in samples:
        if _file.endswith('.tif') and _file.find('cells_') == -1:
            files.append(_file)
            loc = os.path.join(samples_path, _file)
            img = imread(loc)
            for i in range(len(img[0])):
                img[i, :, :] = median(img)
            background = opening(img, selem=selem)
            sans_background = img - background
            blobs = blob_log(sans_background, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
            cell_counts.append(len(blobs[:, 0]))
            cells_img = np.zeros(img.shape, dtype=np.uint16)

            for i in range(img.shape[0]):
                for j in range(len(blobs[:, 0])):
                    if blobs[j, 0] == i:
                        r, c, rad = int(blobs[j, 1]), int(blobs[j, 2]), int(np.round(np.sqrt(blobs[j, 3] * 3)))
                        rr, cc = circle_perimeter(r, c, rad)
                        try:
                            cells_img[i, rr, cc] = 65535
                        except index_error:
                            cells_img[i, r, c] = 65535
            if out_prefix == None:
                file_name = 'cells_' + _file
            else:
                file_name = 'cells_' + out_prefix + '_' + _file
            file_path = os.path.join(samples_path, file_name)
            with tiff_writer(file_path) as tif:
                tif.save(cells_img)
            print('{} cells were identified in {}.\n_the locations of these were saved in {}'.format(len(blobs[:, 0]),
                                                                                                    _file, file_name))
    df = pd.data_frame()
    df['file_name'] = files
    df['cell_counts'] = cell_counts
    output_path = os.path.join(samples_path, log_name)
    df.to_csv(output_path, sep='\t')
    return df
