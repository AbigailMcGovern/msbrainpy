import nrrd
import os
import h5py
import re
from tifffile import imread
import numpy as np
import pandas as pd
from tifffile import TiffWriter
from skimage import io
from skimage import util
from skimage.transform import rescale

# -------------------------------------------- Read tiffs from directory -----------------------------------------------
def imread_i(directory, i):
    files = os.listdir(directory)
    file_ = files[i]
    filepath = os.path.join(directory, file_)
    img = imread(filepath)
    return img


def generate_from_diretory(directory, prefix=None, **kwargs):
    if kwargs:
        prefix = kwargs.get('prefix')
    allfiles = os.listdir(directory)
    files = []
    if prefix is None:
        for _file in allfiles:
            if _file.endswith('.tif'):
                files.append(_file)
    if prefix is not None:
        for _file in allfiles:
            if _file.endswith('.tif') and _file.startswith(prefix):
                files.append(_file)
    for _file in files:
        file_path = os.path.join(directory, _file)
        img = imread(file_path)
        # out = {0: img, 1: _file}
        # print(out[1])
        yield img, _file


# -------------------------------------------------- Atlas related -----------------------------------------------------
def read_annotation(ann_name, crop=None):  # , ori = None):
    annotation, header = nrrd.read(ann_name)
    print("nrrd header for the annotation file:")
    print(header)
    annotation = annotation.astype('int32')
    print('Annotation file has been converted to int32')
    # if ori != None:
    # (0, 1, 2) = z:D-V, y:R-C, x:L-R
    # move indexes to swap, -ve to flip
    # if ori == (1, 2, 0):
    #
    print('the annotation has a shape of {}, {}, {}'.format(annotation.shape[0], annotation.shape[1],
                                                            annotation.shape[2]))
    if crop != None:
        # use with np.s_[:, :, :]
        annotation = annotation[crop]
        print('the annotation was cropped to shape {}, {}, {}'.format(annotation.shape[0], annotation.shape[1],
                                                                      annotation.shape[2]))
    return annotation


# ---------------------------------------------- tera_stitcher output ---------------------------------------------------
def write_cropped_tiff_stack(directory, out_dir, y_ind, x_ind):
    file_list = os.listdir(directory)
    for _file in file_list:
        im_path = os.path.join(directory, _file)
        image = imread(im_path)
        image = image[:, y_ind[0]:y_ind[1], x_ind[0]:x_ind[1]]
        im_path = os.path.join(out_dir, _file)
        print('for {}, the new image shape is {}'.format(_file, image.shape))
        with tiff_writer(im_path, bigtiff=True) as tif:
            for i in range(image.shape[0]):
                tif.save(image[i, :, :])
        print('the image {} was saved at:\n {}'.format(_file, im_path))

# ----------------------------------------------------- H5 data --------------------------------------------------------
def write_hdf5_seq(filename, filedir, dataset_name, img_dir):
    '''
    FUNCTION: write a HDF5 file from a series of 3D images (stacked in Z). Writes hdf5 file and a table describing where
        each tiff lives in the hdf5
    ARGUMENTS:
        filename = name for hdf5 output (str)
        filedir = output directory
        dataset_name = name to which the  data will be assigned within the hdf5 file (str)
        img_dir = directory in which tiffs are sequentially listed (str)
    RETURNS: hdf5 dataset object
    '''
    filepath = os.path.join(filedir, filename)
    files = []
    for _file in os.listdir(img_dir):
        if _file.endswith('.tif'):
            files.append(_file)

    im1path = os.path.join(img_dir, files[0])
    im1 = imread(im1path)
    shape1 = im1.shape
    del im1

    with h5py.File(filepath, 'w') as h5:
        dset = h5.create_dataset(dataset_name, shape1, dtype=np.uint16, maxshape=(None, shape1[1], shape1[2]),
                                 chunks=True)
        started = False
        start_at = 0
        startpos = []
        endpos = []
        for _file in files:
            img_path = os.path.join(img_dir, _file)
            img = imread(img_path)

            if started == False:
                dset[start_at:start_at + img.shape[0], :, :] = img
                started = True
            else:
                dset.resize(dset.shape[0] + img.shape[0], axis=0)
                dset[start_at:start_at + img.shape[0], :, :] = img

            print('{} was saved at Z-indicies {}-{} in {}'.format(_file, start_at, dset.shape[0], filepath))
            startpos.append(start_at)
            endpos.append(dset.shape[0])
            start_at += img.shape[0]

    positions = pd.data_frame()
    positions['Files'] = files
    positions['start_in_h5'] = startpos
    positions['end_in_h5'] = endpos
    positions.to_csv(os.path.join(img_dir, '{}z_table.txt'.format(filename)), sep='\t')
    return dset


def open_h5data(filedir, filename, dataset=True):
    filepath = os.path.join(filedir, filename)
    h5file = h5py.File(filepath, 'r')
    print('A hdf5 file was found at {}. \n_this contains the following datasets:'.format(filepath))
    for key in h5file.keys():
        print(key)
    keys = list(h5file.keys())
    data = h5file[keys[0]]
    print('dataset: {} \nshape: {}, {}, {}'.format(keys[0], data.shape[0], data.shape[1], data.shape[2]))
    return data

# --------------------------------------------- writing from H5 data ---------------------------------------------------
def save_random(data, n, size, prefix, save_dir):
    shape = data.shape
    sample_z = np.random.choice(shape[0] - size[0], n)
    sample_y = np.random.choice(shape[1] - size[1], n)
    sample_x = np.random.choice(shape[2] - size[2], n)
    for i in range(n):
        z_end = sample_z[i] + size[0]
        y_end = sample_y[i] + size[1]
        x_end = sample_x[i] + size[2]
        sample = data[sample_z[i]:z_end, sample_y[i]:y_end, sample_x[i]:x_end]
        file_coord = '_z{}-{}_y{}-{}_x{}-{}.tif'.format(
            sample_z[i], z_end, sample_y[i], y_end, sample_x[i], x_end)
        file_name = prefix + file_coord
        file_path = os.path.join(save_dir, file_name)
        with tiff_writer(file_path) as tif:
            tif.save(sample)


def write_subs_tiff(filedir, filename, data, subsection):
    x_st = subsection.pop('x_start')
    x_en = subsection.pop('x_end')
    y_st = subsection.pop('y_start')
    y_en = subsection.pop('x_end')
    img = data[:, y_st:y_en, x_st:x_en]
    filepath = os.path.join(filedir, filename)
    with tiff_writer(filepath, bigtiff=True) as tiff:
        for i in range(img.shape[0]):
            tiff.save(img[i, :, :])
    return img


# ------------------------------------------------ 2D slice images ----------------------------------------------------

def save_resized_from_directory(directory, image_pattern=r'\d*_pc1_greyscale.tif', scale=0.125, verbose=False,
                                anti_aliasing=None):
    image_list = []
    image_pattern_re = re.compile(image_pattern)
    for file in os.listdir(directory):
        match = image_pattern_re.match(file)
        if match is not None:
            image_list.append(match[0])
    if verbose:
        print('files will be saved in the following location:')
        print(directory)
    for file in image_list:
        image = io.imread(os.path.join(directory, file))
        image = rescale(image, scale, anti_aliasing=anti_aliasing)
        image = util.img_as_ubyte(image)
        new_name = file[:file.find('.tif')] + '_scale-' + str(scale) + '.tif'
        if verbose:
            print('Saving file:', new_name)
        new_name = os.path.join(directory, new_name)
        io.imsave(new_name, image)
