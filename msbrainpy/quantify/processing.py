import os
import re
import numpy as np
from skimage import io
from skimage import exposure
from skimage import filters
from skimage import measure
from skimage import morphology
from skimage import util
from skimage.color import rgb2grey
from skimage.filters import threshold_otsu
from skimage.transform import rescale
from sklearn.decomposition import PCA
# old imports
from skimage.morphology import opening
from skimage.morphology import cube
from skimage.feature import blob_log
from skimage.draw import circle_perimeter
from skimage.filters import median
from tifffile import TiffWriter
from msbrainpy.base import extra_list_nesting


# ------------------------------------------ 2D In Situ Related Functions ----------------------------------------------

def do_gene_series(directory, out_directory, image_name_pattern=r'image_id-\d*.jpeg', find_tissue_mask=False, 
                   save_pc1=True):
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)
    # find the files
    files = []
    image_name_pattern = re.compile(image_name_pattern)
    for file in sorted(os.listdir(directory)):
        match = image_name_pattern.match(file)
        if match is not None:
            files.append(os.path.join(directory, match[0]))
    # iterate though the images
    image_collection = io.imread_collection(files)
    for i in range(len(image_collection)):
        white_balanced = white_balancing(image_collection[i])
        if find_tissue_mask:
            tissue_mask, greyscale = tissue_segmentation(white_balanced)
            mask_name = str(i) + '_tissue_mask.tif'
            mask_path = os.path.join(out_directory, mask_name)
            tissue_name = str(i) + '_slice.tif'
            tissue_path = os.path.join(out_directory, tissue_name)
            io.imsave(mask_path, tissue_mask)
            io.imsave(tissue_path, greyscale)
        segmented, pc1_image = find_expressing_pixels(white_balanced)
        segmented_name = str(i) + '_segmentation.tif'
        segmented_path = os.path.join(out_directory, segmented_name)
        io.imsave(segmented_path, segmented)
        if save_pc1:
            pc1_image = util.img_as_ubyte(pc1_image)
            pc1_image_name = str(i) + '_pc1_greyscale.tif'
            pc1_image_path = os.path.join(out_directory, pc1_image_name)
            io.imsave(pc1_image_path, pc1_image)


def find_expressing_pixels(image_rgb):
    highest_var = rgb_PCA(image_rgb.copy())
    median = filters.median(highest_var)
    edge_enhanced = filters.sobel(median)
    binary = edge_enhanced >= threshold_otsu(edge_enhanced)
    binary = morphology.binary_closing(binary, morphology.disk(3))
    binary = morphology.remove_small_holes(binary)
    return binary, highest_var


def rgb_PCA(rgb_image, n_components=1, get_component=None):
    image_dims = rgb_image.shape
    flattened_image = rgb_image.reshape(-1, 3)
    pca = PCA(n_components=n_components)
    highest_var = pca.fit_transform(flattened_image)
    if get_component is not None:
        highest_var = highest_var[:, get_component]
    highest_var = highest_var.reshape((image_dims[0], image_dims[1]))
    for y in range(highest_var.shape[0]):
        for x in range(highest_var.shape[1]):
            if highest_var[y, x] < 0:
                highest_var[y, x] = 0
    maximum = np.max(highest_var)
    highest_var = highest_var/maximum
    return highest_var


# get masks for tissue location for a series of ISH images
def find_gene_series_masks(directory, out_directory, image_name_pattern=r'image_id-\d*.jpeg'):
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)
    # find the files
    files = []
    image_name_pattern = re.compile(image_name_pattern)
    for file in sorted(os.listdir(directory)):
        match = image_name_pattern.match(file)
        if match is not None:
            files.append(os.path.join(directory, match[0]))
    # iterate though the images
    image_collection = io.imread_collection(files)
    for i in range(len(image_collection)):
        white_balanced = white_balancing(image_collection[i])
        tissue_mask, greyscale = tissue_segmentation(white_balanced)
        mask_name = str(i) + '_tissue_mask.tif'
        mask_path = os.path.join(out_directory, mask_name)
        tissue_name = str(i) + '_slice.tif'
        tissue_path = os.path.join(out_directory, tissue_name)
        io.imsave(mask_path, tissue_mask)
        io.imsave(tissue_path, greyscale)


# Tissue segmentation - this needs to be improved but seems to work for aligning a series into a volume (so far)
def tissue_segmentation(image_rgb, disk_denominator=100, tissue_scale=0.125):
    grey = rgb2grey(image_rgb.copy())
    grey = rescale(grey, tissue_scale)
    grey = exposure.equalize_hist(grey)
    grey_open = morphology.opening(grey, morphology.square(6))
    grey = filters.median(grey)
    threshold = threshold_otsu(grey_open)
    binary = grey_open <= threshold
    binary = morphology.binary_closing(binary, morphology.disk(6))
    labels = measure.label(binary)
    props = measure.regionprops(labels)
    biggest = np.argmax([props[i]['area'] for i in range(len(props))]) + 1
    tissue_mask = labels == biggest
    selm = morphology.disk(np.mean([grey.shape[0], grey.shape[1]]) // disk_denominator)
    tissue_mask = morphology.binary_closing(tissue_mask, selm)
    tissue_mask = morphology.remove_small_holes(tissue_mask, area_threshold=1000)
    grey = util.img_as_ubyte(grey)
    return tissue_mask, grey


def white_balancing(image_rgb):
    image = image_rgb.copy()
    grey = rgb2grey(image_rgb.copy())
    brightest_index = np.unravel_index(np.argmax(grey, axis=None), grey.shape)
    r, g, b = image[:, :, 0], image[:, :, 1], image[:, :, 2]
    brightest_pixel = image[brightest_index[0], brightest_index[1], :]
    wr, wg, wb = brightest_pixel[0], brightest_pixel[1], brightest_pixel[2]
    lum = wr + wg + wb
    r = r * lum / wr
    g = g * lum / wg
    b = b * lum / wb
    return image


# ---------------------------------------- Chain function_dict_list functions --------------------------------------------

def nuclear_detection_list(selem=6, min_sigma=0.25, max_sigma=5, threshold=0.04,
                         write_out=False, directory_path=None, name=None, cube_selem=True):
    """
    FUNCTION: write a list of function dictionaries (as per Chain; chain_method)
    ARGUMENTS:
        selem: structure element for morphological opening. If cube_selem, skimage.morphology.cube() will be applied to
            this arg. In this case, selem should be int.
        min_sigma: standard deviation of the smallest Gaussian kernel to be used in the scale-space representation
            in Laplacian of Gaussian blob detection (skimage.feature.blob_log()). Determines the smallest blobs
            that can be detected.
        max_sigma: standard deviation of the largest Gaussian kernel to be used in the scale-space representation
            in Laplacian of Gaussian blob detection (skimage.feature.blob_log()). Determines the largest blobs
            that can be detected.
        threshold: lower threshold of lo_g local minima that will be detected as a blob. Determines the lowest intensity
            of blob that can be detected
        write_out: should images labeling blob coordinates be written out? (bool)
        directory_path: if write_out, this is the directory to which the blob images will be saved.
        name: if write_out, this name will be saved as part of the blob image file names. Will probably reflect the
            parameters chosen. If write_out and name is None, this will be automatically defined.
        cube_selem: should skimage.morphology.cube() be applied to the inputted selem? (bool)
    DEPENDENCIES: numpy as np (when implemented in chain)
    Native: (1)median_filter3_d()/skimage.filters.median(); (2)remove_background()/skimage.morphology.opening();
        get_blobs__lo_g()/(3)skimage.feature.blob_log(); (4)write_blobs_img()/skimage.draw.circle_perimeter
        1) Median filter sequentially applied to x-y planes along z
        2) Applies morphological opening to generate a background image. This is subtracted from the original and the
            result is returned.
        3) A Laplacian of Gaussian kernel returns an image whereby local minima reflect blobs (~ spherical areas of
            highest local intensity) similar to the size of the kernels sigma. When applied across a scale space, blobs
            of a variety of sizes can be identified. In order to define blobs, some absolute threshold for magnitude of
            minima is chosen. Returns np.ndarray with shape (n, 4) where columns are z, y, x, sigma.
        4) Writes images demarcating the locations of blobs by writing a series of x-y images in which the centre of a
            blob in that plane is encircled by a circle of radius equal to the sigma recorded for the blob.
    RETURNS: list ([{...}, {...}, {...}])
    """
    _selem = selem
    if cube_selem:
        selem = cube(selem)
    med = {'function': median_filter3_d, 'suffix': 'MF', 'write_out': None}
    rem_back = {'function': remove_background, 'suffix': 'sans_bk', 'write_out': None, 'args': [selem]}
    if write_out:
        if name is None:
            name = 'o{}_min{}_max{}_thr{}'.format(_selem, min_sigma, max_sigma, threshold)
        blobs = {'function': get_blobs__lo_g, 'suffix': 'blob_log',
                 'kwargs': {'min_sigma': min_sigma, 'max_sigma': max_sigma, 'threshold': threshold},
                 'write_out': {'function': write_blobs_img, 'prior_link_kw': 'img', 'new_link_kw': 'blobs',
                              'write_info_kw': 'file_', 'kwargs': {'out_dir': directory_path, 'name': name}}}
    else:
        blobs = {'function': get_blobs__lo_g, 'suffix': 'blob_log',
                 'kwargs': {'min_sigma': min_sigma, 'max_sigma': max_sigma, 'threshold': threshold},
                 'write_out': None}
    func_list = [med, rem_back, blobs]
    return func_list


# to be continued...

# ---------------------------------------- Basic image processing functions --------------------------------------------

def get_blobs__lo_g(img, min_sigma, max_sigma, threshold, *args, **kwargs):
    blobs = blob_log(img, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
    print('{} blobs detected'.format(len(blobs)))
    return blobs


def median_filter3_d(img):
    im = img
    for i in range(im.shape[0]):
        im[i, :, :] = median(im[i, :, :])
    return im


def remove_background(data, selem=cube(6), *args, **kwargs):
    img = data
    background = opening(img, selem=selem)
    sans_background = img - background
    return sans_background


# --------------------------------------------- write_out functions ----------------------------------------------------

def write_blobs_img(img, blobs, file_, out_dir, name=None, dtype='auto', **kwargs):
    if dtype == 'auto':
        dtype = type(img[0, 0, 0])
    cells_img = np.zeros(img.shape, dtype=dtype)
    for i in range(img.shape[0]):
        for j in range(len(blobs[:, 0])):
            if blobs[j, 0] == i:
                r, c, rad = int(blobs[j, 1]), int(blobs[j, 2]), int(np.round(blobs[j, 3]))
                rr, cc = circle_perimeter(r, c, rad)
                try:
                    cells_img[i, rr, cc] = 65535
                except index_error:
                    cells_img[i, r, c] = 65535
    if name is None:
        name = ''
    file_name = str(file_[:file_.find('.tif')]) + str(name) + '_blobs.tif'
    file_path = os.path.join(out_dir, file_name)
    with tiff_writer(file_path) as tif:
        tif.save(cells_img)
    print('{} blobs were identified in {}.\n_the locations of these were saved in {}'.format(len(blobs[:, 0]), file_,
                                                                                            file_name))
    return cells_img


# -------------------------------------------- record_method functions --------------------------------------------------

def img_file_count(blobs, file_, **kwargs):
    count = len(blobs)
    return [count, file_]


def correct_for_chunk(coords, chunk, save_dir, prefix, **kwargs):
    z_Pos = chunk['z_start']
    shape0 = chunk['z_end'] - z_Pos
    y_st = chunk['y_start']
    shape1 = chunk['y_end'] - y_st
    x_st = chunk['x_start']
    shape2 = chunk['x_end'] - x_st
    edge = int(np.round(overlap / 2))
    x_lim = shape2 - edge
    y_lim = shape1 - edge
    z_lim = shape0 - edge
    out_of_bounds = []
    for i in range(len(coords[:, 0])):
        if coords[i, 0] > z_lim or coords[i, 1] > y_lim or coords[i, 2] > x_lim:
            out_of_bounds.append(i)
        if z_start != 0:
            if coords[i, 0] < edge:
                out_of_bounds.append(i)
        if y_st != 0:
            if coords[i, 1] < edge:
                out_of_bounds.append(i)
        if x_st != 0:
            if coords[i, 2] < edge:
                out_of_bounds.append(i)
        coords[i, 0] = coords[i, 0] + z_pos
        coords[i, 1] = coords[i, 1] + y_st
        coords[i, 2] = coords[i, 2] + x_st
    coords = np.delete(coords, out_of_bounds, 0)
    if len(coords) != 0:
        save_name = prefix + "_coords_z{}-{}_substack{}".format(z_pos0, chunk['z_end'], chunk['stack_id'])
        save_path = os.path.join(save_dir, save_name)
        np.save(save_path, cells, allow_pickle=False)  # save files throughout as this can be a lengthy process.
        # what if the power disconnects and the computer turns off? Also, need to make sure there is a method for
        # starting part way through z.
        return extra_list_nesting(coords)
    else:
        return None


# ------------------------------------------- extract_method functions --------------------------------------------------

def extract_img_file_counts(record, directory_path, prefix, record_name, **kwargs):
    record = np.array(record)
    df = pd.data_frame({record_name: record[:, 0], 'file': record[:, 1]})
    new_name = prefix + '_' + record_name + '.csv'
    new_path = os.path.join(directory_path, new_name)
    df.to_csv(new_path)
    return df


def extract_chunk_corrected(result, prefix, save_dir, subsection=None):
    result_concat = np.concatenate(result)
    if subsection is not None:
        info = '_{}-x{}:{}-y{}:{}'.format(subsection['stack_id'], subsection['x_start'], subsection['x_end'],
                                                                     subsection['y_start'], subsection['y_end'])
    else:
        info = ''
    name = prefix + info + '_coords.txt'
    out_path = os.path.join(save_dir, name)
    np.savetxt(out_path, result_concat)
    # could add to this. Subsequently remove smaller files perhaps?
    return result_concat
