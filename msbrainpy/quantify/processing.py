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
from msbrainpy.base import extraListNesting


# ------------------------------------------ 2D In Situ Related Functions ----------------------------------------------

def do_gene_series(directory, out_directory, 
                   image_name_pattern=r'image_id-\d*.jpeg', 
                   find_tissue_mask=True, save_pc1=True):
    """
    Process all images in an ish series (i.e., one brain). 
    

    Parameters
    ----------
    directory: str
        path to the series of images
    out_directory: str
        path to which to save the output
    image_name_pattern: r string
        pattern of images to process. Defaults to the naming convention used
        when downloading images from the Allen Brain database 
            see: msbrainpy.query.download_gene_images
    find_tissue_mask: bool
        Should a downsampled, contrast enhanced image + tissue mask be produced
    save_pc1: bool
        Shouls a pc1 signal image be produced. This image represents the 
        RGB combination (~ in situ dye) associated with the greatest shared
        varience in the data
    """
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
        if save_pc1:
            pc1_image = rgb_PCA(white_balanced)
            pc1_image = util.img_as_ubyte(pc1_image)
            pc1_image_name = str(i) + '_pc1_greyscale.tif'
            pc1_image_path = os.path.join(out_directory, pc1_image_name)
            io.imsave(pc1_image_path, pc1_image)
        # add in segmentation function


def ish_expressing_pixels(image):
    """
    function to segment expressing pixles from pre-processed in situ images
    """
    # this needs to be ~ equivalent to method used in allen brain institutes 
    # neuroinformatics/anatomical gene expression pipeline
    # Method must segment ~ cell sized objects and dense structures (e.g., 
    # hippocampus, granular layer cerebellum, etc.)
    # Can't be contrast based as this changes between images and could lead 
    # to false detection in cellularly dense areas, which tend to pick up the
    # dye without appreciable expression.
    pass


def rgb_PCA(rgb_image, n_components=1, get_component=None):
    """
    Get a pca-based intensity image. Designed to be used to extract highest
    shared variance components from RGB colour channels. PCA implemented
    via scikit-learn.
    E.g., extract a component corresponding to the dye colour in an in situ 
    image (purple-blue colour comprised shared R~B variance)

    Parameters
    ----------
    rgb_image: numpy.ndarray
    n_components: int
        number of components to include in the pca model
    get_component: None or int or tupule of int or slice
        The componets that should be retrived from the PCA. 
        pca shape is (num-pixels, n_components)
        get components chooses components --> pixel values 

    Notes
    -----
    Have only tested when retriving one image. 
    Quite slow.
    References:
    [1] The Image Processing Handbook, John C. Russ and F. Brent Neal. 
    CRC Press, Boca Raton, FL, 2015, 1053 pp. <include pages>
    ISBN: 978-1498740265.
    [2] Neal, B. and Russ, J.C., 2004. Principal components analysis of 
    multispectral image data. Microscopy Today, 12(5), pp.36-39.
    """
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
    '''
    Apply tissue_segmentation function to a gene series
    '''
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
    '''
    Segement a tissue mask from a downsampled slice image

    Parameters
    ----------
    image_rgb: np.ndarray
        rgb in situ image (full size)
    disk_denominator: int
        used to choose the selm for binary closing morphological opperation
        Used as a scaling factor for the disk
        bigger >> smaller disk 
        larger >> bigger disk
    tissue_scale: rescaling factor for the image

    Notes
    -----
    Needs some work. There are some holes in the masks where the tissue is 
    not above background. Some of the circular artifacts are covered by the 
    closing opperation. Perfectly sufficient for aligning images when building
    a volume.
    '''
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


# ---------------------------------------- Chain functionDictList functions --------------------------------------------

def nuclearDetectionList(selem=6, min_sigma=0.25, max_sigma=5, threshold=0.04,
                         writeOut=False, directoryPath=None, name=None, cube_selem=True):
    """
    FUNCTION: write a list of function dictionaries (as per Chain; ChainMethod)
    ARGUMENTS:
        selem: structure element for morphological opening. If cube_selem, skimage.morphology.cube() will be applied to
            this arg. In this case, selem should be int.
        min_sigma: standard deviation of the smallest Gaussian kernel to be used in the scale-space representation
            in Laplacian of Gaussian blob detection (skimage.feature.blob_log()). Determines the smallest blobs
            that can be detected.
        max_sigma: standard deviation of the largest Gaussian kernel to be used in the scale-space representation
            in Laplacian of Gaussian blob detection (skimage.feature.blob_log()). Determines the largest blobs
            that can be detected.
        threshold: lower threshold of LoG local minima that will be detected as a blob. Determines the lowest intensity
            of blob that can be detected
        writeOut: should images labeling blob coordinates be written out? (bool)
        directoryPath: if writeOut, this is the directory to which the blob images will be saved.
        name: if writeOut, this name will be saved as part of the blob image file names. Will probably reflect the
            parameters chosen. If writeOut and name is None, this will be automatically defined.
        cube_selem: should skimage.morphology.cube() be applied to the inputted selem? (bool)
    DEPENDENCIES: numpy as np (when implemented in chain)
    Native: (1)medianFilter3D()/skimage.filters.median(); (2)removeBackground()/skimage.morphology.opening();
        getBlobs_LoG()/(3)skimage.feature.blob_log(); (4)writeBlobsImg()/skimage.draw.circle_perimeter
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
    med = {'function': medianFilter3D, 'suffix': 'MF', 'writeOut': None}
    remBack = {'function': removeBackground, 'suffix': 'sansBk', 'writeOut': None, 'args': [selem]}
    if writeOut:
        if name is None:
            name = 'o{}_min{}_max{}_thr{}'.format(_selem, min_sigma, max_sigma, threshold)
        blobs = {'function': getBlobs_LoG, 'suffix': 'blob_log',
                 'kwargs': {'min_sigma': min_sigma, 'max_sigma': max_sigma, 'threshold': threshold},
                 'writeOut': {'function': writeBlobsImg, 'priorLink_kw': 'img', 'newLink_kw': 'blobs',
                              'writeInfo_kw': 'file_', 'kwargs': {'outDir': directoryPath, 'name': name}}}
    else:
        blobs = {'function': getBlobs_LoG, 'suffix': 'blob_log',
                 'kwargs': {'min_sigma': min_sigma, 'max_sigma': max_sigma, 'threshold': threshold},
                 'writeOut': None}
    funcList = [med, remBack, blobs]
    return funcList


# to be continued...

# ---------------------------------------- Basic image processing functions --------------------------------------------

def getBlobs_LoG(img, min_sigma, max_sigma, threshold, *args, **kwargs):
    blobs = blob_log(img, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
    print('{} blobs detected'.format(len(blobs)))
    return blobs


def medianFilter3D(img):
    im = img
    for i in range(im.shape[0]):
        im[i, :, :] = median(im[i, :, :])
    return im


def removeBackground(data, selem=cube(6), *args, **kwargs):
    img = data
    background = opening(img, selem=selem)
    sansBackground = img - background
    return sansBackground


# --------------------------------------------- writeOut functions ----------------------------------------------------

def writeBlobsImg(img, blobs, file_, outDir, name=None, dtype='auto', **kwargs):
    if dtype == 'auto':
        dtype = type(img[0, 0, 0])
    cellsImg = np.zeros(img.shape, dtype=dtype)
    for i in range(img.shape[0]):
        for j in range(len(blobs[:, 0])):
            if blobs[j, 0] == i:
                r, c, rad = int(blobs[j, 1]), int(blobs[j, 2]), int(np.round(blobs[j, 3]))
                rr, cc = circle_perimeter(r, c, rad)
                try:
                    cellsImg[i, rr, cc] = 65535
                except IndexError:
                    cellsImg[i, r, c] = 65535
    if name is None:
        name = ''
    fileName = str(file_[:file_.find('.tif')]) + str(name) + '_blobs.tif'
    filePath = os.path.join(outDir, fileName)
    with TiffWriter(filePath) as tif:
        tif.save(cellsImg)
    print('{} blobs were identified in {}.\nThe locations of these were saved in {}'.format(len(blobs[:, 0]), file_,
                                                                                            fileName))
    return cellsImg


# -------------------------------------------- recordMethod functions --------------------------------------------------

def imgFileCount(blobs, file_, **kwargs):
    count = len(blobs)
    return [count, file_]


def correctForChunk(coords, chunk, saveDir, prefix, **kwargs):
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
    outOfBounds = []
    for i in range(len(coords[:, 0])):
        if coords[i, 0] > z_lim or coords[i, 1] > y_lim or coords[i, 2] > x_lim:
            outOfBounds.append(i)
        if z_start != 0:
            if coords[i, 0] < edge:
                outOfBounds.append(i)
        if y_st != 0:
            if coords[i, 1] < edge:
                outOfBounds.append(i)
        if x_st != 0:
            if coords[i, 2] < edge:
                outOfBounds.append(i)
        coords[i, 0] = coords[i, 0] + zPos
        coords[i, 1] = coords[i, 1] + y_st
        coords[i, 2] = coords[i, 2] + x_st
    coords = np.delete(coords, outOfBounds, 0)
    if len(coords) != 0:
        saveName = prefix + "_coords_z{}-{}_substack{}".format(zPos0, chunk['z_end'], chunk['stackID'])
        savePath = os.path.join(saveDir, saveName)
        np.save(savePath, cells, allow_pickle=False)  # save files throughout as this can be a lengthy process.
        # what if the power disconnects and the computer turns off? Also, need to make sure there is a method for
        # starting part way through z.
        return extraListNesting(coords)
    else:
        return None


# ------------------------------------------- extractMethod functions --------------------------------------------------

def extractImgFileCounts(record, directoryPath, prefix, recordName, **kwargs):
    record = np.array(record)
    df = pd.DataFrame({recordName: record[:, 0], 'file': record[:, 1]})
    newName = prefix + '_' + recordName + '.csv'
    newPath = os.path.join(directoryPath, newName)
    df.to_csv(newPath)
    return df


def extractChunkCorrected(result, prefix, saveDir, subsection=None):
    result_concat = np.concatenate(result)
    if subsection is not None:
        info = '_{}-x{}:{}-y{}:{}'.format(subsection['stackID'], subsection['x_start'], subsection['x_end'],
                                                                     subsection['y_start'], subsection['y_end'])
    else:
        info = ''
    name = prefix + info + '_coords.txt'
    outPath = os.path.join(saveDir, name)
    np.savetxt(outPath, result_concat)
    # could add to this. Subsequently remove smaller files perhaps?
    return result_concat
