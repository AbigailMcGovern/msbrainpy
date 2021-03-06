import os
import re
import numpy as np
import pandas as pd
from skimage import io
from skimage import exposure
from skimage import filters
from skimage import measure
from skimage import morphology
from skimage import util
from skimage.color import rgb2grey
from skimage.filters import threshold_otsu
from skimage.transform import rescale
# old imports
from skimage.morphology import opening
from skimage.morphology import cube
from skimage.feature import blob_log
from skimage.draw import circle_perimeter
from skimage.filters import median
from tifffile import TiffWriter
from msbrainpy.base import extraListNesting


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


def correctForChunk(coords, chunk, saveDir, prefix, overlap, **kwargs):
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
        if z_Pos != 0:
            if coords[i, 0] < edge:
                outOfBounds.append(i)
        if y_st != 0:
            if coords[i, 1] < edge:
                outOfBounds.append(i)
        if x_st != 0:
            if coords[i, 2] < edge:
                outOfBounds.append(i)
        coords[i, 0] = coords[i, 0] + z_Pos
        coords[i, 1] = coords[i, 1] + y_st
        coords[i, 2] = coords[i, 2] + x_st
    coords = np.delete(coords, outOfBounds, 0)
    if len(coords) != 0:
        saveName = prefix + "_coords_z{}-{}_substack{}".format(z_Pos, chunk['z_end'], chunk['stackID'])
        savePath = os.path.join(saveDir, saveName)
        np.save(savePath, coords, allow_pickle=False)  # save files throughout as this can be a lengthy process.
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
