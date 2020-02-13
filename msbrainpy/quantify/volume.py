import os
import numpy as np
import pandas as pd
from skimage.morphology import opening
from skimage.morphology import cube
from skimage.feature import blob_log
from skimage.filters import median
from skimage.draw import circle_perimeter
from msbrainpy.base import extraListNesting, getSubsections
from msbrainpy.chain import Chain, makeChainTemplateDict
from msbrainpy.io import chunkGenerator
from msbrainpy.quantify.processing import nuclearDetectionList, correctForChunk


# note: volume refers to a large volume which must be . Unfortunately confusing choice of wording
# note: many of the below functions have not been tested yet (except for those in Visualising results and
#     Old (hopefully) unnecessary functions).

# ------------------------------------------ Purpose-specific functions ------------------------------------------------
def sequentialNuclearDetection(data, prefix, outDir, no_y_subs=5, no_x_subs=5, overlap=10,
                               selem=6, min_sigma=0.25, max_sigma=5, threshold=0.04, z_size=50):
    """
    FUNCTION: Apply Chain-based nuclear detection to whole data set by dividing data into subsections, which are
        sequentially processed.
    ...
    """
    subsectionList = getSubsections(data.shape, y=no_y_subs, x=no_x_subs, overlap=overlap)
    allCoords = []
    for subsection in subsectionList:
        sub_result, chain = volumeNuclearDetection(data, prefix, outDir, subsection=subsection, selem=selem,
                                                min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold,
                                                z_size=z_size)
        allCoords.append([sub_result])

    allCoords = np.concatenate(allCoords)
    name = prefix+'_allSubsectionCoords.txt'
    savePath = os.path.join(outDir, name)
    np.savetxt(savePath, allCoords)
    return allCoords


def volumeNuclearDetection(data, prefix, saveDir, subsection=None, selem=6, min_sigma=0.25, max_sigma=5, threshold=0.04,
                           z_size=50, alsoReturnChain=True):
    """
    FUNCTION: Apply Chain-based nuclear detection to volume of HDf5 dataset. Particular subsection can be specified
        using subsections produced by the getSubsections() function or a dictionary of the same format.
    """
    functionDictList = nuclearDetectionList(selem=selem, min_sigma=min_sigma, max_sigma=max_sigma,
                                            threshold=threshold)
    result, chain = processVolume(data, functionDictList, subsection, prefix, saveDir, z_size=z_size)
    if alsoReturnChain:
        return result, chain
    else:
        return result


# ------------------------------------------- Chain related functions --------------------------------------------------
# extractChunkCorrected(result, prefix, saveDir, subsection=None)
def processVolume(data, functionDictList, subsection, prefix, saveDir, z_size=50,
                  in_kwargs='chunkGenerator', recordMethod=correctForChunk, record_kwargs='correctForChunk',
                  writeInfo_recordkw='chunk', writeInfoExists=True, extractMethod=extractChunkCorrected,
                  extract_kwargs='extractChunkCorrected', alsoReturnChain=True, noSubsection_overlap=0):
    """
    FUNCTION: process a volume from a hdf5 data set by reading portions of the data (along z) into memory at a time.
        this is designed to work with the dictionaries produced by the getSubsections function (which produces coords
        for portions of the data, depending on how many tiles should be produced in x
    ARGUMENTS: for all possible arguments see Chain.__init__.__doc__ or makeChain.__doc__
         data: a hdf5 dataset from which a part or all of data should be processed
         subsection (chunkGenerator/extractChunkCorrected): if only a subsection of the data should be examined,
            provide a subsection dictionary. Subsection dictionaries are generated by the
            msbrainpy.base.getSubsections() function.
         z_size (chunkGenerator): if using the chunkGenerator function, as is default, the image will be processed in a
            series of overlapping blocks in z. z_size is the voxel depth of these blocks.
         prefix (correctForChunk/extractChunkCorrected): if using correctForChunk() as the record method
            (as is default), the output of each iteration through chain, will be assumed to be an array containing
            coordinates in :, 0-2 (other cols may contain information about objects at these coordinates).
            The function will account for the overlap using info about the chunk (which should be the writeInfo
            in this case - provided by the generator inMethod). The correctForChunk function saves smaller data npy
            files throughout the process (as this may be lengthy). The prefix is used to name these files.
            The names also contain information about the substack number (ID) and the position in z. The prefix is also
            used to name the final concatenated file in the extractChunkCorrected() function.
         saveDir (correctForChunk): The output directory for the smaller files and the final output file.
    DEPENDENCIES: Pkg/mod/etc: depends
    Native: makeProcessVolumeChain() / processVolumeDict() + Chain class, etc.
    RETURNS: if used with default output settings, coordinates, which will be corrected for overlaps in subsections.
        type = np.ndarray
    """
    if subsection is none:
        subsection = {'stackID': 0, 'x_start': 0, 'x_end': data.shape[2], 'y_start': 0, 'y_end': data.shape[1],
                      'overlap': noSubsection_overlap} # one limitation here is that the overlap that will be corrected
        #                                                   for will be this number. This means that objects at the
        #                                                   edges will be eliminated when finding coordinates
    inMethod = chunkGenerator
    if 'chunkGenerator' in in_kwargs:
        in_kwargs = {'subsection': subsection, 'z_size': z_size, 'dowhat': None, 'writeInfoExists': writeInfoExists}
    if 'correctForChunk' in record_kwargs:
        record_kwargs = {'saveDir': saveDir, 'prefix': prefix}
    if 'extractChunkCorrected' in extract_kwargs:
        record_kwargs = {'saveDir': saveDir, 'prefix': prefix, 'subsection': subsection}
    chain = makeProcessVolumeChain(functionDictList, inMethod=inMethod, in_kwargs=in_kwargs, recordMethod=recordMethod,
                                   record_kwargs=record_kwargs, writeInfo_recordkw=writeInfo_recordkw,
                                   writeInfoExists=writeInfoExists, extractMethod=extractMethod,
                                   extract_kwargs=extract_kwargs)
    result = chain.execute(data=data)
    if alsoReturnChain:
        return result, chain
    else:
        return result


# just in case (below).

def makeProcessVolumeChain(functionDictList, inMethod=chunkGenerator, in_kwargs=None, recordMethod=correctForChunk,
                           record_kwargs=None, writeInfo_recordkw='chunk', writeInfoExists=True,
                           extractMethod=extractChunkCorrected, extract_kwargs=None):
    template = processVolumeDict(inMethod=inMethod, in_kwargs=in_kwargs, recordMethod=recordMethod,
                                 record_kwargs=record_kwargs, writeInfo_recordkw=writeInfo_recordkw,
                                 writeInfoExists=writeInfoExists, extractMethod=extractMethod,
                                 extract_kwargs=extract_kwargs)
    template['functionDictList'] = functionDictList
    if in_kwargs is None and inMethod is chunkGenerator:
        print('prior to executing chain, please ensure that in_kwargs is given an appropriate value')
    if record_kwargs is None and recordMethod is chunkGenerator:
        print('prior to executing chain, please ensure that record_kwargs is given an appropriate value')
    if extract_kwargs is None and extractMethod is extractChunkCorrected:
        print('prior to executing chain, please ensure that extract_kwargs is given an appropriate value')
    chain = Chain(**template)
    return chain


def processVolumeDict(inMethod=chunkGenerator, in_kwargs=None, recordMethod=correctForChunk,
                      record_kwargs=None, writeInfo_recordkw='chunk', writeInfoExists=True,
                      extractMethod=np.concatenate, extract_kwargs=None):
    template = makeChainTemplateDict()
    template['inMethod'] = inMethod
    template['in_kwargs'] = in_kwargs
    template['recordMethod'] = recordMethod
    template['record_kwargs'] = record_kwargs
    template['writeInfoExists'] = writeInfoExists
    template['writeInfo_recordkw'] = writeInfo_recordkw
    template['extractMethod'] = extractMethod
    template['extract_kwargs'] = extract_kwargs
    return template


# Sample-to-volume chain conversion ------------------------
def processVolumeFromChain(chain, some, more, arguments):
    pass


def convertChainToVolume(chain):
    chain_ = chain
    pass


# ---------------------------------------------- Visualising results ---------------------------------------------------

def getCellsInPlane(shape, points, order, prefix, index, plsmns):
    cellsImg = np.zeros(shape, dtype=np.uint16)
    blobs = []
    col = order[0]
    for i in range(len(points[:, col])):
        point = int(np.round(points[i, col]))
        if index - plsmns <= point <= index + plsmns:
            blobs.append([points[i, :]])
    blobs = np.concatenate(blobs)
    print(len(blobs))
    for i in range(len(blobs)):
        y = int(np.round(blobs[i, order[1]]))
        x = int(np.round(blobs[i, order[2]]))
        r, c, rad = y, x, 1
        rr, cc = circle_perimeter(r, c, rad)
        try:
            cellsImg[rr, cc] = 65535
        except IndexError:
            print('point {}, {}, out of range'.format(y, x))

    saveName = prefix + '_' + str(col) + '_' + str(index) + '.tif'
    with TiffWriter(saveName) as tif:
        tif.save(cellsImg)
    return cellsImg


# -------------------------------------- Old (hopefully) unnecessary functions -----------------------------------------

def cellsOnly(data, subsection, prefix, outDir, z_size=50,
              min_sigma=0.25, max_sigma=10, threshold=0.04, **kwargs):
    """
    FUNCTION: function to sequentially process chunks to retrive cell coordinates
    ARGUMENTS:
        data = hdf5 data object
        z_size = z size of chunk (int)
        overlap = overlap of chunk (int)
        min_sigma = minimum sigma for LoG blob detection (smaller for detecting smaller blobs)
        max_sigma = maximum sigma for LoG blob detection (larger for detecting larger blobs)
        threshold = intensity threashold for blob detection (lower for fainter cells)
    RETURNS: array with cell coordinates [[z, y, x, sigma], ...]
    """
    sansBackground = chunkGenerator(data, subsection=subsection, z_size=z_size,
                                    dowhat=removeBackground_)
    cellArray = cellsFromBlobChunks(sansBackground, subsection, prefix, outDir, overlap=overlap,
                                    min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)

    return cellArray


def cellsFromBlobChunks(sansBackground, subsection, prefix, outDir, overlap, min_sigma=0.25, max_sigma=10,
                        threshold=0.04):
    cellArray = []
    y_st = subsection['y_start']
    x_st = subsection['x_start']
    stackID = subsection['stackID']
    zPos = 0
    edge = int(np.round(overlap / 2))
    tempName = prefix + "_substack" + str(stackID) + "_cells"
    tempPath = os.path.join(outDir, tempName)
    os.mkdir(tempPath)
    for chunk in sansBackground:
        zPos0 = zPos
        img = chunk
        size = img.shape
        cells = blob_log(img, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
        x_lim = size[2] - edge
        y_lim = size[1] - edge
        z_lim = size[0] - edge
        print("{} cells were initially found at z {}:{} in substack {}".format(len(cells), zPos0, zPos0 + size[0],
                                                                               stackID))
        outOfBounds = []
        for i in range(len(cells[:, 0])):
            if cells[i, 0] > z_lim or cells[i, 1] > y_lim or cells[i, 2] > x_lim:
                outOfBounds.append(i)
            if zPos != 0:
                if cells[i, 0] < edge:
                    outOfBounds.append(i)
            if y_st != 0:
                if cells[i, 1] < edge:
                    outOfBounds.append(i)
            if x_st != 0:
                if cells[i, 2] < edge:
                    outOfBounds.append(i)

            cells[i, 0] = cells[i, 0] + zPos
            cells[i, 1] = cells[i, 1] + y_st
            cells[i, 2] = cells[i, 2] + x_st

        cells = np.delete(cells, outOfBounds, 0)

        if len(cells) != 0:
            cellArray.append(cells)
            saveName = "cells_z{}-{}_substack{}".format(zPos0, zPos0 + size[0], stackID)
            savePath = os.path.join(tempPath, saveName)
            np.save(savePath, cells, allow_pickle=False)

        zPos += size[0] - overlap
        print("{} cells were found at z {}:{} in substack {}".format(len(cells), zPos0, zPos0 + size[0], stackID))

    print('A total of {} cells were found in substack {}'.format(len(cellArray), stackID))
    cellArray = np.concatenate(cellArray)

    cellsFile = tempName + '.npy'
    cellsPath = os.path.join(outDir, cellsFile)
    np.save(cellsPath, cellArray, allow_pickle=False)

    return cellArray


def removeBackground_(data, selem=cube(6)):
    img = data
    for i in range(img.shape[0]):
        img[i, :, :] = median(img[i, :, :])
    background = opening(img)
    sansBackground = img - background
    return sansBackground
