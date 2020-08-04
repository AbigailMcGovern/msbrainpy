import os
import pandas as pd
from tifffile import imread
from skimage.feature import blob_log
from msbrainpy.lightsheet.chain import Chain, makeChainTemplateDict
from msbrainpy.lightsheet.io import generateFromDiretory
from msbrainpy.lightsheet.quantify.processing import nuclearDetectionList, imgFileCount, extractImgFileCounts, \
    extractChunkCorrected, 


# ------------------------------------------ Purpose-specific functions ------------------------------------------------

def sampleNuclearDetection(directoryPath, prefix, name, selem=6, min_sigma=0.25,
                           max_sigma=5, threshold=0.04):
    funcList = nuclearDetectionList(selem=selem, min_sigma=min_sigma, max_sigma=max_sigma,
                                    threshold=threshold, directoryPath=directoryPath, name=name, writeOut=True)
    df = processDirectory(directoryPath, funcList, prefix)
    return df


# to be continued ...
# specifically, may want to extend this process to include segmentation and record other measurements
#   (e.g., size and maybe average intensity) once I am confident in staining and imaging (i.e., that these would be
#   correspond to some biologically meaningful measure - e.g., nuclear size/promoter expression rather than
#   tissue/depth).


# ------------------------------------------- Chain related functions --------------------------------------------------

def processDirectory(directoryPath, functionDictList, prefix, alsoReturnChain=True, recordName='count',
                     inMethod=generateFromDiretory, in_kwargs='generateFromDirectory', recordMethod=imgFileCount,
                     record_kwargs=None, writeInfo_recordkw='file_', writeInfoExists=True,
                     extractMethod=extractImgFileCounts,
                     extract_kwargs='extractImgFileCounts'):
    """
    FUNCTION: A quantification method should be applied to a sample of volumes from a full data set (or multiple
        similar data sets) prior to use on full data sets. Here, a quantification method refers to a series of functions
        with specific parameters. By implementing the chain method with a specific
    ARGUMENTS: see Chain class for more.
        directoryPath: path of the directory that should be processed (str)
        prefix: identifier for the images in the directory. Will be used to find correct files and in saving
            the output (str).
        recordName: a string which can be used as the header for output records in the output data frame in the
            case that extractImgFileCounts is used as the output function.
        alsoReturnChain: Should a chain object be a second returned object.
    DEPENDENCIES: Pkg/mod/etc: pandas as pd, numpy as np
        Native: (1) generateFromDirectory(), (2) imgFileCount(), (3) extractImgFileCounts
                via msbrainpy.chain.Chain()
        1) creates generator object which yields images from a directory of files. Only reads files with the
            designated prefix (arg) and the suffix .tif.
        2) returns a list containing the length of an array and some writeInfo ([record, fileName]).
        3) takes a list of lists of the form [record, fileName] and writes then returns a data frame with
            record and file columns.
    RETURNS: usually a data frame (pd.DataFrame) and chain object (Chain), optionally just a data frame.
    """
    if 'generateFromDirectory' in in_kwargs:
        in_kwargs = {'prefix': prefix}
    if 'extractImgFileCounts' in extract_kwargs:
        extract_kwargs = {'prefix': prefix, 'directoryPath': directoryPath, 'recordName': recordName}
    chain = makeProcessDirChain(functionDictList, inMethod=inMethod, in_kwargs=in_kwargs, recordMethod=recordMethod,
                                record_kwargs=record_kwargs, writeInfo_recordkw=writeInfo_recordkw,
                                writeInfoExists=writeInfoExists, extractMethod=extractMethod,
                                extract_kwargs=extract_kwargs)
    result = chain.execute(data=directoryPath)
    if alsoReturnChain:
        return result, chain
    else:
        return result


def makeProcessDirChain(functionDictList, inMethod=generateFromDiretory, in_kwargs=None,
                        recordMethod=imgFileCount, record_kwargs=None, writeInfo_recordkw='file_', writeInfoExists=True,
                        extractMethod=extractImgFileCounts, extract_kwargs=None):
    template = processDirectoryDict(inMethod=inMethod, in_kwargs=in_kwargs, recordMethod=recordMethod,
                                    record_kwargs=record_kwargs, writeInfo_recordkw=writeInfo_recordkw,
                                    writeInfoExists=writeInfoExists, extractMethod=extractMethod,
                                    extract_kwargs=extract_kwargs)
    template['functionDictList'] = functionDictList
    if in_kwargs is None and inMethod is imgFileCount:
        print('prior to executing chain, please ensure that in_kwargs is given an appropriate value')
    if extract_kwargs is None and extractMethod is extractImgFileCounts:
        print('prior to executing chain, please ensure that extract_kwargs is given an appropriate value')
    chain = Chain(**template)
    return chain


def processDirectoryDict(inMethod=generateFromDiretory, in_kwargs=None, recordMethod=imgFileCount,
                         record_kwargs=None, writeInfo_recordkw='file_', writeInfoExists=True,
                         extractMethod=None, extract_kwargs=None):
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


# ------------------------------------------ Old unnecessary functions -------------------------------------------------

def sampleCellCounts(samplesPath, selem, min_sigma=0.25, max_sigma=5, threshold=0.1,
                     logName='CellCounts_LoG.txt', outPrefix=None):
    samples = os.listdir(samplesPath)
    cellCounts = []
    files = []
    for _file in samples:
        if _file.endswith('.tif') and _file.find('cells_') == -1:
            files.append(_file)
            loc = os.path.join(samplesPath, _file)
            img = imread(loc)
            for i in range(len(img[0])):
                img[i, :, :] = median(img)
            background = opening(img, selem=selem)
            sansBackground = img - background
            blobs = blob_log(sansBackground, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
            cellCounts.append(len(blobs[:, 0]))
            cellsImg = np.zeros(img.shape, dtype=np.uint16)

            for i in range(img.shape[0]):
                for j in range(len(blobs[:, 0])):
                    if blobs[j, 0] == i:
                        r, c, rad = int(blobs[j, 1]), int(blobs[j, 2]), int(np.round(np.sqrt(blobs[j, 3] * 3)))
                        rr, cc = circle_perimeter(r, c, rad)
                        try:
                            cellsImg[i, rr, cc] = 65535
                        except IndexError:
                            cellsImg[i, r, c] = 65535
            if outPrefix == None:
                fileName = 'cells_' + _file
            else:
                fileName = 'cells_' + outPrefix + '_' + _file
            filePath = os.path.join(samplesPath, fileName)
            with TiffWriter(filePath) as tif:
                tif.save(cellsImg)
            print('{} cells were identified in {}.\nThe locations of these were saved in {}'.format(len(blobs[:, 0]),
                                                                                                    _file, fileName))
    df = pd.DataFrame()
    df['FileName'] = files
    df['CellCounts'] = cellCounts
    outputPath = os.path.join(samplesPath, logName)
    df.to_csv(outputPath, sep='\t')
    return df
