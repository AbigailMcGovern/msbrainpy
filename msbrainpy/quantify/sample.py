import os
import pandas as pd
from skimage.morphology import cube
from msbrainpy.base import Chain
from msbrainpy.io import generateFromDiretory
from msbrainpy.quantify.processing import medianFilter3D, removeBackground, getBlobs_LoG, nuclearDetectionList, imgFileCount

# use this to validate The Chain and the processDirectory function. 
# once these work, reconfigure subprocessing to make more dynamic. Consider writing subsection as a Class...
def sampleNuclearDetection(directoryPath, prefix, name,  selem = cube(6), min_sigma = 0.25, max_sigma=5, threshold=0.04):
    funcList = nuclearDetectionList(directoryPath, name, selem = selem, min_sigma = min_sigma, max_sigma=max_sigma, threshold=threshold)
    df = processDirectory(directoryPath, funcList, prefix)
    return df

def processDirectory(directoryPath, functionDictList, prefix, recordMethod = imgFileCount, recordName = 'count', extractMethod = None, writeInfo_recordkw = '_file'):
    imgChain = Chain(directoryPath, inMethod = generateFromDiretory, in_kwargs = {'prefix':prefix})
    imgChain.addScribe(recordMethod = recordMethod, extractMethod = extractMethod, writeInfo_recordkw =writeInfo_recordkw)
    for functionDict in functionDictList:
        if type(functionDict) == dict:
            print(functionDict['function'])
            print(functionDict.get('function')) 
            print(functionDict.get('args'))
            imgChain.addChainMethod(**functionDict)
        else:
            imgChain.addChainMethod(functionDict)
    chainResult = imgChain.execute()
    allfiles = os.listdir(directoryPath)
    files = []
    for _file in allfiles:
        if _file.endswith('.tif') and _file.startswith(prefix):
            files.append(_file)   
    df = pd.DataFrame()
    print(chainResult)
    RN = []
    FN = []
    for row in chainResult:
        RN.append(row[0])
        FN.append(row[1])
    df[recordName] = RN
    df['filename'] = FN
    newName = prefix+recordName+'.txt'
    newPath = os.path.join(directoryPath, newName)
    df.to_csv(newPath, sep='\t')
    return df

# consider this deprecated once chain-based ND works
def sampleCellCounts(samplesPath, selem, min_sigma = 0.25, max_sigma=5, threshold=0.1, logName = 'CellCounts_LoG.txt', outPrefix = None):
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
            blobs = blob_log(sansBackground, min_sigma = min_sigma, max_sigma=max_sigma, threshold=threshold)
            cellCounts.append(len(blobs[:, 0]))
            cellsImg = np.zeros(img.shape, dtype=np.uint16)
        
            for i in range(img.shape[0]):
                for j in range(len(blobs[:, 0])):
                    if blobs[j, 0] == i:
                        r, c, rad = int(blobs[j, 1]), int(blobs[j, 2]), int(np.round(np.sqrt(blobs[j, 3]*3)))
                        rr, cc = circle_perimeter(r, c, rad)
                        try:
                            cellsImg[i, rr, cc] = 65535
                        except IndexError:
                            cellsImg[i, r, c] = 65535
            if outPrefix == None:
                fileName = 'cells_'+_file
            else:
                fileName = 'cells_'+outPrefix+'_'+_file
            filePath = os.path.join(samplesPath, fileName)
            with TiffWriter(filePath) as tif:
                tif.save(cellsImg)
            print('{} cells were identified in {}.\nThe locations of these were saved in {}'.format(len(blobs[:, 0]), _file, fileName))
            
    df = pd.DataFrame()
    df['FileName'] = files
    df['CellCounts'] = cellCounts
    outputPath = os.path.join(samplesPath, logName)
    df.to_csv(outputPath, sep='\t')


