import os
import numpy as np
import re


def resamplePoints(points, origRes, newRes):
    array = points
    for i in range(len(points)):
        array[i, 0] = int(np.ceil(array[i, 0]*origRes[0]/newRes[0]))
        array[i, 1] = int(np.ceil(array[i, 1]*origRes[1]/newRes[1]))
        array[i, 2] = int(np.ceil(array[i, 2]*origRes[2]/newRes[2]))
    print(array[:5, :])
    return array

## To be useful, this would require the transformation to be inverse.
def transformPoints(transformixBin, pointsTXT, parameterFiles, outDir, prefix, outName = 'points_transformed'):
    outName0 = prefix+"_"+outName
    outPath = os.path.join(outDir, outName0)
    outDirNew = os.path.join(outPath, 'outputpoints.txt')
    if os.path.exists(outPath) != True:
        os.mkdir(outPath)
    sourcePath = pointsTXT
    points = ' -def '+sourcePath   
    out = ' -out '+outPath
    param = ' -tp '+parameterFiles
    cmd = transformixBin+points+out+param
    print('Executing the following:')
    print(cmd)
    res = os.system(cmd)
    if res != 0:
        print('An error occured, check the parameters and associated files')
    return outDirNew

def transformImg(transformixBin, inTIFF, parameterFiles, outDir, prefix, outName = 'img_transformed'):
    outName0 = prefix+"_"+outName
    outPath = os.path.join(outDir, outName0)
    outDirNew = os.path.join(outPath, 'result.mhd') 
    if os.path.exists(outPath) != True:
        os.mkdir(outPath)
    sourcePath = inTIFF
    points = ' -in '+sourcePath   
    out = ' -out '+outPath
    param = ' -tp '+parameterFiles
    cmd = transformixBin+points+out+param
    print('Executing the following:')
    print(cmd)
    res = os.system(cmd)
    if res != 0:
        print('An error occured, check the parameters and associated files')
    return outDirNew

def runElastix(elastixBin, movingImage, fixedImage, parameterFiles, outDir):
    if os.path.exists(outDir) != True:
        os.mkdir(outDir)
    threads = ' -threads 16'
    moving = ' -m '+movingImage
    fixed = ' -f '+fixedImage
    if type(parameterFiles) == str:
        param = ' -p '+parameterFiles
    if type(parameterFiles) == list:
        params = []
        for p in parameterFiles:
            par = ' -p '+p
            params.append(par)
        param = ''.join(params)
    out = ' -out '+outDir
    cmd = elastixBin+threads+moving+fixed+param+out
    print('Attempting to execute:')
    print(cmd)
    result = os.system(cmd)
    return outDir

def getAtlasPoints(img, points):
    cellImg = np.zeros(img.shape, dtype=np.uint32)
    indices = np.zeros(points.shape, dtype=int)
    for i in range(len(points)):
        indices[i, :] = [int(np.round(points[i, 2])), int(np.round(points[i, 1])), int(np.round(points[i, 0]))]
    for i in range(len(indices)):
        
        if indices[i, 0] >= cellImg.shape[0]:
            indices[i, 0] = cellImg.shape[0] - 1
        if indices[i, 1] >= cellImg.shape[1]:
            indices[i, 1] = cellImg.shape[1] - 1
        if indices[i, 2] >= cellImg.shape[2]:
            indices[i, 2] = cellImg.shape[2] - 1
        
        cellImg[indices[i, 0], indices[i, 1], indices[i, 2]] += 1
            
    return cellImg

def correctToXYZ(cellArray):
    newArray = np.array([cellArray[:, 2], cellArray[:, 1], cellArray[:, 0]])
    newArray = newArray.T
    return newArray

def writePointsForTransformix(points, prefix, saveDir, res):
    saveName = prefix+'_'+str(res)+'um_coords_xyz_elastixFormat.txt'
    savePath = os.path.join(saveDir, saveName)
    header = 'point\n{}'.format(len(points))
    np.savetxt(savePath, points, header = header, comments='')
    return savePath

def parseTransformixOutput(filePath):
    transPoints = open(filePath, 'r')
    lines = transPoints.readlines()
    array = np.empty([len(lines), 3], dtype = int)
    for i in range(len(lines)):
        string = lines[i]
        match = re.findall(r'OutputPoint\s=\s\[\s\d+\.\d+\s\d+\.\d+\s\d+\.\d+', string)
        if len(match) != 0:
            xyz = re.findall("\d+\.\d+", match[0])
            row = []
            for cood in xyz:
                asfloat = float(cood)
                asint = int(np.round(asfloat))
                row.append(asint)
            array[i, :] = row
    return array

# Density analysis
def getDensityImg(saveDir, saveName, points, shape, order):
    img = getVoxelCounts(points, shape, order)
    img = np.round(img*1000000/(16*16*16))
    img = img.astype(np.uint16)
    savePath = os.path.join(saveDir, saveName)
    with TiffWriter(savePath) as tiff:
        tiff.save(img)
    return img

def getVoxelCounts(points, shape, order):
    '''
    Expects int style points, not float
    '''
    z = order[0]    
    y = order[1]
    x = order[2]
    img = np.zeros(shape, dtype = np.uint16)
    for point in points:
        try:
            ind = (point[z], point[y], point[x])
            img[ind] += 1
        except IndexError:
            print('point {}, {}, {} out of range'.format(point[z], point[y], point[x]))

    return img

### Basic transformation stuff
def affineTransform(parameters, COR, points):
    '''
    Assumes input points in 2d array w/ r: x, y, z
    Assumes input points are np.ndarray
    '''
    mat = np.array([
        [parameters[0], parameters[1], parameters[2]], 
        [parameters[3], parameters[4], parameters[5]], 
        [parameters[6], parameters[7], parameters[8]]])
    trans = np.array([parameters[9], parameters[10], parameters[11]])
    output = []
    for row in points:
        row = row.astype(float)
        row = row - COR
        output_row = mat @ row
        output_row = output_row + trans + COR
        output.append([output_row])
    output = np.concatenate(output)

    return output