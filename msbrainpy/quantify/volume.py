import os
import h5py as h5
import numpy as np
import pandas as pd
from skimage.morphology import opening
from skimage.morphology import cube
from skimage.feature import blob_log
from skimage.feature import peak_local_max
from skimage.filters import median
from skimage.transform import rescale
from skimage.draw import circle_perimeter
from skimage.util import img_as_uint
from shutil import rmtree

def cellsOnly(data, subsection, prefix, outDir, z_size = 50, overlap = 10,  
              min_sigma = 0.25, max_sigma=10, threshold=0.04, **kwargs):
    '''
    FUNCTION: function to sequentially process chunks to retrive cell coordinates
    ARGUMENTS:
        data = hdf5 data object
        z_size = z size of chunk (int)
        overlap = overlap of chunk (int)
        min_sigma = minimum sigma for LoG blob detection (smaller for detecting smaller blobs)
        max_sigma = maximum sigma for LoG blob detection (larger for detecting larger blobs)
        threshold = intensity threashold for blob detection (lower for fainter cells)
    RETURNS: array with cell coordinates [[z, y, x, sigma], ...]
    '''
    if kwargs:
        data = kwargs.get('data')
        subsection = kwargs.get('subsection')
        prefix = kwargs.get('prefix')
        outDir = kwargs.get('outDir')
        if kwargs.get('z_size')!= None:
            z_size = kwargs.get('z_size')
        if kwargs.get('overlap')!= None:
            overlap = kwargs.get('overlap')
        if kwargs.get('min_sigma')!= None:
            min_sigma = kwargs.get('min_sigma')
        if kwargs.get('max_sigma')!= None:
            max_sigma = kwargs.get('max_sigma')
        if kwargs.get('threshold')!= None:
            threshold = kwargs.get('threshold')
        
    sansBackground = chunkGenerator(data, subsection = subsection, z_size = z_size, overlap = overlap, 
                                    dowhat = removeBackground)
    cellArray = cellsFromBlobChunks(sansBackground, subsection, prefix, outDir, overlap = overlap, 
                                   min_sigma = min_sigma, max_sigma=max_sigma, threshold=threshold)
    
    return cellArray

#def cellsFromh5Blobs(data, subsection, z_size, overlap, **kwargs):
    #chunks = chunkGenerator(data, subsection = subsection, z_size = z_size, overlap = overlap, dowhat = None)
    #cellArray = cellsFromBlobChunks(sansBackground, subsection = subsection, overlap = overlap)
    #return cellArray
             
def cellsFromBlobChunks(sansBackground, subsection, prefix, outDir, overlap, min_sigma = 0.25, max_sigma=10, threshold=0.04):
    cellArray = []
    y_st= subsection['y_start']
    x_st= subsection['x_start']
    stackID = subsection['stackID']
    zPos = 0 
    edge = int(np.round(overlap/2))
    tempName = prefix+"_substack"+str(stackID)+"_cells"
    tempPath = os.path.join(outDir, tempName)
    os.mkdir(tempPath)
    for chunk in sansBackground:
        zPos0 = zPos
        img = chunk
        size = img.shape
        cells = blob_log(img, min_sigma = min_sigma, max_sigma=max_sigma, threshold=threshold)
        x_lim = size[2] - edge
        y_lim = size[1] - edge
        z_lim = size[0] - edge
        print("{} cells were initialy found at z {}:{} in substack {}".format(len(cells), zPos0, zPos0+size[0], stackID))
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
            saveName = "cells_z{}-{}_substack{}".format(zPos0, zPos0+size[0], stackID)
            savePath = os.path.join(tempPath, saveName)
            np.save(savePath, cells, allow_pickle=False)
            
        zPos += size[0]-overlap
        print("{} cells were found at z {}:{} in substack {}".format(len(cells), zPos0, zPos0+size[0], stackID))
        
        
    print('A total of {} cells were found in substack {}'.format(len(cellArray), stackID))    
    cellArray = np.concatenate(cellArray)
    
    cellsFile = tempName+'.npy'
    cellsPath = os.path.join(outDir, cellsFile)
    np.save(cellsPath, cellArray, allow_pickle=False)
    
    return cellArray 
    

def removeBackground(data, selem=cube(6)):
    img = data
    for i in range(img.shape[0]):
        img[i, :, :] = median(img[i, :, :])
    background = opening(img)
    sansBackground = img - background
    return sansBackground

