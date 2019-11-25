import os
import numpy as np
from tifffile import TiffWriter
from skimage.morphology import opening
from skimage.morphology import cube
from skimage.feature import blob_dog
from skimage.feature import blob_log
from skimage.draw import circle_perimeter
from skimage.filters import median


def nuclearDetectionList(directoryPath, name, selem = cube(6), min_sigma = 0.25, max_sigma=5, threshold=0.04):
    med = {'function':medianFilter3D, 'suffix':'MF', 'setting':0, 'writeOut':None}
    remBack = {'function':removeBackground, 'suffix':'sansBk', 'setting':1, 'writeOut':None, 'args':[selem]}
    blobs = {'function':getBlobs_LoG, 'suffix':'blob_log', 'setting':2, 
    'kwargs':{'min_sigma':min_sigma, 'max_sigma':max_sigma, 'threshold':threshold}, 
    'writeOut': {'function':writeBlobsImg, 'kwargs':{'outDir':directoryPath, 'name':name}, 'setting':2}}
    funcList = [med, remBack, blobs]
    return funcList

def getBlobs_LoG(img, min_sigma, max_sigma, threshold, **kwargs):
    if kwargs:
        min_sigma = kwargs.get('min_sigma')
        max_sigma = kwargs.get('max_sigma')
        threshold = kwargs.get('threshold')
        outDir = kwargs.get('outDir')
    blobs = blob_log(img, min_sigma = min_sigma, max_sigma=max_sigma, threshold=threshold)
    print('{} blobs detected'.format(len(blobs)))
    return blobs

def medianFilter3D(img):
    im = img
    for i in range(im.shape[0]):
        im[i, :, :] = median(im[i, :, :])
    return im

def removeBackground(data, selem=cube(6), *args):
    # to use with chain --> 1 arg (img sep)
    if args:
        selem = arg[0]
    img = data
    background = opening(img, selem=selem)
    sansBackground = img - background
    return sansBackground 

def writeBlobsImg(img, blobs, outDir, name, _file = None, dtype='auto', **kwargs):
    if kwargs:
        print(kwargs)
        #outDir = kwargs.get('outDir')
        #name = kwargs.get('name')
        if kwargs.get('writeInfo') != None:
            _file = kwargs.get('writeInfo')
            print('successfully obtained filename via chain')
            print(_file)
        #if kwargs.get('writeInfo') == None:
            #_file = kwargs.get('_file')
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
        if name == None:
            fileName = 'cells_'+_file
        if name != None:
            fileName = name +'_cells_'+_file
        print(fileName)
        filePath = os.path.join(outDir, fileName)
        with TiffWriter(filePath) as tif:
            tif.save(cellsImg)
        print('{} blobs were identified in {}.\nThe locations of these were saved in {}'.format(len(blobs[:, 0]), _file, fileName))
        return cellsImg

def imgFileCount(blobs, _file, **kwargs):
    count = len(blobs)
    return [count, _file]

def writeSubsTiff(filedir, filename, data, subsection):
    x_st = subsection.pop('x_start')
    x_en = subsection.pop('x_end')
    y_st = subsection.pop('y_start')
    y_en = subsection.pop('x_end')
    img = data[:, y_st:y_en, x_st:x_en]
    filepath = os.path.join(filedir, filename)
    with TiffWriter(filepath, bigtiff = True) as tiff:
        for i in range(img.shape[0]):
            tiff.save(img[i, :, :])
    return img

#def lableCells(img, blobs, subsection = None):     
  #  cellsImg = np.zeros(img.shape, dtype=np.uint16)
    
   # if subsection != None:
     #   x_st = subsection.pop('x_start')
     #   y_st = subsection.pop('y_start')
      #  blobs[:, 1] = blobs[:, 1] - y_st
      #  blobs[:, 2] = blobs[:, 2] - x_st
        
      #  for i in range(img.shape[0]):
         #   for j in range(len(blobs[:, 0])):
         #       if blobs[j, 0] == i:
          #          r, c, rad = int(blobs[j, 1]), int(blobs[j, 2]), int(np.round(np.sqrt(blobs[j, 3]*3)))
          #          rr, cc = circle_perimeter(r, c, rad)
          #          try:
         #               cellsImg[i, rr, cc] = 65535
           #         except IndexError:
               #        cellsImg[i, r, c] = 65535
       # if outPrefix == None:
      #      fileName = 'cells_'+_file
      #  else:
      #      fileName = 'cells_'+outPrefix+'_'+_file
     #   filePath = os.path.join(samplesPath, fileName)
     #   with TiffWriter(filePath) as tif:
      #      tif.save(cellsImg)
      #  print('{} cells were identified in {}.\nThe locations of these were saved in {}'.format(len(blobs[:, 0])))