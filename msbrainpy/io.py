import nrrd
import os
import h5py as h5
from tifffile import imread
import numpy as np
from tifffile import TiffWriter


def imread_i(directory, i):
    files = os.listdir(directory)
    file_ = files[i]
    filepath = os.path.join(directory, file_)
    img = imread(filepath)
    return img

def generateFromDiretory(directory, prefix = None, **kwargs):
    if kwargs:
        prefix = kwargs.get('prefix')

    allfiles = os.listdir(directory)
    files = []
    if prefix == None:
        for _file in allfiles:
            if _file.endswith('.tif'):
                files.append(_file)
    if prefix != None:
        for _file in allfiles:
            if _file.endswith('.tif') and _file.startswith(prefix):
                files.append(_file)            
    for _file in files:
        filePath = os.path.join(directory, _file)
        img = imread(filePath)
        yield img, _file

def readAnnotation(annName, crop = None): #, ori = None):
    annotation, header = nrrd.read(annName)
    print("nrrd header for the annotation file:")
    print(header)
    annotation = annotation.astype('int32')
    print('Annotation file has been converted to int32')
    #if ori != None:
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

def openH5data(filedir, filename, dataset = True):
    filepath = os.path.join(filedir, filename)
    h5file = h5.File(filepath, 'r')
    print('A hdf5 file was found at {}. \nThis contains the following datasets:'.format(filepath))
    for key in h5file.keys():
        print(key)
    keys = list(h5file.keys())
    data = h5file[keys[0]]
    print('dataset: {} \nshape: {}, {}, {}'.format(keys[0], data.shape[0], data.shape[1], data.shape[2]))

    return data

def saveRandom(data, n, size, prefix, saveDir):
    shape = data.shape
    sample_z = np.random.choice(shape[0]-size[0], n)
    sample_y = np.random.choice(shape[1]-size[1], n)
    sample_x = np.random.choice(shape[2]-size[2], n)
    for i in range(n):
        z_end = sample_z[i]+size[0]
        y_end = sample_y[i]+size[1]
        x_end = sample_x[i]+size[2]
        sample = data[sample_z[i]:z_end, sample_y[i]:y_end, sample_x[i]:x_end]
        fileCoord = '_z{}-{}_y{}-{}_x{}-{}.tif'.format(
            sample_z[i], z_end, sample_y[i], y_end, sample_x[i], x_end)
        fileName = prefix+fileCoord
        filePath = os.path.join(saveDir, fileName)
        with TiffWriter(filePath) as tif:
            tif.save(sample)

def writeCroppedTiffStack(directory, outDir,  y_ind, x_ind):
    fileList = os.listdir(directory)
    for _file in fileList:
        imPath = os.path.join(directory, _file)
        image = imread(imPath)
        image = image[:, y_ind[0]:y_ind[1], x_ind[0]:x_ind[1]]
        imPath = os.path.join(outDir, _file)
        print('for {}, the new image shape is {}'.format(_file, image.shape))
        with TiffWriter(imPath, bigtiff = True) as tif:
            for i in range(image.shape[0]):
                tif.save(image[i, :, :])
        print('the image {} was saved at:\n {}'.format(_file, imPath))

def writeHDF5_seq(filename, filedir, datasetName, imgDir):
    '''
    FUNCTION: write a HDF5 file from a series of 3D images (stacked in Z). Writes hdf5 file and a table describing where 
        each tiff lives in the hdf5
    ARGUMENTS:
        filename = name for hdf5 output (str)
        filedir = output directory
        datasetName = name to which the  data will be assigned within the hdf5 file (str)
        imgDir = directory in which tiffs are sequentially listed (str)
    RETURNS: hdf5 dataset object
    '''
    filepath = os.path.join(filedir, filename)
    files = []
    for _file in os.listdir(imgDir):
        if _file.endswith('.tif'):
            files.append(_file)
            
    im1path = os.path.join(imgDir, files[0])
    im1 = imread(im1path)
    shape1 = im1.shape
    del im1
    
    with h5py.File(filepath, 'w') as h5:
        dset = h5.create_dataset(datasetName, shape1, dtype=np.uint16, maxshape=(None, shape1[1], shape1[2]), chunks=True)

        started = False
        startAt = 0
        startpos = []
        endpos = []
        for _file in files:
            imgPath = os.path.join(imgDir, _file)
            img = imread(imgPath)
    
            if started == False:
                dset[startAt:startAt+img.shape[0], :, :] = img
                started = True
            else:
                dset.resize(dset.shape[0]+img.shape[0], axis = 0)
                dset[startAt:startAt+img.shape[0], :, :] = img
    
            print('{} was saved at Z-indicies {}-{} in {}'.format(_file, startAt, dset.shape[0], filepath))
            startpos.append(startAt)
            endpos.append(dset.shape[0])
            startAt += img.shape[0]
        
    positions = pd.DataFrame()
    positions['Files'] = files
    positions['startInH5'] = startpos
    positions['endInH5'] = endpos
    positions.to_csv(os.path.join(imgDir, '{}zTable.txt'.format(filename)), sep='\t')
    
    return dset