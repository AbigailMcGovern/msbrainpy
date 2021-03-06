import os
import numpy as np
from shutil import rmtree
from multiprocessing import Pool
from skimage.transform import resize
from skimage.util import img_as_uint


# from shutil import rmtree --> eventual removal of accumulated files

# ------------------------------------------ Basic image processing functions ------------------------------------------

def rescaleImg(data, filedir, filename, originalRes, finalRes, verbose=True):
    """
    NOTE: x-y rescaled image must be able to fit in memory
    """
    z_scl = originalRes[0] / finalRes[0]
    y_scl = originalRes[1] / finalRes[1]
    x_scl = originalRes[2] / finalRes[2]
    if verbose:
        print('the image will be scaled by ({}, {}, {})'.format(z_scl, y_scl, x_scl))
        print('the image will be scaled along x & y then along z')
    size0 = (int(np.ceil(y_scl * data.shape[1])), int(np.ceil(x_scl * data.shape[2])))
    shape0 = [data.shape[0], int(np.ceil(y_scl * data.shape[1])), int(np.ceil(x_scl * data.shape[2]))]
    size1 = (int(np.ceil(z_scl * data.shape[0])), int(np.ceil(y_scl * data.shape[1])))
    shape1 = [int(np.ceil(z_scl * data.shape[0])), int(np.ceil(y_scl * data.shape[1])),
              int(np.ceil(x_scl * data.shape[2]))]
    if verbose:
        print(shape0)
    Zblocks = chunkGenerator(data, subsection=None, z_size=50, dowhat=None, zPos=0, **kwargs)
    count = 0
    resampled1 = np.zeros(shape0, dtype=np.float64)
    for block in Zblocks:
        smlBlock = np.zeros([block.shape[0], shape0[1], shape0[2]], dtype=np.float64)
        count_ = count
        for i in range(block.shape[0]):
            plane = block[i, :, :]
            planeYX = resize(plane, size0, order=3, clip=True, anti_aliasing=True, mode='reflect')
            smlBlock[i, :, :] = planeYX
        count += block.shape[0]
        resampled1[count_:count, :, :] = smlBlock
        if verbose:
            print('Values along x & y were resampled at z{}-{}'.format(count_, count))
            print('Original = {}, {}, {}'.format(block.shape[0], block.shape[1], block.shape[2]))
            print('New = {}, {}, {}'.format(block.shape[0], resampled1.shape[1], resampled1.shape[2]))
    filepath = os.path.join(filedir, 'XY_resampled.tif')
    with TiffWriter(filepath, bigtiff=True) as tiff:
        for i in range(resampled1.shape[0]):
            tiff.save(img_as_uint(resampled1[i, :, :]))
    if verbose:
        print('The image was saved at {}'.format(filepath))
    out = np.zeros(shape1, dtype=np.float64)
    for i in range(resampled1.shape[2]):
        plane = resampled1[:, :, i]
        planeZY = resize(plane, size1, order=3, clip=True, anti_aliasing=True, mode='reflect')
        out[:, :, i] = planeZY
    out = img_as_uint(out)
    if verbose:
        print('Values along z were resampled')
        print('The final image has a size of ({}, {}, {})'.format(out.shape[0], out.shape[1], out.shape[2]))
    filepath = os.path.join(filedir, filename)
    with TiffWriter(filepath) as tiff:
        tiff.save(out)
    if verbose:
        print('The image was saved at {}'.format(filepath))
    return out


# replace the second part of the above

def resize_z(resampled, filedir, filename, shape, originalRes, finalRes):
    z_scl = originalRes[0] / finalRes[0]
    y_scl = originalRes[1] / finalRes[1]
    x_scl = originalRes[2] / finalRes[2]
    shape1 = [int(np.ceil(z_scl * shape[0])), int(np.ceil(y_scl * shape[1])), int(np.ceil(x_scl * shape[2]))]
    size1 = (int(np.ceil(z_scl * shape[0])), int(np.ceil(y_scl * shape[1])))
    out = np.zeros(shape1, dtype=np.float64)
    for i in range(resampled.shape[2]):
        plane = resampled[:, :, i]
        planeZY = resize(plane, size1, order=3, clip=True, anti_aliasing=True, mode='reflect')
        out[:, :, i] = planeZY
    out = img_as_uint(out)
    print('Values along z were resampled')
    print('The final image has a size of ({}, {}, {})'.format(out.shape[0], out.shape[1], out.shape[2]))
    filepath = os.path.join(filedir, filename)
    with TiffWriter(filepath) as tiff:
        tiff.save(out)
    print('The image was saved at {}'.format(filepath))


# ----------------------------------------------- Subsection-related  --------------------------------------------------

def getSubsections(shape, y=5, x=5, overlap=10, verbose=False):
    yRange = shape[1]
    xRange = shape[2]
    x_st, x_en = blocks(xRange, x, overlap=overlap)
    y_st, y_en = blocks(yRange, y, overlap=overlap)
    if verbose:
        print('yRange = {}'.format(yRange))
        print('xRange = {}'.format(xRange))
        print('finding x blocks:')
        print('len(x_en) = {}'.format(len(x_en)))
        print('finding y blocks:')
        print('len(y_en) = {}'.format(len(y_en)))
    ends = []
    starts = []
    for i in range(y):
        for j in range(x):
            starts.append([y_st[i], x_st[j]])
    starts = np.array(starts)
    if verbose:
        print('Start index: \nstarting with {}\nending with {}'.format(starts[0], starts[-1]))
    for i in range(y):
        for j in range(x):
            ends.append([y_en[i], x_en[j]])
    ends = np.array(ends)
    if verbose:
        print('End index: \nstarting with {}\nending with {}'.format(ends[0], ends[-1]))
    dictList = []
    for i in range(x * y):
        dictList.append({
            'stackID': i,
            'x_start': starts[i, 1],
            'x_end': ends[i, 1],
            'y_start': starts[i, 0],
            'y_end': ends[i, 0],
            'overlap': overlap
        })
    return dictList


def blocks(axisRange, num, overlap=10):
    cl = int(np.round(axisRange / num))
    print('blocks size = {}'.format(cl))
    count = 0
    st = []
    en = []
    while count < axisRange:
        print('count = {}'.format(count))
        st.append(count)
        if len(en) == num - 1:
            en.append(axisRange)
        else:
            en.append(count + cl)
        count += cl - overlap
    return st, en


# ----------------------------------------------- Data subdivision  ----------------------------------------------------

def chunkGenerator(data, subsection=None, z_size=50, dowhat=None, zPos=0, zOverlap=0, writeInfoExists=False):
    """
    FUNCTION: generator function to process images larger than RAM
    ARGUMENTS:
        data = hdf5 data object
        subsection = appropriate subsection dictionary (dict)
        z_size = z size of chunk (int)
        dowhat = function to execute (func)
        zPos = where to start on z in stack (int)
        zOverlap = if no subsection data is provided and the img chunks should be generated with an overlap, the
            number of voxels overlap should be provided here.
    YEILDS: each iteration yields an ndarray with newly appended data (as per overlap)
    """
    if subsection is not None:
        y_st = subsection['y_start']
        y_en = subsection['y_end']
        x_st = subsection['x_start']
        x_en = subsection['x_end']
        overlap = subsection['overlap']
    else:
        y_st = 0
        y_en = data.shape[1]
        x_st = 0
        x_en = data.shape[2]
        overlap = zOverlap
    if dowhat is None:
        def getResult(data_): return readChunk(data_, zPos, z_size, y_st, y_en, x_st, x_en)
    else:
        def getResult(data_): return dowhat(readChunk(data_, zPos, z_size, y_st, y_en, x_st, x_en))
    while zPos < data.shape[0]:
        out = getResult(data)
        if writeInfoExists:
            chunk = subsection.copy()
            chunk['z_start'] = zPos
            chunk['z_end'] = zPos + z_size
            yield out, chunk
            zPos += z_size - overlap
        else:
            yield out


def readChunk(data, start, z_size, y_st, y_en, x_st, x_en):
    end = start + z_size
    try:
        chunk = data[start:end, y_st:y_en, x_st:x_en]
    except IndexError:
        chunk = data[start:, y_st:y_en, x_st:x_en]
    return chunk


# ----------------------------------------------- Basic helper functions -----------------------------------------------

def extraListNesting(thing):
    thing = [thing]
    return thing


def getProp(props, IDs, key):
    out = []
    for ID in IDs:
        found = False
        for prop in props:
            if prop.label == ID:
                out.append(prop[key])
                found = True
        if found == False:
            out.append('NA')
    return out
