import os
import h5py as h5
import numpy as np
import pandas as pd
from tifffile import TiffWriter
from tifffile import imread
from shutil import rmtree
from multiprocessing import Pool
from skimage.transform import resize
from skimage.util import img_as_uint

def parallelProcess(function, processes = 16, x = 5, y = 5, 
                    concat = np.concatenate, exclude = None, **kwargs):
    overlap = kwargs.get('overlap')
    data = kwargs.get('data')
    
    outDir = kwargs.get('outDir')
    prefix = kwargs.get('prefix')
    saveName = prefix+'cellCoords.npy'
    savePath = os.path.join(outDir, saveName)
    
    subsList = getSubsections(data.shape, x=x, y=y, overlap=overlap)
    
    args = []
    if exclude != None:
        for i in range(x*y):
            subs = subsList[i]
            if subs['stackID'] not in exclude:
                param = kwargs.copy()
                param['subsection'] = subs
                args.append(param)
    if exclude == None:
        for i in range(x*y):
            subs = subsList[i]
            param = kwargs.copy()
            param['subsection'] = subs
            args.append(param)
    done = False
    pos = 0
    output = []
    while done != True:
        if len(args)-pos > processes:
            do = args[pos:pos+processes]
        if len(args)-pos < processes:
            do = args[pos:] 
        with Pool(processes = processes) as pool:
            result = pool.starmap(function, do)
        
    output.append(concat(result))
    output = concat(output)
    np.save(savePath, output, allow_pickle=False)
    return output
    
#

def rescaleImg(data, filedir, filename, originalRes, finalRes): 
    '''
    NOTE:output image must be able to fit in memory 
    '''
    z_scl = originalRes[0]/finalRes[0]
    y_scl = originalRes[1]/finalRes[1]
    x_scl = originalRes[2]/finalRes[2]
    print('the image will be scaled by ({}, {}, {})'.format(z_scl, y_scl, x_scl))
    print('the image will be scaled along x & y then along z')
    size0 = (int(np.ceil(y_scl*data.shape[1])), int(np.ceil(x_scl*data.shape[2])))
    shape0 = [data.shape[0], int(np.ceil(y_scl*data.shape[1])), int(np.ceil(x_scl*data.shape[2]))]
    size1 = (int(np.ceil(z_scl*data.shape[0])), int(np.ceil(y_scl*data.shape[1])))
    shape1 = [int(np.ceil(z_scl*data.shape[0])), int(np.ceil(y_scl*data.shape[1])), int(np.ceil(x_scl*data.shape[2]))]
    print(shape0)
    
    Zblocks = chunkGenerator(data, subsection = None, z_size = 50, overlap = 0, dowhat = None, zPos = 0)
    
    count = 0
    resampled1 = np.zeros(shape0, dtype=np.float64)
    for block in Zblocks:
        smlBlock = np.zeros([block.shape[0], shape0[1], shape0[2]], dtype=np.float64)
        count_ = count
        for i in range(block.shape[0]):
            plane = block[i, :, :]
            planeYX = resize(plane, size0, order = 3, clip=True, anti_aliasing=True, mode='reflect')
            smlBlock[i, :, :] = planeYX
        count += block.shape[0]
        resampled1[count_:count, :, :] = smlBlock
        
        print('Values along x & y were resampled at z{}-{}'.format(count_, count))
        print('Original = {}, {}, {}'.format(block.shape[0], block.shape[1], block.shape[2]))
        print('New = {}, {}, {}'.format(block.shape[0], resampled1.shape[1], resampled1.shape[2]))
    
    filepath = os.path.join(filedir, 'XY_resampled.tif')
    with TiffWriter(filepath, bigtiff=True) as tiff:
        for i in range(resampled1.shape[0]):
            tiff.save(img_as_uint(resampled1[i, :, :]))
    print('The image was saved at {}'.format(filepath))
    
    out = np.zeros(shape1, dtype=np.float64)
    for i in range(resampled1.shape[2]):
        plane = resampled1[:, :, i]
        planeZY = resize(plane, size1, order = 3, clip=True, anti_aliasing=True, mode='reflect')
        out[:, :, i] = planeZY
    
    out = img_as_uint(out)
    print('Values along z were resampled')
    print('The final image has a size of ({}, {}, {})'.format(out.shape[0], out.shape[1], out.shape[2]))    
    
    filepath = os.path.join(filedir, filename)
    with TiffWriter(filepath) as tiff:
        tiff.save(out)
        
    print('The image was saved at {}'.format(filepath))
    
    return out

def resize_z(resampled, filedir, filename, shape, originalRes, finalRes):
    z_scl = originalRes[0]/finalRes[0]
    y_scl = originalRes[1]/finalRes[1]
    x_scl = originalRes[2]/finalRes[2]
    shape1 = [int(np.ceil(z_scl*shape[0])), int(np.ceil(y_scl*shape[1])), int(np.ceil(x_scl*shape[2]))]
    size1 = (int(np.ceil(z_scl*shape[0])), int(np.ceil(y_scl*shape[1])))
    out = np.zeros(shape1, dtype=np.float64)
    for i in range(resampled.shape[2]):
        plane = resampled[:, :, i]
        planeZY = resize(plane, size1, order = 3, clip=True, anti_aliasing=True, mode='reflect')
        out[:, :, i] = planeZY
    
    out = img_as_uint(out)
    print('Values along z were resampled')
    print('The final image has a size of ({}, {}, {})'.format(out.shape[0], out.shape[1], out.shape[2]))    
    
    filepath = os.path.join(filedir, filename)
    with TiffWriter(filepath) as tiff:
        tiff.save(out)
        
    print('The image was saved at {}'.format(filepath))


# GET SUBSECTIONS HAS BEEN VALIDATED
def getSubsections(shape, y = 5, x = 5, overlap = 10):
    yRange = shape[1]
    print('yRange = {}'.format(yRange))
    xRange = shape[2]
    print('xRange = {}'.format(xRange))
    print('finding x blocks:')
    x_st, x_en = blocks(xRange, x, overlap=overlap)
    print('len(x_en) = {}'.format(len(x_en)))
    
    
    print('finding y blocks:')
    y_st, y_en = blocks(yRange, y, overlap=overlap)
    print('len(y_en) = {}'.format(len(y_en)))
    
    ends = []
    starts = []
    
    for i in range(y):
        for j in range(x):
            starts.append([y_st[i], x_st[j]])
    
    starts = np.array(starts)
    print('Start index: \nstarting with {}\nending with {}'.format(starts[0], starts[-1]))
    
    for i in range(y):
        for j in range(x):
            ends.append([y_en[i], x_en[j]])
    ends = np.array(ends)
    print('End index: \nstarting with {}\nending with {}'.format(ends[0], ends[-1]))
    
    dictList = []
    for i in range(x*y):
        dictList.append({
            'stackID' : i, 
            'x_start' : starts[i, 1], 
            'x_end' : ends[i, 1], 
            'y_start' : starts[i, 0], 
            'y_end' : ends[i, 0]
        })
    return dictList

def chunkGenerator(data, subsection = None, z_size = 50, overlap = 10, dowhat = None, zPos = 0):
    '''
    FUNCTION: generator function to process images larger than RAM
    ARGUMENTS:
        data = hdf5 data object
        subsection = appropriate subsection dictionary (dict)
        z_size = z size of chunk (int)
        overlap = overlap of chunk (int)
        dowhat = function to execute (func)
        zPos = where to start on z in stack (int)
    YEILDS: each iteration yields an ndarray with newly appended data (as per overlap)
    '''
    
    if subsection != None:       
        y_st= subsection['y_start']
        y_en = subsection['y_end']
        x_st= subsection['x_start']
        x_en = subsection['x_end']
    if subsection == None:
        y_st= 0
        y_en = data.shape[1]
        x_st= 0
        x_en = data.shape[2]
        
    if dowhat == None:
        while zPos < data.shape[0]:
            img = readChunk(data, zPos, z_size, y_st, y_en, x_st, x_en)
            yield img
            zPos += z_size-overlap
    if dowhat != None:        
        while zPos < data.shape[0]:
            img = readChunk(data, zPos, z_size, y_st, y_en, x_st, x_en)
            img = dowhat(img)
            zPos += z_size-overlap
            yield img
            
# valid function                        
def readChunk(data, start, z_size, y_st, y_en, x_st, x_en):
    end = start+z_size
    try:
        chunk = data[start:end, y_st:y_en, x_st:x_en]
    except IndexError:
        chunk = data[start:, y_st:y_en, x_st:x_en]
    return chunk

def blocks(axisRange, num, overlap=10):
    cl = int(np.round(axisRange/num))
    print('blocks size = {}'.format(cl))
    count = 0 
    st = []
    en = []
    
    while count < axisRange:
        print('count = {}'.format(count))
        st.append(count)
        if len(en) == num-1:
            en.append(axisRange)
        else:
            en.append(count+cl)
        count += cl-overlap
        
    return st, en

## The chain (Rumors reference intended) 

class Chain:
    def __init__(self, preiterable, inMethod = None, in_kwargs = None, writeInfo = True):
        self.preiterable = preiterable # probably a substack block generator or director of files
        self.inMethod = inMethod # function for reading input into chain
        self.in_kwargs = in_kwargs
        self.chain = [] # list of ChainMethods 
        self.writeInfo = writeInfo

    # function for generating/manipulating output over the iteration
    def addScribe(self, recordMethod = None, extractMethod = None, writeInfo_recordkw = None):
        self.scribe = Scribe(recordMethod = recordMethod, extractMethod = extractMethod, writeInfo_recordkw = writeInfo_recordkw)

    def addChainMethod(self, **functionDict):
        newLink = ChainMethod(**functionDict)
        self.chain.append(newLink)

    def updateChainMethod(self, i, paramName, update, keyName = 'function'): 
        self.chain[i].funct[keyName][paramName] = update

    def execute(self):
        if self.inMethod != None:
            if self.in_kwargs == None:
                chainObjs = self.inMethod(self.preiterable)
            if self.in_kwargs != None:
                chainObjs = self.inMethod(self.preiterable, **self.in_kwargs)
            
        else:
            chainObjs = self.preiterable

        if self.writeInfo == True: 
            for chainObj, writeInfo in chainObjs:
                print(writeInfo)
                for chainMethod in self.chain:
                    if chainMethod.funct['writeOut'] != None:
                        chainObj = chainMethod.withWriteOut(chainObj, writeInfo = writeInfo)
                    if chainMethod.funct['writeOut'] == None:
                        chainObj = chainMethod.run(chainObj)
                    self.scribe.addRecord(chainObj, writeInfo = writeInfo)
        if self.writeInfo == False:
            for chainObj in chainObjs:
                for chainMethod in self.chain:
                    if chainMethod.funct['writeOut'] != None:
                        chainMethod.withWriteOut(chainObj)
                    if chainMethod.funct['writeOut'] == None:
                        chainMethod.run(chainObj)
                    self.scribe.addRecord(chainObj)
        output = self.scribe.extract()

        return output

class ChainMethod:
    def __init__(self, **functionDict):
        if functionDict:
            self.funct = functionDict
            self.args = self.funct.get('args')
            self.kwargs = self.funct.get('kwargs')
            self.setting = self.funct.get('setting')
            self.writeOut = self.funct.get('writeOut')

    def update_args(self, update):
        self.args = update

    def update_kwargs(self, update):
        self.kwargs = update

    def update_setting(self, update):
        self.setting = update

    def update(self, attributes, updates): 
        for i in range(len(attributes)):
            if attributes[i] == 'args':
                self.args = updates[i]
            if attributes[i] == 'kwargs':
                self.kwargs = updates[i]
            if attributes[i] == 'setting':
                self.setting = updates[i]

    def withWriteOut(self, chainObj, writeInfo = None):
        write = self.funct['writeOut']['function'] 
        print(self.funct['writeOut']['function'])
        print(write)
        try:
            argsList = self.funct['writeOut']['args']
        except KeyError:
            argsList = None
        try:
            kwargsDict = self.funct['writeOut']['kwargs']
            print(self.funct['writeOut']['kwargs'])
        except KeyError:
            kwargsDict = {}
        
        kwargsDict['writeInfo'] = writeInfo
        print('the kwargs for the write function are')
        print(kwargsDict)

        out = self.run(chainObj)
        print('in withWriteOut method of ChainMethod')
        print(len(out))
        if argsList == None and kwargsDict != None:
            writeOut = write(chainObj, out, **kwargsDict)
        if argsList != None and kwargsDict != None:
            writeOut = write(chainObj, out, *argsList, **kwargsDict)
        #if self.returnSetting == 0:
           # out = functOut
        #if self.returnSetting == 1:
          #  out = writeOut
        return out

    def run(self, chainObj):
        if self.setting == 0:
            chainOut = self.function(chainObj)
        if self.setting == 1:
            chainOut = self.function(chainObj, *self.args)
        if self.setting == 2:
            chainOut = self.function(chainObj, **self.kwargs)
        return chainOut

    def function(self, chainObj, *args, **kwargs):
        funct = self.funct['function']
        return funct(chainObj, *args, **kwargs)

class Scribe:
    def __init__(self, recordMethod = None, extractMethod = None, writeInfo_recordkw = None): 
        self.record = []
        self.recordMethod = recordMethod
        self.extractMethod = extractMethod
        self.writeInfo_recordkw = writeInfo_recordkw
        print(self.writeInfo_recordkw)

    def addRecord(self, chainObj, writeInfo=None):
        if self.recordMethod != None:
            if writeInfo != None:
                kwargs = {self.writeInfo_recordkw: writeInfo}
                rec = self.recordMethod(chainObj, **kwargs)
            else:
                rec = self.recordMethod(chainObj)
        else:
            rec = chainObj
        
        self.record.append(rec)

    def extract(self):
        if self.extractMethod != None:
            output = self.extractMethod(self.record)
        else:
            output = self.record
        return output
