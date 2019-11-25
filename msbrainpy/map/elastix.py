import os
import numpy as np
from tifffile import imread
from tifffile import TiffWriter

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

