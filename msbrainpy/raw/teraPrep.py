import os
import re
import shutil
from collections import defaultdict
from pathlib import Path

def writeTeraFiles(sourceDir, outDir, prefix, leftLS = 'UltraII Filter0000.tif', 
                        rightLS = 'UltraII Filter0001.tif', channels = r"_C[0-9]{2}_xyz", chanName = r"C[0-9]{2}", 
                        leftTiles = [" x 00", " x 01"], rightTiles = [" x 02", " x 03"], xTiles = 4, yTiles = 4, 
                        X = 2.032, Y = 2.032, Z = 3, xWidth = 2160, yLength = 2560, zDepth = 2031, overlap = .2, 
                        tileString = "[0{0} x 0{1}]", zSlice = r"xyz-Table\sZ\d{4}"):

    sourceDir = Path(sourceDir)
    outDir = Path(outDir)
    newName = prefix+'_TS'
    newPath = outDir/newName
    if newPath.exists() != True:
        os.mkdir(newPath)
    
    xSize = X*xWidth
    ySize = Y*yLength
    zSize = Z*zDepth
    invOL = 1-overlap
    
    # get file names for eah channel
    filelist = os.listdir(sourceDir)
    lightsheets = getLightsheets(filelist, leftLS, rightLS)
    channelsDict = getChannels(lightsheets, newPath, channels, chanName, leftTiles, rightTiles)
    writeChannelFiles(channelsDict, sourceDir, xTiles, yTiles, ySize, xSize, Z, invOL, zSlice, tileString)
    
    
def getLightsheets(filelist, leftLS, rightLS):
    leftLightsheet = []
    rightLightsheet = []
    for _file in filelist:
        if _file.endswith(leftLS):
            leftLightsheet.append(_file)
        if _file.endswith(rightLS):
            rightLightsheet.append(_file)
    print('{} and {} files were found for the left and right lightsheets, respectively'.format(len(leftLightsheet), len(rightLightsheet)))
    LR = [leftLightsheet, rightLightsheet]

    return LR

def getChannels(lightsheets, newPath, channels, chanName, leftTiles, rightTiles):
    uniqueChannels = []
    for _file in lightsheets[0]:
        a = re.search(channels, _file)
        if a.group(0) not in uniqueChannels:
            uniqueChannels.append(a.group(0))
    print('Unique channels found: {}'.format(len(uniqueChannels)))
    newName = newPath.stem
    
    channelDict = defaultdict()
    channelPaths = []
    for channel in uniqueChannels:
        b = re.search(chanName, channel)
        b = b.group(0)
        channelName = newName+"_"+b
        channelPath = newPath/channelName
        channelPaths.append(channelPath)
        os.mkdir(channelPath)
        
        leftFiles = lightsheets[0]
        rigthFiles = lightsheets[1]
        
        channelFiles = []
        
        for _file in leftFiles:
            for tile in leftTiles:
                if _file.find(tile) != -1 and _file.find(channel) != -1:
                    channelFiles.append(_file)
        
        for _file in leftFiles:
            for tile in rightTiles:
                if _file.find(tile) != -1 and _file.find(channel) != -1:
                    channelFiles.append(_file)
        print('{} files were found for channel {}'.format(len(channelFiles), b))
        print('Only appropriately sided tiles were added')
        channelDict[b] = defaultdict()
        channelDict[b]['chanString'] = channel
        channelDict[b]['files'] = channelFiles
        channelDict[b]['path'] = channelPath               
        
    return channelDict
        
def writeChannelFiles(channelsDict, sourceDir, xTiles, yTiles, ySize, xSize, Z, invOL, zSlice, tileString):
    for channel in channelsDict.keys():
        channelFiles = channelsDict[channel]['files']
        channelPath = channelsDict[channel]['path']
        for row in range(yTiles):
            micron10thPosition = ySize*invOL*row*10
            micron10thPosition = round(micron10thPosition)
            micron10thPosition = "{:06d}".format(micron10thPosition)
            newPath = channelPath/micron10thPosition
            os.mkdir(newPath)
                        
        rowsInChannel = os.listdir(channelPath)
        for row in rowsInChannel:
            rowDir = channelPath/row
            for tile in range(xTiles):
                micron10thPosition = xSize*invOL*tile*10
                micron10thPosition = round(micron10thPosition)
                micron10thPosition = "{:06d}".format(micron10thPosition)
                newDir = row+"_"+micron10thPosition
                newPath = rowDir/newDir
                os.mkdir(newPath)
                        
        for i in range(yTiles):
            rowDir = channelPath/rowsInChannel[i]
            colsInRow = os.listdir(rowDir)
            for j in range(xTiles):
                colDir = rowDir/colsInRow[j]
                tile = tileString.format(i, j)
        
                filesToCopy = []
                for _file in channelFiles:
                    if _file.find(tile) != -1:
                        filesToCopy.append(_file)
            
                for _file in filesToCopy:
                    fileSource = sourceDir/_file
                    fileDest = colDir/_file
                    shutil.copy2(fileSource, fileDest)
                
                stackFiles = os.listdir(colDir)
                for Zslice in stackFiles:
                    a = re.search(zSlice, Zslice)
                    sliceNo = re.search(r"\d{4}", a.group(0))
                    sliceNo = int(sliceNo.group(0))
                    position = round(sliceNo*Z*10)
                    posString = '{:06d}'.format(position)
                    newName = posString+'.tif'
                    oldAdress = colDir/Zslice
                    newAdress = colDir/newName
                    os.rename(oldAdress, newAdress)
                print('{} files were saved to {}. \nFiles were renamed according to position'.format(len(stackFiles), colDir))

#### Re write this ugly, ugly shit! Generators. Obviously.... not yet done

def makeTeraDirs(channelsDict, yTiles, ySize, invOL):
    for channel in channelsDict.keys():
        channelFiles = channelsDict[channel]['files']
        channelPath = channelsDict[channel]['path']
        for row in range(yTiles):
            micron10thPosition = ySize*invOL*row*10
            micron10thPosition = round(micron10thPosition)
            micron10thPosition = "{:06d}".format(micron10thPosition)
            newPath = channelPath/micron10thPosition
            os.mkdir(newPath)
                        
        rowsInChannel = os.listdir(channelPath)
        for row in rowsInChannel:
            rowDir = channelPath/row
            for tile in range(xTiles):
                micron10thPosition = xSize*invOL*tile*10
                micron10thPosition = round(micron10thPosition)
                micron10thPosition = "{:06d}".format(micron10thPosition)
                newDir = row+micron10thPosition
                newPath = rowDir/newDir
                os.mkdir(newPath)

    return rowsInChannel


def copy2TeraDir(rowsInChannel, yTiles, xTiles, tileString):
    for i in range(yTiles):
        rowDir = channelPath/rowsInChannel[i]
        colsInRow = os.listdir(rowDir)
        for j in range(xTiles):
            colDir = rowDir/colsInRow[j]
            tile = tileString.format(i, j)
        
            filesToCopy = []
            for _file in channelFiles:
                if _file.find(tile) != -1:
                    filesToCopy.append(_file)
            
            for _file in filesToCopy:
                fileSource = sourceDir/_file
                fileDest = colDir/_file
                shutil.copy2(fileSource, fileDest)
                yield colDir

def renameInStackDir(colDir, zSlice, Z):
    stackFiles = os.listdir(colDir)
    for Zslice in stackFiles:
        a = re.search(zSlice, Zslice)
        sliceNo = re.search(r"\d{4}", a.group(0))
        sliceNo = int(sliceNo.group(0))
        position = round(sliceNo*Z*10)
        posString = '{:06d}'.format(position)
        newName = posString+'.tif'
        oldAdress = colDir/Zslice
        newAdress = colDir/newName
        os.rename(oldAdress, newAdress)
    