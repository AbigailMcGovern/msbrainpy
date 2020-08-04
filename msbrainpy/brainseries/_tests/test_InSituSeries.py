import os
import numpy as np
from skimage import io
from msbrainpy.brainseries import InSituSeries


def test_InsituSeries():
    """
    test the current functional capacity of the 
    """
    temp = get_temp_directory()
    make_fake_data(temp, 4000, 5000, 5)
    name_pattern = r'.*\.jpg'
    fake_brain = InSituSeries(temp, name_pattern)
    # test the target brain (registration target) data type and dimensions
    target_brain_tests(fake_brain)
    # get the full image volume, a dask array
    full_brain_tests(fake_brain)
    remove_temp_data(temp)


def get_temp_directory():
    """
    create temp directory
    """
    current = os.getcwd()
    temp = os.path.join(current, "temp")
    os.makedirs(temp, exist_ok=True)
    return temp


def make_fake_data(temp, y, x, num):
    """
    make several jpeg images full of random ints. Images have
    imcreasing dimentions.

    Parameters
    ----------
    temp: str
        temporary directory in which to save the data
    y: int
        minimum size of y axis
    x: int
        minimum size of the x axis
    num: int
        number of data files to produce

    """
    for num in range(num):
        y = y + num * 100
        x = x + num * 100
        image = np.random.randint(0, 256, size=y*x*3, dtype=np.uint8).reshape(y, x, 3)
        name = os.path.join(temp, (str(num) + '.jpg'))
        io.imsave(name, image, plugin='pil')


def remove_temp_data(temp):
    """
    Remove the data in the temp directory (param: temp)
    """
    files = os.listdir(temp)
    for file_ in files:
        path = os.path.join(temp, file_)
        os.remove(path)
    os.rmdir(temp)


def target_brain_tests(fake_brain):
    """
    Test the data type and dimensions of the target brain
    I.e., the brain that will ultimately be used for atlas alignment

    Parameters
    ----------
    fake_brain: InSituSeries 
    """
    # get the volume to be used for registration
    target = fake_brain.target_volume
    # check the type of data for the target brain
    pf = 'InSituSeries: '
    ck = target.dtype
    m = f'target returned incorrect type {ck}'
    assert ck == np.uint8, pf + m
    # check that the brain has 3 dimentions
    ck = target.shape
    m = f'the target volume is of incorrect dimentsions with shape {ck}'
    assert len(ck) == 3, pf + m


def full_brain_tests(fake_brain):
    """
    Test running compute on the full size image volume and assess the 
    shape and type of data.

    Parameters
    ----------
    fake_brain: InSituSeries 

    Notes
    -----
    Should only be around 3 GB and easily fit into memory
    In [1]: from dask.utils import format_bytes
    In [2]: format_bytes(5*4500*5500*3*8)                                                                                                                                                   
    Out[2]: '2.97 GB'
    """
    # get the full image volume, a dask array
    full_image = fake_brain.image_volume
    full_image = full_image.compute()
    pf = 'InSituSeries: '
    # check that this has returned a numpy array
    ck = type(full_image)
    m = f'computing the image volume is not an ndarray and has type {ck}'
    assert ck == np.ndarray, pf + m
    # check that the array has the correct data type
    ck = full_image.dtype
    m = f'computing the image volume has returned an incorrect dtype of {ck}'
    assert ck == np.uint8
    # check that the data has four dimensions
    ck = full_image.shape
    m = f'the image volume is of the wrong dimensionality with shape {ck}'
    # check that last dimension is 3, representing the RGB channel
    m = f'the image volume is not in RGB format with shape {ck}'
    assert ck[-1] == 3, pf + m
