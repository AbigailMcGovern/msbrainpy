import pytest
import numpy as np
from msbrainpy.brainseries._preprocessing import image_PCA, _grey_scaled, _white_balancing

# Fake Data
# ---------
def make_fake_RGB(x, y):
    img = np.random.randint(0, 256, size=x*y*3, dtype=np.uint8).reshape(y, x, 3)
    return img 

image = make_fake_RGB(100, 100)

# Image PCA
# ---------

def test_image_PCA(image=image):
    """
    Test image_PCA function
    """
    pf = 'image_PCA: '
    _image_PCA_PC1(image, pf)
    _image_PCA_PC2(image, pf)
    _image_PCA_ncomp(image, pf)


def _image_PCA_PC1(image, pf):
    """
    Test the image PCA on RGB-like data and retrive a single channel image
    representing the first principle component. 
    """
    pca_im = image_PCA(image, n_components=1, get_component=0)
    # check that the image is of the correct shape
    spf = '1st component: '
    ck = pca_im.shape
    m = f"the image is of incorrect shape with shape {ck}"
    assert ck == (100, 100), pf + spf + m


def _image_PCA_PC2(image, pf):
    """
    Test the image PCA on RGB-like data and retrive a single channel image
    representing the second principle component. 
    """
    pca_im = image_PCA(image, n_components=2, get_component=1)
    # check that the image is of the correct shape
    spf = '2nd component: '
    ck = pca_im.shape
    m = f"the image is of incorrect shape with shape {ck}"
    assert ck == (100, 100), pf + spf + m


def _image_PCA_ncomp(image, pf):
    """
    Test the image PCA on RGB-like data and retrive a 3 channel image
    representing components of the PCA. 
    """
    pca_im = image_PCA(image, n_components=3, get_component=None)
    # check that the image is of the correct shape
    spf = 'All components: '
    ck = pca_im.shape
    m = f"the image is of incorrect shape with shape {ck}"
    assert ck == (100, 100, 3), pf + spf + m


# Greyscale + rescale
# -------------------

def test_grey_scaled(image=image):
    """
    Test the skimage.transform.rescale + skimage.color.rgb2grey wrapper function
    """
    img = _grey_scaled(image, scale=0.75)
    # check that the image is of the correct shape
    pf = '_grey_scaled: '
    ck = img.shape
    m = f"the image is of incorrect shape with shape {ck}"
    assert ck == (75, 75), pf + m


# white balancing
# ---------------

def test_white_balancing(image=image):
    """
    Test the white balancing function
    """
    img = _white_balancing(image)
    pf = '_white_balancing: '
    ck = img.shape
    m = f"the image is of incorrect shape with shape {ck}"
    assert ck == (100, 100, 3), pf + m
  