import numpy as np
from skimage import exposure
from skimage import filters
from skimage import measure
from skimage import morphology
from skimage import util
from ._preprocessing import _grey_scaled, image_PCA


# File contains functions specific to obraining biologically meaningful data
# from in situ hybridisation images
# -----------------------------------------------------------------------------
# Signal qunatification
# ---------------------

def expressing_pixels(image):
    """
    function to segment expressing pixles from pre-processed in situ images
    """
    # this needs to be ~ equivalent to method used in allen brain institutes 
    # neuroinformatics/anatomical gene expression pipeline
    # Method must segment ~ cell sized objects and dense structures (e.g., 
    # hippocampus, granular layer cerebellum, etc.)
    # Can't be local contrast based as this could lead to false detection in 
    # cellularly dense areas, which tend to pick up dye without appreciable 
    # expression and may look high contrast when no real signal is present 
    pass


# -----------------------------------------------------------------------------
# Tissue mask
# -----------

def tissue_segmentation(image_rgb, disk_denominator=100, scale=0.125, 
                        processed_grey=False):
    '''
    Segement a tissue mask from a downsampled slice image

    Parameters
    ----------
    image_rgb: np.ndarray
        rgb in situ image (full size)
    disk_denominator: int
        used to choose the selm for binary closing morphological opperation
        Used as a scaling factor for the disk
        bigger >> smaller disk 
        larger >> bigger disk
    scale: float
        rescaling factor for the image
    processed_grey: bool
        If True return the processed greyscale image from which the mask was
        extracted. 

    Returns
    -------
    tissue_mask or (tissue_mask, grey): np.ndarray of bool or tuple of 
    np.ndarray of bool and np.ndarray of unit8


    Notes
    -----
    Needs some work. There are some holes in the masks where the tissue is 
    not above background. Some of the circular artifacts are covered by the 
    closing opperation. Perfectly sufficient for aligning images when building
    a volume.
    '''
    grey = _grey_scaled(image_rgb, scale=scale, invert=False)
    grey = exposure.equalize_hist(grey)
    grey_open = morphology.opening(grey, morphology.square(6))
    grey = filters.median(grey)
    threshold = filters.threshold_otsu(grey_open)
    binary = grey_open <= threshold
    binary = morphology.binary_closing(binary, morphology.disk(6))
    labels = measure.label(binary)
    props = measure.regionprops(labels)
    biggest = np.argmax([props[i]['area'] for i in range(len(props))]) + 1
    tissue_mask = labels == biggest
    selm = morphology.disk(np.mean([grey.shape[0], grey.shape[1]]) // disk_denominator)
    tissue_mask = morphology.binary_closing(tissue_mask, selm)
    tissue_mask = morphology.remove_small_holes(tissue_mask, area_threshold=1000)
    grey = util.img_as_ubyte(grey)
    if processed_grey:
        return tissue_mask, grey
    else:
        return tissue_mask

