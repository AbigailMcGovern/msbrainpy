import os
import numpy as np
from skimage import exposure
from skimage import filters
from skimage import measure
from skimage import morphology
from skimage import util
from ._BrainSeries import BrainSeries
from ._preprocessing import image_PCA, _grey_scaled, _white_balancing

# File contains functions specific to processing images in a series of in situ
# images that represent a single brain and resolving information about 
# individual samples.

# -----------------------------------------------------------------------------
# InSituSeries
# ------------
class InSituSeries(BrainSeries):
    def __init__(self, brain_directory, name_pattern, mask_scale=0.125, target_scale=0.05):
        """
        Class that contains all 
        
        Several options exist as to the type of volume that should be buit:
         * RGB
         * PC1 signal
         * Greyscale (interted or original)
         * TODO: Segmented signal
         * TODO: aligned volume + template

         TODO: 
        """
        self._get_mask = self._ish_tissue_segmentation
        self._get_target_slices = self._ish_target_slices
        self._get_signal_slices = self._ish_expressing_pixels
        self._mask_scale = mask_scale
        self._target_scale = target_scale

        super().__init__(
            brain_directory=brain_directory,
            name_pattern=name_pattern,
            get_mask=self._get_mask, 
            mask_scale=self._mask_scale, 
            get_target_slices=self._get_target_slices, 
            target_scale=self._target_scale, 
            get_signal_slices=self._get_signal_slices, 
            get_signal_slices_kw={}
        )

    # -----------------------------------------------------------------------------
    # Tissue mask
    # -----------

    def _ish_tissue_segmentation(self, image_rgb, disk_denominator=100, scale=0.125, 
                                 **kwargs):
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
        return tissue_mask

    
    # -------------------------------------------------------------------------
    # target slices
    # -------------

    def _ish_target_slices(self, image):
        grey = _grey_scaled(image, self._target_scale)
        grey = util.img_as_ubyte(grey)
        return grey
    
    # -------------------------------------------------------------------------
    # Signal qunatification
    # ---------------------

    def _ish_expressing_pixels(self, image, **kwargs):
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
# Write out
# ---------
# add optional save formats (zarr, tiff). Set zarr to default?
def save_brain_series(directory, out_directory, 
                      image_name_pattern=r'image_id-\d*.jpeg', 
                      find_tissue_mask=True, save_pc1=True):
    """
    Process all images in an ish series (i.e., one brain). 
    

    Parameters
    ----------
    directory: str
        path to the series of images
    out_directory: str
        path to which to save the output
    image_name_pattern: r string
        pattern of images to process. Defaults to the naming convention used
        when downloading images from the Allen Brain database 
            see: msbrainpy.query.download_gene_images
    find_tissue_mask: bool
        Should a downsampled, contrast enhanced image + tissue mask be produced
    save_pc1: bool
        Shouls a pc1 signal image be produced. This image represents the 
        RGB combination (~ in situ dye) associated with the greatest shared
        varience in the data
    """
    #os.makedirs(out_directory, exist_ok=True)
    # find the files
    #files = []
    #image_name_pattern = re.compile(image_name_pattern)
    #for file in sorted(os.listdir(directory)):
        #match = image_name_pattern.match(file)
        #if match is not None:
            #files.append(os.path.join(directory, match[0]))
    # iterate though the images
    #image_collection = io.imread_collection(files)
    #for i in range(len(image_collection)):
        #if find_tissue_mask:
         #   save_series_masks(directory, out_directory, 
         #   image_name_pattern=image_name_pattern)
        #if save_pc1:
          #  white_balanced = _white_balancing()
          #  pc1_image = image_PCA(white_balanced)
          #  pc1_image = util.img_as_ubyte(pc1_image)
          #  pc1_image_name = str(i) + '_pc1_greyscale.tif'
          #  pc1_image_path = os.path.join(out_directory, pc1_image_name)
          #  io.imsave(pc1_image_path, pc1_image)
    #return out_directory
    pass


def save_series_masks(directory, out_directory, 
                       image_name_pattern=r'image_id-\d*.jpeg'):
    '''
    Apply tissue_segmentation function to a gene series
    '''
    #os.makedirs(out_directory, exist_ok=True)
    # find the files
    #files = []
    #image_name_pattern = re.compile(image_name_pattern)
    #for file in sorted(os.listdir(directory)):
        #match = image_name_pattern.match(file)
        #if match is not None:
            #files.append(os.path.join(directory, match[0]))
    # iterate though the images
    #image_collection = io.imread_collection(files)
    #for i in range(len(image_collection)):
        #white_balanced = white_balancing(image_collection[i])
        #tissue_mask, greyscale = tissue_segmentation(white_balanced)
        #mask_name = str(i) + '_tissue_mask.tif'
        #mask_path = os.path.join(out_directory, mask_name)
        #tissue_name = str(i) + '_slice.tif'
        #tissue_path = os.path.join(out_directory, tissue_name)
        #io.imsave(mask_path, tissue_mask)
        #io.imsave(tissue_path, greyscale)
    #return out_directory
    pass


def save_pc1_signal():
    pass