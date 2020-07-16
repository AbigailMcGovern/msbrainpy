import os
from 
from ._volume_assemble import VolumeAssembler
# File contains functions specific to processing images in a series of in situ
# images that represent a single brain and resolving information about 
# individual samples.
# -----------------------------------------------------------------------------
# Signal qunatification
# ---------------------
def atlas_quantification(...):
    pass


def _segmented_signal(...):
    pass


# -----------------------------------------------------------------------------
# Volume construction
# -------------------
# get the tissue masks with which to align the course volume
# note that this function could be changed to accomodate simply obtaining info
# about bounding box & centroid with an option for saving .... ? 

# OLD - 
# required for aligning images to the Allen Brain Atlas
# ensure files for a single brain are stored in a directory
# neither of these functions will be relevant

def volume_assembler(tissue_mask_directory, scale_factor=8, 
                         tissue_mask_pattern=r'\d*_tissue_mask.tif'):
    dictionary = get_volume_assembler_dict(tissue_mask_directory, 
                                           scale_factor=scale_factor,
                                           tissue_mask_pattern=tissue_mask_pattern)
    volume_assembler = VolumeAssembler(dictionary)
    return volume_assembler


def volume_from_series(volume_assembler, series_pattern, directory, 
                           scale_factor=1):
    files = []
    pattern = re.compile(series_pattern)
    for file in os.listdir(directory):
        match = pattern.match(file)
        if match is not None:
            files.append(match[0])
    volume = volume_assembler.generate_volume(directory, files)
    return volume


# -----------------------------------------------------------------------------
# Atlas allignment
# ----------------


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
    os.makedirs(out_directory, exist_ok=True)
    # find the files
    files = []
    image_name_pattern = re.compile(image_name_pattern)
    for file in sorted(os.listdir(directory)):
        match = image_name_pattern.match(file)
        if match is not None:
            files.append(os.path.join(directory, match[0]))
    # iterate though the images
    image_collection = io.imread_collection(files)
    for i in range(len(image_collection)):
        if find_tissue_mask:
            save_series_masks(directory, out_directory, 
            image_name_pattern=image_name_pattern)
        if save_pc1:
            pc1_image = rgb_PCA(white_balanced)
            pc1_image = util.img_as_ubyte(pc1_image)
            pc1_image_name = str(i) + '_pc1_greyscale.tif'
            pc1_image_path = os.path.join(out_directory, pc1_image_name)
            io.imsave(pc1_image_path, pc1_image)
    return out_directory


def save_series_masks(directory, out_directory, 
                       image_name_pattern=r'image_id-\d*.jpeg'):
    '''
    Apply tissue_segmentation function to a gene series
    '''
    os.makedirs(out_directory, exist_ok=True)
    # find the files
    files = []
    image_name_pattern = re.compile(image_name_pattern)
    for file in sorted(os.listdir(directory)):
        match = image_name_pattern.match(file)
        if match is not None:
            files.append(os.path.join(directory, match[0]))
    # iterate though the images
    image_collection = io.imread_collection(files)
    for i in range(len(image_collection)):
        white_balanced = white_balancing(image_collection[i])
        tissue_mask, greyscale = tissue_segmentation(white_balanced)
        mask_name = str(i) + '_tissue_mask.tif'
        mask_path = os.path.join(out_directory, mask_name)
        tissue_name = str(i) + '_slice.tif'
        tissue_path = os.path.join(out_directory, tissue_name)
        io.imsave(mask_path, tissue_mask)
        io.imsave(tissue_path, greyscale)
    return out_directory


def save_pc1_signal():
    pass