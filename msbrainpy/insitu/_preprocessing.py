import numpy as np
from skimage.color import rgb2grey
from skimage import util
from skimage.transform import rescale
from sklearn.decomposition import PCA


# File contains functions related to processing in situ hybridisation images
# such that they can be used to resolve biologically meaningful information
# -----------------------------------------------------------------------------
# RGB preprocessing
# -----------------
# preprocessing for signal quantification
def image_PCA(rgb_image, n_components=1, get_component=None, clip_neg=True):
    """
    Get a pca-based intensity image. Designed to be used to extract highest
    shared variance components from RGB colour channels. PCA implemented
    via scikit-learn.
    E.g., extract a component corresponding to the dye colour in an in situ 
    image (purple-blue colour comprised shared R~B variance)

    Parameters
    ----------
    rgb_image: numpy.ndarray
    n_components: int
        number of components to include in the pca model
    get_component: None or int or tupule of int or slice, optional
        The componets that should be retrived from the PCA. 
        pca shape is (num-pixels, n_components)
        get components chooses components --> pixel values
    clip_neg: bool, optional
        should the negative PCA values be cliped out to reduced noise

    Returns
    -------
    image_pca: np.ndarray
        image whose pixels are represented by a vector of PCA components
        or a single chosen componet.

    Notes
    -----
    Have only tested when retriving one component. 
    Quite slow. Haven't tested, but would be the for loop.
    This can be fixed, i'm sure.

    References:
    [1] The Image Processing Handbook, John C. Russ and F. Brent Neal. 
    CRC Press, Boca Raton, FL, 2015, 1053 pp. <include pages>
    ISBN: 978-1498740265.
    [2] Neal, B. and Russ, J.C., 2004. Principal components analysis of 
    multispectral image data. Microscopy Today, 12(5), pp.36-39.
    """
    rgb_image = _white_balancing(rgb_image)
    image_dims = rgb_image.shape
    flattened_image = rgb_image.reshape(-1, 3)
    pca = PCA(n_components=n_components)
    image_pca = pca.fit_transform(flattened_image)
    if get_component is not None:
        image_pca = image_pca[:, get_component]
    image_pca = image_pca.reshape((image_dims[0], image_dims[1]))
    if clip_neg:
        # replace with a faster method. This could be unusably slow
        # for many images
        for y in range(image_pca.shape[0]):
            for x in range(image_pca.shape[1]):
                if image_pca[y, x] < 0:
                    image_pca[y, x] = 0
    # scale to [0,1]
    maximum = np.max(highest_var)
    highest_var = highest_var/maximum
    return highest_var


# preprocessing for both tissue segmentaiton and image PCA
def _white_balancing(rgb_image):
    '''
    White balancing of RGB images. Finds brightest pixel in greyscale image and
    scales RGB values so as to set this pixel to white (should correct for hue)

    Parameters
    ----------
    image_rgb: np.ndarray
        Shape is (y, x, 3)

    Returns
    -------
    image: np.ndarray
        White balanced RGB image
    '''
    image = rgb_image.copy()
    grey = rgb2grey(rgb_image.copy())
    brightest_index = np.unravel_index(np.argmax(grey, axis=None), grey.shape)
    r, g, b = image[:, :, 0], image[:, :, 1], image[:, :, 2]
    brightest_pixel = image[brightest_index[0], brightest_index[1], :]
    wr, wg, wb = brightest_pixel[0], brightest_pixel[1], brightest_pixel[2]
    lum = wr + wg + wb
    r = r * lum / wr
    g = g * lum / wg
    b = b * lum / wb
    return image

# -----------------------------------------------------------------------------
# Alignment preprocessing
# -----------------------
# preprocessing for signal quantification
def _grey_scaled(image_rgb, scale=0.125, invert=True):
    """
    Helper function to white balance RGB image, convert to greyscale and invert
    if required. 

    Parameters
    ----------
    image_rgb: np.ndarray of int 
        RGB image of shape (y, x, 3)
    scale: float, int
        factor by which to scale the image
    invert: bool
        should the inverted image be returned
    
    Returns
    -------
    grey: np.ndarray of float (check this)
        greyscale image 

    Notes
    -----
    Uses bicubic interpolation for rescaling operation.
    """
    image_rgb = _white_balancing(image_rgb)
    grey = rgb2grey(image_rgb.copy())
    grey = rescale(grey, scale, order=3)
    if invert:
        grey = util.invert(grey)
    return grey
