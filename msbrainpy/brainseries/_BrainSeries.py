import os
import re
import dask.array as da
import numpy as np
from dask import delayed
from dask.delayed import Delayed
from dask import distributed
from skimage import measure
from skimage.io import imread
# from toolz import curry
# need to get the allen brain atlas template

# -----------------------------------------------------------------------------
# BrainSeries Class
# -----------------------------------------------------------------------------
class BrainSeries:
    def __init__(self, brain_directory, name_pattern, get_mask=None, 
                 mask_scale=0.125, get_target_slices=None, 
                 target_scale=0.05, get_signal_slices=None, 
                 get_signal_slices_kw={}, # registration_3d=True
                 ): 
        """
        BrainSereies class: Automated processing of a series of brain slices.

        Scientists can easily produce sereies of images (typically slices in 
        either the coronal or sagittal planes) that represent a single brain 
        sample. These images can be difficult to understand in the context of 
        the original 3D biologial sample without (1) computational assembly 
        such that the sample can be visualised as a whole and (2) without 
        easy unbiased alignment to an atlas.

        This class must be provided functions to do the following:
        (will depend upon type of data. Functions may exist in a subclass)
            [1] a means to find a segmenation of the tissue vs background
                (used to find initial coordinates of the image in the volume)
            [2] a means to obtain images which will be used to build the 
                template for registration.
            [3] a means to obtain a representation of the signal to be 
                quantified against the atlas. This representaiton could have
                the same shape as the original image or 

        This class contains (or may oneday) methods for the following:
            * lazy loading of images (dask delayed skimage.io.imread)
            * lazy generation of tissue masks from which to compute initial
                position of image in volume (based on provided function [1])
            * TODO: calculation of volume positions from tissue masks
            * TODO: affine correction of slice alignment (via scikit-image)
            * TODO: the actual atlas alignment (via scikit-image)
            * TODO: execute the signal quantification 

        Parameters
        ----------
        brain_directory: str
            Directory in which to find the series of brain images.
            (preferably absolute file path)
        name_pattern: str
            Regular expression with which to find images (r-string)
        get_mask: function
            Function with which to process original images to obtain a
            binary mask depicting the location of the tissue.
        get_mask_kw: dict
            Keyword arguments for get_mask function. Default {}.
        get_target_slices: function
            Function with which to process original images to obtain 
            target slices with which to register brain to atlas.
        get_template_slices_kws: dict
            Keyword arguments for get_target_slices function.
        get_signal_slices: function
            Function with which to obtain signal to be quantified against atlas
        get_signal_slices_kw: dict
            Keyword arguments for get_signal_slices function.
        # registration_3d: bool
            # Not yet relevant. If a 2D-3D registration method becomes available
            # when False, will perform 2D-3D registration but will still build 
            # volume for visualisation perposes.

        Attributes
        ----------
        brain_directory: str
            The absolute file path to the directory (hopefully... :) )
        name_pattern: str
            An r-string with which to find the brain images. Note that the images
            must be named such that the sorted list of names will be in order.
            To be fair I've not yet decided which order, hasn't mattered yet. 
                (i.e., rostral-caudal [probably] or caudal-rostral) 
        image_paths: list of str
            absolute file paths to images of the brain in question
        volume_shape: tuple
            The shape of the volume. For some reason I haven't included the 
            RGB channel ... not sure why??? (TODO?)
        image_volume: dask array from delayed
            Mearly the idea of how to put the images into an appropriate order
            and place. 
            TODO: Note that the volume generation process needs to be improved 
            to put a better value than zero in the padding. In situ images 
            have white backgrounds. Perhaps the average value of the edge
            values in the first image?
        target_volume: np.ndarray


        Methods
        -------
        ...
        """

        m0 = "BrainSeries class cannot be instantiated without get_mask"
        m1 = " function: \nsee documentation" 
        m = m0 + m1
        # check that the function to obtain the tissue mask exists
        if get_mask is not None:
            self._get_mask = get_mask
        else:
            raise NotImplementedError(m)
        # check that the function to obtain the obtain target slices exists
        if get_target_slices is not None:
            self._get_target_slices = get_target_slices
        else:
            raise NotImplementedError(m)
        # check that the function to obtain the signal quantification exists
        if get_signal_slices is not None:
            self._get_signal_slices = get_signal_slices
        else:
            raise NotImplementedError(m)
        
        # BASIC ATTRIBUTES
        # ----------------
        self._brain_directory = brain_directory
        self._name_pattern = re.compile(name_pattern)
        self._image_paths = self._get_image_paths()
        self._n = len(self._image_paths)
        # mask scale (i.e., how much smaller should the tissue masks be)
        self._mask_scale = mask_scale 
        # target volume scale for registration
        self._target_scale = target_scale 
        # tissue masks list, which will first be computed when building the 
        # volume plan. It will be persisted in memory until it is no longer
        # needed, i.e., when the volume is built in self._build_mask_volume()

        # USEFUL LISTS
        # ------------
        # get the image shapes, the list of mask images, and the list of
        # target images
        the_lists = self._get_the_lists()
        self._dims_list, self._masks_list, self._targets_list = the_lists
        del the_lists

        # VOLUME ATTRIBUTES
        # -----------------
        # calculate initial positions in volume
        self._volume_info = {'coordinates_list' : []}
        self._get_volume_plan() 
        # the z, y, x dims (does not include RGB)
        self._shape = self._volume_info['shape']
        self._centroid = self._volume_info['local_centroid']    

        # VOLUMES
        # -------
        # move the images from lists into the appropriate course position in 
        # the volume.
        # delayed dask array
        self._image_volume = self._lazy_image_volume()
        # mask volume (np.ndarray)
        # self._mask_volume = self._build_mask_volume()
        # nd.array registration volume (needs to be small, not a memory waste)
        self._target_volume = self._build_target_volume()

    # -------------------------------------------------------------------------
    # BASIC
    # -----

    @property
    def brain_directory(self):
        """
        Absolute file path to the directory
        """
        return self._brain_directory


    @property
    def name_pattern(self):
        """
        """
        return self._name_pattern


    @property
    def shape(self):
        """
        """
        return self._shape


    @property
    def image_paths(self):
        return self._image_paths


    def _get_image_paths(self):
        """
        """
        image_paths = []
        for file_ in sorted(os.listdir(self._brain_directory)):
            match = self._name_pattern.match(file_)
            if match is not None:
                image_paths.append(os.path.join(self._brain_directory, 
                                                      match[0]))
        return image_paths


    def _lazy_images(self):
        """
        """
        imread_lazy = delayed(imread, pure=True)
        lazy_images = [imread_lazy(path) for path in self._image_paths]
        return lazy_images


    # -----------------------
    # Apply supplied function
    # -----------------------
    # 
    def _get_the_lists(self):
        """
        """
        lazy_images = self._lazy_images()
        masks = []
        targets = []
        shapes = []
        for image in lazy_images:
            image = image.compute()
            shape = image.shape
            shapes.append(shape)
            mask = self._get_mask_function(image) 
            masks.append(mask)
            target = self._get_target_function(image)
            targets.append(target)
        return shapes, masks, targets


    def _get_mask_function(self, image):
        """
        """
        if self._get_mask is not None:
            return self._get_mask(image)

 
    def _get_target_function(self, image):
        """
        """
        if self._get_target_slices is not None:
            return self._get_target_slices(image)

    
    # -------------------------------------------------------------------------
    # VOLUME
    # ------
    # Volumes
    # -------

    # image volume pyramid (comprised of all in situ images)
    @property
    def image_volume(self):
        """
        Volumetric dask array comprised of the original RGB images.
        Primarily designed for viewing the data set as a whole. 
        """
        # return self._image_volume 
        return self._image_volume


    def save_image_volume(self, directory, format='zarr'):
        """
        Save the image volume to file as a pyramid.

        PARAMETERS
        ----------
        directory: str
            absolute file path describing where the image should be saved
        format: str ('zarr' or 'tiff')
            File format. Defaults to Zarr.
        """
        if format is 'tiff':
            pass
        if format is 'zarr':
            pass                         

    
    # tissue mask volumetric image
    @property
    def mask_volume(self):
        """
        Dask array comprised of tissue masks. 
        For the visualisation of 
        """
        masks = self.mask_volume.compute()
        # return masks
        return masks


    def save_mask_volume(self, directory, format='tiff'):
        """
        """
        if format is 'tiff':
            pass
        if format is 'zarr':
            pass   

    
    # registration template volumetric image
    @property
    def target_volume(self):
        """
        """
        # return self._signal_volume
        return self._target_volume


    def save_template_volume(self, directory, format='tiff'):
        """
        """
        if format is 'tiff':
            pass
        if format is 'zarr':
            pass   

    
    # ---------------
    # Volume Building
    # --------------- 

    def _lazy_image_volume(self):
        """
        """
        # original images with which to build the volume
        lazy_images = self._lazy_images()
        image_volume = self._generate_volume(lazy_images)
        return image_volume
        

    def _build_mask_volume(self):
        """
        """
        scale = self._mask_scale
        mask_volume = self._generate_volume(self._masks_list, scale=scale)
        mask_volume = mask_volume.compute()
        return mask_volume


    def _build_target_volume(self):
        """
        """
        scale = self._target_scale
        target_volume = self._generate_volume(self._targets_list, scale=scale)
        target_volume = target_volume.compute()
        return target_volume


    def _generate_volume(self, image_list, scale=1):
        """
        """
        image_info = self._volume_info['coordinates_list']
        # get the shape of each plane in the x y axis
        shape = self._shape
        shape = np.round(np.array(shape) * scale).astype(int)
        shape = (shape[1], shape[2])
        # find out if the image is RGB shaped
        if isinstance(image_list[0], Delayed):
            image = image_list[0].compute()
        else:
            image = image_list[0]
        # get the data type
        dtype = image.dtype
        # add the RGB dim if necessary
        if image.shape[-1] == 3:
            shape = (shape[0], shape[1], 3) 
        del image
        # get a list of delayed arrays representing padded images
        arrays = [da.from_delayed(self._padded_image(z, image_list, 
                                  image_info, shape, scale), 
                                  shape, dtype=dtype) 
                  for z in range(self._n)]
        # get dask array representing image volume
        volume = da.stack(arrays, axis=0)
        return volume


    @delayed
    def _padded_image(self, z, image_list, image_info, shape, scale):
        """
        """
        image = image_list[z]
        # initialise the full size image plane
        new_image = np.zeros(shape=shape, dtype=np.uint8)
        # find the bounding box for the tissue in the image
        bounding_box = np.array(image_info[z]['bounding_box'])
        # scale according to the factor difference from full scale
        bounding_box = np.floor(bounding_box * scale).astype(int)
        min_y, min_x, max_y, max_x = bounding_box
        # get only the bounded tissue
        image = image[min_y:max_y, min_x:max_x]
        # find the coordinate to move the tissue into in the full plane
        coordinates = np.array(image_info[z]['volume_coordinates'])
        coordinates = np.floor(coordinates * scale).astype(int)
        v_min_y, v_min_x, v_max_y, v_max_x = coordinates
        try:
            new_image[v_min_y:v_max_y, v_min_x:v_max_x] = image
        except ValueError:
            wrong_shape = np.array(new_image[v_min_y:v_max_y, v_min_x:v_max_x]
                                   .shape[:2])
            img_shape = np.array(image.shape[:2])
            new_coords = self._correct_coords_for_shape(img_shape, 
                                                        wrong_shape, 
                                                        (v_min_y, v_max_y,
                                                        v_min_x, v_max_x))
            v_min_y, v_max_y, v_min_x, v_max_x = new_coords
            new_image[v_min_y:v_max_y, v_min_x:v_max_x] = image
        return new_image


    def _correct_coords_for_shape(self, image_shape, wrong_shape, 
                                  coords):
        """
        Probably poor attempt at function to force image into the padded image.
        Corrects the coordinates so that they fit the image to be placed.
        """
        diff = image_shape - wrong_shape
        sign = lambda x: x / abs(x)
        half = lambda x: int(abs(x) // 2)
        new_coords = []
        counter = 0
        for c in diff:
            if c != 0:
                h = half(c)
                s = sign(c)
                corr = h * s
                mod = abs(c) % 2
                if mod:
                    min_ = coords[counter] - corr
                    max_ = coords[counter + 1] + corr + s
                if not mod:
                    min_ = coords[counter] - corr
                    max_ = coords[counter + 1] + corr
                new_coords.append(int(min_))
                new_coords.append(int(max_))
            else:
                new_coords.append(coords[counter])
                new_coords.append(coords[counter + 1])
            counter += 2
        return tuple(new_coords)


    # -------------------
    # Volume Calculations
    # -------------------

    def _get_volume_plan(self):
        """
        """
        # find the information about the location of the tissue slice in each image
        self._get_mask_info()
        self._calculate_positions()
        

    def _get_mask_info(self):
        """
        """
        # get the list of tissue masks to compute
        #lazy_masks = self._tissue_masks_lazy()
        # initialise the dictionary which will be filled with information about 
        volume_plan = self._volume_info
        # get the number with which to scale coordinates to full size
        scale = np.divide(1, self._mask_scale)
        # find the important properties
        for i in range(self._n):
            mask = self._masks_list[i]
            mask = mask.astype(int)
            # get region properties
            props = measure.regionprops(mask)[0]
            # find the bounding box
            bounding_box = np.round(np.array(props['bbox']) * scale)
            bounding_box = bounding_box.astype(int)
            # find the centroid within the bounding box
            centroid = np.round(np.array(props['local_centroid']) * scale)
            centroid = centroid.astype(int)
            volume_plan['coordinates_list'].append({})
            volume_plan['coordinates_list'][i]['index'] = i
            volume_plan['coordinates_list'][i]['bounding_box'] = bounding_box
            volume_plan['coordinates_list'][i]['local_centroid'] = centroid
            self._volume_info = volume_plan


    def _calculate_positions(self):
        """
        """
        # the information 
        image_info  = self._volume_info['coordinates_list']
        # get the number of pixels greater than the centroid in y-axis within
        # the bounding box for every image.
            # max y - min y (size of bounding box) - centroid y
        y_upper = [img['bounding_box'][2] - img['bounding_box'][0] \
                    - img['local_centroid'][0] for img in image_info]
        # get the number of pixels between the edge of bounding box and 
        # the centroid in the y-axis.
            # centroid y
        y_lower = [img['local_centroid'][0] for img in image_info]
        # get the number of pixels greater than the centroid in y-axis within
        # the bounding box for every image.
            # max x - min x (size of bounding box) - centroid x
        x_upper = [img['bounding_box'][3] - img['bounding_box'][1] \
                     - img['local_centroid'][1] for img in image_info]
        # get the number of pixels between the edge of bounding box and 
        # the centroid in the x-axis.
            # centroid y
        x_lower = [image['local_centroid'][1] for image in image_info]
        # the y dimention: based on maximum distances from the y centre point
        y_range = np.max(y_upper) + np.max(y_lower) + 2
        # the y centroid for the volume
        y_cent = np.max(y_lower)
        # the x dimention: based on maximum distances from the y centre point
        x_range = np.max(x_upper) + np.max(x_lower) + 2
        # the x centroid for the volume
        x_cent = np.max(x_lower)
        # the z dimention (i.e., num images)
        z_range = self._n
        # volume shape
        shape = (z_range, y_range, x_range)
        # x-y centroid
        centroid = (y_cent, x_cent)
        self._volume_info['shape'] = shape
        self._volume_info['local_centroid'] = centroid
        for i in range(len(image_info)):
            y_max = centroid[0] + y_upper[i] # y
            y_min = centroid[0] - y_lower[i]
            x_max = centroid[1] + x_upper[i]
            x_min = centroid[1] - x_lower[i]
            c = (y_min, x_min, y_max, x_max)
            self._volume_info['coordinates_list'][i]['volume_coordinates'] = c
        return shape, centroid


    # -------------------------------------------------------------------------
    # REGISTRATION?
    # -------------


    # -------------------------------------------------------------------------
    # QUANTIFICATION METHODS
    # ----------------------

    # -------------------------------------------------------------------------
    # VISUALISATION METHODS
    # ---------------------
        