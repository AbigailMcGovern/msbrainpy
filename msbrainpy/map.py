import os
import numpy as np
import re
from skimage import io
from skimage import measure
from msbrainpy.quantify.processing import find_gene_series_masks


# ------------------------------------------ 3D Mapping Related Functions ----------------------------------------------

# ---------------------------------------------- Get transformations ---------------------------------------------------

def run_elastix(elastix_bin, moving_image, fixed_image, parameter_files, out_dir):
    if os.path.exists(out_dir) != True:
        os.mkdir(out_dir)
    threads = ' -threads 16'
    moving = ' -m ' + moving_image
    fixed = ' -f ' + fixed_image
    if type(parameter_files) == str:
        param = ' -p ' + parameter_files
    if type(parameter_files) == list:
        params = []
        for p in parameter_files:
            par = ' -p ' + p
            params.append(par)
        param = ''.join(params)
    out = ' -out ' + out_dir
    cmd = elastix_bin + threads + moving + fixed + param + out
    print('Attempting to execute:')
    print(cmd)
    result = os.system(cmd)
    return out_dir


# ------------------------------------------- Convert coord resolution  ------------------------------------------------

def resample_points(points, orig_res, new_res):
    array = points
    for i in range(len(points)):
        array[i, 0] = int(np.ceil(array[i, 0] * orig_res[0] / new_res[0]))
        array[i, 1] = int(np.ceil(array[i, 1] * orig_res[1] / new_res[1]))
        array[i, 2] = int(np.ceil(array[i, 2] * orig_res[2] / new_res[2]))
    print(array[:5, :])
    return array


def correct_to_xyz(cell_array):
    new_array = np.array([cell_array[:, 2], cell_array[:, 1], cell_array[:, 0]])
    new_array = new_array.T
    return new_array


# ------------------------------------------------ Density analysis  ---------------------------------------------------

def get_density_img(save_dir, save_name, points, shape, order):
    img = get_voxel_counts(points, shape, order)
    img = np.round(img * 1000000 / (16 * 16 * 16))
    img = img.astype(np.uint16)
    save_path = os.path.join(save_dir, save_name)
    with tiff_writer(save_path) as tiff:
        tiff.save(img)
    return img


def get_voxel_counts(points, shape, order):
    """
    Expects int style points, not float
    """
    z = order[0]
    y = order[1]
    x = order[2]
    img = np.zeros(shape, dtype=np.uint16)
    for point in points:
        try:
            ind = (point[z], point[y], point[x])
            img[ind] += 1
        except index_error:
            print('point {}, {}, {} out of range'.format(point[z], point[y], point[x]))
    return img


def transform_img(transformix_bin, in_tiff, parameter_files, out_dir, prefix, out_name='img_transformed'):
    out_name0 = prefix + "_" + out_name
    out_path = os.path.join(out_dir, out_name0)
    out_dir_new = os.path.join(out_path, 'result.mhd')
    if os.path.exists(out_path) != True:
        os.mkdir(out_path)
    source_path = in_tiff
    points = ' -in ' + source_path
    out = ' -out ' + out_path
    param = ' -tp ' + parameter_files
    cmd = transformix_bin + points + out + param
    print('Executing the following:')
    print(cmd)
    res = os.system(cmd)
    if res != 0:
        print('An error occured, check the parameters and associated files')
    return out_dir_new


# ------------------------------------------ Coordinates transformation  -----------------------------------------------

def write_points_for_transformix(points, prefix, save_dir, res):
    save_name = prefix + '_' + str(res) + 'um_coords_xyz_elastix_format.txt'
    save_path = os.path.join(save_dir, save_name)
    header = 'point\n{}'.format(len(points))
    np.savetxt(save_path, points, header=header, comments='')
    return save_path


# To be useful, this would require the transformation to be inverse.
def transform_points(transformix_bin, points_txt, parameter_files, out_dir, prefix, out_name='points_transformed'):
    out_name0 = prefix + "_" + out_name
    out_path = os.path.join(out_dir, out_name0)
    out_dir_new = os.path.join(out_path, 'outputpoints.txt')
    if os.path.exists(out_path) != True:
        os.mkdir(out_path)
    source_path = points_txt
    points = ' -def ' + source_path
    out = ' -out ' + out_path
    param = ' -tp ' + parameter_files
    cmd = transformix_bin + points + out + param
    print('Executing the following:')
    print(cmd)
    res = os.system(cmd)
    if res != 0:
        print('An error occured, check the parameters and associated files')
    return out_dir_new


def parse_transformix_output(file_path):
    trans_points = open(file_path, 'r')
    lines = trans_points.readlines()
    array = np.empty([len(lines), 3], dtype=int)
    for i in range(len(lines)):
        string = lines[i]
        match = re.findall(r'output_point\s=\s\[\s\d+\.\d+\s\d+\.\d+\s\d+\.\d+', string)
        if len(match) != 0:
            xyz = re.findall("\d+\.\d+", match[0])
            row = []
            for cood in xyz:
                asfloat = float(cood)
                asint = int(np.round(asfloat))
                row.append(asint)
            array[i, :] = row
    return array

# ------------------------------------------------------ Other ---------------------------------------------------------

def get_atlas_points(img, points):
    cell_img = np.zeros(img.shape, dtype=np.uint32)
    indices = np.zeros(points.shape, dtype=int)
    for i in range(len(points)):
        indices[i, :] = [int(np.round(points[i, 2])), int(np.round(points[i, 1])), int(np.round(points[i, 0]))]
    for i in range(len(indices)):
        if indices[i, 0] >= cell_img.shape[0]:
            indices[i, 0] = cell_img.shape[0] - 1
        if indices[i, 1] >= cell_img.shape[1]:
            indices[i, 1] = cell_img.shape[1] - 1
        if indices[i, 2] >= cell_img.shape[2]:
            indices[i, 2] = cell_img.shape[2] - 1
        cell_img[indices[i, 0], indices[i, 1], indices[i, 2]] += 1
    return cell_img


def affine_transform(parameters, COR, points):
    """
    Assumes input points in 2d array w/ r: x, y, z
    Assumes input points are np.ndarray
    """
    mat = np.array([
        [parameters[0], parameters[1], parameters[2]],
        [parameters[3], parameters[4], parameters[5]],
        [parameters[6], parameters[7], parameters[8]]])
    trans = np.array([parameters[9], parameters[10], parameters[11]])
    output = []
    for row in points:
        row = row.astype(float)
        row = row - COR
        output_row = mat @ row
        output_row = output_row + trans + COR
        output.append([output_row])
    output = np.concatenate(output)
    return output


# ------------------------------------------ 2D Mapping Related Functions ----------------------------------------------

# ---------------------------------------- Find approximate in situ volume ---------------------------------------------
# required for aligning images to the Allen Brain Atlas
# ensure files for a single brain are stored in a directory
def find_volume_coordinates(directory, out_directory, image_name_pattern=r'image_id-\d*.jpeg', scale_factor=8):
    """
    directory: directory in which to find rgb images
    out_directory: directory to which tissue mask images should be saved
    """
    find_gene_series_masks(directory, out_directory, image_name_pattern=image_name_pattern)
    volume_assembler = get_volume_assembler(out_directory, scale_factor=scale_factor)
    return volume_assembler


def get_volume_assembler(tissue_mask_directory, scale_factor=8, tissue_mask_pattern=r'\d*_tissue_mask.tif'):
    dictionary = get_volume_assembler_dict(tissue_mask_directory, scale_factor=scale_factor,
                                           tissue_mask_pattern=tissue_mask_pattern)
    volume_assembler = volume_assembler(dictionary)
    return volume_assembler


def get_volume_assembler_dict(tissue_mask_directory, scale_factor=8, tissue_mask_pattern=r'\d*_tissue_mask.tif'):
    # initialise the dictionary which will be filled with information about volume assembly (for volume_assembler object)
    transformation_dictionary = {'coordinates_list' : []}
    # find the files containing tissue masks
    tissue_mask_files = []
    tissue_mask_pattern = re.compile(tissue_mask_pattern)
    for file in sorted(os.listdir(tissue_mask_directory)):
        match = tissue_mask_pattern.match(file)
        if match is not None:
            tissue_mask_files.append(match[0])
    # find the important parameters
    for i in range(len(tissue_mask_files)):
        file = tissue_mask_files[i]
        file_path = os.path.join(tissue_mask_directory, file)
        index = str(re.findall(r'\d*', file)[0])
        mask_image = io.imread(file_path)
        mask_image = mask_image.astype(int)
        props = measure.regionprops(mask_image)[0]
        bounding_box = (np.round(props['bbox'][0]*scale_factor).astype(int),
                        np.round(props['bbox'][1]*scale_factor).astype(int),
                        np.round(props['bbox'][2]*scale_factor).astype(int),
                        np.round(props['bbox'][3]*scale_factor).astype(int))
        centroid = (np.round(props['local_centroid'][0]*scale_factor).astype(int),
                    np.round(props['local_centroid'][1]*scale_factor).astype(int) )
        transformation_dictionary['coordinates_list'].append({})
        transformation_dictionary['coordinates_list'][i]['index'] = index
        transformation_dictionary['coordinates_list'][i]['bounding_box'] = bounding_box
        transformation_dictionary['coordinates_list'][i]['local_centroid'] = centroid
    return transformation_dictionary


def get_volume_from_series(volume_assembler, series_pattern, directory, scale_factor=1):
    files = []
    pattern = re.compile(series_pattern)
    for file in os.listdir(directory):
        match = pattern.match(file)
        if match is not None:
            files.append(match[0])
    volume = volume_assembler.generate_volume(directory, files)
    return volume

# -------------------------------------------- volume_assembler Class ---------------------------------------------------
class VolumeAssembler:
    def __init__(self, volume_dictionary):
        """
        :param volume_dictionary: only requires coordinates list with the local centroid, bounding_box, indexes
        {'shape': tuple, 'coordinates_list': [
                                             {'index' : int, 'bounding_box' : tuple, 'centroid': , 'area':},
                                              'volume_coordinates' : ]}
        """
        self.assembly_info = volume_dictionary
        # self.shape = volume_dictionary['shape']
        self.image_info = volume_dictionary['coordinates_list']
        self.shape, self.centroid = self.calculate_positions()


    def calculate_positions(self):
        y_upper = [image['bounding_box'][2] - image['bounding_box'][0] - image['local_centroid'][0]
                   for image in self.image_info]
        y_lower = [image['local_centroid'][0] for image in self.image_info]
        x_upper = [image['bounding_box'][3] - image['bounding_box'][1] - image['local_centroid'][1]
                   for image in self.image_info]
        x_lower = [image['local_centroid'][1] for image in self.image_info]
        y_range = np.max(y_upper) + np.max(y_lower) + 1
        y_cent = np.max(y_lower)
        x_range = np.max(x_upper) + np.max(x_lower) + 1
        x_cent = np.max(x_lower)
        z_range = len(self.image_info)
        shape = (z_range, y_range, x_range)
        centroid = (y_cent, x_cent)
        self.assembly_info['shape'] = shape
        self.assembly_info['local_centroid'] = centroid
        for i in range(len(self.image_info)):
            row_max = centroid[0] + y_upper[i]
            row_min = centroid[0] - y_lower[i]
            col_max = centroid[1] + x_upper[i]
            col_min = centroid[1] - x_lower[i]
            self.assembly_info['coordinates_list'][i]['volume_coordinates'] = (row_min, col_min, row_max, col_max)
        return shape, centroid


    def generate_volume(self, directory, out_name, image_file_pattern=r'\d*_pc1_greyscale.tif',
                        verbose=False, scale_factor=1):
        # check that each index has a corresponding file
        clean_list = []
        image_file_pattern_re = re.compile(image_file_pattern)
        for file in sorted(os.listdir(directory)):
            match = image_file_pattern_re.match(file)
            if match is not None:
                clean_list.append(match[0])
        # generate an empty volume with the correct shape
        if verbose:
            print('initialising an array of shape:', self.shape)
        volume = np.zeros(self.shape, dtype=np.uint8)
        for i in range(len(self.image_info)):
            z = int(self.image_info[i]['index'])
            if verbose:
                print('adding the data at z = {}'.format(z))
            # get the bounding box from each image
            image_file_0 = image_file_pattern[:image_file_pattern.find('\d*')]
            image_file_1 = image_file_pattern[image_file_pattern.find('\d*')+3:]
            image_file = image_file_0+str(z)+image_file_1
            image_file = os.path.join(directory, image_file)
            image = io.imread(image_file)
            min_row = np.round(self.image_info[i]['bounding_box'][0]*scale_factor).astype(int)
            min_col = np.round(self.image_info[i]['bounding_box'][1]*scale_factor).astype(int)
            max_row = np.round(self.image_info[i]['bounding_box'][2]*scale_factor).astype(int)
            max_col = np.round(self.image_info[i]['bounding_box'][3]*scale_factor).astype(int)
            image = image[min_row:max_row, min_col:max_col]
            # place the bounded image into the correct location in the volume
            min_row = np.round(self.image_info[i]['volume_coordinates'][0]*scale_factor).astype(int)
            min_col = np.round(self.image_info[i]['volume_coordinates'][1]*scale_factor).astype(int)
            max_row = np.round(self.image_info[i]['volume_coordinates'][2]*scale_factor).astype(int)
            max_col = np.round(self.image_info[i]['volume_coordinates'][3]*scale_factor).astype(int)
            if verbose:
                print('Data is being added at indices [{}:{}, {}:{}]'.format(min_row, max_row, min_col, max_col))
            volume[z, min_row:max_row, min_col:max_col] = image
        io.imsave(os.path.join(directory, out_name), volume, bigtiff=True)
        if verbose:
            print('An image was saved at:')
            print(os.path.join(directory, out_name))
        return image


# Some Example Code --------------------------------------------
# out_name = 'expression_quantification'
# out_directory = '/Users/amcg0011/Data/in_situ_data/entrez_id_12064_Bdnf/plane_of_section-1/age_id-15_id-79587720'
# out_directory = os.path.join(out_directory, out_name)
# save_resized_from_directory(out_directory, image_pattern=r'\d*_pc1_greyscale.tif', scale=0.125)
#
# FIND 1/8th ORIGINAL SIZE VOLUME (scale factor refers to how to scale coordinates from the tissue masks, which are 1/8)
# transform_dict = get_volume_assembler_dict(out_directory, scale_factor=1, tissue_mask_pattern=r'\d*_tissue_mask.tif') 
# volume_assembler = volume_assembler(transform_dict)
# volume = volume_assembler.generate_volume(out_directory, 'age_id-15_id-79587720_pc1.tif', 
#                                           image_file_pattern=r'\d*_pc1_greyscale_scale-0.125.tif', verbose=True)
