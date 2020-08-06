import os
from msbrainpy.data import download_gene_images, get_image_addresses
from msbrainpy._test_helpers import get_temp_directory, remove_temp_data, \
    timeout_function
from skimage.io import imread


grin3a = 242443 # entrez id for gene
timeout = 120 # number of seconds to allow download to continue


def test_download_gene_images(gene=grin3a, timeout=timeout):
    """
    Test that both the download gene image function and the function
    to find directories and images on file are working.

    Parameters
    ----------
    gene: int
        entrez id of the gene for which to collect the data
    timeout: int
        number of seconds to wait before cancelling the operation
    """
    temp = get_temp_directory()
    _download_images_wrapper(temp, gene, timeout)
    _get_image_test(temp)
    remove_temp_data(temp)


def _download_images_wrapper(temp, gene, timeout):
    """
    A wrapper function to feed download_gene_images and the key word arguments
    to the function that cancels the operation after the specified time
    """
    kw = {'entrez_id' : gene, 'directory_path' : temp}
    timeout_function(download_gene_images, kw, timeout)


def _get_image_test(temp):
    """
    This tests get_image_adresses to ensure that the appropriate directory
    structure can be found and that at least one image has been downloaded
    """
    data_tree = get_image_addresses(temp)
    a_gene = list(data_tree.keys())[0]
    a_plane = list(data_tree[a_gene].keys())[0]
    a_brain = list(data_tree[a_gene][a_plane].keys())[0]
    brain_directory = os.path.join(temp, a_gene, a_plane, a_brain)
    a_slice = data_tree[a_gene][a_plane][a_brain][0]
    pf = 'Download Allen ISH images: '
    m = "the path to the first brain directory does not exist"
    ck = os.path.exists(brain_directory)
    assert ck, pf + m
    m = f'no image was found at {a_slice}'
    ck = len(imread(a_slice).shape) == 3
    assert ck, pf + ck


test_download_gene_images(gene=grin3a, timeout=timeout)

# def func(line):
    # with open('a.txt', 'w') as f:
        # f.write(line)
    # time.sleep(10)
    # with open('b.txt', 'w') as f:
        # f.write(line)
#_timeout_function(func, {'line': 'a line'}, 2)


