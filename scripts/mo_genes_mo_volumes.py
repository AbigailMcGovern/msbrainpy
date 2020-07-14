import os
from msbrainpy.atlas import query, wrangle
from msbrainpy.quantify import processing
from msbrainpy.io import save_resized_from_directory
from msbrainpy import map

# where will the data be stored
gene_directory = '/Users/amcg0011/Data/InSituData'

# get in situ data for the gene Grin3a (ID: 242443)
#   WHY: developmentally important NDMA-R component thought to be involved in delaying synaptic maturation
query.download_gene_images(entrez_id=242443, directory_path=gene_directory)

# get the in situ data for the gene Mc4r (ID: 17202)
#   WHY: very different (more selecive) expression pattern to Grin3a and BDNF expression (quite ubiquitous)
query.download_gene_images(entrez_id=17202, directory_path=gene_directory)

# get the in situ data for gene Grin2b (ID: 14812)
#   WHY: have aquired an equivalent image series from gene paint
query.download_gene_images(entrez_id=14812, directory_path=gene_directory)

# get the image paths for genes avaliable in the gene data directory
gene_data_tree = wrangle.get_image_addresses(gene_directory, age_id=None)

# get segmented images, PC1 images, and tissue masks for one of the samples in this gene data set
genes = list(gene_data_tree.keys())
# Out[12]:
# ['entrez_id_12064_Bdnf',
#  'entrez_id_14812_Grin2b',
#  'entrez_id_17202_Mc4r',
#  'entrez_id_242443_Grin3a']

# ------------------------------------------------- Grin2b -------------------------------------------------------------
grin2b = genes[1]
grin2b_coronal_brains = list(gene_data_tree[grin2b]['plane_of_section-1'].keys())
grin2b_brain_directory = os.path.join(gene_directory, grin2b, 'plane_of_section-1', grin2b_coronal_brains[0])
os.path.exists(grin2b_brain_directory)
# True
grin2b_out_directory = os.path.join(grin2b_brain_directory, 'pipeline_work')
os.path.exists(grin2b_brain_directory)
# False, do_gene_series will make this directory
processing.do_gene_series(directory=grin2b_brain_directory, out_directory=grin2b_out_directory, find_tissue_mask=True)

# get a 1/8 volume ISH reconstruction for the gene series
save_resized_from_directory(grin2b_out_directory, image_pattern=r'\d*_pc1_greyscale.tif', scale=0.125)
volume_assembler = map.get_volume_assembler(grin2b_out_directory, scale_factor=1,
                                            tissue_mask_pattern=r'\d*_tissue_mask.tif')
name = grin2b_coronal_brains[0] + '_pc1.tif'
volume_assembler.generate_volume(grin2b_out_directory, name, scale_factor=1,
                                 image_file_pattern=r'\d*_pc1_greyscale_scale-0.125.tif')

# The GenePaint coronal Grin2b brain
gp_grin2b_brain_directory = os.path.join(gene_directory, grin2b, 'plane_of_section-1', 'GP-P56-HB605')
os.path.exists(gp_grin2b_brain_directory)
gp_grin2b_out_directory = os.path.join(gp_grin2b_brain_directory, 'pipeline_work')
processing.do_gene_series(directory=gp_grin2b_brain_directory, out_directory=gp_grin2b_out_directory,
                          find_tissue_mask=True, image_name_pattern=r'HB\d*_\d*B.jpg')

save_resized_from_directory(gp_grin2b_out_directory, image_pattern=r'\d*_pc1_greyscale.tif', scale=0.125)
gp_grin2b_volume_assembler = map.get_volume_assembler(gp_grin2b_out_directory, scale_factor=1,
                                                      tissue_mask_pattern=r'\d*_tissue_mask.tif')
name_0 = 'GP-P56-HB605_pc1.tif'
gp_grin2b_volume_assembler.generate_volume(gp_grin2b_out_directory, name_0, scale_factor=1,
                                           image_file_pattern=r'\d*_pc1_greyscale_scale-0.125.tif')