import os
from msbrainpy.atlas import query, wrangle

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
#  'entrez_id_242443_Grin3a'
