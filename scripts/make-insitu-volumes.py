import os
from msbrainpy.brainseries import InSituSeries
from msbrainpy.atlas import wrangle # this takes a while because AllenSDK structure tree is imported ... hmmm 
from dask.distributed import Client
from tifffile import TiffWriter

client = Client() # pretty sure this is unnecessary
# The dask docs indicate that a client will be created and closed
# where required, but was unsure. 

gene_directory = gene_directory = '/Users/amcg0011/Data/InSituData'
gene_data_tree = wrangle.get_image_addresses(gene_directory, age_id=None)
genes = list(gene_data_tree.keys())
a_gene = genes[1]
print(f'results from gene lucky dip: {a_gene}')
# results from gene lucky dip: entrez_id_14812_Grin2b
coronal_brains = list(gene_data_tree[a_gene]['plane_of_section-1'].keys())
a_brain = coronal_brains[0]
print(f'Your brain: {a_brain}')
# Your brain: age_id-15_id-74988710
brain_directory = os.path.join(gene_directory, a_gene, 'plane_of_section-1', a_brain)
print(f'Your brain directory: {brain_directory}')
# Your brain directory: /Users/amcg0011/Data/InSituData/plane_of_section-1/entrez_id_14812_Grin2b/age_id-15_id-74988710
print(os.path.exists(brain_directory))
# True
name_pattern = r'image_id-\d*.jpeg'
# the amoount of work required to get a brain directory is too much. 
# This needs to be simplified. Should be as simple as:
# data_tree = get_data_tree(base_directory) 
# data_tree.genes[0].planes[0].brains[0] --> absolute dir path
# or somthing equivalently easy
my_brain = InSituSeries(brain_directory, name_pattern)
target = my_brain._target_volume
save_at = os.path.join(gene_directory, 'InSituVolumes', 'Grin2b-age_id-15_id-74988710.tif')
with TiffWriter(save_at) as tiff:
    for i in range(target.shape[0]):
        tiff.save(target[i, :, :])

a_gene = genes[2]
print(f'results from gene lucky dip: {a_gene}')
# results from gene lucky dip: entrez_id_17202_Mc4r
coronal_brains = list(gene_data_tree[a_gene]['plane_of_section-1'].keys())
a_brain = coronal_brains[0]
print(f'Your brain: {a_brain}')
# Your brain: age_id-15_id-79556630
brain_directory = os.path.join(gene_directory, a_gene, 'plane_of_section-1', a_brain)
print(f'Your brain directory: {brain_directory}')
# Your brain directory: /Users/amcg0011/Data/InSituData/entrez_id_17202_Mc4r/plane_of_section-1/age_id-15_id-79556630
print(os.path.exists(brain_directory))
my_brain = InSituSeries(brain_directory, name_pattern)
target = my_brain._target_volume
save_at = os.path.join(gene_directory, 'InSituVolumes', 'Mcr4-age_id-15_id-74988710.tif')
with TiffWriter(save_at) as tiff:
    for i in range(target.shape[0]):
        tiff.save(target[i, :, :])
client.close()