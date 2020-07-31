import os
from msbrainpy.brainseries import InSituSeries
from msbrainpy.atlas import wrangle # once I have general wrapper functions for the functions in this 

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
brain_directory = os.path.join(gene_directory, 'plane_of_section-1', a_gene, a_brain)
print(f'Your brain directory: {brain_directory}')
# Your brain directory: /Users/amcg0011/Data/InSituData/plane_of_section-1/entrez_id_14812_Grin2b/age_id-15_id-74988710