import pandas as pd


read_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/ChromFiles/cosmic_chrom_X.tsv'

write_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/chrom_X_positions.tsv'

positions = pd.read_csv(read_path, sep='\t', usecols=['MAPPING_POS'])

positions.to_csv(write_path, sep='\t', index=False)

