import pandas as pd

for chrom in range(1, 23):
    read_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/ChromFiles/cosmic_chrom_' + \
                str(chrom) + '.tsv'

    write_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/chrom_' + str(chrom) + '_positions.tsv'

    positions = pd.read_csv(read_path, sep='\t', usecols=['MAPPING_POS'])

    positions.to_csv(write_path, sep='\t', index=False)

