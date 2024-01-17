import pandas as pd

frame = pd.read_csv('/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_breast_pos_unfiltered.tsv', sep='\t')

frame['CHROM'] = frame[' MUTATION_GENOME_POSITION'].str[:2].str.replace(':', '')

frame['START_POS'] = frame[' MUTATION_GENOME_POSITION'].str.extract(r'(?<=\:)(.*?)(?=\-)')
frame['END_POS'] = frame[' MUTATION_GENOME_POSITION'].str.extract(r'\W(\d+)$')

frame['MAPPING_POS'] = (frame['START_POS'].astype('Int64') + frame['END_POS'].astype('Int64')) / 2
frame['MAPPING_POS'] = frame['MAPPING_POS'].astype('Int64')

frame.drop([' MUTATION_GENOME_POSITION'], axis=1, inplace=True)

frame.sort_values(by=['CHROM'], inplace=True)
frame.to_csv('/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_breast_parsed.tsv', sep='\t', index=False)
