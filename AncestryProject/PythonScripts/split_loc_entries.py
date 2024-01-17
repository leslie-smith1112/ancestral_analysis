import pandas as pd

frame = pd.read_csv('/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_breast_parsed.tsv', sep='\t')

no_locs = frame[frame['START_POS'].isna()]
locs = frame.dropna(subset=['START_POS'])

locs['START_POS'] = locs['START_POS'].astype('Int64')
locs['END_POS'] = locs['END_POS'].astype('Int64')
locs['MAPPING_POS'] = locs['MAPPING_POS'].astype('Int64')

no_locs.drop(columns=['START_POS', 'END_POS', 'MAPPING_POS'], inplace=True)

no_locs.to_csv('/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_breast_no_locs.tsv', sep='\t', index=False)
locs.to_csv('/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_breast_original_locs.tsv', sep='\t', index=False)

