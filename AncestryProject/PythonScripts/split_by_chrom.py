import pandas as pd

super_frame = pd.read_csv(
    '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_breast_original_locs.tsv', sep='\t')

frames = [y for x, y in super_frame.groupby('CHROM', as_index=False)]

for frame in frames:
    chrom = str(int(frame['CHROM'].values[0]))
    path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/ChromFiles/cosmic_chrom_' + chrom + '.tsv'
    frame.to_csv(path, sep='\t', index=False)

