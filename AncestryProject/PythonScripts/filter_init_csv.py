import pandas as pd

frame = pd.read_csv('/blue/kgraim/lucaspereira/ancestry_project/Cosmic_Breast_Mutations.csv',
                    sep=',')

test = frame.columns

frame.drop([' ACCESSION_NUMBER', ' GENE_CDS_LENGTH', ' HGNC_ID', ' SAMPLE_NAME', ' ID_SAMPLE', ' ID_TUMOUR',
            ' SITE_SUBTYPE_1', ' SITE_SUBTYPE_2', ' SITE_SUBTYPE_3', ' GENOME_WIDE_SCREEN', ' LEGACY_MUTATION_ID',
            ' MUTATION_CDS', ' MUTATION_AA', ' MUTATION_ZYGOSITY', ' LOH', ' GRCH', ' MUTATION_STRAND',
            ' RESISTANCE_MUTATION',
            ' PUBMED_PMID', ' ID_STUDY', ' SAMPLE_TYPE', ' TUMOUR_ORIGIN', ' HGVSP', ' HGVSC', ' HGVSG'],
           axis=1, inplace=True)


frame.to_csv('/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_breast_pos_unfiltered.tsv',
             sep='\t', index=False)
