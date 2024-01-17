This is pretty messy, a lot of the filepaths are hardcoded etc. so sorry in advance 

The general idea is to:
1. filter for the snps you want to map (filter_init_csv.py), this will need to be custom-made for whatever database you're working with though
2. split up the target file into a file for each chromosome (split_by_chrom.py)
       -  if you want to use the gnomAD files they should already be split up at /blue/kgraim/lucaspereira/ancestry_project/enhanced_exomeAF_split_files
3. Get the locations of each SNP (get_chrom_locations.py)
      - look at /blue/kgraim/lucaspereira/ancestry_project/ExomePositions for the gnomAD positions
4. Find the line in the gnomAD exome position file that corresponds to the closest SNP match (cosmic_ind_mapping.py)
5. Use these matched line numbers to make a file with ancestral allele frequencies attatched to the SNPs in your target database (cosmic_to_exome.py)

If you have any questions lmk!
