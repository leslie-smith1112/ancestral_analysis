import os


fin_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_with_exome.tsv'
fin_file = open(fin_path, 'a')

log_file = open('/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_to_exome_custom_log.txt', 'a')

c_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/ChromFiles/cosmic_chrom_X.tsv'
e_path = '/blue/kgraim/lucaspereira/ancestry_project/exomeAF_split_files/slice_chrom_X.tsv'
i_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/MappedIndices/chrom_X_indices.tsv'

c_file = open(c_path, 'r')
c_file.readline()
c_line = c_file.readline()

e_file = open(e_path, 'r')
e_obj = e_file.readlines()

i_file = open(i_path, 'r')
i_file.readline()
i_line = i_file.readline()

while c_line and i_line:
    ind = i_line[:i_line.find('\n')]
    ind = int(float(ind)) + 1

    log_string = str(ind) + ' ' + 'X\n'
    log_file.write(log_string)
    log_file.flush()
    os.fsync(log_file.fileno())

    exome_string = e_obj[ind]

    log_file.write('exome String: ')
    log_file.write(exome_string)
    log_file.write('\n')
    log_file.flush()
    os.fsync(log_file.fileno())
        
    newline = c_line[:c_line.find('\n')] + '\t' + exome_string

    fin_file.write(newline)

    c_line = c_file.readline()
    i_line = i_file.readline()

c_file.close()
e_file.close()
i_file.close()

log_file.close()

