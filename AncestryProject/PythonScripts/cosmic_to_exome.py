import os


fin_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_with_exome.tsv'
fin_file = open(fin_path, 'w')

temp_e_file = open('/blue/kgraim/lucaspereira/ancestry_project/exomeAF_split_files/slice_chrom_1.tsv', 'r')
temp_c_file = open('/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/ChromFiles/cosmic_chrom_1.tsv')

c_init = temp_c_file.readline()
init = c_init[:c_init.find('\n')] + '\t' + temp_e_file.readline()
fin_file.write(init)

temp_e_file.close()
temp_c_file.close()

log_file = open('/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/cosmic_to_exome_custom_log.txt', 'w')

for chrom in range(1, 23):
    c_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/ChromFiles/cosmic_chrom_' + str(chrom) + \
             '.tsv'
    e_path = '/blue/kgraim/lucaspereira/ancestry_project/exomeAF_split_files/slice_chrom_' + str(chrom) + '.tsv'
    i_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/MappedIndices/chrom_' + str(chrom) + \
             '_indices.tsv'

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

        log_string = str(ind) + ' ' + str(chrom) + '\n'
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

