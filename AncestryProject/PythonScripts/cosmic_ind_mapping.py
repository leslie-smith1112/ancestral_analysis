import numpy as np


def get_arrays(chrom):
    cosmic_path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/PositionArrays/chrom_' + chrom + '_positions.tsv'
    exome_path = '/blue/kgraim/lucaspereira/ancestry_project/ExomePositions/chrom_' + chrom + '_positions.tsv'

    cosmic_f = open(cosmic_path, 'r')
    exome_f = open(exome_path, 'r')

    cosmic_l = []
    exome_l = []

    exome_f.readline()
    e_line = exome_f.readline()
    while e_line:
        e_line = e_line[e_line.find('\t') + 1:e_line.find('\n')] 
        exome_l.append(float(e_line))
        e_line = exome_f.readline()

    cosmic_f.readline()
    c_line = cosmic_f.readline()
    while c_line:
        c_line = c_line[:c_line.find('\n')]
        cosmic_l.append(float(c_line))
        c_line = cosmic_f.readline()

    cosmic_l = np.array(cosmic_l)
    exome_l = np.array(exome_l)

    return cosmic_l, exome_l


def get_inds(t, chrom):
    cosmic_a = t[0]
    exome_a = t[1]
    
    a_list = np.array_split(cosmic_a, 150)
    path = '/blue/kgraim/lucaspereira/ancestry_project/CosmicBreastMapping/MappedIndices/chrom_' + chrom + '_indices.tsv'
    _file = open(path, 'w')
    _file.write('exome indices for cosmic file\n')
    _file.close()

    for i in a_list:
        ind = np.abs(exome_a - i[:, None]).argmin(1)
        ind  = ind.astype(int)
        with open(path, 'a') as f:
            np.savetxt(f, ind)


if __name__ == '__main__':
    for i in range(1, 23):
        tup = get_arrays(str(i))
        get_inds(tup, str(i))

