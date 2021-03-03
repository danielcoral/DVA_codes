# Indexing UK Biobank BGEN files

# Libraries
import os
import hail as hl

dir = '../files/ukb_index'

# BGEN files to index
bgen_files = [
    '/gpfs/gpfs0/Raw_Data_Archive/UKBB/imp/ukb_imp_chr{}_v3.bgen'.format(x) for x in range(1,23)
    ]

# Mapping each BGEN file to its index
file_map = {'/gpfs/gpfs0/Raw_Data_Archive/UKBB/imp/ukb_imp_chr{}_v3.bgen'.format(x):
            args.dir + '/chr{}.idx2'.format(x) for x in range(1,23)}

# Chromosome (aka contig) formats
contigs = {'0{}'.format(x):str(x) for x in range(1, 10)}

# Indexing
hl.index_bgen(bgen_files,
              contig_recoding = contigs,
              reference_genome = None,
              index_file_map = file_map)
