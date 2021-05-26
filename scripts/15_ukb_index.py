# Indexing UK Biobank BGEN files

# Libraries
import os
import hail as hl

dir = '/ludc/Home/daniel_c/dva/files/ukb_index'

# BGEN files to index
bgen_files = [
    '/ludc/Raw_Data_Archive/UKBB/imp/ukb_imp_chr{}_v3.bgen'.format(x) for x in range(1,23)
    ]

# Mapping each BGEN file to its index
file_map = {'/ludc/Raw_Data_Archive/UKBB/imp/ukb_imp_chr{}_v3.bgen'.format(x):
            dir + '/chr{}.idx2'.format(x) for x in range(1,23)}

# Chromosome (aka contig) formats
contigs = {'0{}'.format(x):str(x) for x in range(1, 10)}

# Indexing
hl.index_bgen(bgen_files,
              contig_recoding = contigs,
              reference_genome = 'GRCh37',
              index_file_map = file_map,
              skip_invalid_loci = True)
