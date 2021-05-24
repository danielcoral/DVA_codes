## Building PRS for each trait identified

import argparse
import pandas as pd
import hail as hl
import os
import glob

parser = argparse.ArgumentParser(description = 'Multi PRS UKB')
parser.add_argument('ss', help = 'Summary stats of traits of interest')
args = parser.parse_args()

coltypes = {'trait':object, 'rsid':object, 'chr':object, 'hg19_pos':object, 'ea':object}

ss_tab = pd.read_table(args.ss,
                      usecols = ['trait', 'chr', 'hg19_pos', 'ea'],
                      dtype = coltypes)

traitname = ss_tab.trait.unique()[0]

ss_tab = hl.Table.from_pandas(ss_tab)

ss_tab = ss_tab.annotate(
    locus = hl.struct(contig = ss_tab.chr, position = hl.int32(ss_tab.hg19_pos)),
    interval = hl.interval(
        hl.struct(contig = ss_tab.chr, position = hl.int32(ss_tab.hg19_pos)),
        hl.struct(contig = ss_tab.chr, position = hl.int32(ss_tab.hg19_pos)),
        includes_end=True
        )
    ).key_by('locus')

# BGEN files
bgen_files = [
    '/gpfs/gpfs0/Raw_Data_Archive/UKBB/imp/ukb_imp_chr{}_v3.bgen'.format(x) for x in range(1,23)
    ]

# Mapping each BGEN file to its index
file_map = {'/gpfs/gpfs0/Raw_Data_Archive/UKBB/imp/ukb_imp_chr{}_v3.bgen'.format(x):
            '../files/ukb_index/chr{}.idx2'.format(x) for x in range(1,23)}

# Creating MatrixTable
mt = hl.import_bgen(bgen_files,
                    entry_fields = ['GT'],
                    sample_file = '../files/ukb_57232/ukb57232_imp_chr1_v3_s487296.sample',
                    index_file_map = file_map)

# Extracting SNPs of interest
mt = hl.filter_intervals(mt, ss_tab.interval.collect()).repartition(500)

# Adding information of SNPs to MatrixTable
mt = mt.annotate_rows(**ss_tab[mt.locus])

# Harmonizing
flip = hl.case().when(
    mt.ea == mt.alleles[0], True
).when(
    mt.ea == mt.alleles[1], False
).or_missing()

mt = mt.annotate_rows(flip = flip)

# Calculating scores
mt = mt.annotate_cols(
    prs = hl.agg.sum(hl.if_else(mt.flip,
                                2 - mt.GT.n_alt_alleles(),
                                mt.GT.n_alt_alleles()))
    )

filename = '../files/multiprs/' + traitname + '_prs.tsv'

mt.key_cols_by().cols().export(filename)

## Removing log file
for x in glob.glob('*.log'):
    os.remove(x)
