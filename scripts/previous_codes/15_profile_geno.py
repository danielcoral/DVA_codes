# Hail script for scores in UKBB

import os
import glob
import pandas as pd
import hail as hl

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

# Concordant and discordant SNPs
mix = pd.read_table('../files/mix.tsv', usecols = ['chr', 'hg19_pos', 'ea', 'disc'],
                    dtype = {'chr':object, 'hg19_pos':object, 'ea':object, 'disc':int})

conc_snps = hl.Table.from_pandas(mix[mix.disc == 0])

disc_snps = hl.Table.from_pandas(mix[mix.disc == 1])

profiles = [conc_snps, disc_snps]

prs_res = []

for profile in profiles:
    # Preparing loci to be extracted
    profile_tab = profile.annotate(
        locus = hl.struct(contig = profile.chr, position = hl.int32(profile.hg19_pos)),
        interval = hl.interval(hl.struct(contig = profile.chr,
                                         position = hl.int32(profile.hg19_pos)),
                               hl.struct(contig = profile.chr,
                                         position = hl.int32(profile.hg19_pos)),
                               includes_end=True)
    ).key_by('locus')
    # Extracting SNPs of interest
    mtprofile = hl.filter_intervals(mt, profile_tab.interval.collect()).repartition(500)
    # Adding information of SNPs to MatrixTable
    mtprofile = mtprofile.annotate_rows(**profile_tab[mtprofile.locus])
    # Harmonizing
    flip = hl.case().when(
        mtprofile.ea == mtprofile.alleles[0], True
    ).when(
        mtprofile.ea == mtprofile.alleles[1], False
    ).or_missing()
    mtprofile = mtprofile.annotate_rows(flip = flip)
    # Calculating scores
    res = mtprofile.annotate_cols(
        prs = hl.agg.sum(hl.if_else(mtprofile.flip, 2 - mtprofile.GT.n_alt_alleles(),
                                    mtprofile.GT.n_alt_alleles())))
    prs_res.append(res)

# Saving results
filenames = ['../files/{}_prs.tsv'.format(x) for x in ['conc', 'disc']]

for prs, filename in zip(prs_res, filenames):
    prs.key_cols_by().cols().export(filename)

## Removing log file
logfile = glob.glob('*.log')
os.remove(logfile[0])
