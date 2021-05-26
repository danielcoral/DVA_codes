# Hail script for scores in UKBB

import os
import glob
import pandas as pd
import hail as hl

# Concordant and discordant SNPs
mix = pd.read_table('/ludc/Home/daniel_c/dva/files/mix.tsv',
                    header = 0, usecols = [0, 1, 3],
                    names = ['chrom', 'hg19_pos', 'ea'],
                    dtype = {'chrom':object, 'hg19_pos':int, 'ea':object})

# UKB Sample file
ukb_sf = '/ludc/Home/daniel_c/dva/files/ukb_57232/ukb57232_imp_chr1_v3_s487266.sample'

# Extracting information per chromosome
for ch in mix['chrom'].unique():
    print('Analyzing chromosome ' + ch)
    chromdat = mix[mix['chrom'] == ch]
    posranges = [str(p) + '-' + str(p + 1) for p in chromdat.hg19_pos]
    loci = [ch + ':' + p for p in  posranges]
    ploci = [hl.parse_locus_interval(l, reference_genome='GRCh37') for l in loci]
    # BGEN file
    bgen_file = '/ludc/Raw_Data_Archive/UKBB/imp/ukb_imp_chr{}_v3.bgen'.format(ch)
    # Mapping BGEN file to index file
    file_map = {'/ludc/Raw_Data_Archive/UKBB/imp/ukb_imp_chr{}_v3.bgen'.format(ch):
                '/ludc/Home/daniel_c/dva/files/ukb_index/chr{}.idx2'.format(ch)}
    # Creating MatrixTable
    mt = hl.import_bgen(bgen_file, entry_fields = ['GT'], sample_file = ukb_sf,
                        index_file_map = file_map, _row_fields = ['rsid'])
    # Extracting SNPs of interest
    mt_f = hl.filter_intervals(mt, ploci)
    mt_f = hl.variant_qc(mt_f)
    chromdat['chrompos'] = chromdat['chrom'] + ':' + chromdat['hg19_pos'].astype(str)
    chromdat_hl = hl.Table.from_pandas(chromdat)
    chromdat_hl = chromdat_hl.annotate(
        locus = hl.parse_locus(chromdat_hl.chrompos,
                               reference_genome='GRCh37')
        )
    chromdat_hl = chromdat_hl.key_by('locus')
    mt_f = mt_f.annotate_rows(**chromdat_hl[mt_f.locus])
    flip = hl.case().when(mt_f.ea == mt_f.alleles[0], True).when(
        mt_f.ea == mt_f.alleles[1], False).or_missing()
    mt_f = mt_f.annotate_rows(flip = flip)
    mt_f = mt_f.annotate_rows(
        prior = 2 * hl.if_else(mt_f.flip, mt_f.variant_qc.AF[0],
                               mt_f.variant_qc.AF[1]))
    mt_f = mt_f.select_entries(
        G = hl.coalesce(hl.if_else(mt_f.flip,
                                   2 - mt_f.GT.n_alt_alleles(),
                                   mt_f.GT.n_alt_alleles()),
                        mt_f.prior)
        )
    ## Exporting result
    output = '/ludc/Home/daniel_c/dva/files/ukbgeno/chrom{}.vcf.bgz'.format(ch)
    hl.export_vcf(mt_f, output)

## Removing log file
logfile = glob.glob('*.log')
os.remove(logfile[0])
