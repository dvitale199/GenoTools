from CNV.cnv import idat_snp_metrics, call_cnvs

parser = argparse.ArgumentParser(description='Arguments for Running CNV Pipeline.')    
parser.add_argument('--idat_path', type=str, default='Nope.', help='Path to directory containing idats. Assumes idats are stored by SentrixBarcode_A i.e. /path/to/123456789, which contains all idats under that barcode such as 123456789_R01C01_Red.idat, etc.')
parser.add_argument('--out_path', type=str, default='Nope.', help='Path to directory for outputs. right now, will output a large number of files. This will change in a future release.')
parser.add_argument('--bpm', type=str, default='Nope.', help='Path to Illumina manifest .bpm file.')
parser.add_argument('--bpm_csv', type=str, default='Nope.', help='Path to Illumina manifest .csv file.')
parser.add_argument('--egt', type=str, default='Nope.', help='Path to Illumina clusterfile .egt.')
parser.add_argument('--ref_fasta', type=str, default='Nope.', help='Path to Reference FASTA File.')
parser.add_argument('--iaap', type=str, default='Nope.', help='Path to Illumina iaap executable.')
# parser.add_argument('--intervals', type=str, default="glist_hg38_intervals.csv", help='Gene or other feature intervals to analyze. Header is [NAME,CHR,START,STOP], one line per itnervals.Autosomes only.')
# parser.add_argument('--min_variants', type=int, default=10, help='Minimum number of variants to run the CNV algorithm on per gene.')
# parser.add_argument('--kb_window', type=int, default=100, help='Kilobase window around each interval, a value of 100 would mean +/- 100kb.')

idat_path = args.idat_path
bpm = args.bpm
bpm_csv = args.bpm_csv
egt = args.egt
ref_fasta = args.ref_fasta
iaap = args.iaap
out_path = args.out_path
# intervals = args.intervals
# min_variants = args.min_variants
# kb_window = args.kb_window

# idat_path = '/data/vitaled2/cnv_test/206046180074'
# test_out = '/data/vitaled2/cnv_test'
snp_metrics_out = idat_snp_metrics(idat_path=idat_path, bpm=bpm, bpm_csv=bpm_csv, egt=egt, ref_fasta=ref_fasta, out_path=out_path, iaap=iaap)

for mfile in snp_metrics_out:
    out_prefix = mfile.replace('.csv','').replace('snp_metrics_', '')
    call_cnvs(snp_metrics_file=mfile, out_path=out_prefix, intervals_file=intervals, min_variants=10, kb_window=100)