from CNV.cnv import idat_snp_metrics, call_cnvs
import argparse

parser = argparse.ArgumentParser(description='Arguments for Running CNV Pipeline.')    
parser.add_argument('--metrics', type=str, default='Nope.', help='Path to SNP metrics file (output from run_snp_metrics_pipeline.py')
parser.add_argument('--out_path', type=str, default='Nope.', help='Path to output CNV report.')
parser.add_argument('--intervals', type=str, default="glist_hg38_intervals.csv", help='Gene or other feature intervals to analyze. Header is [NAME,CHR,START,STOP], one line per itnervals.Autosomes only.')
parser.add_argument('--min_variants', type=int, default=10, help='Minimum number of variants to run the CNV algorithm on per gene.')
parser.add_argument('--kb_window', type=int, default=100, help='Kilobase window around each interval, a value of 100 would mean +/- 100kb.')

args = parser.parse_args()

metrics = args.metrics
out_path = args.out_path
intervals = args.intervals
min_variants = args.min_variants
kb_window = args.kb_window

call_cnvs(snp_metrics_file=metrics, out_path=out_path, intervals_file=intervals, min_variants=min_variants, kb_window=kb_window)