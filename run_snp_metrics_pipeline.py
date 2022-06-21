from CNV.cnv import idat_snp_metrics
import argparse

parser = argparse.ArgumentParser(description='Arguments for Running CNV Pipeline.')    
parser.add_argument('--idat_path', type=str, default='Nope.', help='Path to directory containing idats. Assumes idats are stored by SentrixBarcode_A i.e. /path/to/123456789, which contains all idats under that barcode such as 123456789_R01C01_Red.idat, etc.')
parser.add_argument('--out_path', type=str, default='Nope.', help='Path to directory for outputs. right now, will output a large number of files. This will change in a future release.')
parser.add_argument('--bpm', type=str, default='Nope.', help='Path to Illumina manifest .bpm file.')
parser.add_argument('--bpm_csv', type=str, default='Nope.', help='Path to Illumina manifest .csv file.')
parser.add_argument('--egt', type=str, default='Nope.', help='Path to Illumina clusterfile .egt.')
parser.add_argument('--ref_fasta', type=str, default='Nope.', help='Path to Reference FASTA File.')
parser.add_argument('--iaap', type=str, default='Nope.', help='Path to Illumina iaap executable.')
parser.add_argument('--bcftools_plugins_path', type=str, default='/data/vitaled2/bin', help='Path to bcftools plugins for gtc2vcf')
args = parser.parse_args()

idat_path = args.idat_path
bpm = args.bpm
bpm_csv = args.bpm_csv
egt = args.egt
ref_fasta = args.ref_fasta
iaap = args.iaap
out_path = args.out_path
bcftools_plugins_path = args.bcftools_plugins_path

snp_metrics_out = idat_snp_metrics(idat_path, bpm, bpm_csv, egt, ref_fasta, out_path, iaap=iaap, bcftools_plugins_path=bcftools_plugins_path)
