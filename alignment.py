import subprocess
import pandas as pd

def create_index(fasta, star_index_output, processors):
	 p = str(processors)

	 subprocess.run(['STAR', '--runMode', 'genomeGenerate', '--runThreadN', p,
		 	 '--genomeDir', star_index_output, '--genomeFastaFiles', fasta,
			 '--genomeSAindexNbases', '6'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, universal_newlines=True)

def align_fastqs(star_index_output, star_output_prefix, fastqr1, fastqr2, processors):
	p = str(processors)

	subprocess.run(['STAR', '--runThreadN', p, '--genomeDir', star_index_output, '--readFilesIn',
			fastqr1, fastqr2,
			'--outFileNamePrefix', star_output_prefix, '--alignEndsType', 'Extend5pOfRead1',
			'--outSAMmultNmax', '-1', '--outReadsUnmapped', 'Fastx', '--alignEndsProtrude', '50',
			'ConcordantPair', '--seedMultimapNmax', '1000000'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, universal_newlines=True)

def get_stats(star_output_prefix):
	df = pd.read_csv(star_output_prefix + 'Log.final.out', sep = '\t', header = None)

	total = df.iloc[4][1]
	mapped = df.iloc[7][1]
	perc_mapped = str(float(df.iloc[8][1][:-1]) + float(df.iloc[23][1][:-1])) + '%'
	perc_mapped_uniq = df.iloc[8][1]
	perc_mapped_multi = df.iloc[23][1]
	perc_unmapped = df.iloc[30][1]

	return total, mapped, perc_mapped, perc_mapped_uniq, perc_mapped_multi, perc_unmapped

def convert_and_sort_bam(star_output_prefix, processors):
	p = str(processors)

	sam_name = star_output_prefix + 'Aligned.out.sam'
	bam_name = star_output_prefix + '.bam'
	sorted_bam_name = star_output_prefix + 'sorted.bam'

	subprocess.run(['samtools', 'view', '-@', p, '-b', '-o', bam_name, sam_name], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, universal_newlines=True)

	subprocess.run(['samtools', 'sort', bam_name, '-@', p, '-o', sorted_bam_name], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, universal_newlines=True)

	subprocess.run(['samtools', 'index', sorted_bam_name], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, universal_newlines=True)

	subprocess.run(['rm', bam_name], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, universal_newlines=True)
