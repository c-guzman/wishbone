import pandas as pd
import pysam

def read_df(tsv):
	df = pd.read_csv(tsv, sep = '\t', index_col = None, header = 0)

	insert_names, rna, dna = df['inserts'], df['rna'], df['dna']

	return df, insert_names, rna, dna

def read_bam(bam):
	bam = pysam.AlignmentFile(bam, "rb")

	return bam
		
def read_bed(bed):
	df = pd.read_csv(bed, sep = '\t', index_col = None, header = None)

	insert_names = df.loc[:, 3]

	return df, insert_names
