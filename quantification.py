import argparse

import pysam
import pandas as pd
import numpy as np
import itertools
import csv

def get_insert_names_from_bam(bam):
	insert_names = bam.references

	return list(insert_names)

def create_count_df(insert_names, insert_size):
	l = len(insert_names) # row size

	count_df = pd.DataFrame(index = list(range(0, l)))
	count_df = pd.concat([count_df, pd.DataFrame(columns = list(range(1, insert_size + 1)))])

	for col in count_df.columns:
		count_df[col].values[:] = 0

	count_df['inserts'] = insert_names
	count_df['rna'] = 0
	count_df['dna'] = 0

	return count_df

def count_csrna_rna(bam, bed, count_df, shift):
	df = count_df.copy()

	for r in range(bed.shape[0]):
		chrom = str(bed.iloc[r][0])
		start = int(bed.iloc[r][1])
		end = int(bed.iloc[r][2])
		strand = str(bed.iloc[r][5])

		reads = bam.fetch(chrom, start, end)

		if strand == '+':
			for read in reads:
				if not read.is_secondary and not read.is_unmapped and not read.is_reverse and read.reference_start > start:
					tags = read.tags[0]
					if not tags[1] > 1:
						if shift == 0:
							pos = read.reference_start - start + shift
						elif shift == 1:
							pos = read.reference_start - start + shift
						else:
							pos = read.reference_start - start + (shift + 1)

						df.at[r, pos] += 1
						df.at[r, 'rna'] += 1

		elif strand == '-':
			for read in reads:
				if not read.is_secondary and not read.is_unmapped and read.is_reverse and read.reference_end < end:
					tags = read.tags[0]
					if not tags[1] > 1:
						pos = end - read.reference_end + shift

						df.at[r, pos] += 1
						df.at[r, 'rna'] += 1

	return df

def count_mpra_rna(bam, insert_names, count_df):

	for i, insert in enumerate(insert_names):

		reads = bam.fetch(insert)

		for read in reads:
			if read.is_proper_pair and not read.is_supplementary and not read.is_secondary and not read.is_unmapped and read.is_read1 and read.reference_start != 0:
				tags = read.tags[0]
				if not tags[1] > 1:
					count_df.at[i, read.reference_start] += 1

			elif read.is_proper_pair and not read.is_supplementary and not read.is_secondary and not read.is_unmapped and read.is_read2:
				tags = read.tags[0]
				if not tags[1] > 1:
					count_df.at[i, 'rna'] += 1

	return count_df

def count_mpra_dna(bam, insert_names, count_df):

	for i, insert in enumerate(insert_names):

		reads = bam.fetch(insert)

		for read in reads:
			if read.is_proper_pair and not read.is_supplementary and not read.is_secondary and not read.is_unmapped and read.is_read2:
				tags = read.tags[0]
				if not tags[1] > 1:
					count_df.at[i, 'dna'] += 1

	return count_df

def create_bdg(count_df, insert_names, insert_size):
	count_df = count_df.drop([x for x in ['rna', 'dna', 'inserts'] if x in count_df.columns], axis = 1)

	inserts = list(itertools.chain.from_iterable(itertools.repeat(x, insert_size) for x in insert_names))
	positions = list(range(1, insert_size + 1))
	start_end_list = positions * count_df.shape[0]
	coverage_list = []

	for i in range(count_df.shape[0]):
		coverage_list = np.concatenate((coverage_list, list(count_df.iloc[i, :])))

	bdg_df = pd.DataFrame(columns = ['inserts', 'start', 'end', 'coverage'])

	bdg_df['inserts'], bdg_df['start'], bdg_df['end'], bdg_df['coverage'] = inserts, start_end_list, start_end_list, coverage_list

	return bdg_df
