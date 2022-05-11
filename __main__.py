import argparse
import pandas as pd
import numpy as np
import csv

from reading import read_df, read_bam, read_bed
from quantification import get_insert_names_from_bam, create_count_df, count_csrna_rna, count_mpra_rna, count_mpra_dna, create_bdg # log_norm
from normalization import sum2one
from distance import tss_dist, transform_arr_to_df
from alignment import create_index, align_fastqs, get_stats, convert_and_sort_bam
from transformation import create_count_matrix, CountMatrix, tpm_transform

def main():
	parser = argparse.ArgumentParser(prog = 'wishbone', description = 'Tool to investigate, analyze, and compare TSS-MPRA data')

	setup_parsers(parser)

	args = parser.parse_args()

	if args.command == 'quant':
		run_quant(args)

	elif args.command == 'norm':
		run_norm(args)

	elif args.command == 'dist':
		run_dist(args)

	elif args.command == 'align':
		run_align(args)

	elif args.command == 'transform':
		run_transform(args)

def setup_parsers(parser):

	subparsers = parser.add_subparsers(dest='command')

	align_parser = subparsers.add_parser('align', help='Align unzipped reporter fastq files')

	quant_parser = subparsers.add_parser('quant', help="Quantify 5' reads and barcode frequencies of a sorted and indexed BAM file")

	transform_parser = subparsers.add_parser('transform', help="Transform raw quantified read counts into normalized counts")

	norm_parser = subparsers.add_parser('norm', help="Normalize transformed read counts for distance calculations")

	distance_parser = subparsers.add_parser('dist', help="Calculate distance metrics of a quant file")

	### START QUANT ARGUMENTS ###
	quant_parser.add_argument('-b', '--bam', help='A sorted and indexed BAM file', action='store', type=str, required=True)

	quant_parser.add_argument('-o', '--output_file', help='Name of the resulting output file', action='store', type=str, default='output.counts.txt')

	quant_parser.add_argument('-size', '--insert_size', help='Size of the insert pool. All inserts must be the same size. (Default: 200)', action='store', type=int, default=200, required=True)

	quant_parser.add_argument('-bed', '--bed', help='A bed file of genomic positions to quantify if using csRNA-seq mode', action='store', type=str)

#	quant_parser.add_argument('-norm', '--norm', help='Turn on normalization', action='store_true')

	quant_parser.add_argument('-dna', '--dna', help='A sorted and indexed BAM file', action='store', type=str)

	quant_parser.add_argument('-bdg', '--bedgraph', help='Create a bedgraph from count information', action='store_true')

	quant_parser.add_argument('-csrna', '--csrna', help='Quantification is being performed on csRNA-seq mapped data', action='store_true')

	quant_parser.add_argument('-s', '--shift', help='Number of bp to shift csRNA-seq positions to align with inserts. (Default: 0)', action='store', type=int, default=0)
	### END QUANT ARGUMENTS

	### START NORM ARGUMENTS ###
	norm_parser.add_argument('-i', '--infile', help='TSV file of counts output by quant', action='store', type=str, required=True)

	norm_parser.add_argument('-o', '--output_file', help='Name of the resulting output file', action='store', type=str, default='output.normed_counts.txt')

	norm_parser.add_argument('-size', '--insert_size', help='Size of the insert pool. All inserts must be the same size. (Default: 200)', action='store', type=int, default=200, required=True)
	### END NORM ARGUMENTS ###

	### START DISTANCE ARGUMENTS ###
	distance_parser.add_argument('-i', '--infiles', help='TSV file(s) of counts output by quant', action='store', type=str, required=True, nargs='+')

	distance_parser.add_argument('-o', '--output_file', help='Name of the resulting output file', action='store', type=str, default='output.distances.txt')

	distance_parser.add_argument('-k', '--ksize', help='Window sliding length for distance calculation (Default: 5)', action='store', type=int, default=5)
	### END DISTANCE ARGUMENTS ###

	### START ALIGNMENT ARGUMENTS ###
	align_parser.add_argument('-r1', '--read1', help='Read 1 unzipped fastq file', action='store', type=str, required=True)

	align_parser.add_argument('-r2', '--read2', help='Read 2 unzipped fastq file', action='store', type=str, required=True)

	align_parser.add_argument('-fa', '--fasta', help='The insert pool fasta file.', action='store', type=str, required=True)

	align_parser.add_argument('-o', '--output', help='Prefix name for output file. (Default: output.)', type=str, default='output.')

	align_parser.add_argument('-oi', '--output_index', help='The name of the output directory for the STAR index. (Default: index/)', type=str, default='index/')

	align_parser.add_argument('-p', '--processors', help='The number of processors to use when possible. (Default: 1)', type=int, default=1)

	### START TRANSFORM ARGUMENTS ###
	transform_parser.add_argument('-i', '--infiles', help='TSV file(s) of counts output by quant', action='store', type=str, required=True, nargs='+')

	transform_parser.add_argument('-o', '--output_file', help='Name of the resulting output file', action='store', type=str, default='output.transformed_counts.txt')

	transform_parser.add_argument('-t', '--transformation', help='The type of transformation normalization to apply to raw counts (Default: tpm)', action='store', type=str, choices=['vst', 'rlog', 'tmm', 'tpm'], default='tpm')

	transform_parser.add_argument('-size', '--insert_size', help='Size of the insert pool. All inserts must be the same size. (Default: 200)', action='store', type=int, default=200, required=True)
	### END TRANSFORM ARGUMENTS ###

def run_quant(args):
	print(f'\nRunning RNA quantification on file: {args.bam}.\n')
	print(f'Input data is csRNA-seq: {args.csrna}\n')
#	print(f'Normalize Counts: {args.norm}\n')
	print(f'Creating a bedgraph of counts: {args.bedgraph}\n')

	if not args.csrna:
		rna_bam = read_bam(args.bam)
		insert_names = get_insert_names_from_bam(rna_bam)
		count_df = create_count_df(insert_names, args.insert_size)
		
		raw_count_df = count_mpra_rna(rna_bam, insert_names, count_df)

		rna_bam.close()

	else:
		rna_bam = read_bam(args.bam)
		bed, insert_names = read_bed(args.bed)
		count_df = create_count_df(insert_names, args.insert_size)

		raw_count_df = count_csrna_rna(rna_bam, bed, count_df, args.shift)

		rna_bam.close()

	if args.dna is not None:
		print(f'Running quantification on DNA file: {args.dna}.\n')

		dna_bam = read_bam(args.dna)

		raw_count_df = count_mpra_dna(dna_bam, insert_names, raw_count_df)

		dna_bam.close()

	else:
		print(f'No matched DNA file used.\n')

#	if args.norm:
#		norm_count_df = log_norm(raw_count_df, args.insert_size, scaling_factor = 10000)
#		norm_count_df.to_csv(args.output_file, sep='\t', index = False, quoting = csv.QUOTE_NONE)
#
#		if args.bedgraph:
#			bdg_df = create_bdg(norm_count_df, insert_names, args.insert_size)
#			bdg_out_name = args.output_file + '.bdg'
#			bdg_df.to_csv(bdg_out_name, sep='\t', index = False, quoting = csv.QUOTE_NONE, header = False)
#
#	else:
	raw_count_df.to_csv(args.output_file, sep='\t', index = False, quoting = csv.QUOTE_NONE)

	if args.bedgraph:
		bdg_df = create_bdg(raw_count_df, insert_names, args.insert_size)
		bdg_out_name = args.output_file + '.bdg'
		bdg_df.to_csv(bdg_out_name, sep='\t', index = False, quoting = csv.QUOTE_NONE, header = False)

def run_norm(args):
	print(f'\nPerforming normalization on file: {args.infile}\n')

	df, insert_names, rna, dna = read_df(args.infile)

	sum2one_df = sum2one(df, args.insert_size)

	sum2one_df.to_csv(args.output_file, sep = '\t', index = False, quoting = csv.QUOTE_NONE)

def run_dist(args):
	print(f'\nWill merge the following count files for pairwise tss distance calculations: {args.infiles}\n')

	dfs = []
	insert_names_list = []
	rna_list = []
	dna_list = []

	for tsv in args.infiles:
		df, insert_names, rna, dna = read_df(tsv)
		dfs.append(df)
		insert_names_list.append(insert_names)
		rna_list.append(rna)
		dna_list.append(dna)

	count_df = pd.concat(dfs, axis = 0).reset_index(drop=True)
	count_df = count_df.drop([x for x in ['rna', 'dna', 'inserts'] if x in count_df.columns], axis = 1)
	print(f'Concatenated DataFrame has {count_df.shape[0]} rows and {count_df.shape[1]} columns.\n')

	insert_names = [i for sublist in insert_names_list for i in sublist]
	rna_counts = [i for sublist in rna_list for i in sublist]
	dna_counts = [i for sublist in dna_list for i in sublist]

	print(f'Calculating distances ...\n')
	kmer_arr = tss_dist(count_df.values, args.ksize)

	print(f'Done with distance calculation, transforming matrix into DataFrame for output ...\n')
	kmer_df = transform_arr_to_df(kmer_arr, insert_names, rna_counts, dna_counts)

	kmer_df.to_csv(args.output_file, sep = '\t', index = False, quoting = csv.QUOTE_NONE)

	kmer_df.to_feather(args.output_file + '.feather')

def run_align(args):
	print(f'\nCreating index file {args.output_index} using fasta file {args.fasta}\n')
	create_index(args.fasta, args.output_index, args.processors)

	print(f'Aligning reads {args.read1} and {args.read2} using genome directory {args.output_index} with {args.processors} processors.\n')
	print(f'Output files prefix will be: {args.output}\n')
	align_fastqs(args.output_index, args.output, args.read1, args.read2, args.processors)

	print(f'Done mapping ... Calculating mapping stats\n')
	total, mapped, perc_mapped, perc_mapped_uniq, perc_mapped_multi, perc_unmapped = get_stats(args.output)

	print(f'Total number of reads attempted to map: {total}')
	print(f'Number of reads that mapped uniquely: {mapped}')
	print(f'Percentage of all reads that mapped (uniquely + multimapped): {perc_mapped}')
	print(f'Percentage of reads that mapped uniquely only: {perc_mapped_uniq}')
	print(f'Percentage of reads that are multimappers: {perc_mapped_multi}')
	print(f'Percentage of reads that were unmapped: {perc_unmapped}\n')

	print(f'Converting SAM to sorted BAM and generating index\n')
	convert_and_sort_bam(args.output, args.processors)

def run_transform(args):
	print(f'\nWill merge the following count files for transformation: {args.infiles}\n')

	dfs = []
	insert_names_list = []
	rna_list = []
	dna_list = []

	for tsv in args.infiles:
		df, insert_names, rna, dna = read_df(tsv)
		dfs.append(df)
		insert_names_list.append(insert_names)
		rna_list.append(rna)
		dna_list.append(dna)

	count_df = create_count_matrix(dfs, args.insert_size)
	print(f'Concatenated DataFrame has {count_df.shape[0]} rows and {count_df.shape[1]} columns.\n')

	countTransformer = CountMatrix(count_df)

	print(f'Applying {args.transformation} transformation to counts ...\n')
	if args.transformation == 'vst':
		transformed_count_df = countTransformer.vst_transform()

	elif args.transformation == 'rlog':
		transformed_count_df = countTransformer.rlog_transform()

	elif args.transformation == 'tmm':
		transformed_count_df = countTransformer.tmm_transform()

	elif args.transformation == 'tpm':
		transformed_count_df = tpm_transform(count_df, args.insert_size)

	print(f'Finished transformation, reshaping transformed counts into DataFrame for output ...\n')
	final_count_df = pd.DataFrame(transformed_count_df.T.values.reshape(-1, args.insert_size), columns = list(range(1, args.insert_size + 1)))

	final_count_df['inserts'], final_count_df['rna'], final_count_df['dna'] = pd.concat(insert_names_list, ignore_index = True), pd.concat(rna_list, ignore_index = True), pd.concat(dna_list, ignore_index = True)

	final_count_df.to_csv(args.output_file, sep = '\t', index = False, quoting = csv.QUOTE_NONE)

if __name__ == '__main__':
	main()
