import pandas as pd
import numpy as np

def sum2one(count_df, insert_size):
	'''
	Input:
		Pandas DataFrame with raw count, insert, rna, and dna columns

	Returns:
		Pandas DataFrame where count columns are normalized so that they sum up to 1
	'''

	norm_df = count_df[list(map(str, range(1, insert_size + 1)))].div(count_df[list(map(str, range(1, insert_size + 1)))].sum(axis = 1), axis = 0)
	norm_df.fillna(0, inplace=True)

	norm_df['inserts'], norm_df['rna'], norm_df['dna'] = count_df['inserts'], count_df['rna'], count_df['dna']

	return norm_df

def log_norm(count_df, insert_size, scaling_factor = 10000):
	'''
	Input:
		Pandas DataFrame with raw count, insert, rna, and dna columns

	Returns:
		Pandas DataFrame where count columns are log-normalized as in Seurat
		
	Normalization:
		1. Divide each cell by the total number of molecules measured in the cell
		2. Multiply that number by a scaling factor (default = 10,000)
		3. Add 1 and take the natural log
	'''

	summed_count = df[list(map(str, range(1, insert_size + 1)))].to_numpy().sum()

	norm_df = df[list(map(str, range(1, insert_size + 1)))].apply(lambda x: np.log(((x / summed_count) * scaling_factor) + 1))

	norm_df['inserts'], norm_df['rna'], norm_df['dna'] = count_df['inserts'], count_df['rna'], count_df['dna']

	return norm_df

