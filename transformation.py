import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, numpy2ri, Formula
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import numpy as np

deseq = importr('DESeq2')
edger = importr('edgeR')
summarized_experiment = importr('SummarizedExperiment')

to_dataframe = robjects.r('function(x) data.frame(x)')

def create_count_matrix(dfs, insert_size):
	'''
	Input:
		A list of dataframes

	Output:
		A pandas DataFrame where each row is a insert/gene nucleotide and column is a separate quant file
	'''

	mat = pd.DataFrame(columns = list(range(len(dfs))),
			   index = [i + '_nt_' + str(j) for i in dfs[0]['inserts'] for j in list(range(1, insert_size + 1))])

	for i, df in enumerate(dfs):
		mat[i] = df[list(map(str, range(1, insert_size + 1)))].values.reshape(df.shape[0] * insert_size, -1)

	return mat

def tpm_transform(count_df, insert_size):
	df_norm = count_df.divide(insert_size / 1000)
	scaling_factors = count_df.sum() / 1000000
	df_norm = count_df.divide(scaling_factors)

	return df_norm

class CountMatrix:

	def __init__(self, count_matrix, threshold = 1):
		self.deseq_dds = None
		self.edger_dds = None
		self.vst = None
		self.rlog = None
		self.tmm = None
		self.original_count_matrix = count_matrix.reset_index(drop = True)
		self.filt_count_matrix = self.original_count_matrix[(self.original_count_matrix > threshold).any(1)]
		print(f'Concatenated DataFrame has {self.filt_count_matrix.shape[0]} rows and {self.filt_count_matrix.shape[1]} columns after filtering.\n')
		with localconverter(robjects.default_converter + pandas2ri.converter):
			self.count_matrix = robjects.conversion.py2rpy(self.filt_count_matrix)
			self.deseq_design_matrix = robjects.conversion.py2rpy(pd.DataFrame({'treatment':['ctrl' for i in range(self.original_count_matrix.shape[1])]}))
		self.deseq_design_formula = Formula('~ 1')

	def tmm_transform(self):
		self.edger_dds = edger.DGEList(counts = self.count_matrix)
		self.edger_dds = edger.calcNormFactors(self.edger_dds, method="TMM")
		self.tmm = edger.cpm(self.edger_dds)
		self.tmm = to_dataframe(self.tmm)
		with localconverter(robjects.default_converter + pandas2ri.converter):
			self.tmm = robjects.conversion.rpy2py(self.tmm)
		self.tmm.columns = self.original_count_matrix.columns.values
		self.tmm.index = map(int, self.tmm.index.values)
		self.original_count_matrix.update(self.tmm)
		return self.original_count_matrix

	def rlog_transform(self):
		self.deseq_dds = deseq.DESeqDataSetFromMatrix(countData=self.count_matrix, colData=self.deseq_design_matrix, design=self.deseq_design_formula)
		self.rlog = deseq.rlog(self.deseq_dds, blind=True)
		self.rlog = summarized_experiment.assay(self.rlog)
		self.rlog = to_dataframe(self.rlog)
		with localconverter(robjects.default_converter + pandas2ri.converter):
			self.rlog = robjects.conversion.rpy2py(self.rlog)
		self.rlog.columns = self.original_count_matrix.columns.values
		self.rlog.index = map(int, self.rlog.index.values)
		self.original_count_matrix.update(self.rlog)
		return self.original_count_matrix

	def vst_transform(self):
		self.deseq_dds = deseq.DESeqDataSetFromMatrix(countData=self.count_matrix, colData=self.deseq_design_matrix, design=self.deseq_design_formula)
		#self.vst = deseq.vst(self.deseq_dds, blind=True)
		self.vst = deseq.varianceStabilizingTransformation(self.deseq_dds, blind=True)
		self.vst = summarized_experiment.assay(self.vst)
		self.vst = to_dataframe(self.vst)
		with localconverter(robjects.default_converter + pandas2ri.converter):
			self.vst = robjects.conversion.rpy2py(self.vst)
		self.vst.columns = self.original_count_matrix.columns.values
		self.vst.index = map(int, self.vst.index.values)
		self.original_count_matrix.update(self.vst)
		return self.original_count_matrix
