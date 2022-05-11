import numpy as np
import numba as nb
import pandas as pd

@nb.njit
def tss_dist(arr, k):

	h, w = arr.shape[0], arr.shape[1]
	diffs = np.zeros((h, h), dtype=arr.dtype)
	max_k = k + 1

	for n in range(1, max_k):
		arr1 = np.empty(shape=(h, w-n+1), dtype=arr.dtype)

		for i in range(h):
			# Efficient sliding window algorithm
			assert w >= n
			s = np.sum(arr[i, 0:n])
			arr1[i, 0] = s

			for j in range(n, w):
				s -= arr[i, j-n]
				s += arr[i, j]
				arr1[i, j-n+1] = s

		# Efficient distance matrix computation
		for i in range(h):
			for j in range(i+1, h):
				s = 0
				for k in range(w-n+1):
					s += np.abs(arr1[i,k] - arr1[j,k])
				diffs[i, j] += s * n

	# Fill in lower triangular part
	for i in range(h):
		for j in range(i):
			diffs[i, j] = diffs[j, i]

	return diffs

def transform_arr_to_df(arr, insert_names, rna, dna):
	indexes_rows = list(range(arr.shape[0]))
	indexes_columns = list(map(str, range(arr.shape[0])))

	diffs_df = pd.DataFrame(arr, index = indexes_rows, columns = indexes_columns)

	diffs_df['inserts'], diffs_df['rna'], diffs_df['dna'] = insert_names, rna, dna

	return diffs_df
