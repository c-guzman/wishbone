Wishbone is a tool to investigate, analyze, and compare TSS-MPRA data.

To install simply `git clone` this repository and then create a `conda` or `mamba` environment using the environment.yml file available.

# DOCUMENTATION
Wishbone contains five major commands that are meant to be run one after another.

1. align - create a reference index from a fasta file and then aligns TSS-MPRA to the generated reference index using STAR. Generated SAM files are converted into sorted and indexed BAM files.
2.quant - quantifies the 5' read positions and barcode frequencies from a sorteda nd indexed BAM file.
3. transform - transforms raw quantified read counts (but not barcode frequencies) into sequencing depth normalized counts (by default we use TPM)
4. norm - normalizes transformed read counts using a sum2one normalization where the sum of all values in the sequence sum to one. This is required for WIP distance calculations.
5. dist - calculates WIP dissimilarity scores from a transformed-normalized quant file.

# TYPICAL WORKFLOW

You can always get information about the available commands and parameters for each command by doing:

`python3 __main__.py -h`
`python3 __main__.py <command> -h`

1. Align read 1 and read 2 fastq files (ensure they are not compressed) to a custom fasta file that contains sequence information about your oligo pool

`python3 __main__.py align -r1 <read1.fq> -r2 <read2.fq> -fa <fasta.fa> -o <output_file_name>`

Output is a sorted and indexed BAM file

2. Quantify the number of reads at each nucleotide in a sequence as well as the barcode frequencies of each sequence from a sorted and indexed BAM file.

`python3 __main__.py quant -b <file.sorted.bam> -o <output_file_name> -size <size_of_sequences_in_pool> -dna <dna.sorted.bam> -bdg`

Output is a text file that contains columns 1-x where x represents the integar used in the -size parameter followed by the name of the sequence from the fasta file and the total RNA and DNA barcode counts

3. Transform raw TSS read counts to account for sequencing depth.

`python3 __main__.py transform -i <quant_files.txt> -o <output_file_name> -size <size_of_sequences_in_pool>`

Output is same as quant command except the number of reads at each nucleotide in a sequence are transformed for sequencing depth. While wishbone offers several options including vst, rlog, tmm and tpm, all of the work done in our paper using TPM transformed data.

RNA and DNA barcode frequencies are NOT transformed.

You can include multiple quant.txt files at this step to combine different oligo pool libraries.

4. Normalize transformed TSS read so that all counts sum up to one.

`python3 __main__.py norm -i <transformed.quant.txt> -o <output_file_name> -size <size_of_sequences_in_pool>`

Output is same as transform command except the number of reads at each nucleotide in a sequence have been normalized so that they sum to one.

5. Calculate all pairwise WIP dissimilarity scores from a transformed-normalized quant file. 

`python3 __main__.py dist -i <transformed.normalized.quant.txt> -o <output_file_name> -k <max_window_size>`

Output is a matrix where each sequence is compared against every other sequence in the quant file given and where each number represents the WIP dissimilarity score between those two sequences.

The -k parameter is the maximum window size that the algorithm will use when running. For our paper our -k parameter was always set to 5, but we have seen that using different parameters have interesting outputs (for example setting -k to 20 seems very useful in identifying highly differently expressed genes in csRNA-seq or CAGE data). 
