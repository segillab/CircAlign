# CircAlign
#### Written By: Litao Tao and [Frankie James]
---
 This script will align reads (in fastq format) to the mitochondrial genome or any other circular DNA.

Required Packages:
- [BWA]
- [Samtools]
- [Bedtools]
- [bedGraphToBigWig]

To install Packages using Conda, run:
```
bash install_dependencies.sh
```
If you do not have Conda installed, instructions are [Here].

Note: If you would like to manually install dependencies, comment out the line `source activate circAlign` in CircAlign.sh.

Input:
1) read1 and read2 (if paired-end) fastq file (can be gzipped or not)
2) shifting length (int, should be slightly bigger than read length)
3) reference genome sequence in fasta format (one single sequence)

Output:
1) Coverage bedGraph and bigwig file
2) BedGraph and bigwig files for transposition (for ATACseq) or Cutting (for DNase) on forward and reverse stand

To Run (PE):
```
./CircAlign -f /path/to/genome.fasta -1 /path/to/read1.fastq -2 /path/to/read2.fastq.gz -r 55
```

Usage:
```
Circular Aligner
----------------
	Required Arguments:
		-f [FASTA_FILE] Fasta file with Chromosome M in one sequence
		-1 [FASTQ_FILE] Fastq file for read (Read 1 if Paired End)
		-r [READ_LENGTH] Length of reads found in [FASTQ_FILE]
	Optional Arguments:
		-2 [FASTQ_FILE] Fastq file for read 2 (Paired end data)
		-o [OUTPUT_DIRECTORY] Directory to place output results into (./ is default)
		-t [THREAD_#] Number of threads to use
		-d This will run the data as DNase (ATACseq is default)
		-s This will force the BigWigs generated to not be scaled (Default is to scale by read number aligned to [FASTA_FILE])
```

[bwa]: <https://github.com/lh3/bwa>
[samtools]: <https://github.com/samtools/samtools>
[bedtools]: <http://bedtools.readthedocs.io/en/latest/>
[bedgraphtobigwig]: <http://hgdownload.soe.ucsc.edu/admin/exe/>
[frankie james]: <http://github.com/fjames003>
[Here]: <https://conda.io/docs/user-guide/install/index.html#regular-installation>
