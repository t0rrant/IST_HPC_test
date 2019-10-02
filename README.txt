flowcraft build -t "reads_download downsample_fastq integrity_coverage fastqc_trimmomatic check_coverage fastqc spades assembly_mapping pilon mlst" -o simple_innuca.nf


Query in SRA (https://www.ncbi.nlm.nih.gov/sra) : ("illumina"[Platform]) AND "Streptococcus pneumoniae"[orgn:__txid1313] 

Run with slurm+shifter:
nextflow run -profile slurm_shifter simple_innuca.nf

Monitor progress:
flowcraft inspect -m broadcast
