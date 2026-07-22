# the script assumes that it will be executed from the original project directory

mkdir -p ../../results/fastqc_results
fastqc -t 2 -o ../../results/fastqc_results/ ../../data/*.fastq.gz