# Index BAM files.
samtools index $1 #1 input.sorted.bam

#Calculate the RNA-seq reads coverage over gene body.
geneBody_coverage.py -r $3 -i $1 -o $2 #2 Prefix of output files. #3 Reference gene model in bed fomat.