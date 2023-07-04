# Install seqtk with github.
git clone https://github.com/lh3/seqtk.git;

#Subsample 10000 read pairs from two large paired FASTQ files:
seqtk sample -s100 read1.fq 10000 > sub1-1k.fq
seqtk sample -s100 read2.fq 10000 > sub2-1k.fq
