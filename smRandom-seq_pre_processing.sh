
# STAR Creating a genome index
STAR --runThreadN N --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles $genomeFastaFile --sjdbGTFfile $sjdbGTFfile --sjdbOverhang 122 --sjdbGTFfeatureExon gene --sjdbGTFtagExonParentTranscript gene_id --genomeSAindexNbases 10

# STAR read alignment.
# The output filetype was BAM.
STAR --runThreadN N --runMode alignReads --genomeDir $starIndex --readFilesIn ${sample}_2.fq.extracted  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 41143265264 --outFileNamePrefix $sample.

# Exon based gene expression matrix. 
featureCounts -t exon -g gene_id -a $gtf -o gene_assigned_gene -R BAM $sample.Aligned.sortedByCoord.out.bam -T N # FeatureCounts is used to count single bacteria RNA-seq reads for genomic features in BAM files. Alignments are provided in BAM format. Annotations for gene regions are provided in GTF format.
mv $sample.Aligned.sortedByCoord.out.bam.featureCounts.bam $sample.Aligned.sortedByCoord.featureCountsUniq.bam # Rename the BAM file.

# Sort alignments 
samtools sort -@ N $sample.Aligned.sortedByCoord.featureCountsUniq.bam -o $sample.Aligned.sortedByCoord.featureCountsUniq.sorted.bam

# Index the coordinate-sorted BGZIP-compressed BAM files for fast random access.
samtools index $sample.Aligned.sortedByCoord.featureCountsUniq.sorted.bam

# Count reads per gene from BAM using UMIs and mapping coordinates
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I $sample.Aligned.sortedByCoord.featureCountsUniq.sorted.bam -S $sample.counts.tsv

# Select cell barcodes.
python selectResult.py $sample.counts.tsv $thres $sample.selected #sample.counts.tsv: Expression matrix. #thres: Cell number. #sample.selected:outPrefix
