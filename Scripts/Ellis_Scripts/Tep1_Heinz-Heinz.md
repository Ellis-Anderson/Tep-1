# Mapping Heinz to Heinz

## Files

Files include the reference genome  `S_lycopersicum_chromosome.2.50.fa` as well as the reads from the Lavelle_tie1 experiment `SRR404081.fastq.gz` 

## Reference Index
An indexed reference file for `S_lycopersicum_chromosome.2.50.fa` was made by running `bwa index S_lycopersicum_chromosome.2.50.fa`

## Read Mapping
Heinz reads were subsequently mapped to the reference genome using 


`bwa mem -R '@RG\tID:Heinz\tSM:Heinz' S_lycopersicum_chromosomes.2.50.fa SRR404081.fastq.gz | samtools view -Sbu | samtools rmdup - HeinzG_tep1_rmdup_unsorted.bam`

`samtools sort -m 8000000000 HeinzG_tep1_rmdup_unsorted.bam HeinzG_tep1_rmdup_sort.bam`

These reads were merged with Leslie's `tep_M82_merge_1.bam` and `tep_M82_merge_2.bam` files with the following commands

	picard-tools AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT TMP_DIR=./ INPUT=tep_M82_merge_1.bam OUTPUT=tep_M82_merge_rg_1.bam RGSM=heinz RGLB=NA RGPL=Illumina RGPU=NA RGID=3

	picard-tools AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT TMP_DIR=./ INPUT=tep_M82_merge_2.bam OUTPUT=tep_M82_merge_rg_2.bam RGSM=heinz RGLB=NA RGPL=Illumina RGPU=NA RGID=3

	picard-tools MergeSamFiles VALIDATION_STRINGENCY=LENIENT TMP_DIR=./ OUTPUT=tep_M82_heinz_merge_1.bam INPUT=tep_M82_merge_rg_1.bam INPUT=HeinzG_tep1_rmdup.bam CREATE_INDEX=true

	picard-tools MergeSamFiles VALIDATION_STRINGENCY=LENIENT TMP_DIR=./ OUTPUT=tep_M82_heinz_merge_1.bam INPUT=tep_M82_merge_rg_2.bam INPUT=HeinzG_tep1_rmdup.bam CREATE_INDEX=true

### SNP Calling

SNP calling was done with Freebayes using the same flags and settings that leslie used. 

	freebayes -f S_lycopersicum_chromosomes.2.50.fa -X -u -m 30 -F 0.4 tep_M82_heinz_merge_1.bam > tep_M82_heinz_1.vcf
	
	freebayes -f S_lycopersicum_chromosomes.2.50.fa -X -u -m 30 -F 0.4 tep_M82_heinz_merge_2.bam > tep_M82_heinz_2.vcf








