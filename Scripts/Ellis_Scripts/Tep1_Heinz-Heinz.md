# Mapping Heinz to Heinz

## Files

Files include the reference genome  `S_lycopersicum_chromosome.2.50.fa` as well as the reads from the Lavelle_tie1 experiment `SRR404081.fastq.gz`

## Reference Index
An indexed reference file for `S_lycopersicum_chromosome.2.50.fa` was made by running `bwa index S_lycopersicum_chromosome.2.50.fa`

## Read Mapping
Heinz reads were subsequently mapped to the reference genome using


`bwa mem -R '@RG\tID:Heinz\tSM:Heinz' S_lycopersicum_chromosomes.2.50.fa SRR404081.fastq.gz | samtools view -Sbu | samtools rmdup - HeinzG_tep1_rmdup_unsorted.bam`

`samtools sort -m 8000000000 HeinzG_tep1_rmdup_unsorted.bam HeinzG_tep1_rmdup.bam`

These reads were merged with Leslie's `tep_M82_merge_1.bam` and `tep_M82_merge_2.bam` files with the following commands

	picard-tools AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT TMP_DIR=./ INPUT=HeinzG_tep1_rmdup.bam OUTPUT=HeinzG_tep1_rmdup_rg.bam RGSM=heinz RGLB=NA RGPL=Illumina RGPU=NA RGID=3

	picard-tools MergeSamFiles VALIDATION_STRINGENCY=LENIENT TMP_DIR=./ OUTPUT=tep_M82_heinz_merge_2.bam INPUT=tep_M82_merge_2.bam INPUT=HeinzG_tep1_rmdup_rg.bam CREATE_INDEX=true

	freebayes -f S_lycopersicum_chromosomes.2.50.fa -X -u -m 30 -F 0.4 tep_M82_heinz_merge_2.bam > heinz_tep_M.vcf

# Determining loss of depth

Start with bam files without duplicates removed.

	samtools depth tep_paired_2.bam tep_unpaired_2.bam > tep_2.depth

	cat tep_2.depth | awk '{sum+=$3} END { print "Average = ",sum/NR}'

Returned a mean depth of 12.6641

Next check coverage of files with duplicates removed.

	samtools depth tep_merge_rmdup_2.bam > tep_rmdup_2.depth

	cat tep_rmdup_2.depth | awk '{sum+=$3} END { print "Average = ",sum/NR}'

Returned a mean depth of 11.9954

After SNP calling mean depth of tep SNPs was 12.71 and after filtering the mean depth dropped to 9.90. On chromosome 9 (our chromosome of interest) mean read depth before filtering was 14.96 and 8.077 after. These numbers were determined using only the tep-1 read depths in R.

# What genomic area are we looking at?

	#determine genomic region where runmean.15 = 100
	area.int <- subset(tep_snp_data, tep.runmean.15 == 100)
	min(area.int$POS)
	max(area.int$POS)

returns values of 67664897 and 67942531 respectively

# Mapping for deletions

reads were mapped using subread-1.5.0-p2

Create fasta file for chr09

	samtools faidx S_lycopersicum_chromosomes.2.50.fa SL2.50ch09 > S_lycopersicum_ch09.2.50.fa

Create index for subread to use

	subread-buildindex -o S_lycopersicum_ch09.2.50.fa.index S_lycopersicum_ch09.2.50.fa

####Map tep, M82 and heinz to ch09.2.50.

__Tep-1:__

	subread-align -t 1 -T 3 -i S_lycopersicum_ch09.2.50.fa.index -I 200 -r Day1_R1.trim.paired.fq.gz -R Day1_R2.trim.paired.fq.gz -o Day1_SR.bam

	subread-align -t 1 -T 3 -i S_lycopersicum_ch09.2.50.fa.index -I 200 -r Day2_R1.trim.paired.fq.gz -R Day2_R2.trim.paired.fq.gz -o Day2_SR.bam

	subread-align -t 1 -T 3 -i S_lycopersicum_ch09.2.50.fa.index -I 200 -r Day3_R1.trim.paired.fq.gz -R Day3_R2.trim.paired.fq.gz -o Day3_SR.bam

	subread-align -t 1 -T 3 -i S_lycopersicum_ch09.2.50.fa.index -I 200 -r Day4_R1.trim.paired.fq.gz -R Day4_R2.trim.paired.fq.gz -o Day4_SR.bam

	subread-align -t 1 -T 3 -i S_lycopersicum_ch09.2.50.fa.index -I 200 -r Day1_all.trim.unpaired.fq.gz  -o Day1_SR_unpaired.bam

	subread-align -t 1 -T 3 -i S_lycopersicum_ch09.2.50.fa.index -I 200 -r Day2_all.trim.unpaired.fq.gz  -o Day2_SR_unpaired.bam

	subread-align -t 1 -T 3 -i S_lycopersicum_ch09.2.50.fa.index -I 200 -r Day3_all.trim.unpaired.fq.gz  -o Day3_SR_unpaired.bam

	subread-align -t 1 -T 3 -i S_lycopersicum_ch09.2.50.fa.index -I 200 -r Day4_all.trim.unpaired.fq.gz  -o Day4_SR_unpaired.bam
