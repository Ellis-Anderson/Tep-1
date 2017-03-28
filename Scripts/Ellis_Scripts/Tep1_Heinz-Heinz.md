# Mapping Heinz to Heinz

## Files

Files include the reference genome  `S_lycopersicum_chromosome.2.50.fa` as well as the reads from the Lavelle_tie1 experiment `SRR404081.fastq.gz`

## Reference Index
An indexed reference file for `S_lycopersicum_chromosome.2.50.fa` was made by running `bwa index S_lycopersicum_chromosome.2.50.fa`

## Read Mapping
Heinz reads were subsequently mapped to the reference genome using


`bwa mem -R '@RG\tID:Heinz\tSM:Heinz' S_lycopersicum_chromosomes.2.50.fa SRR404081.fastq.gz | samtools view -Sbu | samtools rmdup - HeinzG_tep1_rmdup_unsorted.bam`

`samtools sort -m 8000000000 HeinzG_tep1_rmdup_unsorted.bam HeinzG_tep1_rmdup.bam`

These reads were merged with Leslie's `tep_M82_merge_2.bam` files with the following commands

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

#### determine genomic region where runmean.15 = 100
	area.int <- subset(tep_snp_data, tep.runmean.15 == 100)
	min(area.int$POS)
	max(area.int$POS)

returns values of 67664897 and 67942531 respectively

__What genes are in our region of interest?__
- [Solyc09g075990.2](https://solgenomics.net/feature/17936457/details) - 5&apos-nucleotidase surE (AHRD V1 ***- D6YVT2_WADCW); contains Interpro domain(s) IPR002828 Survival protein SurE-like phosphatase/nucleotidase
- [Solyc09g076000.2](https://solgenomics.net/feature/17936476/details) - Serine/threonine protein kinase (AHRD V1 *-** D0NNA3_PHYIN); contains Interpro domain(s) IPR002290 Serine/threonine protein kinase
- [Solyc09g076010.2](https://solgenomics.net/feature/17936491/details) - PHD finger family protein (AHRD V1 *--- D7MI23_ARALY); contains Interpro domain(s) IPR019787 Zinc finger, PHD-finger
- [Solyc09g076020.2](https://solgenomics.net/feature/17936508/details) - Imidazoleglycerol-phosphate dehydratase (AHRD V1 **** B9RTS0_RICCO); contains Interpro domain(s) IPR000807 Imidazole glycerol-phosphate dehydratase
- [Solyc09g076040.2](https://solgenomics.net/feature/17936536/details) - Protein SUPPRESSOR OF GENE SILENCING 3 homolog (AHRD V1 *--- SGN-S3 contains Interpro domain(s) IPR005380 Region of unknown function XS
- [Solyc09g076050.2](https://solgenomics.net/feature/17936543/details) - FRIGIDA (Fragment) (AHRD V1 *--- B0Z034_ARALP); contains Interpro domain(s) IPR012474 Frigida-like
- [Solyc09g082050.2](https://solgenomics.net/feature/17936550/details) - Histone-lysine N-methyltransferase (AHRD V1 *-** B6K768_SCHJY); contains Interpro domain(s) IPR003105 SRA-YDG
- [Solyc09g082060.2](https://solgenomics.net/feature/17936557/details) - Cysteine synthase (AHRD V1 **** Q3LAG5_TOBAC); contains Interpro domain(s) IPR005859 Cysteine synthase A
- [Solyc09g082090.1](https://solgenomics.net/feature/17936581/details) - Unknown Protein (AHRD V1)


```{unix}
wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG2.4_release/ITAG2.4_proteins_full_desc.fasta
grep "^>" ITAG2.4_proteins_full_desc.fasta
cat ITAG2.4_fasta_header.txt | sed ':a;s/^\(\([^"]*"[^"]*"[^"]*\)*[^"]*"[^"]*\) /\1_/;ta' | sed 's/evidence_code:[A-Za-Z0-9]*\ //' | sed 's/go_terms:[G,O,:,0-9]*\ //' | sed 's/GO:[0-9][0-9][0-9][0-9][0-9]*\ //' | tr " " "," > itag2.4_fasta_header.csv
```


# Mapping for deletions

Leslie's raw trimmed and barcode-split files from her Day_1 Day_2 Day_3 and Day_4 subdirectories were aligned to SL2.50
using STAR.

The Genome Index Files were generated with the following commands

	STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ./S_lycopersicum_chromosomes.2.50.fa --sjdbGTFfile ../Annotations/ITAG2.4_gene_models.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100

Reads were mapped in two subsequent rounds; one for paired end reads and one for unpaired reads.

	STAR --runThreadN 8 --genomeDir ../Genome/ --readFilesIn ../Reads/Day1_R1.trim.paired.fq.gz,../Reads/Day2_R1.trim.paired.fq.gz,../Reads/Day3_R1.trim.paired.fq.gz,../Reads/Day4_R1.trim.paired.fq.gz ../Reads/Day1_R2.trim.paired.fq.gz,../Reads/Day2_R2.trim.paired.fq.gz,../Reads/Day3_R2.trim.paired.fq.gz,../Reads/Day4_R2.trim.paired.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

	STAR --runThreadN 12 --genomeDir ../Genome/ --readFilesIn ../Reads/Day1_all.trim.unpaired.fq.gz,../Reads/Day2_all.trim.unpaired.fq.gz,../Reads/Day3_all.trim.unpaired.fq.gz,../Reads/Day4_all.trim.unpaired.fq.gz  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

These were saved as all_days.pair.ann.sorted.bam and all_days.pair.ann.sorted.bam respectively.
These two files were merged using samtools merged.

	samtools merge all_days.ann.sorted.bam ./Star_GFF/all_days.pair.ann.sorted.bam ./Star_GFF_unpaired/all_days.unpair.ann.sorted.bam

The combined BAM file was used to generate a BedGraph file to be analyzed in IGV.

	STAR --runThreadN 10 --runMode inputAlignmentsFromBAM --inputBAMfile all_days.ann.sorted.bam --outWigType bedGraph --outWigStrand Unstranded

The gff3 file containing gene annotations for sl2.50 was downloaded from solgenomics.net using

	wget ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/annotation/ITAG2.4_release/ITAG2.4_gene_models.gff3

#### Use snpEff to predict possible outcomes of mutations

heinz_tep_M.vcf file was subset to include the original header and chromosome 9 using

	grep '^#' heinz_tep_M.vcf >> heinz_tep_M_ch09.vcf
	grep '^SL2.50ch09' heinz_tep_M.vcf >> heinz_tep_M_ch09.vcf

This vcf file was analyzed using snpEff and snpEffs database for SL2.50

	java -Xmx4g -jar ~/snpEff/snpEff.jar -v SL2.50 heinz_tep_M_ch09.vcf
	java -Xmx4g -jar ~/snpEff/snpEff.jar -v SL2.50 heinz_tep_M_ch09.vcf > heinz_tep_M_ch09_ann.vcf

This was analyzed in R studio with the following commands.
The script can be found in `~/Documents/Maloof_lab/Tep-1/Scripts/Ellis_Scripts/des_tep_snpeff.rmd`
Nothing of particular interest was found.
