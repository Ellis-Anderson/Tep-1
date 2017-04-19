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


``` {unix}
wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG2.4_release/ITAG2.4_proteins_full_desc.fasta
grep "^>" ITAG2.4_proteins_full_desc.fasta
cat ITAG2.4_fasta_header.txt | sed ':a;s/^\(\([^"]*"[^"]*"[^"]*\)*[^"]*"[^"]*\) /\1_/;ta' | sed 's/evidence_code:[A-Za-Z0-9]*\ //' | sed 's/go_terms:[G,O,:,0-9]*\ //' | sed 's/GO:[0-9][0-9][0-9][0-9][0-9]*\ //' | tr " " "," > itag2.4_fasta_header.csv
```
```

# Mapping for deletions

#### Mapping with standard STAR protocol
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

#### Analyzing in IGV

I searched the genomic area determined to be 100% Tep-1 for deletions in exonic regions of genic DNA.
I found two suggested deletions in exonic DNA.
One in [Solyc09g076020.2](https://solgenomics.net/feature/17936508/details) - Imidazoleglycerol-phosphate dehydratase (AHRD V1 **** B9RTS0_RICCO); contains Interpro domain(s) IPR000807 Imidazole glycerol-phosphate dehydratase
and another in [Solyc09g082090.1](https://solgenomics.net/feature/17936581/details) - Unknown Protein (AHRD V1)

#### Mapping Using Alexander Dobin's Suggestions for genomic DNA

The first mapping attempt split reads over huge distances, I want to try again without using annotations, instead using a 2 pass method.
I also want to try mapping for genomic data using the suggestions Alexander Dobin, author of STAR, suggests.

	STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ./S_lycopersicum_chromosomes.2.50.fa

	STAR --runThreadN 12 --genomeDir ../Genome_No_Ann/ --readFilesIn ../Reads/Day1_R1.trim.paired.fq.gz,../Reads/Day2_R1.trim.paired.fq.gz,../Reads/Day3_R1.trim.paired.fq.gz,../Reads/Day4_R1.trim.paired.fq.gz ../Reads/Day1_R2.trim.paired.fq.gz,../Reads/Day2_R2.trim.paired.fq.gz,../Reads/Day3_R2.trim.paired.fq.gz,../Reads/Day4_R2.trim.paired.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --twopassMode Basic

	STAR --runThreadN 12 --genomeDir ../Genome_No_Ann/ --readFilesIn ../Reads/Day1_all.trim.unpaired.fq.gz,../Reads/Day2_all.trim.unpaired.fq.gz,../Reads/Day3_all.trim.unpaired.fq.gz,../Reads/Day4_all.trim.unpaired.fq.gz  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --twopassMode Basic

	samtools merge all_days.sorted.bam ./Star_2pass_pair/all_days.pair.sorted.bam ./Star_2pass_unpair/all_days.unpair.sorted.bam



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



## Re-barcoding and Trimming tep reads

### Barcoding with Autobarcode

### Barcoding with Python
As auto-barcode disturbs the order of the reads in the fastq files it returns, further downstream analysis with BWA became challenging.
Overall quality of SNPs was also unimportant in this stage as we were looking for deletions within the region of interest on chromosome 9 and not performing any variant calling.
Therefore instead of using auto-barcode to seperate the reads by day, to later be recombined, I wrote a simple python script to remove the first 6 bases(the barcodes) of all the reads in the JMAS010 fastq files.
This can be found in my Tep-1 github repo under `Tep-1/Data/Barcoding/trim-fq.py`.
I called the program with the following commands:

	python trim-fq.py JMAS010.1.fq
	python trim-fq.py JMAS010.2.fq

These calls produced the `bar.JMAS010.1.fq` and `bar.JMAS010.2.fq` files in the barcoding folder.
According to some simple math and visual inspection this program seems to have done its job correctly.

### Trimming
After barcoding the reads still needed to be trimmed for quality.
Visual inspection of the original fastq files showed some reads with zero usable information.
The `bar.JMAS010.*.fq` files were trimmed using the following command:

	trimmomatic PE -threads 10 -phred33 ../Barcoding/bar.JMAS010.1.fq ../Barcoding/bar.JMAS010.2.fq bar.JMAS010.1.paired.fq bar.JMAS010.1.unpaired.fq bar.JMAS010.2.paired.fq  bar.JMAS010.2.unpaired.fq ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

__Output:__

	Input Read Pairs: 206012240 Both Surviving: 118942857 (57.74%) Forward Only Surviving: 11943836 (5.80%) Reverse Only Surviving: 8409328 (4.08%) Dropped: 66716219 (32.38%)

### Mapping with BWA (Again)
After trimming the read pairs I used BWA mem to align them to the SL2.50 genome.
With the extra coverage provided by the now included reads, we hoped to get a higher fidelity picture of deletions in our region of interest.

	bwa mem -t 8 S_lycopersicum_chromosomes.2.50.fa  ../Trimming/bar.JMAS010.1.paired.fq ../Trimming/bar.JMAS010.2.paired.fq | samtools view  -u - | samtools sort -m 2000000000 - | samtools rmdup -  tep1.PE.sort.rmdup.bam

Upon inspection of the newly generated bam file, the potential deletion in [Solyc09g076020.2](https://solgenomics.net/feature/17936508/details) now had coverage, but the potential deletion in [Solyc09g082090.1](https://solgenomics.net/feature/17936581/details) continued to have zero coverage.












# Failed Code

### Auto-Barcode
We found potential spots of interest so now we want to include more data(if possible) to confirm our speculative findings before proceeding.
Using Mike Covingtons autobarcode to remove barcodes

	auto_barcode/barcode_split_trim.pl --id tep1 -m 1 -l -b epi1_barcode.list -o R1.split JMAS010.1.fq
	auto_barcode/barcode_split_trim.pl --id tep1 -m 1 -l -b epi1_barcode.list -o R2.split JMAS010.2.fq

This produced several log files:
	log_barcode_counts.fq_JMAS010.1.bar_epi1_barcode.list
	log_barcodes_observed.fq_JMAS010.1.bar_epi1_barcode.list
	log_barcode_counts.fq_JMAS010.2.bar_epi1_barcode.list
	log_barcodes_observed.fq_JMAS010.2.bar_epi1_barcode.list
These can be found in the R1.split and R2.split directories respectively

### Trimmomatic
Split Reads from Days 1-4 were trimmed using the following commands

##### Day 1

	trimmomatic PE -threads 10 -phred33 ../Barcoding/R1.split/Day1.AAGAC.fq ../Barcoding/R2.split/Day1.AAGAC.fq R1.D1.trim.paired.fq R1.D1.trim.unpaired.fq R2.D1.trim.paired.fq  R2.D1.trim.unpaired.fq ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

	__output:__ Input Read Pairs: 38244690 Both Surviving: 20710085 (54.15%) Forward Only Surviving: 8075309 (21.11%) Reverse Only Surviving: 6663326 (17.42%) Dropped: 2795970 (7.31%)

unpaired reads were merged with the following

	cat R1.D1.trim.unpaired.fq R2.D1.trim.unpaired.fq > all.D1.trim.unpaired.fq

##### Day 2

	trimmomatic PE -threads 10 -phred33 ../Barcoding/R1.split/Day2.CATTA.fq ../Barcoding/R2.split/Day2.CATTA.fq Day2/R1.D2.trim.paired.fq Day2/R1.D2.trim.unpaired.fq Day2/R2.D2.trim.paired.fq  Day2/R2.D2.trim.unpaired.fq ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

	__Output:__ Input Read Pairs: 41372544 Both Surviving: 22267370 (53.82%) Forward Only Surviving: 8618299 (20.83%) Reverse Only Surviving: 7486798 (18.10%) Dropped: 3000077 (7.25%)

unpaired reads were merged with the following

	cat R1.D2.trim.unpaired.fq R2.D2.trim.unpaired.fq > all.D2.trim.unpaired.fq

##### Day 3

	trimmomatic PE -threads 10 -phred33 ../Barcoding/R1.split/Day3.GTAGG.fq ../Barcoding/R2.split/Day3.GTAGG.fq Day3/R1.D3.trim.paired.fq Day3/R1.D3.trim.unpaired.fq Day3/R2.D3.trim.paired.fq  Day3/R2.D3.trim.unpaired.fq ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

	__Output:__ Input Read Pairs: 34362137 Both Surviving: 9341134 (27.18%) Forward Only Surviving: 6986574 (20.33%) Reverse Only Surviving: 10229269 (29.77%) Dropped: 7805160 (22.71%)

unpaired reads were merged with the following

	cat R1.D3.trim.unpaired.fq R2.D3.trim.unpaired.fq > all.D3.trim.unpaired.fq

#### Day 4

	trimmomatic PE -threads 10 -phred33 ../Barcoding/R1.split/Day4.TGCCT.fq ../Barcoding/R2.split/Day4.TGCCT.fq Day4/R1.D4.trim.paired.fq Day4/R1.D4.trim.unpaired.fq Day4/R2.D4.trim.paired.fq  Day4/R2.D4.trim.unpaired.fq ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

	__Output:__ Input Read Pairs: 41629059 Both Surviving: 15829602 (38.03%) Forward Only Surviving: 9815036 (23.58%) Reverse Only Surviving: 9839762 (23.64%) Dropped: 6144659 (14.76%)

unpaired reads were merged with the following

	cat R1.D4.trim.unpaired.fq R2.D4.trim.unpaired.fq > all.D4.trim.unpaired.fq

## Mapping with BWA (Again)
To simplify the number of bam files that require merging, fastq files were merged with the following commands.

	cat Day1/R1.D1.trim.paired.fq Day2/R1.D2.trim.paired.fq Day3/R1.D3.trim.paired.fq Day4/R1.D4.trim.paired.fq > All_Days/R1.trim.paired.fq

	cat Day1/R2.D1.trim.paired.fq Day2/R2.D2.trim.paired.fq Day3/R2.D3.trim.paired.fq Day4/R2.D4.trim.paired.fq > All_Days/R2.trim.paired.fq

	cat Day1/all.D1.trim.unpaired.fq Day2/all.D2.trim.unpaired.fq Day3/all.D3.trim.unpaired.fq Day4/all.D4.trim.unpaired.fq > All_Days/all.trim.unpaired.fq

These concatenated files were then mapped with the following commands.

	bwa mem -t 8 S_lycopersicum_chromosomes.2.50.fa  ../Trimming/All_Days/R1.trim.paired.fq ../Trimming/All_Days/R2.trim.paired.fq | samtools view  -u - | samtools sort -m 2000000000 - | samtools rmdup -  tep1.PE.sort.rmdup.bam

	bwa mem -t 8 S_lycopersicum_chromosomes.2.50.fa  ../Trimming/Day1/R1.D1.trim.paired.fq ../Trimming/Day1/R2.D1.trim.paired.fq | samtools view  -u - | samtools sort -m 2000000000 -o D1.PE.sort.bam

	bwa mem -t 8 S_lycopersicum_chromosomes.2.50.fa  ../Trimming/Day2/R1.D2.trim.paired.fq ../Trimming/Day2/R2.D2.trim.paired.fq | samtools view  -u - | samtools sort -m 2000000000 -o D2.PE.sort.bam

	bwa mem -t 8 S_lycopersicum_chromosomes.2.50.fa  ../Trimming/All_Days/all.trim.unpaired.fq| samtools view  -bu - | samtools sort -m 2000000000 -o tep1.unpaired -

	`bwa mem -R '@RG\tID:Heinz\tSM:Heinz' S_lycopersicum_chromosomes.2.50.fa SRR404081.fastq.gz | samtools view -Sbu | samtools rmdup - HeinzG_tep1_rmdup_unsorted.bam`

	`samtools sort -m 8000000000 HeinzG_tep1_rmdup_unsorted.bam HeinzG_tep1_rmdup.bam`

	`cat R1.D1.trim.paired.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > Sorted/R1.D1.PE.trim.sort.fq`

# WTF

##### Day 1

	trimmomatic PE -threads 10 -phred33 ../Barcoding/R1.split/Day1.fq ../Barcoding/R2.split/Day1.fq R1.Day1.trim.paired.fq R1.Day1.trim.unpaired.fq R2.Day1.trim.paired.fq  R2.Day1.trim.unpaired.fq ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

	__output:__ Input Read Pairs: 38244690 Both Surviving: 20710085 (54.15%) Forward Only Surviving: 8075309 (21.11%) Reverse Only Surviving: 6663326 (17.42%) Dropped: 2795970 (7.31%)

	bwa mem -t 8 S_lycopersicum_chromosomes.2.50.fa  ../Trimming/R1.Day1.trim.paired.fq ../Trimming/R2.Day1.trim.paired.fq | samtools view  -u - | samtools sort -m 2000000000 - | samtools rmdup -  tep1.D1.PE.sort.rmdup.bam

	bwa mem -t 8 S_lycopersicum_chromosomes.2.50.fa  ../Barcoding/R1.split/Day1.fq ../Barcoding/R2.split/Day1.fq | samtools view  -u - | samtools sort -m 2000000000 - | samtools rmdup -  tep1.D1.PE.sort.rmdup.bam

##### Day 1 - 3rd time?

	auto_barcode/barcode_split_trim.pl --id tep1 -l -b epi1_barcode.list -o for.split JMAS010.1.fq
	auto_barcode/barcode_split_trim.pl --id tep1 -l -b epi1_barcode.list -o rev.split JMAS010.2.fq

	bwa mem -t 8 S_lycopersicum_chromosomes.2.50.fa  ../Barcoding/for.split/Day1.fq ../Barcoding/rev.split/Day1.fq | samtools view  -u - | samtools sort -m 2000000000 - | samtools rmdup -  tep1.D1.PE.sort.rmdup.bam

	~/bbmap/repair.sh in=../Barcoding/R1.split/Day1.fq in2=../Barcoding/R1.split/Day1.fq out=d1.r1.matched out2=d2.r2.matched

	bwa mem -t 8 S_lycopersicum_chromosomes.2.50.fa  ../Barcoding/bar.JMAS010.1.fq ../Barcoding/bar.JMAS010.2.fq | samtools view  -u - | samtools sort -m 2000000000 - | samtools rmdup -  tep1.PE.sort.rmdup.bam

	trimmomatic PE -threads 10 -phred33 ../Barcoding/bar.JMAS010.1.fq ../Barcoding/bar.JMAS010.2.fq bar.JMAS010.1.paired.fq bar.JMAS010.1.unpaired.fq bar.JMAS010.2.paired.fq  bar.JMAS010.2.unpaired.fq ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

	Input Read Pairs: 206012240 Both Surviving: 118942857 (57.74%) Forward Only Surviving: 11943836 (5.80%) Reverse Only Surviving: 8409328 (4.08%) Dropped: 66716219 (32.38%)

	bwa mem -t 8 S_lycopersicum_chromosomes.2.50.fa  ../Trimming/bar.JMAS010.1.paired.fq ../Trimming/bar.JMAS010.2.paired.fq | samtools view  -u - | samtools sort -m 2000000000 - | samtools rmdup -  tep1.PE.sort.rmdup.bam
