tomato elongated petiole (tep-1): Shade Avoidance Mutant
========================================================
author: Ellis Anderson
date: 5/30/2016
Background slides courtesy of Leslie Herrera, May 2nd 2015



What is Shade Avoidance?
========================================================
incremental: true

- Low Red:Far-Red light is indicative of shade
- Plants respond via
  - Internode and petiole elongation
  - Increased apical dominance
  - This is theorized to affect yield
  
***

![Shade Avoidance](./images/photo1)

Research Question
========================================================
type: sub-section

What genes in domesticated tomato are involved in shade avoidance?

Elongated Mutant Screen
========================================================

Wild-type: Standard for comparison

![Wild Type Plants](./images/photo2)

***

Mutant of interest: Consitutive elongation

![Mutant Plants](./images/photo3)

tomato elongated petiole (tep-1)
========================================================

![tep-1](./images/photo4.png)

Phenotype Characteristics
========================================================
type: sub-section

Extended Petioles
========================================================
Mutant of interest (tep-1) has petioles that are longer than the control (M82) by the 5th week of development

![Petiole](./images/photo5)

***
![M82 vs tep-1 - Petiole](./images/photo6)

Reduced Hypocotyl Length
========================================================
tep-1 mutant shows less hypocotyl elongation in response to shade

![M82 vs tep-1 - Hypocotyl](./images/photo7)

***
![Hypocotyl](./images/photo8.png)

Finding the Gene: Bulk Segregant Analysis
========================================================
type: sub-section

A bioinformatics approach to find the mutation responsible for the tep-1 phenotype

Summary
========================================================
![BSA1](./images/photo9)

Summary
========================================================
![BSA2](./images/photo10)

Tissue Collection From Elongated F2 Plants
========================================================
left: 70%

![tep and heinz F2s](./images/photo11)

***

![25% homozygous F2s](./images/photo12)

Sequencing and Computational Analysis
========================================================
left: 30%

- Filter low quality reads
- Map reads to reference genome
- Generate graphical output

***

![illumina seq](./images/photo13)

Graphical Output to Locate Gene of Interest
=========================================================

- Take advantage of different parental backgrounds to find region of interest
- Should ideally produce one spot with 100% M82 background

Initial Results
=========================================================

![Leslie's Results](./images/photo14)

My Work on tep-1
=========================================================
type: section
incremental: true

- Check over Leslie's code and reproduce her results
- Map Heinz genome to Heinz genome and add those results to the analysis
- Look for deletions from fast neutron mutagenesis in region of interest

Recreating Leslie's Results
=========================================================
Run through Leslie's code step by step

Results:
![Results 1](./images/photo15)


Proofreading R Script
=========================================================
Filtering based on genotype

```r
# Filter to keep rows where M82 = "1/1"
tep_snp_data <-subset(tep_data, tep.gt != "0/0")

tep_snp_data<-tep_snp_data[grep("TYPE=snp", tep_snp_data$INFO),]

tep_snp_data <- tep_snp_data[!is.na(tep_snp_data$tep.read.depth),]
```

Proofreading R Script
=========================================================
Filtering based on genotype

```r
# Filter to keep rows where M82 = "1/1"
tep_snp_data <-subset(tep_data, tep.gt != "0/0")

tep_snp_data <- tep_snp_data[grepl("1/1",tep_snp_data$M82.gt),]

tep_snp_data<-tep_snp_data[grep("TYPE=snp", tep_snp_data$INFO),]

tep_snp_data <- tep_snp_data[!is.na(tep_snp_data$tep.read.depth),]
```

New Results
=========================================================
![Results 2](./images/photo16)

Mapping Heinz to Heinz
=========================================================
- Used heinz reads from Lavelle tie1 experiment: `SRR404081.fastq.gz`
- Mapped these to `S_lycopersicum_chromosomes.2.50.fa`
- Reads were mapped using `bwa mem`

Merging Heinz file to Leslie's
=========================================================
- Used `picard-tools` functions to add read groups and combine all reads
- Used `freebayes` to produce a VCF file for analysis in R

R analysis
==========================================================
type: sub-section
Added additional filtering criteria to remove any Heinz-Heinz polymorphisms


```r
tep_snp_data <- tep_snp_data[grepl("0/0",tep_snp_data$Heinz.gt),]
```

Final Results:
==========================================================
Chromosome 9:
![Ch 09 - Final](./images/photo18)

Final Results:
==========================================================
All Chromosomes:
![All Chromosomes - Final](./images/photo17)

Next Steps
==========================================================
- Area of about 300kb identified as 100% M82
- Fast neutron mutagenesis often leads to large deletions
- Map reads again, but this time look for large deletions in the area of interest

Further Work
==========================================================
Use molecular biology techniques to ensure identified gene results in tep-1 phenotype

Special Thanks
==========================================================
type: sub-section
- Julin Maloof
- Cody Markelz
- Leslie Herrera & Amanda Schrager Lavelle
- Zamir Lab
- Maloof and Harmer Lab members


